//go:build soma

package soma

import (
	"fmt"
	"math"
	"os"
	"sync"

	tiledb "github.com/TileDB-Inc/TileDB-Go"
)

// Reader provides minimal SOMA reads via TileDB arrays.
type Reader struct {
	experimentURI string
	ctx           *tiledb.Context

	geneOnce sync.Once
	geneMap  map[string]int64 // gene_id -> gene soma_joinid
	geneErr  error

	obsIdxMu    sync.Mutex
	obsIdxCache map[string]map[string][]int64 // column -> value -> []cell_joinid
}

func NewReader(somaPath string) (*Reader, error) {
	uri, err := ResolveExperimentURI(somaPath)
	if err != nil {
		return nil, err
	}
	if _, statErr := os.Stat(uri); statErr != nil {
		return nil, fmt.Errorf("soma experiment not found at %s: %w", uri, statErr)
	}

	ctx, err := tiledb.NewContext(nil)
	if err != nil {
		return nil, fmt.Errorf("failed to create TileDB context: %w", err)
	}

	return &Reader{
		experimentURI: uri,
		ctx:           ctx,
	}, nil
}

func (r *Reader) Supported() bool { return true }

func (r *Reader) ExperimentURI() string { return r.experimentURI }

func (r *Reader) GeneJoinID(gene string) (int64, error) {
	r.geneOnce.Do(func() { r.geneErr = r.loadGeneMap() })
	if r.geneErr != nil {
		return 0, r.geneErr
	}
	id, ok := r.geneMap[gene]
	if !ok {
		return 0, fmt.Errorf("gene not found in SOMA var: %s", gene)
	}
	return id, nil
}

// ExpressionByCellJoinID reads expression values for one gene at given cell joinids.
// Returns a sparse map: cell_joinid -> value (only non-zero entries).
func (r *Reader) ExpressionByCellJoinID(gene string, cellJoinIDs []int64) (geneJoinID int64, values map[int64]float32, err error) {
	geneJoinID, err = r.GeneJoinID(gene)
	if err != nil {
		return 0, nil, err
	}
	if len(cellJoinIDs) == 0 {
		return geneJoinID, map[int64]float32{}, nil
	}

	xURI := r.experimentURI + "/ms/RNA/X/data"
	arr, err := tiledb.NewArray(r.ctx, xURI)
	if err != nil {
		return 0, nil, fmt.Errorf("failed to open X array (%s): %w", xURI, err)
	}
	defer arr.Free()
	if err := arr.Open(tiledb.TILEDB_READ); err != nil {
		return 0, nil, fmt.Errorf("failed to open X array for read: %w", err)
	}
	defer arr.Close()

	sub, err := arr.NewSubarray()
	if err != nil {
		return 0, nil, fmt.Errorf("failed to create subarray: %w", err)
	}
	defer sub.Free()

	// Filter to selected cells and one gene.
	// NOTE: TileDB subarray supports multiple ranges per dimension; for point queries we add [id,id] ranges.
	for _, cid := range cellJoinIDs {
		if err := sub.AddRangeByName("soma_dim_0", tiledb.MakeRange[int64](cid, cid)); err != nil {
			return 0, nil, fmt.Errorf("failed to add cell range: %w", err)
		}
	}
	if err := sub.AddRangeByName("soma_dim_1", tiledb.MakeRange[int64](geneJoinID, geneJoinID)); err != nil {
		return 0, nil, fmt.Errorf("failed to add gene range: %w", err)
	}

	q, err := tiledb.NewQuery(r.ctx, arr)
	if err != nil {
		return 0, nil, fmt.Errorf("failed to create query: %w", err)
	}
	defer q.Free()

	if err := q.SetSubarray(sub); err != nil {
		return 0, nil, fmt.Errorf("failed to set subarray: %w", err)
	}
	// For sparse reads, unordered is generally fine.
	_ = q.SetLayout(tiledb.TILEDB_UNORDERED)

	// Worst-case nnz for (cells subset) x (one gene) is len(cells).
	n := len(cellJoinIDs)
	outCell := make([]int64, n)
	outGene := make([]int64, n)
	outVal := make([]float32, n)
	valNullable, err := attributeNullable(arr, "soma_data")
	if err != nil {
		return 0, nil, fmt.Errorf("failed to inspect soma_data nullable: %w", err)
	}
	var outValValid []uint8
	if valNullable {
		outValValid = make([]uint8, n)
	}

	if _, err := q.SetDataBuffer("soma_dim_0", outCell); err != nil {
		return 0, nil, fmt.Errorf("failed to set buffer soma_dim_0: %w", err)
	}
	if _, err := q.SetDataBuffer("soma_dim_1", outGene); err != nil {
		return 0, nil, fmt.Errorf("failed to set buffer soma_dim_1: %w", err)
	}
	if _, err := q.SetDataBuffer("soma_data", outVal); err != nil {
		return 0, nil, fmt.Errorf("failed to set buffer soma_data: %w", err)
	}
	// If soma_data is nullable, validity buffer is required.
	if valNullable {
		if _, err := q.SetValidityBuffer("soma_data", outValValid); err != nil {
			return 0, nil, fmt.Errorf("failed to set validity buffer soma_data: %w", err)
		}
	}

	if err := q.Submit(); err != nil {
		return 0, nil, fmt.Errorf("query submit failed: %w", err)
	}
	status, err := q.Status()
	if err != nil {
		return 0, nil, fmt.Errorf("query status failed: %w", err)
	}
	if status != tiledb.TILEDB_COMPLETED && status != tiledb.TILEDB_INCOMPLETE {
		return 0, nil, fmt.Errorf("unexpected query status: %v", status)
	}

	elems, err := q.ResultBufferElements()
	if err != nil {
		return 0, nil, fmt.Errorf("failed to get result buffer elements: %w", err)
	}
	got := int(elems["soma_data"][1])
	if got > len(outVal) {
		got = len(outVal)
	}
	gotValid := 0
	if valNullable {
		gotValid = int(elems["soma_data"][2])
		if gotValid > len(outValValid) {
			gotValid = len(outValValid)
		}
	}

	values = make(map[int64]float32, got)
	for i := 0; i < got; i++ {
		// If validity was returned and marks null, skip.
		if valNullable && i < gotValid && outValValid[i] == 0 {
			continue
		}
		// Expect outGene[i] == geneJoinID; but we don't require it.
		values[outCell[i]] = outVal[i]
	}
	return geneJoinID, values, nil
}

func (r *Reader) loadGeneMap() error {
	varURI := r.experimentURI + "/ms/RNA/var"
	arr, err := tiledb.NewArray(r.ctx, varURI)
	if err != nil {
		return fmt.Errorf("failed to open var array (%s): %w", varURI, err)
	}
	defer arr.Free()
	if err := arr.Open(tiledb.TILEDB_READ); err != nil {
		return fmt.Errorf("failed to open var array for read: %w", err)
	}
	defer arr.Close()

	// Use non-empty domain to avoid relying on potentially unbounded dimension domains.
	ned, isEmpty, err := arr.NonEmptyDomainFromName("soma_joinid")
	if err != nil {
		return fmt.Errorf("failed to get var non-empty domain: %w", err)
	}
	if isEmpty || ned == nil {
		r.geneMap = map[string]int64{}
		return nil
	}
	minID, maxID, err := boundsMinMaxInt64(ned.Bounds)
	if err != nil {
		return fmt.Errorf("failed to parse var non-empty domain bounds: %w", err)
	}

	sub, err := arr.NewSubarray()
	if err != nil {
		return fmt.Errorf("failed to create var subarray: %w", err)
	}
	defer sub.Free()
	if err := sub.AddRangeByName("soma_joinid", tiledb.MakeRange[int64](minID, maxID)); err != nil {
		return fmt.Errorf("failed to set var range: %w", err)
	}

	q, err := tiledb.NewQuery(r.ctx, arr)
	if err != nil {
		return fmt.Errorf("failed to create var query: %w", err)
	}
	defer q.Free()
	if err := q.SetSubarray(sub); err != nil {
		return fmt.Errorf("failed to set var subarray: %w", err)
	}
	if err := q.SetLayout(tiledb.TILEDB_ROW_MAJOR); err != nil {
		return fmt.Errorf("failed to set var query layout: %w", err)
	}

	// Stream in chunks to avoid huge allocations and to handle unbounded domains safely.
	const chunkRows = 4096
	joinIDs := make([]int64, chunkRows)
	offsets := make([]uint64, chunkRows)
	geneNullable, err := attributeNullable(arr, "gene_id")
	if err != nil {
		return fmt.Errorf("failed to inspect gene_id nullable: %w", err)
	}
	var validity []uint8
	if geneNullable {
		validity = make([]uint8, chunkRows)
	}
	dataBytes := make([]byte, 1024*1024) // 1MB for var-length gene_id bytes

	m := make(map[string]int64, 32768)
	for {
		// Reset buffers each submit so TileDB sees full capacities (buffer sizes are in/out params).
		if _, err := q.SetDataBuffer("soma_joinid", joinIDs); err != nil {
			return fmt.Errorf("failed to set buffer soma_joinid: %w", err)
		}
		if _, err := q.SetOffsetsBuffer("gene_id", offsets); err != nil {
			return fmt.Errorf("failed to set offsets buffer gene_id: %w", err)
		}
		if _, err := q.SetDataBuffer("gene_id", dataBytes); err != nil {
			return fmt.Errorf("failed to set data buffer gene_id: %w", err)
		}
		if geneNullable {
			if _, err := q.SetValidityBuffer("gene_id", validity); err != nil {
				return fmt.Errorf("failed to set validity buffer gene_id: %w", err)
			}
		}

		if err := q.Submit(); err != nil {
			return fmt.Errorf("var query submit failed: %w", err)
		}
		status, err := q.Status()
		if err != nil {
			return fmt.Errorf("var query status failed: %w", err)
		}
		elems, err := q.ResultBufferElements()
		if err != nil {
			return fmt.Errorf("var query ResultBufferElements failed: %w", err)
		}

		usedJoin := int(elems["soma_joinid"][1])
		usedOffsets := int(elems["gene_id"][0])
		usedBytes := int(elems["gene_id"][1])
		usedValid := 0
		if geneNullable {
			usedValid = int(elems["gene_id"][2])
		}
		if usedJoin > len(joinIDs) {
			usedJoin = len(joinIDs)
		}
		if usedOffsets > len(offsets) {
			usedOffsets = len(offsets)
		}
		if usedBytes > len(dataBytes) {
			usedBytes = len(dataBytes)
		}
		if geneNullable {
			if usedValid > len(validity) {
				usedValid = len(validity)
			}
		}

		// If buffers are too small to return even a single row, grow and retry.
		if status == tiledb.TILEDB_INCOMPLETE && usedOffsets == 0 && usedBytes == 0 && usedJoin == 0 {
			if len(dataBytes) < 64*1024*1024 {
				dataBytes = make([]byte, len(dataBytes)*2)
				continue
			}
			return fmt.Errorf("var query buffers too small (gene_id); grew to %d bytes and still no progress", len(dataBytes))
		}

		join := joinIDs[:usedJoin]
		off := offsets[:usedOffsets]
		data := dataBytes[:usedBytes]
		var val []uint8
		if geneNullable {
			val = validity[:usedValid]
		}

		lim := usedJoin
		if usedOffsets < lim {
			lim = usedOffsets
		}
		if geneNullable && usedValid > 0 && usedValid < lim {
			lim = usedValid
		}
		for i := 0; i < lim; i++ {
			if geneNullable && usedValid > 0 && val[i] == 0 {
				continue
			}
			start := int(off[i])
			end := len(data)
			if i+1 < usedOffsets {
				end = int(off[i+1])
			}
			if start < 0 || end < start || end > len(data) {
				continue
			}
			g := string(data[start:end])
			if g != "" {
				m[g] = join[i]
			}
		}

		if status == tiledb.TILEDB_COMPLETED {
			r.geneMap = m
			return nil
		}
		if status != tiledb.TILEDB_INCOMPLETE {
			return fmt.Errorf("unexpected TileDB query status for var: %v", status)
		}
	}
}

func boundsMinMaxInt64(bounds interface{}) (int64, int64, error) {
	switch v := bounds.(type) {
	case []int64:
		if len(v) >= 2 {
			return v[0], v[1], nil
		}
	case []int32:
		if len(v) >= 2 {
			return int64(v[0]), int64(v[1]), nil
		}
	case []uint64:
		if len(v) >= 2 {
			if v[0] > math.MaxInt64 || v[1] > math.MaxInt64 {
				return 0, 0, fmt.Errorf("uint64 bounds exceed int64 range")
			}
			return int64(v[0]), int64(v[1]), nil
		}
	case []uint32:
		if len(v) >= 2 {
			return int64(v[0]), int64(v[1]), nil
		}
	}
	return 0, 0, fmt.Errorf("unsupported bounds type for non-empty domain")
}

func attributeNullable(arr *tiledb.Array, name string) (bool, error) {
	schema, err := arr.Schema()
	if err != nil {
		return false, err
	}
	defer schema.Free()
	attr, err := schema.AttributeFromName(name)
	if err != nil {
		return false, err
	}
	defer attr.Free()
	return attr.Nullable()
}

// GeneRange returns the min and max soma_joinid of genes in the var DataFrame.
func (r *Reader) GeneRange() (minID, maxID int64, err error) {
	varURI := r.experimentURI + "/ms/RNA/var"
	arr, err := tiledb.NewArray(r.ctx, varURI)
	if err != nil {
		return 0, 0, fmt.Errorf("failed to open var array: %w", err)
	}
	defer arr.Free()
	if err := arr.Open(tiledb.TILEDB_READ); err != nil {
		return 0, 0, fmt.Errorf("failed to open var array for read: %w", err)
	}
	defer arr.Close()

	ned, isEmpty, err := arr.NonEmptyDomainFromName("soma_joinid")
	if err != nil {
		return 0, 0, fmt.Errorf("failed to get var non-empty domain: %w", err)
	}
	if isEmpty || ned == nil {
		return 0, 0, nil
	}
	return boundsMinMaxInt64(ned.Bounds)
}

// AllGenes returns a map of gene_id -> soma_joinid for all genes.
func (r *Reader) AllGenes() (map[string]int64, error) {
	r.geneOnce.Do(func() { r.geneErr = r.loadGeneMap() })
	if r.geneErr != nil {
		return nil, r.geneErr
	}
	return r.geneMap, nil
}

// ScanXForCells streams through ms/RNA/X/data for the given cell joinids (all genes).
// For each non-zero entry, calls onRow(cellJoinID, geneJoinID, value).
// This processes in batches to handle large datasets.
func (r *Reader) ScanXForCells(cellJoinIDs []int64, onRow func(cell, gene int64, val float32)) error {
	if len(cellJoinIDs) == 0 {
		return nil
	}

	xURI := r.experimentURI + "/ms/RNA/X/data"
	arr, err := tiledb.NewArray(r.ctx, xURI)
	if err != nil {
		return fmt.Errorf("failed to open X array: %w", err)
	}
	defer arr.Free()
	if err := arr.Open(tiledb.TILEDB_READ); err != nil {
		return fmt.Errorf("failed to open X array for read: %w", err)
	}
	defer arr.Close()

	// Get gene range
	geneMin, geneMax, err := r.GeneRange()
	if err != nil {
		return fmt.Errorf("failed to get gene range: %w", err)
	}

	sub, err := arr.NewSubarray()
	if err != nil {
		return fmt.Errorf("failed to create X subarray: %w", err)
	}
	defer sub.Free()

	// Add cell ranges
	for _, cid := range cellJoinIDs {
		if err := sub.AddRangeByName("soma_dim_0", tiledb.MakeRange[int64](cid, cid)); err != nil {
			return fmt.Errorf("failed to add cell range: %w", err)
		}
	}
	// Add gene range (all genes)
	if err := sub.AddRangeByName("soma_dim_1", tiledb.MakeRange[int64](geneMin, geneMax)); err != nil {
		return fmt.Errorf("failed to add gene range: %w", err)
	}

	q, err := tiledb.NewQuery(r.ctx, arr)
	if err != nil {
		return fmt.Errorf("failed to create X query: %w", err)
	}
	defer q.Free()

	if err := q.SetSubarray(sub); err != nil {
		return fmt.Errorf("failed to set X subarray: %w", err)
	}
	_ = q.SetLayout(tiledb.TILEDB_UNORDERED)

	// Buffer size: expect sparse data, allocate reasonable buffer
	const bufSize = 1024 * 1024 // 1M entries
	outCell := make([]int64, bufSize)
	outGene := make([]int64, bufSize)
	outVal := make([]float32, bufSize)
	valNullable, err := attributeNullable(arr, "soma_data")
	if err != nil {
		return fmt.Errorf("failed to inspect soma_data nullable: %w", err)
	}
	var outValValid []uint8
	if valNullable {
		outValValid = make([]uint8, bufSize)
	}

	for {
		if _, err := q.SetDataBuffer("soma_dim_0", outCell); err != nil {
			return fmt.Errorf("failed to set buffer soma_dim_0: %w", err)
		}
		if _, err := q.SetDataBuffer("soma_dim_1", outGene); err != nil {
			return fmt.Errorf("failed to set buffer soma_dim_1: %w", err)
		}
		if _, err := q.SetDataBuffer("soma_data", outVal); err != nil {
			return fmt.Errorf("failed to set buffer soma_data: %w", err)
		}
		if valNullable {
			if _, err := q.SetValidityBuffer("soma_data", outValValid); err != nil {
				return fmt.Errorf("failed to set validity buffer soma_data: %w", err)
			}
		}

		if err := q.Submit(); err != nil {
			return fmt.Errorf("X query submit failed: %w", err)
		}
		status, err := q.Status()
		if err != nil {
			return fmt.Errorf("X query status failed: %w", err)
		}

		elems, err := q.ResultBufferElements()
		if err != nil {
			return fmt.Errorf("X query ResultBufferElements failed: %w", err)
		}
		got := int(elems["soma_data"][1])
		if got > len(outVal) {
			got = len(outVal)
		}
		gotValid := 0
		if valNullable {
			gotValid = int(elems["soma_data"][2])
			if gotValid > len(outValValid) {
				gotValid = len(outValValid)
			}
		}

		for i := 0; i < got; i++ {
			if valNullable && i < gotValid && outValValid[i] == 0 {
				continue
			}
			onRow(outCell[i], outGene[i], outVal[i])
		}

		if status == tiledb.TILEDB_COMPLETED {
			return nil
		}
		if status != tiledb.TILEDB_INCOMPLETE {
			return fmt.Errorf("unexpected X query status: %v", status)
		}
	}
}

// ObsColumns returns the list of attribute names in the obs DataFrame.
func (r *Reader) ObsColumns() ([]string, error) {
	obsURI := r.experimentURI + "/obs"
	arr, err := tiledb.NewArray(r.ctx, obsURI)
	if err != nil {
		return nil, fmt.Errorf("failed to open obs array: %w", err)
	}
	defer arr.Free()
	if err := arr.Open(tiledb.TILEDB_READ); err != nil {
		return nil, fmt.Errorf("failed to open obs array for read: %w", err)
	}
	defer arr.Close()

	schema, err := arr.Schema()
	if err != nil {
		return nil, fmt.Errorf("failed to get obs schema: %w", err)
	}
	defer schema.Free()

	nattrs, err := schema.AttributeNum()
	if err != nil {
		return nil, fmt.Errorf("failed to get attribute count: %w", err)
	}

	var columns []string
	for i := uint(0); i < nattrs; i++ {
		attr, err := schema.AttributeFromIndex(i)
		if err != nil {
			continue
		}
		name, err := attr.Name()
		attr.Free()
		if err != nil {
			continue
		}
		// Skip internal columns
		if name == "soma_joinid" {
			continue
		}
		columns = append(columns, name)
	}
	return columns, nil
}

// ObsColumnValues returns unique values for a string column in obs.
func (r *Reader) ObsColumnValues(column string) ([]string, error) {
	idx, err := r.ObsGroupIndex(column)
	if err != nil {
		return nil, err
	}
	values := make([]string, 0, len(idx))
	for v := range idx {
		values = append(values, v)
	}
	return values, nil
}

// ObsGroupIndex returns a map of column value -> cell joinids for a string column.
// Results are cached per column.
func (r *Reader) ObsGroupIndex(column string) (map[string][]int64, error) {
	r.obsIdxMu.Lock()
	defer r.obsIdxMu.Unlock()

	if r.obsIdxCache == nil {
		r.obsIdxCache = make(map[string]map[string][]int64)
	}
	if cached, ok := r.obsIdxCache[column]; ok {
		return cached, nil
	}

	idx, err := r.loadObsGroupIndex(column)
	if err != nil {
		return nil, err
	}
	r.obsIdxCache[column] = idx
	return idx, nil
}

func (r *Reader) loadObsGroupIndex(column string) (map[string][]int64, error) {
	obsURI := r.experimentURI + "/obs"
	arr, err := tiledb.NewArray(r.ctx, obsURI)
	if err != nil {
		return nil, fmt.Errorf("failed to open obs array: %w", err)
	}
	defer arr.Free()
	if err := arr.Open(tiledb.TILEDB_READ); err != nil {
		return nil, fmt.Errorf("failed to open obs array for read: %w", err)
	}
	defer arr.Close()

	// Check if column is varlen string
	schema, err := arr.Schema()
	if err != nil {
		return nil, fmt.Errorf("failed to get obs schema: %w", err)
	}
	defer schema.Free()

	attr, err := schema.AttributeFromName(column)
	if err != nil {
		return nil, fmt.Errorf("column not found in obs: %s", column)
	}
	defer attr.Free()

	// Get non-empty domain for soma_joinid
	ned, isEmpty, err := arr.NonEmptyDomainFromName("soma_joinid")
	if err != nil {
		return nil, fmt.Errorf("failed to get obs non-empty domain: %w", err)
	}
	if isEmpty || ned == nil {
		return map[string][]int64{}, nil
	}
	minID, maxID, err := boundsMinMaxInt64(ned.Bounds)
	if err != nil {
		return nil, fmt.Errorf("failed to parse obs non-empty domain: %w", err)
	}

	sub, err := arr.NewSubarray()
	if err != nil {
		return nil, fmt.Errorf("failed to create obs subarray: %w", err)
	}
	defer sub.Free()
	if err := sub.AddRangeByName("soma_joinid", tiledb.MakeRange[int64](minID, maxID)); err != nil {
		return nil, fmt.Errorf("failed to set obs range: %w", err)
	}

	q, err := tiledb.NewQuery(r.ctx, arr)
	if err != nil {
		return nil, fmt.Errorf("failed to create obs query: %w", err)
	}
	defer q.Free()
	if err := q.SetSubarray(sub); err != nil {
		return nil, fmt.Errorf("failed to set obs subarray: %w", err)
	}
	if err := q.SetLayout(tiledb.TILEDB_ROW_MAJOR); err != nil {
		return nil, fmt.Errorf("failed to set obs query layout: %w", err)
	}

	// Stream in chunks
	const chunkRows = 8192
	joinIDs := make([]int64, chunkRows)
	offsets := make([]uint64, chunkRows)
	colNullable, _ := attributeNullable(arr, column)
	var validity []uint8
	if colNullable {
		validity = make([]uint8, chunkRows)
	}
	dataBytes := make([]byte, 2*1024*1024) // 2MB for var-length column bytes

	result := make(map[string][]int64)
	for {
		if _, err := q.SetDataBuffer("soma_joinid", joinIDs); err != nil {
			return nil, fmt.Errorf("failed to set buffer soma_joinid: %w", err)
		}
		if _, err := q.SetOffsetsBuffer(column, offsets); err != nil {
			return nil, fmt.Errorf("failed to set offsets buffer %s: %w", column, err)
		}
		if _, err := q.SetDataBuffer(column, dataBytes); err != nil {
			return nil, fmt.Errorf("failed to set data buffer %s: %w", column, err)
		}
		if colNullable {
			if _, err := q.SetValidityBuffer(column, validity); err != nil {
				return nil, fmt.Errorf("failed to set validity buffer %s: %w", column, err)
			}
		}

		if err := q.Submit(); err != nil {
			return nil, fmt.Errorf("obs query submit failed: %w", err)
		}
		status, err := q.Status()
		if err != nil {
			return nil, fmt.Errorf("obs query status failed: %w", err)
		}
		elems, err := q.ResultBufferElements()
		if err != nil {
			return nil, fmt.Errorf("obs query ResultBufferElements failed: %w", err)
		}

		usedJoin := int(elems["soma_joinid"][1])
		usedOffsets := int(elems[column][0])
		usedBytes := int(elems[column][1])
		usedValid := 0
		if colNullable {
			usedValid = int(elems[column][2])
		}
		if usedJoin > len(joinIDs) {
			usedJoin = len(joinIDs)
		}
		if usedOffsets > len(offsets) {
			usedOffsets = len(offsets)
		}
		if usedBytes > len(dataBytes) {
			usedBytes = len(dataBytes)
		}
		if colNullable && usedValid > len(validity) {
			usedValid = len(validity)
		}

		// Grow buffer if needed
		if status == tiledb.TILEDB_INCOMPLETE && usedOffsets == 0 && usedBytes == 0 && usedJoin == 0 {
			if len(dataBytes) < 64*1024*1024 {
				dataBytes = make([]byte, len(dataBytes)*2)
				continue
			}
			return nil, fmt.Errorf("obs query buffers too small for column %s", column)
		}

		join := joinIDs[:usedJoin]
		off := offsets[:usedOffsets]
		data := dataBytes[:usedBytes]
		var val []uint8
		if colNullable {
			val = validity[:usedValid]
		}

		lim := usedJoin
		if usedOffsets < lim {
			lim = usedOffsets
		}
		if colNullable && usedValid > 0 && usedValid < lim {
			lim = usedValid
		}
		for i := 0; i < lim; i++ {
			if colNullable && usedValid > 0 && val[i] == 0 {
				continue
			}
			start := int(off[i])
			end := len(data)
			if i+1 < usedOffsets {
				end = int(off[i+1])
			}
			if start < 0 || end < start || end > len(data) {
				continue
			}
			v := string(data[start:end])
			if v != "" {
				result[v] = append(result[v], join[i])
			}
		}

		if status == tiledb.TILEDB_COMPLETED {
			return result, nil
		}
		if status != tiledb.TILEDB_INCOMPLETE {
			return nil, fmt.Errorf("unexpected TileDB query status for obs: %v", status)
		}
	}
}


