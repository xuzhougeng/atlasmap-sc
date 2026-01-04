package zarr

import (
	"os"
	"path/filepath"
	"testing"
)

func testZarrPath(t *testing.T) string {
	t.Helper()

	wd, err := os.Getwd()
	if err != nil {
		t.Fatalf("failed to get working dir: %v", err)
	}

	// Package directory: server/internal/data/zarr
	// Repo root: ../../../../
	p := filepath.Clean(filepath.Join(wd, "../../../../data/preprocessed/zarr/bins.zarr"))
	if _, err := os.Stat(p); os.IsNotExist(err) {
		t.Skipf("test data not found at %s, skipping", p)
	}
	return p
}

func TestReader_GetGeneStats_Zoom0(t *testing.T) {
	r, err := NewReader(testZarrPath(t))
	if err != nil {
		t.Fatalf("failed to create reader: %v", err)
	}
	defer r.Close()

	genes := r.Metadata().Genes
	if len(genes) == 0 {
		t.Skip("no genes in metadata")
	}
	gene := genes[0]
	if _, ok := r.Metadata().GeneIndex["ENSG00000167286"]; ok {
		gene = "ENSG00000167286"
	}

	stats, err := r.GetGeneStats(gene, 0)
	if err != nil {
		t.Fatalf("GetGeneStats(%q, 0) error: %v", gene, err)
	}
	if stats.Gene != gene {
		t.Fatalf("unexpected gene: got %q want %q", stats.Gene, gene)
	}
	if stats.TotalBins <= 0 {
		t.Fatalf("expected positive total_bins, got %d", stats.TotalBins)
	}
}

func TestReader_GetExpression_MultiChunk(t *testing.T) {
	r, err := NewReader(testZarrPath(t))
	if err != nil {
		t.Fatalf("failed to create reader: %v", err)
	}
	defer r.Close()

	genes := r.Metadata().Genes
	if len(genes) == 0 {
		t.Skip("no genes in metadata")
	}
	gene := genes[0]
	geneIdx, ok := r.Metadata().GeneIndex[gene]
	if !ok {
		t.Fatalf("gene %q missing from gene_index", gene)
	}

	// Pick a zoom level that is chunked in the bundled dataset.
	zoom := r.Metadata().ZoomLevels - 1
	if zoom < 0 {
		t.Skip("no zoom levels")
	}

	expr, err := r.GetExpression(zoom, geneIdx, "mean")
	if err != nil {
		t.Fatalf("GetExpression(zoom=%d,gene=%q) error: %v", zoom, gene, err)
	}
	if len(expr) == 0 {
		t.Fatalf("expected non-empty expression vector at zoom %d", zoom)
	}
}
