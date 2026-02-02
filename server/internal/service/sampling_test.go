// Package service provides business logic for the tile server.
package service

import (
	"testing"
)

func TestCellHash(t *testing.T) {
	// Test that hash is deterministic
	seed := int64(42)
	joinID := int64(12345)

	hash1 := cellHash(seed, joinID)
	hash2 := cellHash(seed, joinID)

	if hash1 != hash2 {
		t.Errorf("cellHash is not deterministic: %d != %d", hash1, hash2)
	}

	// Test that different seeds produce different hashes
	hash3 := cellHash(43, joinID)
	if hash1 == hash3 {
		t.Errorf("Different seeds should produce different hashes: %d == %d", hash1, hash3)
	}

	// Test that different joinIDs produce different hashes
	hash4 := cellHash(seed, 12346)
	if hash1 == hash4 {
		t.Errorf("Different joinIDs should produce different hashes: %d == %d", hash1, hash4)
	}
}

func TestDeterministicSample(t *testing.T) {
	// Create test cells
	cells := make([]CellInfo, 100)
	for i := 0; i < 100; i++ {
		cells[i] = CellInfo{
			JoinID: int64(i),
			X:      float32(i),
			Y:      float32(i * 2),
		}
	}

	// Test that sampling is deterministic
	seed := int64(42)
	k := 10

	sample1 := deterministicSample(cells, k, seed)
	sample2 := deterministicSample(cells, k, seed)

	if len(sample1) != len(sample2) {
		t.Errorf("Sample lengths differ: %d != %d", len(sample1), len(sample2))
	}

	for i := range sample1 {
		if sample1[i].JoinID != sample2[i].JoinID {
			t.Errorf("Samples differ at index %d: %d != %d", i, sample1[i].JoinID, sample2[i].JoinID)
		}
	}

	// Test that different seeds produce different samples
	sample3 := deterministicSample(cells, k, 43)
	sameCount := 0
	for i := range sample1 {
		if sample1[i].JoinID == sample3[i].JoinID {
			sameCount++
		}
	}
	// With k=10 and different seeds, it's extremely unlikely all 10 would be the same
	if sameCount == k {
		t.Errorf("Different seeds should produce different samples")
	}

	// Test that sample size is correct
	if len(sample1) != k {
		t.Errorf("Sample size should be %d, got %d", k, len(sample1))
	}

	// Test with k > len(cells) - should return all cells
	largeSample := deterministicSample(cells, 200, seed)
	if len(largeSample) != len(cells) {
		t.Errorf("Sample with k > len(cells) should return all cells: %d != %d", len(largeSample), len(cells))
	}

	// Test with k = 0 - should return empty
	emptySample := deterministicSample(cells, 0, seed)
	if len(emptySample) != 0 {
		t.Errorf("Sample with k=0 should return empty slice, got %d", len(emptySample))
	}
}

func TestDeterministicSampleStability(t *testing.T) {
	// Test that the same bbox+seed+limit produces stable results across calls
	cells := make([]CellInfo, 1000)
	for i := 0; i < 1000; i++ {
		cells[i] = CellInfo{
			JoinID:   int64(i * 7 + 13), // Non-sequential IDs
			X:        float32(i) * 0.1,
			Y:        float32(i) * 0.2,
			Category: "type_" + string(rune('A'+i%5)),
		}
	}

	seed := int64(0) // Default seed
	k := 50

	// Run multiple times
	var firstSample []CellInfo
	for run := 0; run < 5; run++ {
		sample := deterministicSample(cells, k, seed)
		if firstSample == nil {
			firstSample = sample
		} else {
			for i := range sample {
				if sample[i].JoinID != firstSample[i].JoinID {
					t.Errorf("Run %d: sample differs at index %d: %d != %d",
						run, i, sample[i].JoinID, firstSample[i].JoinID)
				}
			}
		}
	}
}

func BenchmarkDeterministicSample(b *testing.B) {
	// Create a large dataset similar to what we'd see in production
	cells := make([]CellInfo, 50000)
	for i := 0; i < 50000; i++ {
		cells[i] = CellInfo{
			JoinID: int64(i),
			X:      float32(i) * 0.01,
			Y:      float32(i) * 0.02,
		}
	}

	seed := int64(0)
	k := 5000

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = deterministicSample(cells, k, seed)
	}
}
