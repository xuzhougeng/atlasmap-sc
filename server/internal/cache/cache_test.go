package cache

import "testing"

func TestCategoryTileKey(t *testing.T) {
	base := "cat:0/1/2:cell_type:ps=1.000"

	t.Run("nilFilters", func(t *testing.T) {
		got := CategoryTileKey(0, 1, 2, "cell_type", nil, 1.0)
		if got != base {
			t.Fatalf("expected %q, got %q", base, got)
		}
	})

	t.Run("emptyFilters", func(t *testing.T) {
		got := CategoryTileKey(0, 1, 2, "cell_type", []string{}, 1.0)
		want := base + ":none"
		if got != want {
			t.Fatalf("expected %q, got %q", want, got)
		}
	})

	t.Run("sortedFilters", func(t *testing.T) {
		key1 := CategoryTileKey(0, 1, 2, "cell_type", []string{"B", "A"}, 1.0)
		key2 := CategoryTileKey(0, 1, 2, "cell_type", []string{"A", "B"}, 1.0)
		if key1 != key2 {
			t.Fatalf("expected stable key, got %q vs %q", key1, key2)
		}
		if key1 == base || key1 == base+":none" {
			t.Fatalf("expected filtered key to differ from base, got %q", key1)
		}
	})
}
