package colormap

import (
	"image/color"
	"testing"
)

func TestSeuratColormapEndpoints(t *testing.T) {
	t.Parallel()

	c0, ok := Seurat.At(0).(color.RGBA)
	if !ok {
		t.Fatalf("expected color.RGBA at t=0")
	}
	if c0 != (color.RGBA{R: 211, G: 211, B: 211, A: 255}) {
		t.Fatalf("unexpected Seurat.At(0): %#v", c0)
	}

	c1, ok := Seurat.At(1).(color.RGBA)
	if !ok {
		t.Fatalf("expected color.RGBA at t=1")
	}
	if c1 != (color.RGBA{R: 255, G: 0, B: 0, A: 255}) {
		t.Fatalf("unexpected Seurat.At(1): %#v", c1)
	}
}

