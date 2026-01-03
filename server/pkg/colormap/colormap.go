// Package colormap provides color schemes for visualization.
package colormap

import (
	"image/color"
)

// Colormap maps normalized values [0, 1] to colors.
type Colormap interface {
	At(t float64) color.Color
	AtIndex(i int) color.Color
}

// LinearColormap is a linear interpolation colormap.
type LinearColormap struct {
	colors []color.RGBA
}

// At returns the color at position t (0-1).
func (c LinearColormap) At(t float64) color.Color {
	if t <= 0 {
		return c.colors[0]
	}
	if t >= 1 {
		return c.colors[len(c.colors)-1]
	}

	idx := t * float64(len(c.colors)-1)
	lower := int(idx)
	upper := lower + 1
	if upper >= len(c.colors) {
		upper = len(c.colors) - 1
	}

	frac := idx - float64(lower)
	return interpolate(c.colors[lower], c.colors[upper], frac)
}

// AtIndex returns color at index i (wraps around).
func (c LinearColormap) AtIndex(i int) color.Color {
	return c.colors[i%len(c.colors)]
}

func interpolate(c1, c2 color.RGBA, t float64) color.RGBA {
	return color.RGBA{
		R: uint8(float64(c1.R) + t*(float64(c2.R)-float64(c1.R))),
		G: uint8(float64(c1.G) + t*(float64(c2.G)-float64(c1.G))),
		B: uint8(float64(c1.B) + t*(float64(c2.B)-float64(c1.B))),
		A: 255,
	}
}

// Viridis colormap (matplotlib viridis)
var Viridis = LinearColormap{
	colors: []color.RGBA{
		{68, 1, 84, 255},
		{72, 35, 116, 255},
		{64, 67, 135, 255},
		{52, 94, 141, 255},
		{41, 120, 142, 255},
		{32, 144, 140, 255},
		{34, 167, 132, 255},
		{68, 190, 112, 255},
		{121, 209, 81, 255},
		{189, 222, 38, 255},
		{253, 231, 37, 255},
	},
}

// Plasma colormap
var Plasma = LinearColormap{
	colors: []color.RGBA{
		{13, 8, 135, 255},
		{75, 3, 161, 255},
		{125, 3, 168, 255},
		{168, 34, 150, 255},
		{203, 70, 121, 255},
		{229, 107, 93, 255},
		{248, 148, 65, 255},
		{253, 195, 40, 255},
		{240, 249, 33, 255},
	},
}

// Inferno colormap
var Inferno = LinearColormap{
	colors: []color.RGBA{
		{0, 0, 4, 255},
		{40, 11, 84, 255},
		{101, 21, 110, 255},
		{159, 42, 99, 255},
		{212, 72, 66, 255},
		{245, 125, 21, 255},
		{250, 193, 39, 255},
		{252, 255, 164, 255},
	},
}

// Magma colormap
var Magma = LinearColormap{
	colors: []color.RGBA{
		{0, 0, 4, 255},
		{28, 16, 68, 255},
		{79, 18, 123, 255},
		{129, 37, 129, 255},
		{181, 54, 122, 255},
		{229, 80, 100, 255},
		{251, 135, 97, 255},
		{254, 194, 135, 255},
		{252, 253, 191, 255},
	},
}

// CategoricalColormap provides distinct colors for categories.
type CategoricalColormap struct {
	colors []color.RGBA
}

// At returns color at position t.
func (c CategoricalColormap) At(t float64) color.Color {
	idx := int(t * float64(len(c.colors)))
	if idx >= len(c.colors) {
		idx = len(c.colors) - 1
	}
	return c.colors[idx]
}

// AtIndex returns color at index.
func (c CategoricalColormap) AtIndex(i int) color.Color {
	return c.colors[i%len(c.colors)]
}

// Categorical colormap with 20 distinct colors
var Categorical = CategoricalColormap{
	colors: []color.RGBA{
		{31, 119, 180, 255},   // Blue
		{255, 127, 14, 255},   // Orange
		{44, 160, 44, 255},    // Green
		{214, 39, 40, 255},    // Red
		{148, 103, 189, 255},  // Purple
		{140, 86, 75, 255},    // Brown
		{227, 119, 194, 255},  // Pink
		{127, 127, 127, 255},  // Gray
		{188, 189, 34, 255},   // Olive
		{23, 190, 207, 255},   // Cyan
		{174, 199, 232, 255},  // Light blue
		{255, 187, 120, 255},  // Light orange
		{152, 223, 138, 255},  // Light green
		{255, 152, 150, 255},  // Light red
		{197, 176, 213, 255},  // Light purple
		{196, 156, 148, 255},  // Light brown
		{247, 182, 210, 255},  // Light pink
		{199, 199, 199, 255},  // Light gray
		{219, 219, 141, 255},  // Light olive
		{158, 218, 229, 255},  // Light cyan
	},
}
