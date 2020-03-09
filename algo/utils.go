package algo

import "math"

func Sign(x float64) float64 {
	return x / math.Abs(x)
}

func Min(xs ...float64) (float64, int) {
	minX := xs[0]
	minIndex := 0
	for i := 1; i < len(xs); i++ {
		if xs[i] < minX {
			minX = xs[i]
			minIndex = i
		}
	}
	return minX, minIndex
}

func Det(coords [] Coord) float64 {
	if len(coords) != 3 {
		panic("Not supported!")
	}
	return 2 *
		((coords[1].X-coords[2].X)*coords[0].Y +
			(coords[2].X-coords[0].X)*coords[1].Y +
			(coords[0].X-coords[1].X)*coords[2].Y)
}

