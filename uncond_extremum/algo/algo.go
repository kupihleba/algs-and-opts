package algo

import (
	"math"
)

func SvannMethod(xStart float64, step float64, f func(float64) float64) Interval {
	fun_1 := f(xStart - step)
	fun_2 := f(xStart)
	fun_3 := f(xStart + step)
	interval := Interval{
		xStart - step,
		xStart,
	}
	xs := []float64{xStart}
	if fun_1 >= fun_2 && fun_2 <= fun_3 {
		return interval
	} else if fun_1 <= fun_2 && fun_2 >= fun_3 {
		panic("Interval invalid!")
	}
	var Δ float64
	k := 1
	if fun_1 >= fun_2 && fun_2 >= fun_3 {
		Δ = step
		interval.Start = xs[0]
		xs = append(xs, xStart+step)
	} else if fun_1 <= fun_2 && fun_2 <= fun_3 {
		Δ = -step
		interval.End = xs[0]
		xs = append(xs, xStart-step)
	}
	for {
		for len(xs) <= k+1 {
			xs = append(xs, 0.0)
		}
		xs[k+1] = xs[k] + math.Pow(2.0, float64(k))*Δ
		if f(xs[k+1]) >= f(xs[k]) {
			if Δ > 0 {
				interval.End = xs[k+1]
			} else if Δ < 0 {
				interval.Start = xs[k+1]
			}
		} else {
			if Δ > 0 {
				interval.Start = xs[k]
			} else if Δ < 0 {
				interval.End = xs[k]
			}
		}
		if f(xs[k+1]) >= f(xs[k]) {
			break
		}
		k++
	}
	return interval
}

func BisectionMethod(ε float64, interval Interval, f func(float64) float64) (float64, Info) {
	xMid := interval.Center()
	it := 0
	for true {
		it++
		leftMid := interval.Start + interval.Length()/4
		rightMid := interval.End - interval.Length()/4

		if f(leftMid) < f(xMid) {
			interval.End = xMid
			xMid = leftMid
		} else if f(rightMid) < f(xMid) {
			interval.Start = xMid
			xMid = rightMid
		} else {
			interval.Start = leftMid
			interval.End = rightMid
		}
		if interval.Length() <= ε {
			break
		}
	}
	return xMid, Info{it}
}

func GoldenSectionMethod(ε float64, interval Interval, f func(float64) float64) (float64, Info) {
	phi := (1 + math.Sqrt(5.0)) / 2
	it := 0
	for interval.Length() > ε {
		it++
		a := interval.End - interval.Length()/phi
		b := interval.Start + interval.Length()/phi
		if f(b) <= f(a) {
			interval.Start = a
		} else {
			interval.End = b
		}
	}
	return interval.Center(), Info{it}
}

func FibonacciMethod(ε float64, interval Interval, f func(float64) float64) (float64, Info) {
	it := 0
	fib := []float64{1.0, 1.0}
	for lst := 1.0; lst < interval.Length()/ε; lst = fib[len(fib)-1] {
		fib = append(fib, lst+fib[len(fib)-2])
	}

	n := len(fib)
	for i := 1; i < n-4; i++ {
		a := interval.Start + fib[n-i-2]/fib[n-i]*interval.Length()
		b := interval.Start + fib[n-i-1]/fib[n-i]*interval.Length()
		it++
		if f(a) <= f(b) {
			interval.End = b
		} else {
			interval.Start = a
		}
	}
	return interval.Center(), Info{it}
}

func QuadraticInterpolation(
	ε float64,
	Δ float64,
	xStart float64,
	step float64,
	f func(float64) float64,
) (float64, Info) {
	pts := make([]Coord, 3)
	pts[0].X = xStart
	it := 0
	for true {
		pts[1].X = pts[0].X + step
		pts[0].Y, pts[1].Y = f(pts[0].X), f(pts[1].X)
		pts[2].X = pts[0].X + 2*step*Sign(pts[0].Y-pts[1].Y)
		for true {
			it++
			pts[2].Y = f(pts[2].X)
			_, i := Min(pts[0].Y, pts[1].Y, pts[2].Y)
			minCoord := Coord{
				X: pts[i].X,
				Y: pts[i].Y,
			}
			detRes := Det(pts)
			if detRes == 0 {
				pts[0].X = minCoord.X
			} else {
				sqrXs := make([]Coord, 3)
				copy(sqrXs, pts)
				for i := 0; i < len(sqrXs); i++ {
					sqrXs[i].X = math.Pow(sqrXs[i].X, 2) // To the power of 2
				}
				x := Det(sqrXs) / detRes / 2
				if math.Abs((minCoord.Y-f(x))/f(x)) < ε && math.Abs((minCoord.X-x)/x) < Δ {
					return pts[0].X, Info{it}
				}
				if pts[0].X <= x && x <= pts[2].X {
					if x < pts[1].X {
						pts[2].X, pts[1].X = pts[1].X, x
					} else {
						pts[0].X, pts[1].X = pts[1].X, x
					}
				} else {
					pts[0].X = x
					break
				}
			}
		}
	}
	panic("Function ended unexpectedly")
}

func CubicInterpolation(
	ε float64,
	Δ float64,
	xStart float64,
	step float64,
	f func(float64) float64,
	derivF func(float64) float64,
) (float64, Info) {
	pts := make([]Coord, 2)
	pts[0].X, pts[1].X = xStart, xStart
	it := 0
	m := 0.0
	df := derivF(xStart)
	var x float64
	for {
		pts[0] = pts[1]
		pts[1].X += math.Pow(m, 2) * step * Sign(-df)
		m += 1.0
		if derivF(pts[0].X) * derivF(pts[1].X) <= 0.0 {
			break
		}
	}

	for {
		it++
		pts[0].Y = f(pts[0].X)
		pts[1].Y = f(pts[1].X)
		df0 := derivF(pts[0].X)
		df1 := derivF(pts[1].X)

		a := 3 * (pts[0].Y - pts[1].Y) / (pts[1].X - pts[0].X) + df0 + df1
		b := math.Sqrt(a * a - df0 * df1)* Sign(pts[1].X - pts[0].X)
		μ := (df1 + b - a) / (df1 - df0 + 2 * b)
		if μ < 0 {
			x = pts[1].X
		} else if 0 < μ && μ <= 1.0 {
			x = pts[1].X - μ* (pts[1].X - pts[0].X)
		} else {
			x = pts[0].X
		}

		for f(x) > f(pts[0].X) && math.Abs((x - pts[0].X) / x) > Δ {
			x -= (x - pts[0].X) / 2
		}

		if math.Abs(derivF(x)) <= ε && math.Abs((x - pts[0].X)/x) <= Δ {
			return x, Info{it}
		} else {
			if derivF(x) * derivF(pts[0].X) <= 0 {
				pts[1].X = pts[0].X
			}
			pts[0].X = x
		}

	}
}