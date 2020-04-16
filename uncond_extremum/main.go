package main

import (
	"fmt"
	"math"
	"modeling/lab3_mth_opt/algo"
)

func f5(x float64) float64 {
	//return 10*math.Pow(math.E, -x*x) + x*x + math.Abs(x-5)
	return math.Pow(x-1, 2) + math.Pow(x+5, 2)
}

func f1(x float64) float64 {
	return -12*x*x - 2*x + math.Abs(x-3)
}

func df1(x float64) float64 {
	return (x-3)/math.Abs(x-3) - 24*x - 2
}

func df5(x float64) float64 {
	//return -20*math.Pow(math.E, -x*x) + (x-5)/math.Abs(x-5) + 2*x
	return 4 * (x + 2)
}

func main() {
	start := -10.0 // f5
	//start := -3.0 // f1
	eps := 0.01
	step := 0.01
	interval := algo.SvannMethod(start, step, f5)
	bis, inf := algo.BisectionMethod(eps, interval, f5)
	fmt.Printf("Bisection Method:\t%v\t(%d)\n", bis, inf.Iterations)

	gold, inf := algo.GoldenSectionMethod(eps, interval, f5)
	fmt.Printf("Golden Section Method:\t%v\t(%d)\n", gold, inf.Iterations)

	fib, inf := algo.FibonacciMethod(eps, interval, f5)
	fmt.Printf("Fibonacci Method:\t%v\t(%d)\n", fib, inf.Iterations)

	delta := 0.01
	qad, inf := algo.QuadraticInterpolation(eps, delta, -start, step, f5)
	fmt.Printf("Quadratic Interpolation:\t%v\t(%d)\n", qad, inf.Iterations)

	cub, inf := algo.CubicInterpolation(eps, delta, -start, step, f5, df5)
	fmt.Printf("Cubic Interpolation:\t%v\t(%d)\n", cub, inf.Iterations)
}
