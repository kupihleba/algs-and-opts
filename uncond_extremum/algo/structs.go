package algo

type Interval struct {
	Start float64
	End   float64
}

type Info struct {
	Iterations int
}

type Coord struct {
	X float64
	Y float64
}

func (i *Interval) Length() float64 {
	return i.End - i.Start
}
func (i *Interval) Center() float64 {
	return (i.Start + i.End) / 2
}

func (i *Interval) Contains(x float64) bool {
	return x > i.Start && x < i.End
}
