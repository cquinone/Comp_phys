package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"
	"os"
	"strconv"
)


func main() {
	//T, _ := strconv.ParseFloat(os.Args[2], 64)
	N, _ := strconv.Atoi(os.Args[1])

	for k := 1; k < 40; k++ {
		rand.Seed(time.Now().UTC().UnixNano())
		var J float64 = 1.0
		var h float64 = 0.0
		t := float64(k)/10.0
		var num_steps int = 2000
		initial_state := two_d(N)
		// randomly set spins of initial configuration
		var r float64
		for i := 0; i < N; i++ {
			for j := 0; j < N; j ++ {
				r = rand.Float64()
				if r <= .5 {
					initial_state[i][j] = 1
				} else {
					initial_state[i][j] = -1
				}
			}
		}
		//fmt.Println("INITIAL E: ", hamil(J,h,initial_state))
		energies := evolve(initial_state,J,h,t,num_steps)
		fmt.Println(t, mean(energies[1000:])/float64(N*N))
	}
}


// actually evolve matrix
func evolve(state [][]float64, J float64, h float64, T float64, num_steps int) []float64 {
	// loop over "MC STEPS", i.e. a full NXN sweep of possible flips	
	energies := make([]float64, num_steps)
	e_old := hamil(J,h, state)
	N := len(state[0])
	var e_new float64
	var r float64
	test_state := two_d(N)
	for step := 0; step < num_steps; step++ {
		// one full sweep per loop
		for i := 0; i < N; i ++ {
			for j := 0; j < N; j++ {
				test_state = set_matrix(test_state, state)
				test_state[i][j] = -1*test_state[i][j]
				e_new = hamil(J,h, test_state)	
				if e_new - e_old > 0 {   // if the energy change is positive (not towards minimization)
					r = rand.Float64()
					if r < 1.0/(math.Exp((e_new-e_old)/T)) {
						state = set_matrix(state, test_state)
					}
				} else {
					state = set_matrix(state, test_state)   // accept if lowers energy?
				}
				e_old = hamil(J,h, state)
			}
		}
		energies[step] = hamil(J,h, state)
		//fmt.Printf("\rON STEP, E: %d, %f", step, hamil(J,h, state))
	}
	return energies
}

// calculate PBC energy
func hamil(J float64, h float64, state [][]float64) float64 {
	total_E := 0.0
	var right int
	var below int
	N := len(state[0])
	for i := 0; i < N; i ++ {
		for j := 0; j < N; j++ {
			// grab one to right and one below, i is y (going down), j is x, going right
			// grabbing just these means never double counting? --> because PBC
			if i == N-1 {
				below = 0
			} else if  i < N-1 {
				below = i+1
			}
			if j == N-1 {
				right = 0
			} else if j < N-1 {
				right = j+1
			}
			total_E = total_E + -J*state[i][j]*(state[below][j]+state[i][right]) - h*state[i][j]
		}
	}
	return total_E
}


// Go has no math.mean() :(
func mean(data []float64) float64 {
	var avg float64
	var N int = len(data)
	for i := 0; i < N; i ++ {
		avg = avg + data[i]
	}
	return avg/float64(N)
}


// sets up 2-d matrix of floats
func two_d(N int) [][]float64 {
	matrix := make([][]float64, N)
	for i := range matrix {
		matrix[i] = make([]float64, N)
	}
	return matrix
}


// print matrix nicely
func print_m(m [][]float64) {
	length := len(m[0])
	for i := 0; i < length; i ++ {
		fmt.Println(m[i])
	}
}


// helper function to set values of 2d matrix equal to other
func set_matrix(mat_1 [][]float64, mat_2 [][]float64) [][]float64 {
	for i := 0; i < len(mat_1); i++ {
		for j := 0; j < len(mat_1); j++ {
			mat_1[i][j] = mat_2[i][j]
		}
	}
	return mat_1
}