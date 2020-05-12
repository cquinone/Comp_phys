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
	// parameters
	var num_steps int = 10000
	N, _ := strconv.Atoi(os.Args[1])

	// store calculated data
	energies := make([]float64, num_steps)
	mags := make([]float64, num_steps)
	g := make([][]float64, num_steps/2)

	//open file to write spin(t) to (ovito simulation stuff)
	ovito_file, err := os.Create("spins"+"_"+strconv.Itoa(N*N)+".txt")
	if err != nil {
		fmt.Println(err)
	}
	defer ovito_file.Close()

	//open file for mag / E writing
	mag_file , err2 := os.Create("results_"+strconv.Itoa(N)+"_"+strconv.Itoa(num_steps)+".txt")
	if err2 != nil {
		fmt.Println(err2)
	}
	defer mag_file.Close()

	var t float64
	for k :=23; k < 24; k++ {
		rand.Seed(time.Now().UTC().UnixNano())
		var J float64 = 1.0
		var h float64 = 0.0
		t = float64(k)/10.0
		fmt.Println("ON T = ", t)
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

		// write initial ovito
		write_ovito(initial_state, ovito_file)

		// return list of energies for each step of evolution
		energies,mags,g = evolve(initial_state,J,h,t,num_steps,ovito_file)
		avg_g := list_mean(g)
		
		avg_mag := mean(mags[num_steps/2:],false)
		avg_e := mean(energies[num_steps/2:],false)
		mag_fluct := (mean(mags[num_steps/2:],true) - avg_mag*avg_mag)/t
		heat_capac := (mean(energies[num_steps/2:],true) - avg_e*avg_e)/(t*t)

		fmt.Println("AVG E per site = ", avg_e/float64(N*N))
		fmt.Println("AVG MAG = ", avg_mag/float64(N*N))
		fmt.Println("MAG SUCEPT = ", mag_fluct)
		fmt.Println("HEAT CAPAC = ", heat_capac)
	}
}


// actually evolve matrix
func evolve(state [][]float64, J float64, h float64, T float64, num_steps int, file *os.File) ([]float64, []float64, [][]float64){
	// loop over "MC STEPS", i.e. a full NXN sweep of possible flips	
	energies := make([]float64, num_steps)
	mags := make([]float64, num_steps)
	g := make([][]float64, num_steps/2)
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
		mags[step] = mag(state)
		
		// print current avg mag for eq check
		fmt.Println(step, mean(mags[:step+1],false))
		
		// write current state to ovito file
		write_ovito(state, file)

		// calculate g(dist) at this point, return list and set to g of step, doing only second half to do only "eq section"
		if step >= num_steps/2 {
			g[step-(num_steps/2)] = spin_correlate(state)
		}
		if math.Mod(float64(step),1000.0) <= 0.01 {
			fmt.Println("ON STEP = ", step)
		}

	}
	return energies, mags, g
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
			// grabbing just these means never double counting --> because PBC
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


// total magnetization
func mag(state [][]float64) float64 {
	total_mag := 0.0
	N := len(state[0])
	for i := 0; i < N; i ++ {
		for j := 0; j < N; j++ {
			total_mag = total_mag + state[i][j]
		}
	}
	return total_mag
}


// return list of correlations between time separate magnetizations
func time_correlate(mags []float64) []float64 {
	l_max := 500
	correlate := make([]float64, l_max)
	for l := 0; l < l_max; l++ {
		avg_m := mean(mags[:len(mags)-l], false)
		avg_m_2 := mean(mags[:len(mags)-l], true)
		total := 0.0
		for i := 0; i < len(mags)-l; i ++ {
			total = total + mags[i]*mags[i+l]
		}
		avg_m_l := total/float64(len(mags)-l)
		numerator := avg_m_l - (avg_m*avg_m)
		denominator := avg_m_2 - (avg_m*avg_m)
		correlate[l] = numerator/denominator
	}
	return correlate
}


// other time correlation function
func alternate_correlate(mags []float64) []float64 {
	t_max := len(mags)
	l_max := 500
	correlate := make([]float64, l_max)
	for l := 0; l < l_max; l++ {
		combined_sum := 0.0
		for i := 0; i < len(mags)-l; i ++ {
			combined_sum = combined_sum + mags[i]*mags[i+l]
		}
		combined_sum = combined_sum/float64(t_max-l)
		split_sum_m := 0.0
		split_sum_m_l := 0.0
		for i := 0; i < len(mags)-l; i ++ {
			split_sum_m = split_sum_m + mags[i]
		}
		split_sum_m = split_sum_m/float64(t_max-l)
		for i := 0; i < len(mags)-l; i ++ {
			split_sum_m_l = split_sum_m_l + mags[i+l]
		}
		split_sum_m_l = split_sum_m_l/float64(t_max-l)
		correlate[l] = combined_sum - (split_sum_m*split_sum_m_l)
	}
	return correlate
}


// returns spin-spin correlation at some step in ising simulation
func spin_correlate(state [][]float64) []float64 {
	N := len(state[0])
	g := make([]float64, N-1) // length of steps in distance.. must be ints so only this many -->  at each int value will do all columns and rows though
	var total float64
	var count float64
	for r := 1; r < N; r ++ { 
		total = 0
		count = 0
		for i := 0; i < N; i ++ {
			for j := 0; j< N; j ++{
				// so here we can reference a particular cell.. could get both pairs for this one --> row and col pair?
				// (both pairs as in doing left and down, ignore rest as they  will get counted by check if doesnt exist based on dist / location of cell
				if r == 1 && i == N - 1 {
					total = total + state[i][j]*state[0][j]
					count = count + 1
				}
				if r == 1 && j == N -1 {
					total = total + state[i][j]*state[i][0]
					count = count + 1
				}
				if i + r <= N-1 {
					total = total + state[i][j]*state[i+r][j]
					count = count + 1
				}
				if j + r <= N-1 {
					total = total + state[i][j]*state[i][j+r]
					count = count + 1
				}
			}
		}
		g[r-1] = total/count
	}
	return g
}


// averages list of lists, returning one list of values averaged at that index for all the sublists
func list_mean(g [][]float64) []float64 {
	max_dist := len(g[0])
	num_lists := len(g)
	avg_g := make([]float64, max_dist)
	var total float64
	// go over length of all lists in g, then over ever list, for that index grab val and add to curr total
	// after through lists average by num lists, toss in g_avg, go to next index
	for r := 0; r < max_dist ; r ++ {
		total = 0
		for i := 0; i < num_lists; i ++ {
			total = total + g[i][r]
		}
		avg_g[r] = total/float64(num_lists)  // though we are getting the dist val, which starts 1, it is refrenced from 0 instead! (so just assign like g[r])
	}
	return avg_g
}


// Go has no math.mean() :(
func mean(data []float64, squared bool) float64 {
	var avg float64
	var N int = len(data)
	if !(squared) {
		for i := 0; i < N; i ++ {
			avg = avg + data[i]
		}
	} else {
		for i := 0; i < N; i ++ {
			avg = avg + (	data[i]*data[i])
		}	
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


// specific ovito format writing
func write_ovito(state [][]float64, file *os.File) {
	var string_type string
	var string_x string
	var string_y string
	string_N := strconv.Itoa(len(state[0])*len(state[0]))
	_,err := file.WriteString(string_N+"\n")
	_,err2 := file.WriteString("type	x	y"+"\n")
	if err != nil {
		fmt.Println("Error writing num sites")
	}
	if err2 != nil {
		fmt.Println("Error writing type	 x	 y")
	}
	// while at a current timestep, write all the positions of all the spins + type
	for i := 0; i < len(state[0]); i++ {
		for j := 0; j < len(state[0]); j++ {
			if state[i][j] == 1 {
				string_type = "1"
			} else {
				// otherwise must be -1
				string_type = "2"
			}
			string_x = strconv.Itoa(i)
			string_y = strconv.Itoa(j)
			_,err := file.WriteString(string_type+"	"+string_x+"	"+string_y+"\n")
			if err != nil {
				fmt.Println("Error writing spin data")
			}
		}
	}
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