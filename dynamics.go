package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)


type Particle struct {
	x_pos float64
	y_pos float64
	x_real float64
	y_real float64
	x_v   float64
	y_v   float64
} 


func main() {
	// opening and read particle file
	var data = "data.txt"
	file, err := os.Open(data)
	if err != nil {
		fmt.Println("Open failure for data.txt")
		os.Exit(3)
	}
	var particle_data []string = make([]string, 0)
	scanner := bufio.NewScanner(file)
	var N int
	var L float64
	for scanner.Scan() {
		// append it to the particle data slice
		if scanner.Text()[0:1] != "#" {
			particle_data = append(particle_data, scanner.Text())
		} else {
			// grab particle number and box size
			tempdata := strings.Split(scanner.Text(), "	")
			N, _ = strconv.Atoi(tempdata[0][1:len(tempdata[0])])
			L,_  = strconv.ParseFloat(tempdata[1], 64)
		}
	}
	if scanner.Err() != nil {
		fmt.Println("Sorry: there was an error during the file reading")
		os.Exit(3)
	}
	file.Close()

	// process particle data
	parlist := process(particle_data)

	// run actual dynamics, recast L and pass for now
	p_e, k_e, times := dynamics(parlist, N, L)

	// get temperatures too
	//var temps []float64
    //temps = make([]float64, len(k_e))
    //for i := range k_e {
    //	temps[i] = k_e[i]/(3*float64(N))   // say k_b = 1 for now
    //}
	
	fmt.Println("PE CHECK: ", p_e[0])
	fmt.Println("KE CHECK: ", k_e[0])
	fmt.Println("TIME CHECK: ", times[0])
	// write pe, ke, temps, and times to file
	//file_name := "01_100000_12_"+strconv.Itoa(N)+"TEST.txt"
	//write_file(times, k_e, p_e, temps, file_name)
}


// run dynamics algorithm, return energy data after so many steps
func dynamics(parlist []Particle, N int, L float64) ([]float64, []float64, []float64) {
	// set dt, total_steps
	dt := .01
	total_steps := 20000
	sigma := 1.0
	var times []float64 = make([]float64,0)
	
	// store energy vals here
	var k_e []float64 = make([]float64, 0)
	var p_e []float64 = make([]float64, 0)
	
	// calc initial forces (before starting dynamics loop)
	x_forces, y_forces := force_matrices(parlist, N, L, sigma)

	// open ovito file
	periodic_file, err := os.Create("ovito_2d_periodic.xyz")
	if err != nil {
		fmt.Println(err)
	}
	defer periodic_file.Close()
	real_file, err := os.Create("ovito_2d_real.xyz")
	if err != nil {
		fmt.Println(err)
	}
	defer real_file.Close()
	
	// write initial positons
	write_ovito(N, parlist, periodic_file, "periodic")
	write_ovito(N, parlist, real_file, "real")

	// now do main loop, stepping dynamics until time is up
	for step := 0; step <= total_steps; step++ {
		// for each particle, update x and y positions
		for i := 0; i < N; i++ {
			parlist[i].x_pos,parlist[i].x_real = update_pos(parlist[i], x_forces, L, dt, "x", i)
			parlist[i].y_pos,parlist[i].y_real = update_pos(parlist[i], y_forces, L, dt, "y", i)
		}
		
		// now update forces, will use new positions
		x_forces_new, y_forces_new := force_matrices(parlist, N, L, sigma)
		
		// now we can update velocities
		for i := 0; i < N; i++ {
			parlist[i].x_v = update_veloc(parlist[i], x_forces, x_forces_new, dt, "x", i)
			parlist[i].y_v = update_veloc(parlist[i], y_forces, y_forces_new, dt, "y", i)
		}
		
		// update force matrices to new ones
		x_forces = set_matrix(x_forces, x_forces_new)
		y_forces = set_matrix(y_forces, y_forces_new)
		
		// calc and save pe and ke
		p_e = append(p_e, calc_potential(parlist, L))
		k_e = append(k_e, calc_kinetic(parlist))
		
		// update time
		times = append(times, dt*float64(step))

		// write ovito file
		if step >= 10000 && math.Mod(float64(step),100) == 0 {
			write_ovito(N, parlist, periodic_file, "periodic")
			write_ovito(N, parlist, real_file, "real")
		}

		// show progress
		if math.Mod(dt*float64(step), 10) == 0 {
			time_string := fmt.Sprintf("%f", dt*float64(step))
			fmt.Println("AT TIME = "+time_string)
		}
	}
	return p_e, k_e, times
}


// total potential energy of all particles
func calc_potential(parlist []Particle, L float64) float64{
	var pe float64
	for i := range parlist {
		for j := range parlist {
			p1 := parlist[i]
			p2 := parlist[j]
			if (p1.x_pos - p2.x_pos) < -L/2 {
				p2.x_pos = p2.x_pos - L
			} else if (p1.x_pos - p2.x_pos) > L/2 {
				p2.x_pos = p2.x_pos + L
			}
			if (p1.y_pos - p2.y_pos) < -L/2 {
				p2.y_pos = p2.y_pos - L
			} else if (p1.y_pos - p2.y_pos) > L/2 {
				p2.y_pos = p2.y_pos + L
			}
			radius := distance(p1,p2)
			if i != j {
				pe = pe + 4*((1/math.Pow(radius,12)) - (1/math.Pow(radius,6)))
			}
		}
	}
	// cover double counting with .5 !
	return .5*pe
}


// total kinetic energy of all particles
func calc_kinetic(parlist []Particle) float64{
	var ke float64
	for i := range parlist {
		ke = ke + math.Pow(parlist[i].x_v,2) + math.Pow(parlist[i].y_v,2)
	}
	return .5*ke
}


// velocity verlet to update positions
func update_pos(curr_p Particle, forces [][]float64, L float64, dt float64, axis string, index int) (float64, float64) {
	var new_pos float64
	var real_pos float64
	if axis == "x" {
		new_pos = curr_p.x_pos + curr_p.x_v*dt + row_sum(forces, index)*dt*dt*.5
		real_pos = curr_p.x_real + curr_p.x_v*dt + row_sum(forces, index)*dt*dt*.5
	}
	if axis == "y" {
		new_pos = curr_p.y_pos + curr_p.y_v*dt + row_sum(forces, index)*dt*dt*.5
		real_pos = curr_p.y_real + curr_p.y_v*dt + row_sum(forces, index)*dt*dt*.5
	}
	// return periodic modulo of position, and "real" position
	return math.Mod((2*L)+new_pos, L), real_pos
}


// velocity verlet to update velocities
func update_veloc(curr_p Particle, forces [][]float64, forces_new [][]float64, dt float64, axis string, index int) float64 {
	var new_veloc float64
	if axis == "x" {
		new_veloc = curr_p.x_v + row_sum(add_matrix(forces, forces_new), index)*dt*.5
	} else {
		new_veloc = curr_p.y_v + row_sum(add_matrix(forces, forces_new), index)*dt*.5
	}
	return new_veloc
}


// construct matrix of forces between all particles
func force_matrices(parlist []Particle, N int, L float64, sigma float64) ([][]float64, [][]float64) {
	x_forces := two_d(N)
	y_forces := two_d(N)
	// calc all forces on upper diag
	for i := 0; i < N; i++ {
		for j := i + 1; j < N; j++ {
			p1 := parlist[i]
			p2 := parlist[j]
			x_forces[i][j] = force_between(p1, p2, "x", L, sigma)
			y_forces[i][j] = force_between(p1, p2, "y", L, sigma)
		}
	}
	// now reflect negative values on lower diag
	for i := 1; i < N; i++ {
		for j := 0; j < i; j++ {
			x_forces[i][j] = -x_forces[j][i]
			y_forces[i][j] = -y_forces[j][i]
		}
	}
	return x_forces, y_forces
}


// lennard jones force between two particles
func force_between(p1 Particle, p2 Particle, axis string, L float64, sigma float64) float64 {
	// first adjust for min image
	if (p1.x_pos - p2.x_pos) < -L/2 {
		p2.x_pos = p2.x_pos - L
	} else if (p1.x_pos - p2.x_pos) > L/2 {
		p2.x_pos = p2.x_pos + L
	}
	if (p1.y_pos - p2.y_pos) < -L/2 {
		p2.y_pos = p2.y_pos - L
	} else if (p1.y_pos - p2.y_pos) > L/2 {
		p2.y_pos = p2.y_pos + L
	}
	radius := distance(p1, p2)
	// check if past cutoff distance
	if radius >= 3*sigma {
		return 0
	}
	temp_f := (48 / math.Pow(radius, 2)) * (1/(math.Pow(radius, 12)) - .5*(1/math.Pow(radius, 6)))
	var force float64
	if axis == "x" {              
		force = temp_f * (p1.x_pos - p2.x_pos)
	} else {
		// y force between two particles
		force = temp_f * (p1.y_pos - p2.y_pos)
	}
	return force
}


// distance between two particles
func distance(p1 Particle, p2 Particle) float64 {
	return math.Sqrt(math.Pow((p1.x_pos-p2.x_pos), 2) + math.Pow((p1.y_pos-p2.y_pos), 2))
}


// adds up floats in row of matrix
func row_sum(matrix [][]float64, index int) float64 {
	var sum float64
	for i := 0; i < len(matrix); i++ {
		sum = sum + matrix[index][i]
	}
	return sum
}


// sets up 2-d matrix of floats
func two_d(N int) [][]float64 {
	matrix := make([][]float64, N)
	for i := range matrix {
		matrix[i] = make([]float64, N)
	}
	return matrix
}


// helper function to add two matrices
func add_matrix(mat_1 [][]float64, mat_2 [][]float64) [][]float64 {
	mat_3 := two_d(len(mat_1))
	for i := 0; i < len(mat_1); i++ {
		for j := 0; j < len(mat_1); j++ {
			mat_3[i][j] = mat_1[i][j] + mat_2[i][j]
		}
	}
	return mat_3
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


// write energy and time data as tab separated columns
func write_file(times []float64, k_e []float64, p_e []float64, temps []float64, file_name string) {
	file, err := os.Create(file_name)
	if err != nil {
    	fmt.Println(err)
    	return
    }
    defer file.Close()
    for i := range k_e {
    	time := fmt.Sprintf("%f", times[i])
    	kinetic := fmt.Sprintf("%f", k_e[i])
    	potential := fmt.Sprintf("%f", p_e[i])
    	temp := fmt.Sprintf("%f", temps[i])
    	_,err := file.WriteString(time+"	"+kinetic+"	"+potential+"	"+temp+"\n")
    	if err != nil {
    		fmt.Println(err)
    		return
    	}
    }
}


// specific ovito format writing
func write_ovito(N int, parlist []Particle, file *os.File, kind string) {
	string_N := strconv.Itoa(N)
	_,err := file.WriteString(string_N+"\n")
	_,err2 := file.WriteString("x	y"+"\n")
	if err != nil {
		fmt.Println("Error writing num particles")
	}
	if err2 != nil {
		fmt.Println("Error writing x	y")
	}
	// we are at a current timestep, write all the positions of all the particles
	for i := 0; i < len(parlist); i ++ {
		p := parlist[i]
		// assume in "real" case, otherwise revalued if want periodic values
		string_x := fmt.Sprintf("%f", p.x_real)
		string_y := fmt.Sprintf("%f", p.y_real)
		if kind == "periodic" {
			string_x = fmt.Sprintf("%f", p.x_pos)
			string_y = fmt.Sprintf("%f", p.y_pos)
		}
		_,err := file.WriteString(string_x+"	"+string_y+"\n")
		if err != nil {
		fmt.Println("Error writing particle data")
		}
	}
}


// read particle data file and make list of particles
func process(particle_data []string) []Particle {
	var parlist []Particle
	parlist = make([]Particle, len(particle_data))
	for i := 0; i <= len(particle_data)-1; i++ {
		// split into 4 values
		tempdata := strings.Split(particle_data[i], "	")
		// grab float64 conv of each
		parlist[i].x_pos, _ = strconv.ParseFloat(tempdata[0], 64)
		parlist[i].y_pos, _ = strconv.ParseFloat(tempdata[1], 64)
		parlist[i].x_real = parlist[i].x_pos
		parlist[i].y_real = parlist[i].y_pos
		parlist[i].x_v, _ = strconv.ParseFloat(tempdata[2], 64)
		parlist[i].y_v, _ = strconv.ParseFloat(tempdata[3], 64)
	}
	fmt.Println("made it past process!!")
	return parlist
}
