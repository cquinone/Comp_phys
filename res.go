package main

import (
	"bufio"
	"fmt"
	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/vg"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

type Particle struct {
	name    string
	charge  int
	barynum int
	flag    int
	mass    float64
}

func main() {

	// opening and reading mass file
	var masses = "masses.lst"
	file, err := os.Open(masses)
	if err != nil {
		fmt.Println("Open failure for Masses.lst")
		os.Exit(3)
	}
	var lines []string = make([]string, 0)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// append it to the lines slice
		lines = append(lines, scanner.Text())
	}
	if scanner.Err() != nil {
		fmt.Println("Sorry: there was an error during the file reading")
		os.Exit(3)
	}
	file.Close()
	// file data now stored as slice, "lines"
	//particle info processing
	parlist := process(lines)

	//set up intro
	year, month, day := time.Now().Date()
	fmt.Println("")
	fmt.Println("Kinematics Version 3.0 ; 1 + 2 --> 3 + 4")
	fmt.Println("")
	fmt.Println("For one incident momentum (1) and a range of laboratory angles (3)")
	fmt.Println("")
	fmt.Print("              ")
	fmt.Print(year, month, day)
	fmt.Println("")
	fmt.Println("        Carnegie Mellon University   ")
	fmt.Println(" Relativistic 2-body Kinematics, Version 3.0")
	fmt.Println("")
	fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")
	input(parlist, filename)

}

func input(parlist []Particle, filename string) {
	//opening file
	outfile, _ := os.Create(filename)
	defer outfile.Close()
	canvas := CreateNewCanvas(1000,1000)
	canvas.SetLineWidth(2)
	canvas.SetStrokeColor(MakeColor(0, 0, 0))
	canvas.SetFillColor(MakeColor(0, 0, 0))
	//creating plot
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = "Particle 3 Momentum vs Lab Angle"
	p.X.Label.Text = "Lab Angle  (Deg)"
	p.Y.Label.Text = "P3 momentum  (MeV/c)"

	var particles []Particle
	particles = make([]Particle, 4)
	var a, b, c int
	for true {
		var check bool = false
		fmt.Println("Enter TWOBODY reaction (e.g. P(K-,PI-)S+):")
		var i string
		fmt.Scan(&i)
		//checking here for the right form of input
		for index := 0; index < len(i); index++ {
			if string(i[index]) == "(" {
				a = index
			}
			if string(i[index]) == "," {
				b = index
			}
			if string(i[index]) == ")" {
				c = index
			}
		}
		if a < b && b < c {
			check = true
		}	fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")
	input(parlist, filename)
		if check {
			var junk []string
			junk = make([]string, 2)
			//split along the comma
			junk = strings.Split(i, ",")
			var junk2 []string
			junk2 = make([]string, 2)

			// now split along the two paranthese to make 4 total strings
			junk2 = strings.Split(junk[0], "(")
			var junk3 []string
			junk3 = make([]string, 2)
			junk3 = strings.Split(junk[1], ")")

			//force capit	fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")
	input(parlist, filename)alization of particle names
			particles[0].name = strings.ToUpper(junk2[1])
			particles[1].name = strings.ToUpper(junk2[0])
			particles[2].name = strings.ToUpper(junk3[0])
			particles[3].name = strings.ToUpper(junk3[1])
			//search for particle's name and assign the correct data
			var repeat bool = false
			particles, repeat = assign(particles, parlist, repeat)
			// if there was an error
			if repeat == false {
				break	
			}
		} else {
		fmt.Println("Incorrect format")
		fmt.Println("")
		}
	}
	//add date, time, title to file we're writing too
	year, month1, day := time.Now().Date()
	month := month1.String()
	fmt.Fprintln(outfile, "Carnegie Mellon Univeristy")
	fmt.Fprint(outfile, "      ", "-", year, " ", month, " ", day, "-")
	fmt.Fprintln(outfile, "")
	fmt.Fprintln(outfile, "")
	//display particle data after assignment
	fmt.Println("")
	fmt.Println("Particle	Charge	Baryon	Mass")
	fmt.Println("-------------------------------------------")
	for _, i := range particles {
		fmt.Print(i.name)
		fmt.Print("		")
		fmt.Print(i.charge)
		fmt.Print("	")
		fmt.Print(i.barynum)
		fmt.Print("	")
		fmt.Print(i.mass)
		fmt.Println("")
	}
	//putting table in file
	fmt.Fprintln(outfile, "REACTION DATA:")
	fmt.Fprintln(outfile, "Particle	Charge	Baryon	Mass")
	fmt.Fprintf(outfile, "%s", "-------------------------------------------")
	fmt.Fprintln(outfile, "")
	for _, i := range particles {
		fmt.Fprint(outfile, i.name)
		fmt.Fprint(outfile, "		")
		fmt.Fprint(outfile, i.charge)
		fmt.Fprint(outfile, "		")
		fmt.Fprint(outfile, i.barynum)
		fmt.Fprint(outfile, "		")
		fmt.Fprint(outfile, i.mass)
		fmt.Fprintln(outfile,	fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")
	input(parlist, filename) "")
	}
	fmt.Fprintln(outfile, "")

	//check conservation
	fmt.Println("")
	if (particles[0].barynum + particles[1].barynum) != (particles[2].barynum + particles[3].barynum) {
		fmt.Println("-Baryon number is not conserved")
	} else {
		fmt.Println("-Baryon number is conserved")
	}

	if (particles[0].charge + particles[1].charge) != (particles[2].charge + particles[3].charge) {
		fmt.Println("-Charge is not consevred")
	} else {
		fmt.Println("-Charge is conserved")
	}
	// ask for excitation energy
	fmt.Println("")
	fmt.Println("Enter excitation energy of recoiling particle in MeV:")
	var exc string
	fmt.Scan(&exc)
	fmt.Println("")
	fmt.Fprintln(outfile, "The excitation energy is:", exc, " MeV")
	excfloat, _ := strconv.ParseFloat(exc, 64)
	particles[3].mass = particles[3].mass + excfloat
	// ask for momentum
	p1 := threshold(particles, excfloat, outfile)
	//getting angle range
	var thtmn, thtmx, thts float64
	fmt.Println("Enter THETAMIN THETAMAX THETASTEP :")
	fmt.Scanln(&thtmn, &thtmx, &thts) 
	fmt.Println("")
	thtmin := thtmn
	thtmax := thtmx
	thtstp := thts
	//constant for c
	light := float64(.299792458)
	p1mass := particles[0].mass
	p2mass := particles[1].mass	fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")
	input(parlist, filename)
	p3mass := particles[2].mass
	p4mass := particles[3].mass
	var p3 float64
	e1 := math.Pow((math.Pow(p1mass, 2) + math.Pow(p1, 2)), .5)
	beta1 := p1 / e1

	// recip velocity
	var rv1 float64
	if beta1 > 0 {
		rv1 = 1 / (light * beta1)
	} else {
		rv1 = 999
	}
	etot := e1 + p2mass
	t1 := e1 - p1mass
	etotcm := math.Pow((p1mass*p1mass)+(p2mass*p2mass)+2*e1*p2mass, .5)
	e3cm := (etotcm*etotcm + p3mass*p3mass - p4mass*p4mass) / (2 * etotcm)
	p34cm := math.Pow(e3cm*e3cm-p3mass*p3mass, .5)
	//	beta3cm := p3cm/e3cm
	betacm := p1 / etot
	gamma := 1 / (math.Pow(1-betacm*betacm, .5))
	p12cm := gamma * (p1 - betacm*e1)
	//	e1cm  := math.pow((p12cm*p12cm + p1mass*p1mass), .5)
	//	e2cm  := math.pow((p12cm*p12cm + p2mass*p2mass), .5)
	//	e4cm  := e1cm + e2cm - e3cm
	//print out some info
	fmt.Println("")
	fmt.Print("(Lab) ")
	fmt.Printf("%s%.3f%s%s%s%.3f%s%s%.3f%s", "P1: ", p1, " MeV/c", "	", "Beta1: ", beta1, "	 ", "1/v1: ", rv1*1e9, " ns/m")
	fmt.Println("")
	fmt.Println("")
	fmt.Print("(CoM) ")
	fmt.Printf("%s%.3f%s%s%s%.2f%s\n%s%s%.3f%s", "P12: ", p12cm, " MeV/c", "	", "P34: ", p34cm, " MeV/c", "	", "Total Relativistic Energy: ", etotcm/1000, " GeV/c")
	//also write this to file
	fmt.Fprintln(outfile, "")
	fmt.Fprint(outfile, "(Lab) ")
	fmt.Fprintf(outfile, "%s%.3f%s%s%s%.3f%s%s%.3f%s", "P1: ", p1, " MeV/c", "	", "Beta1: ", beta1, "	 ", "1/v1: ", rv1*1e9, " ns/m")
	fmt.Fprintln(outfile, "")
	fmt.Fprintln(outfile, "")
	fmt.Fprint(outfile, "(CoM) ")
	fmt.Fprintf(outfile, "%s%.3f%s%s%s%.2f%s%s%s%.3f%s", "P12: ", p12cm, " MeV/c", "	", "P34: ", p34cm, " MeV/c", "	", "Total Relativistic Energy: ", etotcm/1000, " GeV/c")
	fmt.Fprintln(outfile, "")
	fmt.Fprintln(outfile, "")
	//writing header to file
	fmt.Fprintln(outfile, "Angle3	  Angle3	   P3	         T3	Beta3	        1/v3        Jac.      Angle4	 Angle4        T4      q**2	1/lambda	   -t	         -u")
	fmt.Fprintln(outfile, "(Lab)       (cm)       (Lab)      (Lab)                (ns/m)                 (Lab)    (cm)        (Lab)   (GeV/c)2  (1/Fermi)    (GeV/c)2   (GeV/c)2")
	fmt.Fprintln(outfile, "------------------------------------------------------------------------------------------------------------------------------------------------------")
	fmt.Println("")
	fmt.Println("")
	fmt.Println("Angle3	Angle3	P3	T3	Beta3	1/v3	Jac.	Angle4	Angle4	T4	q**2	1/lambda	-t	-u")
	fmt.Println("(Lab)       (cm)       (Lab)      (Lab)                (ns/m)                 (Lab)    (cm)        (Lab)   (GeV/c)2  (1/Fermi)    (GeV/c)2   (GeV/c)2")
	fmt.Println( "------------------------------------------------------------------------------------------------------------------------------------------------------")
	//slice of points for plot
	pts := make(plotter.XYs, (int(thtmax/thtstp) + 1))
	j := 0

	//main loop of calculation
	for i := thtmin; i <= thtmax; i = i + thtstp {
		dg := i * (math.Pi / 180)
		sin := math.Sin(dg)
		cos := math.Cos(dg)
		p3p, p3m, stop := kinetic(t1, p1mass, p2mass, p3mass, p4mass, dg)
		if stop == 1 {
			fmt.Println("")
			fmt.Println("")
			fmt.Println("End of possible range")
			fmt.Println("")
			break
		}
		// cover all possible solutions
		solutions := make([]float64, 2)
		solutions[0] = p3p	fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")
	input(parlist, filename)
		solutions[1] = p3m
		counter := 0
		for _, moment := range solutions {
			counter++
			if moment >= 0 {
				p3 = moment
				e3 := math.Pow(p3*p3+p3mass*p3mass, .5)	
				t3 := e3 - p3mass
				beta3 := p3 / e3
				r := betacm / beta3	
				thtcm := (180 / math.Pi) * math.Atan2(sin, (gamma*(cos-r)))
				if thtcm < 0 {
					thtcm = thtcm + 180
				}	
				tht4cm := thtcm - 180
				perv := 1 / (beta3 * light)
				if perv <= .001 {
					perv = 0
				}
				if perv >= 999.999 {
					perv = 999.99
				}
				p13 := p1 * p3 * cos
				p4 := math.Pow((p1*p1 + p3*p3 - 2*p13), .5)
				t4 := math.Pow((p4*p4+p4mass*p4mass), .5) - p4mass
				var tht41 float64 = 0
				if p4 != 0 {
					tht41 = (180 / math.Pi) * math.Asin(-(p3/p4)*sin)
				} else {
					tht41 = -90
				}
				qsqrcm := (p12cm*p12cm + p34cm*p34cm - 2*p12cm*p34cm*math.Cos(thtcm*(math.Pi/180))) / 1e6
				perlcm := math.Pow(math.Abs(qsqrcm*1000000), .5) / 197.327053
				rmtot := (p1mass*p1mass + p3mass*p3mass + p2mass*p2mass + p4mass*p4mass) / 1e6
				tchan := (2*(e1*e3-p13) - (p1mass*p1mass + p3mass*p3mass)) / 1e6
				schan := (etotcm * etotcm) / 1e6
				uchan := -(rmtot - schan + tchan)
				factor := math.Pow(1-betacm*betacm*sin*sin-r*(2*cos-r), .5)
				g := factor * (1 - r*cos) * gamma * gamma
				if g < .00001 {
					g = 0
				}
				if g > 99.999 {
					g = 99.999
				}

				//give current point correct data
				pts[j].X = float64(i)
				pts[j].Y = p3
				j++
				// points for half circle, then draw small squares as points
				x := (p3)/10*math.Cos(math.Pi/180 * i) 
				y := (p3)/10*math.Sin(math.Pi/180 * i)
				Square(canvas, x, y, 2)
				//write data to table in file
				if counter <= 1 {


					fmt.Println(float64(i), "	", int(thtcm), "	", int(p3), "	", int(t3), "  ", float32(beta3), "	", int(perv), "	", float32(g), "	", int(tht41), "	", int(tht4cm), "	", int(t4), "	", float32(qsqrcm), "	", float32(perlcm), "	", float32(tchan), "	", float32(uchan))
					fmt.Fprintf(outfile, "%7.1f%s%7.1f%s%7.1f%s%7.1f%s%7.3f%s%7.3f%s%7.3f%s%7.1f%s%7.1f%s%7.1f%s%7.3f%s%7.3f%s%7.3f%s%7.3f\n", float64(i), "	", thtcm, "	", p3, "	", t3, "  ", beta3, "	", perv, "	", g, "	", tht41, "	", tht4cm, "	", t4, "	", qsqrcm, "	", perlcm, "	", tchan, "	", uchan)

				} else {


					fmt.Println("		", int(thtcm), "	", int(p3), "	", int(t3), "	", float32(beta3), "	", int(perv), "	", float32(g), "	", int(tht41), "	", int(tht4cm), "	", int(t4), "	", float32(qsqrcm), "	", float32(perlcm), "	", float32(tchan), "	", float32(uchan))
					fmt.Fprintf(outfile, "%s%7.1f%s%7.1f%s%7.1f%s%7.3f%s%7.3f%s%7.3f%s%7.1f%s%7.1f%s%7.1f%s%7.3f%s%7.3f%s%7.3f%s%7.3f\n", "		", thtcm, "	", p3, "	", t3, "	", beta3, "	", perv, "	", g, "	", tht41, "	", tht4cm, "	", t4, "	", qsqrcm, "	", perlcm, "	", tchan, "	", uchan)
				
				}

			}
		}
	}

	// make and add lines to plot pointsf[x_] := x^2
	line, err := plotter.NewLine(pts)
	if err != nil {
		log.Panic(err)
	}
	p.Add(line)
	p.Add(plotter.NewGrid())	fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")
	input(parlist, filename)
	p.Y.Padding = 0 * vg.Inch
	p.X.Padding = 0 * vg.Inch
	//save in png file
	err = p.Save(6*vg.Inch, 6*vg.Inch, "collisiondata.png")
	if err != nil {
		log.Panic(err)
	}
	fmt.Println("")
	var radius string
	fmt.Println("Name the file (type .png) for radius:")
	fmt.Scan(&radius)
	canvas.SaveToPNG(radius)
	fmt.Println("")
	fmt.Println("DONE")
	os.Exit(3)
}

func assign(particles []Particle, parlist []Particle, repeat bool) ([]Particle, bool) {
	// count if all 4 are possible particle
	convertcount := 0
	for _, i := range parlist {
		if particles[0].name == i.name {
			convertcount++
			particles[0] = i
		}
		if particles[1].name == i.name {
			convertcount++
			particles[1] = i
		}
		if particles[2].name == i.name {
			convertcount++
			particles[2] = i
		}
		if particles[3].name == i.name {
			convertcount++
			particles[3] = i
		}	fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")
	input(parlist, filename)
	}
	if convertcount != 4 {
		fmt.Println("")
		fmt.Println("At least one particle was not found in masses.lst")
		fmt.Println("")
		repeat = true
	}
	return particles, repeat
}

func threshold(particles []Particle, excfloat float64, outfile *os.File) float64 {
	w := particles[2].mass + particles[3].mass
	w2 := math.Pow(w, 2)
	var e1 float64
	var pt float64
	var t1 float64

	if particles[2].mass > 0 {
		e1 = (w2 - math.Pow((particles[0].mass), 2) - math.Pow((particles[1].mass), 2)) / (2 * particles[1].mass)
	} else {
		e1 = 0
	}
	if e1 > particles[0].mass {
		pt = math.Pow((math.Pow(e1, 2) - math.Pow((particles[0].mass), 2)), .5)
	} else {
		pt = 0
	}
	t1 = e1 - particles[0].mass
	if t1 >= 0 {
		fmt.Println("The reaction threshold is", t1, "Mev, or", pt, "Mev/c")
		fmt.Fprintln(outfile, "The reaction threshold is", t1, "Mev, or", pt, "Mev/c")
	} else {
		fmt.Println("The reaction is exothermic by", t1, "Mev")
		fmt.Fprintln(outfile, "The reaction is exothermic by", t1, "Mev")
	}
	p1 := momentum(pt)
	fmt.Println("")
	return p1

}

func momentum(pt float64) float64 {
	fmt.Println("")
	var mom_s string
	var mom float64
	var fake string
	for true {
		fmt.Scanln(&fake)
		fmt.Print(fake)
		fmt.Println("Enter the momentum of incident particle:")
		fmt.Scanln(&mom_s)
		mom, _ = strconv.ParseFloat(mom_s, 64)
		if mom < pt {
			fmt.Println("the given momentum, ", mom, " Mev/c, is below the threshold momentum,", pt)
		}
		if mom >= pt {
			break
		}
	}
	return mom
}

func kinetic(t1, p1mass, p2mass, p3mass, p4mass, i float64) (float64, float64, int) {
	stop := 0
	e1 := t1 + p1mass
	p1 := math.Pow(e1*e1-p1mass*p1mass, .5)
	rm := p1mass*p1mass + p2mass*p2mass + p3mass*p3mass - p4mass*p4mass
	etot := e1 + p2mass
	a := p1 * math.Cos(i)
	b := rm + 2*e1*p2mass
	c := etot*etot - a*a
	sqt := b*b - 4*p3mass*p3mass*c
	if sqt > 0 {
		sqt = math.Pow(sqt, .5)
	} else {
		stop = 1
		p3p := 0.0
		p3m := 0.0
		return p3p, p3m, stop
	}
	p3p := (a*b + etot*sqt) / (2 * c)
	p3m := (a*b - etot*sqt) / (2 * c)
	return p3p, p3m, stop
}

func process(lines []string) []Particle {
	// creates []string that holds the names of the particles (to be added to []Particle cells)
	var parlist []Particle
	parlist = make([]Particle, len(lines))
	//variables for the particle info at some line i
	var tempname string
	var tempbar int
	var tempch int
	var tempflag int
	var tempmass float64
	var linecount int
	for i := 27; i <= len(lines)-1; i++ {
		// turn offset into current number of accrued zero length lines
		var counter int = 0
		if len(lines[i]) != 0 {
			tempstr := strings.Fields(lines[i])
			for i := 0; i < len(tempstr[0]); i++ {
				if string(tempstr[0][i]) == "'" {
					counter = i
				}
			fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")	fmt.Print("Name the output file, or <k> for ")
	fmt.Printf("%s\n", "\"kin.out\":")
	var filename string
	fmt.Scan(&filename)
	if filename == "k" {
	fmt.Println("")
	fmt.Print("The output data will be written in the default file ")
	fmt.Printf("%s\n", "\"kin.out\".")
	}
	fmt.Println("") 
	s := []string{filename, ".txt"}
	filename = strings.Join(s,"")
	input(parlist, filename)
	input(parlist, filename)	}
			linecount++
			tempname = string(tempstr[0][1:counter])
			tempch, _ = strconv.Atoi(string(tempstr[1]))
			tempbar, _ = strconv.Atoi(string(tempstr[2]))
			tempflag, _ = strconv.Atoi(string(tempstr[3]))
			tempmass, _ = strconv.ParseFloat(string(tempstr[4]), 64)
			//set current bary_num, charge, and name with correct index, if on a non empty line
			parlist[linecount].name = tempname
			parlist[linecount].charge = tempch
			parlist[linecount].barynum = tempbar
			parlist[linecount].mass = tempmass
			parlist[linecount].flag = tempflag
			if tempflag == 0 {
				parlist[linecount].mass = tempmass + float64(tempbar)*(931.5016) - (.51100034)*(float64(tempch))
			}
		}
	}
	return parlist
}

func Square(canvas Canvas, x float64, y float64, b float64) Canvas{
	canvas.SetStrokeColor(MakeColor(0, 0, 0))
	canvas.SetFillColor(MakeColor(0, 0, 0))
	//make square by drawing lines between 4 points
		canvas.MoveTo(x,y)
		canvas.LineTo(x,y-b)
		canvas.MoveTo(x,y-b)
		canvas.LineTo(x+b,y-b)
		canvas.MoveTo(x+b,y-b)
		canvas.LineTo(x+b,y)
		canvas.MoveTo(x+b,y)
		canvas.LineTo(x,y)
	canvas.FillStroke()
	canvas.Stroke()
	return canvas
}
