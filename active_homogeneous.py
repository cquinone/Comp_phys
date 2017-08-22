import numpy as np
import scipy.optimize
import scipy.stats
import math
import matplotlib.pyplot as plt
from time import gmtime, strftime

#-------------------------------------------------------------------------------------------------------
#Code to solve nonlinear equations for dynamics of Q, for given theta, range of gamma dot, and activity. 
#Adjust act_hat for activity, g for gammadot, b for theta.
#Written in python 2.7, uses scipy, matplotlib, and numpy
#-------------------------------------------------------------------------------------------------------

np.set_printoptions(precision=5,linewidth=200, suppress=True)

def same(eig, eigtot):
	#determine if an eigen value is already on the totallist
	for testeig in eigtot:
		diff1 = (testeig[0]-eig[0])**2
		diff2 = (testeig[1]-eig[1])**2
		diff3 = (testeig[2]-eig[2])**2
		pct = ((diff1+diff2+diff3)**.5)/(mag(testeig))
		if pct <= .05:
			#if within 5%, consider the one testing against as the same
			return True

	return False

def mag(vect):
	#calculate the magnitude of a vector
	mag = 0
	for comp in vect:
		mag = mag + (comp)**2
	mag = float(mag)
	return (mag**.5)

def funcs2(guess):
	#take in a list, use the list items as list[0] = q11, list[1] = q12, list[2] = q22
	x=guess[0]
	y=guess[1]
	z=guess[2]
	#calculate H + S functions, f1 = the 11 element, f2 = 12, f3 = 22. 
	f1 = -(1 - (b)/3.0)*x - (b)*x*(2*x**2 +2*y**2 + 2*z**2 + 2*x*z) + (1/3.0)*(b)*(x**2 + y**2) - (1/3.0)*(b)*(2*z**2 + 2*x*z)      + (g)*y*((e + 1) - 2*e*(x + (1/3.0)))
	f2 = -(1 - (b)/3.0)*y - (b)*y*(2*x**2 + 2*y**2 + 2*z**2 + 2*x*z) + (b)*(x*y + z*y)                                              + ((z + (1/3.0))*((g)/2.0)*(e + 1) + (x + (1/3.0))*((g)/2.0)*(e - 1) - 2*y**2*(g)*e)
	f3 = -(1 - (b)/3.0)*z - (b)*z*(2*x**2 + 2*y**2 + 2*z**2 + 2*x*z) + (1/3.0)*(b)*(z**2 + y**2) - (1/3.0)*(b)*(2*x**2 + 2*x*z)     + (g)*y*((e - 1) - 2*e*(z + (1/3.0)))
	f = np.array([(f1),(f2),(f3)]) #make it an numpy array
	f.shape = (3,1) #now make it a column vector
	return f

def search():
#solve system of nonlinear eqs, using scipy package
#use newton's method, search from -1/3 to 1/3 as this is only possible range

	x0 = -((1/3.0)+.01)
	solarray = []
	for i in range(10):
		#look for solutions in range -1/3 to 1/3 for each varaible (3 unique vars)
		z0 = -((1/3.0)+.01)
		y0 = -((1/3.0)+.01)
		x0 = x0 + .15
		for j in range(10):
			y0 = -((1/3.0)+.01)
			z0 = z0 + .15
			for k in range(10):
				y0 = y0 + .15
				x_init = [x0,y0,z0]
				try:
					#use newtons method to check that point
					x = scipy.optimize.newton_krylov(funcs2,x_init, maxiter=8)
				except:
					pass
				else:
					#convert to regular float data type
					xreg = [float(x[0]),float(x[1]),float(x[2])]
					#append that solution to total list
					solarray.append(xreg)
	#the rest of this function looks through the total list and gets rid of doubles
	#"pct" is checked, if within 5%, but not the same solution, then reject it as a repeat
	fixsol = solarray
	for testsol in solarray:
		for sol in fixsol:
			if abs(sol[0]+sol[1]+sol[2]) <= .0001:
				fixsol.remove(sol)
				continue
			diff1 = (testsol[0]-sol[0])**2
			diff2 = (testsol[1]-sol[1])**2
			diff3 = (testsol[2]-sol[2])**2
			pct = ((diff1+diff2+diff3)**.5)/(mag(testsol))
			if pct <= .05 and pct != 0:
				fixsol.remove(sol)
				continue

	fixsol2 = fixsol
	for testsol in fixsol:
		for sol in fixsol2:
			if abs(sol[0]+sol[1]+sol[2]) <= .0001:
				fixsol2.remove(sol)                                                                                                                      
				continue
			diff1 = (testsol[0]-sol[0])**2
			diff2 = (testsol[1]-sol[1])**2
			diff3 = (testsol[2]-sol[2])**2
			pct = ((diff1+diff2+diff3)**.5)/(mag(testsol))
			if pct <= .05 and pct != 0:
				fixsol2.remove(sol)
				continue
	#return the cleaned list
	return fixsol2

#define units/constants 
t_hat 	= 1.0		#a_0*GAMMA, not actually included in this code, but technically could be
e 		= .7		#tumbling
b 		= 2.66 		#gamma temp, or theta
n_hat 	= .0567 	#viscosity
g0		= -.5		#gamma dot
act_hat0 = 0.9	#acitivty

#file to write to
filename = str(b) + "_" + str(act_hat0) + "_" + "txt" 
print "filename: ", filename
#how many gamma dots points
max = 100
file = open(filename, "a")

for j in range(1):
	stress 	= []
	stable  = []
	gs		= []
	g = g0
	print "theta: ", b
	#change above range for multiple runs of different activities
	if j == 0:
		color = "red"
		shape="o"
		act_hat = act_hat0
	if j == 1:
		color = "red"
		shape="o"
		act_hat = -.1
	if j ==2:
		color = "yellow"
		shape="o"
		act_hat = 0
	print "activity: ", act_hat
	for i in range(max):
	#for the number of steps
		stress.append([])
		stable.append([])
		totsol = []
		uniqeig = []
		eigvects = []
		print "GAMMA  : ", g
		file.write("gamma: " + str(g)+"\n")
		#find solutions for this gamma dot
		solarray = search()
		count = 0
		#for all the solutions found
		for sol in solarray:
			qtest = sol
			qtot = np.zeros((3,3))
			qtot[0][0] = qtest[0]
			qtot[0][1] = qtest[1]
			qtot[1][0] = qtot[0][1]
			qtot[1][1] = qtest[2]
			qtot[2][2] = -(qtest[0] + qtest[2])
			#construct Q ^
			eig, eigv = np.linalg.eig(qtot)  
			#and grab the eigenvalues / vectors
			if count == 0:
				e1 = float(eig[0])
				e2 = float(eig[1])
				e3 = float(eig[2])
				eig = [e1,e2,e3]
				uniqeig.append(eig)
				eigvects.append(eigv)
				totsol.append(sol)
			else:
				#check if you seen a solution with eigenval already
				if not (same(eig, uniqeig)):
					e1 = float(eig[0])
					e2 = float(eig[1])
					e3 = float(eig[2])
					eig = [e1,e2,e3]
					#append since this is a new solution
					uniqeig.append(eig)
					eigvects.append(eigv)
					totsol.append(sol)
			count = count + 1

		f = 0
		for sol in totsol:
			#for the solutions deemed unique
			print "unique sol: ", sol
			file.write("unique sol: " + str(sol)+"\n")
			print "eigenval: ", uniqeig[f]
			print "eigenvect: ", eigvects[f]
			f = f + 1 
			q = np.zeros((3,3))
			q[0][0] = sol[0]
			q[0][1] = sol[1]
			q[1][0] = sol[1]
			q[1][1] = sol[2]
			q[2][2] = -(sol[0]+sol[2])
			#calculate stability
			stab = np.zeros((3,3))
			stab[0][0] = (b/3.0 -1)-b*(6*sol[0]**2 + 4*sol[0]*sol[2] + 2*sol[1]**2 + 2*sol[2]**2) + (b/3.0)*(2*sol[0])-(b/3.0)*(2*sol[2]) - 2*g*e*sol[1]
			stab[0][1] = -b*(4*sol[1]*sol[0]) + (b/3.0)*(2*sol[1]) + g*((e+1) - 2*e*(sol[0] + 1/3.0)) #got rid of last sol[1] ??
			stab[0][2] = -b*(4*sol[0]*sol[2] + 2*sol[0]**2) + (-b/3.0)*(4*sol[2] + 2*sol[0])
			stab[1][0] = -b*(4*sol[0]*sol[1]+ 2*sol[2]*sol[1]) + b*(sol[1]) + (g/2.0)*(e-1)
			stab[1][1] = (b/3.0 -1)-b*(2*sol[0]**2 + 2*sol[2]**2 + 6*sol[1]**2 + 2*sol[2]*sol[0]) + b*(sol[0]+ sol[2]) - 4*g*e*(sol[1]) #add initial stuff
			stab[1][2] = -b*(4*sol[2]*sol[1] + 2*sol[0]*sol[1]) + b*(sol[1]) + (g/2.0)*(e+1)
			stab[2][0] = -b*(4*sol[0]*sol[2] + 2*sol[2]**2) - (b/3.0)*(4*sol[0] + 2*sol[2])
			stab[2][1] = -b*(4*sol[1])*sol[2] + (b/3.0)*(2*sol[1]) + g*((e-1) - 2*e*(sol[2] + 1/3.0)) #change to q22 in last bit
			stab[2][2] = (b/3.0 -1) - b*(6*sol[2]**2 + 4*sol[0]*sol[2] + 2*sol[0]**2 + 2*sol[1]**2) + (b/3.0)*(2*sol[2]) - (b/3.0)*(2*sol[0]) - g*(2*e*sol[1]) #on 4 change to q11
			eig, eigv = np.linalg.eig(stab)
			print "stability eig: ", eig
			ecount = 0
			for etest in eig:
				#check if this is a stable solution, if so, ecount will increment up to 3
				if isinstance(etest, complex):
					if etest.real < 0:
						ecount = ecount + 1
				else:
					if etest < 0:
						ecount = ecount + 1

			file.write("stability eig: " + str(eig)+"\n")
			print ""

			#some intial values needed to calculate the passive stress
			q2 = np.dot(q,q)
			trq2 = np.trace(q2)
			h = -1*((1 - b/3.0) + b*trq2)*q
			h = np.add(h, b*q2)
			trI = b*trq2/3.0 * np.identity(3)
			h = np.add(h, -1*trI)
			
			#now the overall passive calculation
			q_h = np.dot(q,h)
			h_q = np.dot(h,q)
			q3  = np.add(q, (1/3.0)*np.identity(3) )
			st_pass = np.dot(2*e*q3, q_h) - np.dot(e*h, q3) - np.dot(e*q3,h)+q_h-h_q
			
			#active now
			st_act = act_hat*q
			#visc now
			st_visc = np.zeros((3,3))
			st_visc[0][1] = n_hat*g
			st_tot = st_pass + st_act + st_visc
			stress[-1].append(st_tot[0][1])
			if ecount == 3:
				#if stable
				stable[-1].append(st_tot[0][1])
				if g not in gs:
					gs.append(g)
			print "stress xy: ", st_tot[0][1]
			file.write("stress xy: " + str(st_tot[0][1])+"\n")
			file.write(""+"\n")

		g = g + .01
		print "---------------------------------------------------"
		file.write("---------------------------------------------------"+"\n")

	xs=[]
	for i in range(max):
		xs.append(g0+.01*i)

	print "gamma dots : ", xs
	print "stresses : ", stress
	print ""
	print "gamma dots for stability : ", gs
	time = strftime("%Y-%m-%d %H:%M:%S", gmtime())
	for xe, se in zip(xs, stress):
		plt.scatter([xe] * len(se), se, s=6, c=color, marker=shape,  edgecolors=color)

	for xe, se in zip(gs, stable):
	   plt.scatter([xe] * len(se), se, s=6, c="blue", marker="o",  edgecolors="blue")
	
	title  = "b=" + str(b)+"g="+str(g)+"act="+str(act_hat)
	#plt.title(title)
file.close()
plt.savefig(title+".png")
plt.show()