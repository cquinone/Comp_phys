import numpy as np
import math
from math import cos
from numpy.linalg import inv
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------------------------------------
#Code to perform crank-nicholson finite difference calculation of inhomogenous active fluid
#Takes in initial points, activity, theta (temperature coefficient), viscosity, initial conditions, etc.
#Assumes boundary condition of grad(Q) = 0
#written in Python 2.7, uses numpy and matplotlib
#--------------------------------------------------------------------------------------------------------


np.set_printoptions(precision=5,linewidth=200, suppress=True)
def g_update(Q,stepcount):
	#update gama dot along mesh, based on Q
	for i in range(3*N):
		if i%3 == 0:
			points = [Q[i][0], Q[i+1][0], Q[i+2][0]]
			stress_curr = stress_minus(points)
			g_y[i/3] = ((n_hat*g0)+ sum(stress_sum)/len(stress_sum) - stress_curr)/(n_hat)
	
	return None

def stress_minus(points):
	#calculate stress MINUS the viscous component
	#needed for gamma dot calculation
	p1 = points[0]
	p2 = points[1]
	p3 = points[2]
	q = np.zeros((3,3))
	q[0][0] = p1
	q[0][1] = p2
	q[1][0] = p2
	q[1][1] = p3
	q[2][2] = -(p1+p3)
	
	q2 = np.dot(q,q)
	trq2 = np.trace(q2)
	h = -1*((1 - b/3.0) + b*trq2)*q
	h = np.add(h, b*q2)
	trI = b*trq2/3.0 * np.identity(3)
	h = np.add(h, -1*trI)		
	#passive now
	q_h = np.dot(q,h)
	h_q = np.dot(h,q)
	q3  = np.add(q, (1/3.0)*np.identity(3))
	st_pass = np.dot(2*e*q3, q_h) - np.dot(e*h, q3) - np.dot(e*q3,h)+q_h-h_q
	#active now
	st_act = act_hat*q
	st_tot = st_pass + st_act

	return st_tot[0][1]

def visualize(Q, t_tot, stepcount):
	#plot whatever you want
	g = g0
	xs = []
	stress = []
	for i in range(3*N):
		if i%3 == 0:
			#calculate the stress for each space point
			g	 = g_y[i/3]
			xs.append(i)
			(p1) = Q[i][0]
			(p2) = Q[i+1][0]
			(p3) = Q[i+2][0]
			q = np.zeros((3,3))
			q[0][0] = p1
			q[0][1] = p2
			q[1][0] = p2
			q[1][1] = p3
			q[2][2] = -(p1+p3)
	
			q2 = np.dot(q,q)
			trq2 = np.trace(q2)
			h = -1*((1 - b/3.0) + b*trq2)*q
			h = np.add(h, b*q2)
			trI = b*trq2/3.0 * np.identity(3)
			h = np.add(h, -1*trI)		
			#passive now
			q_h = np.dot(q,h)
			h_q = np.dot(h,q)
			q3  = np.add(q, (1/3.0)*np.identity(3))
			st_pass = np.dot(2*e*q3, q_h) - np.dot(e*h, q3) - np.dot(e*q3,h)+q_h-h_q
			#active now
			st_act = act_hat*q
			#visc now
			st_visc = np.zeros((3,3))
			st_visc[0][1] = n_hat*g
			st_tot = st_pass + st_act + st_visc
			#print "stress: ", st_tot[0][1]
			stress.append(st_tot[0][1])
	
	#see how far in time you are, but jump a bit ahead (because you plot so much)
	tpct = 3*(t_tot/tmax)
	if tpct > 1.0:
		tpct = 1
	#make color redder based on how far you are	
	color = '#%02x%02x%02x' % (tpct*(255), 0, 0)
	
	#plt.scatter(xs, stress, s=5, c="red", marker="o", edgecolors="red")
	#if t_tot == t0:
	#	xs = []
	#	for i in range(len(Q)):
	#		if i%3 == 1:
	#				xs.append(i)
	#
	#	ys = []
	#	for i in range(len(Q)):
	#		if i%3 == 1:
	#			ys.append(Q[i][0])
	#	plt.scatter(xs,ys, s=5, c="black", marker="o", edgecolors="black")
	

	if stepcount%1 == 0:
		#every 100 times, plot
		xs = []
		for i in range(len(Q)):
			if i%3 == 1:
					xs.append(i/3)

		ys = []
		for i in range(len(Q)):
			if i%3 == 1:
				ys.append(Q[i][0])

		#plt.scatter(xs,ys, s=5, c="orange", marker="o", edgecolors="orange")
		plt.xlim([-.5,25])
		#plot gamma dot profile
		plt.scatter(xs,g_y, s=5, c=color, marker="o", edgecolors=color)
	return None

def jacob(Q):
	#calculate df/dQ, without boundaries
	dF = np.zeros((3*N -6,3*N -6))
	for i in range(N-2):
		if i%3 == 0:
			g	 = g_y[i + 1]
			(q1) = Q[i][0]
			(q2) = Q[i+1][0]
			(q3) = Q[i+2][0]

			#set via df/dQ equations
			dF[i][i]		= (b/3.0 -1) - b*(6*((q1)**2) + 2*((q3)**2) + 2*((q2)**2) + 4*(q1)*(q3)) - (b/3.0)*(2*(q3)) - g*(2*e*(q3))
			dF[i][i+1]		= -b*(4*(q2)*(q1)) + 2*b*((q2)/3.0) + g*((e+1) - 2*e*((q1) + 1/3.0))
			dF[i][i+2]		= -b*(4*(q3)*(q1) + 2*((q1)**2)) - (b/3.0)*(4*(q3) + 2*(q1))
			dF[i+1][i]		= -b*(4*(q1)*(q2) + 2*(q3)*(q2)) + b*((q2)) + (g/2.0)*(e-1)
			dF[i+1][i+1]	= (b/3.0 -1) - b*(2*((q1)**2) + 2*((q3)**2) + 6*((q2))**2 + 2*(q1)*(q3)) + b*((q1) + (q3)) - 2*g*e*(2*(q2))
			dF[i+1][i+2]	= -b*(4*(q3)*(q2) + 2*(q1)*(q2)) + b*((q2)) + (g/2.0)*(e+1)
			dF[i+2][i]		= -b*(4*(q1)*(q3) + 2*((q3))**2) - (b/3.0)*(4*(q1) + 2*(q3))
			dF[i+2][i+1]	= -b*(4*(q2)*(q3)) + ((2*b)/3.0)*((q2)) + g*((e-1) - 2*e*((q1) + 1/3.0)) 
			dF[i+2][i+2]	= (b/3.0 -1) - b*(6*((q3)**2) + 2*((q1)**2) + 2*((q2)**2) + 4*(q1)*(q3)) + (2*b/3.0)*((q3)) - (b/3.0)*(2*(q1)) - g*(2*e*(q2))

	return dF

def funcs(x_init,g):
	#take in a list, use the list items as list[0] = q11, list[1] = q12, list[2] = q22
	x=x_init[0]
	y=x_init[1]
	z=x_init[2]
	#calculate H + S functions, f1 = the 11 element, f2 = 12, f3 = 22. 
	f1 = -(1 - (b)/3.0)*x - (b)*x*(2*x**2 +2*y**2 + 2*z**2 + 2*x*z) + (1/3.0)*(b)*(x**2 + y**2) - (1/3.0)*(b)*(2*z**2 + 2*x*z)      + (g)*y*((e + 1) - 2*e*(x + (1/3.0)))
	f2 = -(1 - (b)/3.0)*y - (b)*y*(2*x**2 + 2*y**2 + 2*z**2 + 2*x*z) + (b)*(x*y + z*y)                                              + ((z + (1/3.0))*((g)/2.0)*(e + 1) + (x + (1/3.0))*((g)/2.0)*(e - 1) - 2*y**2*(g)*e)
	f3 = -(1 - (b)/3.0)*z - (b)*z*(2*x**2 + 2*y**2 + 2*z**2 + 2*x*z) + (1/3.0)*(b)*(z**2 + y**2) - (1/3.0)*(b)*(2*x**2 + 2*x*z)     + (g)*y*((e - 1) - 2*e*(z + (1/3.0)))
	return f1,f2,f3

def initial(p1,p2):
	#Q_0 and Q have 3Nx1 shape, each 3 elements are q11, q12, q22 at some space point (of which there are N)
	Q_0		  = np.zeros((3*N,1))
	F_0		  = np.zeros((3*N,1))

	#p1 is left inital point
	Q_0[0][0] = p1[0] 
	Q_0[1][0] = p1[1]
	Q_0[2][0] = p1[2]

	#and p2 the right
	Q_0[3*N -1][0] = p2[2]
	Q_0[3*N -2][0] = p2[1]
	Q_0[3*N -3][0] = p2[0]

	for k in range(3):
		#constants for interpolation
		a = (Q_0[k][0] + Q_0[3*N - (3-k)][0])/2.0
		b = (Q_0[k][0] - Q_0[3*N - (3-k)][0])/2.0
		step = (pi)/(N - 1)
		#half cosine interpolation
		for i in range(1,N -1):
			Q_0[3*i + k][0] = a + b*cos(i*step) 
	
	#initialize stress_sum - needed to update g_y
	for i in range(3*N):
		if i%3 == 0:
			points= [Q_0[i][0],Q_0[i+1][0],Q_0[i+2][0]]
			stress_sum.append(stress_minus(points))

	#initialize g_y
	for i in range(3*N):
		if i%3 == 0:
			points = [Q_0[i][0], Q_0[i+1][0], Q_0[i+2][0]]
			stress_curr = stress_minus(points)
			g_y.append(((n_hat*g0)+sum(stress_sum)/len(stress_sum) - stress_curr)/(n_hat))

	return Q_0

def timestep(Q, stepcount):
	#construct slice of Q, leaving out boundaries
	Qint = np.zeros(((3*N)-6,1))
	for i in range(3,(3*N)-3):
		Qint[i-3][0] = Q[i][0]
	#now you have intermediate to move forward
	#identity matrix I
	I		= np.identity(3*N -6)
	#df/dQ
	dF		= jacob(Qint)
	#two intermediate values
	temp	= np.add(I,(-dt/2.0)*dF)
	temp2	= np.add(temp, (-dt/2.0)*M)
	inverse	= inv(temp2)
	#construct f(q) -- functions in qq dynamics equation that dont include grad^2 term
	F = np.zeros((3*N-6,1))
	for i in range(3*N - 7):
		if i%3 == 0:
			#grab associated gamma dot from g_y
			g = g_y[i/3 + 1]
			#and use correct Q values
			f1,f2,f3 = funcs([Qint[i][0],Qint[i+1][0],Qint[i+2][0]],g)
			F[i][0]		= f1
		if i%3 == 1:
			F[i][0]		= f2
		if i%3 == 2:
			F[i+2][0]	= f3
	#some more values for the calculation
	M_Q		= np.dot(M,Qint)
	dF_Q	= np.dot(dF,Qint)
	temp	= np.add((dt/2.0)*M_Q, -(dt/2.0)*dF_Q)
	temp2	= np.add(temp, Qint)
	temp3	= np.add(temp2, dt*(F))
	#put them together to get new Q (minus boundaries)
	Q_new	= np.dot(inverse, temp3)

	#now update Q (except boundaries)
	for i in range(3*N - 6):
		Q[i+3][0] = Q_new[i][0]

	#now update boundaries of Q
	Q[0][0] = Q_new[0][0]
	Q[1][0] = Q_new[1][0]
	Q[2][0] = Q_new[2][0]
	Q[(3*N)-3][0] = Q_new[(3*N)-9][0]
	Q[(3*N)-2][0] = Q_new[(3*N)-8][0]
	Q[(3*N)-1][0] = Q_new[(3*N)-7][0]

	#update stress_sum
	for i in range(3*N -3):
		if i%3 == 0:
			points= [Q[i][0],Q[i+1][0],Q[i+2][0]]
			stress_sum[i/3] = (stress_minus(points))
	
	#update g_y
	g_update(Q,stepcount)

	return Q

def main():
	t_tot = 0
	#two initial values of Q, from stable points
	#.2 and .3, 0 act, b=2.9
	p1 = [-0.1438040355342654, 0.020053120001894055, -0.15012127862328947]
	p2 = [0.33681739557643914, 0.10683582262212758, -0.16113780342545897]

	#.2 and .4, -.01 act, b= 2.9
	p1 = [-0.14380403553555873, 0.020053120002893846, -0.15012127863317348]
	p2 = [0.34768181087469335, 0.09848309791341625, -0.1689842555475872]

	#.4? p1 = [0.34768181087546546, 0.09848309791603506, -0.1689842555471204]
	#.48? p2 = [0.3544527249259323, 0.0924578582394749, -0.17388925901507707]
	
	print "Initializing"
	print ""
	#give Q some inital value, Q is now a 3N long column vector
	Q_0 = initial(p1,p2)
	print "Q inital: "
	print Q_0
	stepcount = 1
	#continue until tmax is reached
	while t_tot <= tmax:
		if t_tot == t0:
			#step Q in time
			Q = timestep(Q_0, stepcount)
			stepcount = stepcount + 1
			#plot some stuff
			visualize(Q, t_tot, stepcount)
		else:
			#step Q in time
			Q = timestep(Q, stepcount)
			stepcount = stepcount + 1
		#plot some stuff
		#visualize(Q, t_tot, stepcount)
		t_tot = t_tot+dt

	visualize(Q, t_tot, stepcount)
	#finally, save the plot
	plt.savefig("oactcrankcolsmall"+".png")
	return None

#Some constants / variables
#----------------------------------------------------------------#
pi		= math.pi
e 		= .7	#tumbling
b 		= 2.9 	#gamma temp
n_hat 	= .0567 #viscosity
g0		= 0.25  #gama dot between two initial points (average)
act_hat = 0.0		#activity
N		= 5	#how many points too discretize for mesh
t0		= 0				
dy      = 1.0/(float(N))		#length, 1/N so length of mesh is 1
dt    	= 1000*(.01*(dy**2))	#time step, chosen so certain quantities are small enough
steps   = 10000					#number of steps
tmax    = steps*dt				#max time to be reached
psi		= .01					#replaces K constant
stress_sum  = []
#make a list to hold gamma dot for some point in space
g_y			= []

#define matrix that does finite difference derivatives (double derivatives in space)
#3N because of Q size, 3N-6 as we leave out the boundaries, and -1, -1, -1 at boundaries to respect grad(Q) = 0.
M = np.zeros(((3*N)-6,(3*N)-6))
M[0][0] 	= -1.0
M[1][1] 	= -1.0
M[2][2]		= -1.0
M[0][3]		= 1.0
M[1][4]		= 1.0
M[2][5]		= 1.0
M[(3*N)-7][(3*N)-7]	 = -1.0
M[(3*N)-8][(3*N)-8]	 = -1.0
M[(3*N)-9][(3*N)-9]	 = -1.0
M[(3*N)-7][(3*N)-10] = 1.0
M[(3*N)-8][(3*N)-11] = 1.0
M[(3*N)-9][(3*N)-12] = 1.0
for i in range(3,(3*N)-10):
	M[i][i-3]		= 1
	M[i][i]			= -2
	M[i][i+3]		= 1
#scale the matrix by psi
M = (psi/(dy**2))*M
print "M:"
print M
print "dt,N,dy,tmax: ", dt, N, dy, tmax
main()
#----------------------------------------------------------------#