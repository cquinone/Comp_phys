import numpy as np
from scipy.linalg import solve
import math
#finite element code

#accept functions for A, P
def r(x):
	return 25*(2-(x/L)**3)
	#define this as whatever needs to be tested
def A(x):
	return math.pi*(r(x)**2)
	#define this as whatever needs to be tested
def P(x):
	return 2*math.pi*r(x)
	#define this as whatever needs to be tested

#define a trapezoidal method to approximate the intergals we need
def approx(a,b,f, dx):
	n = int((b-a)/dx)
	area = 0
	for i in range(1, n+1):
		x1 = a+(i-1)*dx
		x2 = a+i*dx
		darea = dx*(f(x1) + f(x2))/2
		area = area + darea

	return area

#make numpy arrays easy to see
np.set_printoptions(precision=3,linewidth=200, suppress=True)
#accept values for k, T_infinity, T_0, h, L
T_inf = float(raw_input("T_infinity: "))
k = float(raw_input("k: "))
T_0 = float(raw_input("T_0: "))
h = float(raw_input("h: "))
L = float(raw_input("L: "))
#specify the number of elements
N = int(raw_input("N: "))
#define delta J, size of elements is uniform
delt_j = float(L/N)
#all elements will be of the same size

#say that local coordinate dn_j is .1 (for all j), so that integration can be approximated:
dnj = .1
#now loop over j-2,3,4...N 
#make arrays for global stiffness elements
K1 = []
K2 = []
#global load
f = []
K1.append(0)
K2.append(0)
f.append(0)

#j from 2 to n+2 means j=1,2,3,4,5.....N+1
for j in range(1,N+2):
	#P1 , P2, P3. P4 functions are defined to capture the integrals that need to be done, using 
	#what P(x) has been defined as above
	def P1(x):
		return ((1-x)**2)*P(x)

	def P2(x):
		return ((1-x)*x)*P(x)
	
	def P3(x):
		return (x)*P(x)

	def P4(x):
		return ((x)**2)*P(x)

	#K1 and K2 and k and kprime from the derivation
	if j != N+1:
		K1.append(float((1/(delt_j))*k*approx(0,1,A,dnj) + delt_j*h*approx(0,1,P1,dnj))) 
		K2.append(float(-(1/(delt_j))*k*approx(0,1,A,dnj) + delt_j*h*approx(0,1,P2,dnj))) 
		f.append(delt_j*h*T_inf*approx(0,1,P3, dnj))
	
	#covering boundary condition
	if j == N+1:
		K1.append(float((1/(delt_j))*k*approx(0,1,A,dnj) + delt_j*h*approx(0,1,P4,dnj)) + h*A(L))
		K2.append(float(-(1/(delt_j))*k*approx(0,1,A,dnj) + delt_j*h*approx(0,1,P2,dnj))) 
		f.append(float(delt_j*h*T_inf*approx(0,1,P3, dnj)) + h*A(L)*T_inf)	

kfinal = np.zeros((N+2,N+2))
f_final = np.zeros((N+2,1))
#boundary
f_final[0,0] = T_0
kfinal[0,0] = 1
#set up rest of matrices
i=1
for q in range(1,N+2):
	if q <= N-1:
		kfinal[q,q]=K1[i]
		kfinal[q,q+1]=K2[i] + K1[i]
		kfinal[q,q+2]=K2[i]
	if q == N:
		kfinal[q,q]=K1[i]
		kfinal[q,q+1]=K2[i] + K1[i]
	if q == N+1:
		kfinal[q,q]=K1[i]

	i=i+1

#set up f_final matrix
i=1
for q in range(1,N+2):
	if q <= N:
		f_final[q,0] = f[i] + f[i+1] 
	if q == N+1:
		f_final[q,0] = f[i]
	i = i+1
#finally, we solve for T
#use imported solve function, which solves matrices of form Ax = b
#here we have k_final(T) = f_final
T = solve(kfinal, f_final)
print T


