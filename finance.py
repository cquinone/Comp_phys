import random
import math 
import matplotlib.pyplot as plt

r = float(raw_input("risk free rate: "))
sig = float(raw_input("volatility: "))
T = float(raw_input("Time: "))
S = float(raw_input("Stock original price: "))
n = int(raw_input("number of paths: "))
N = int(raw_input("N for one path: "))
X = float(raw_input("Strike Price: "))

#S=45.66
dt = T/N
s=[]
call=[]
put=[]
#do n paths
for i in range(n):
#do N steps
	s=[]
	s.append(S)
	for j in range(1,N):
		e = random.normalvariate(0,1)
		s.append(s[j-1]*math.exp((r-(sig**2)/2)*dt + sig*e*(dt**.5)))
	c=math.exp(-r*T) * max((s[-1]-X),0)
	p=math.exp(-r*T) * max((X-s[-1]),0)
	put.append(p)
	call.append(c)
	fig = plt.figure()
	plt.plot(range(N), s, 'o', ms=3)
	plt.show()

p = sum(put) / n
c = sum(call) / n
print "call option price: ", c
print "put option price	: ", p

