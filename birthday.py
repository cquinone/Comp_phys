import math
#Calculate probabilities for birthday chances

def binom(n, k):
	ret = 0
	if n == k:
		ret = 1
	elif k == 1:
		ret = n
	else:
		a = math.factorial(n)
		b = math.factorial(k)
		c = math.factorial(n - k)
		ret = a // (b * c)
	return ret

p_b = (364/365.0)
pi = math.pi
print "Table without approximation"
print "N	PN"
print "----------------"
for N in range(1,51):
	p = 1-(p_b)**(.5*(N**2 - N))
	output = "%s	%s" % (N,p)
	print output

print "other"
print "N	PN"
print "----------------"
people = 50
for N in range(1,people+1):
	tot = 0
	for r in range(2,N+1):	
		#b = binom(people,2)
		#print "b: ", b
		#print "sum from 1 to", N+1
		#print "so (pb)^("
		tot = tot + (binom(N,r))
		tot  = (p_b)**tot
		tot = 1 - tot

	output = "%s	%s" % (N,tot)
	print output