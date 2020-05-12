import math
import numpy as np
import time
import sys
import matplotlib.pyplot as plt
import random as rand
import statistics as stats
from collections import defaultdict


class Graph():
	def __init__(self):
		self.graph = defaultdict(list) #sets kind to add (if no key exists) as list
	
	def addEdge(self,u,v): 
		self.graph[u].append(v)

	def connections(self,coord):
		return self.graph[coord]

	def make(self,u):
		self.graph[u] = []

  
def gen_bonds(grid,T,J):
	# generate graph representing bonds that define spin clusters
	# if neighbor states are same, assign bond based on prob distro
	N = len(grid[0])
	g = Graph()
	visited = np.full((N,N),False,dtype=bool)
	for i in range(N):
		for j in range(N):
			# generate bond probability
			bond_prob = 1.0 - math.exp(-2.0*J*(1.0/T))
			# generate PBC coordinates for neighbors
			down = i + 1
			right = j + 1
			up = i - 1
			left = j - 1
			if i == N - 1:
				down = 0
			if j == N - 1:
				right = 0
			if i == 0:
				up = N - 1
			if j == 0:
				left = N - 1
			# actual check, rand.uniform is repeated as each case is separate chance?
			if grid[i,j] == grid[down,j]:
				r = rand.uniform(0, 1)
				if r < bond_prob:
					g.addEdge((i,j),(down,j))
					g.addEdge((down,j),(i,j))
			if grid[i,j] == grid[i,right]:
				r = rand.uniform(0, 1)
				if r < bond_prob:
					g.addEdge((i,j),(i,right))
					g.addEdge((i,right),(i,j))
			if grid[i,j] == grid[i,left]:
				r = rand.uniform(0, 1)
				if r < bond_prob:
					g.addEdge((i,j),(i,left))
					g.addEdge((i,left),(i,j))
			if grid[i,j] == grid[up,j]:
				r = rand.uniform(0, 1)
				if r < bond_prob:
					g.addEdge((i,j),(up,j))
					g.addEdge((up,j),(i,j))
			if grid[i,j] != grid[i,right] and grid[i,j] != grid[down,j] and grid[i,j] != grid[i,left] and grid[i,j] != grid[up,j]:
				g.make((i,j))

	return g


def search_adj(i,j,visited,bond_graph,coords):
	#use i,j and depth first search for region from given site, return visited updated from search
	# also grab coordinates of sites you reach
	visited[i,j] = True #we've visited this point now
	for point in bond_graph.connections((i,j)):
		k = point[0]
		l = point[1]
		if visited[k,l] == False:
			coords.append((k,l))
			search_adj(k,l,visited,bond_graph,coords)
	return visited,coords


def update_clusters(bond_graph,grid,N): 
	# check first cell, find all adjacent similar cells, mark each on other matrix as false
	# continue, skipping any cell that is True as visited in other matrix --> this is DFS
	visited = np.full((N,N),False,dtype=bool)
	curr_clust_coords = []
	for i in range(N):
		for j in range(N):
			if visited[i,j] == False:  #haven't been here yet
				curr_clust_coords.append((i,j))				
				visited,coords = search_adj(i,j,visited,bond_graph,[]) #search here through all connected ones, update visited
				for coord in coords:
					curr_clust_coords.append(coord)
				r = rand.uniform(0,1)
				if r < .5:
					for coord in curr_clust_coords:
						grid[coord[0]][coord[1]] = -1*grid[coord[0]][coord[1]]
				curr_clust_coords = []
	return grid


def time_correlate(mags):
	correlate = []
	# go over range of lengths first ... 1, 2, 3, 4, 5, .. these are lengths in time
	for l in range(500):
		# grab regular average for this l
		avg_m = stats.mean(mags[:len(mags)-l])
		avg_m_2 = mean_squared(mags[:len(mags)-l])
		# calc split avg for this l
		total = 0
		for i in range(0,len(mags)-l):
			total = total + mags[i]*mags[i+l]
		avg_m_l = total/(len(mags)-l)
		numerator = avg_m_l - avg_m**2
		denominator = avg_m_2 - avg_m**2
		correlate.append(numerator/denominator)
	return correlate


def alternate_correlate(mags):
	t_max = len(mags)
	l_max = 500
	correlate = []
	for l in range(l_max):
		combined_sum = 0.0
		for i in range(len(mags)-l):
			combined_sum = combined_sum + mags[i]*mags[i+l]
		combined_sum = combined_sum/float(t_max-l)
		split_sum_m = 0.0
		split_sum_m_l = 0.0
		for i in range(len(mags)-l):
			split_sum_m = split_sum_m + mags[i]
		split_sum_m = split_sum_m/float(t_max-l)
		for i in range(len(mags)-l):
			split_sum_m_l = split_sum_m_l + mags[i+l]
		split_sum_m_l = split_sum_m_l/float(t_max-l)
		correlate.append(combined_sum - (split_sum_m*split_sum_m_l))
	return correlate


def hamil(grid,J):
	# calculate pair - pair energy from given state
	N = len(grid[0])
	total_E = 0
	for i in range(N):
		for j in range(N):
			if i == N-1:
				below = 0
			elif  i < N-1:
				below = i+1
			if j == N-1:
				right = 0
			elif j < N-1:
				right = j+1
			total_E = total_E + -J*grid[i][j]*(grid[below][j]+grid[i][right])

	return total_E


def mag(grid):
	# add up all spins on lattice
	total_mag = 0
	N = len(grid[0])
	for i in range(N):
		for j in range(N):
			total_mag = total_mag + grid[i][j]
	return total_mag


def mean_squared(data):
	total = 0
	for datum in data:
		total = total + datum**2
	return total / len(data)


np.set_printoptions(threshold=sys.maxsize)
J = 1.0
T = 2.24
N = 35
num_steps = 10000
grid = np.zeros((N,N))
for i in range(N):
	for j in range(N):
		r = rand.randint(0,1)
		if r == 0:
			r = -1
		grid[i,j] = r

# now generate graph for this initial grid
bond_graph = gen_bonds(grid,T,J)

for i in range(22,23):
	k = i/10.0
	print("TEMP: ", k)
	energies = []
	mags = []
	for step in range(num_steps):
		# go over bonds, update spin clusters on grid appropriately -> 1/2 prob of all flip for each cluster
		# clear bonds, make new bond_graph, and loop
		grid = update_clusters(bond_graph,grid,N)
		bond_graph = gen_bonds(grid,k,J)
		energies.append(hamil(grid,J))
		mags.append(mag(grid))
	#save_title = "T_"+str(int(k*10))+".png"
	#plt.title("T="+str(k))
	#plt.imshow(grid, interpolation='nearest')
	#plt.savefig(save_title)
	
	print(k, stats.mean(energies[num_steps/2:])/(float(N)*float(N)))
	print(k, stats.mean(mags[num_steps/2:]))

# now calculate correlation function using mag array
other_correlate = alternate_correlate(mags)
for i in range(len(other_correlate)):
	print(i, other_correlate[i])
