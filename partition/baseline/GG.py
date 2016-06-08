import random as rn 
import sys, time
from sys import argv
### Parameters ###
global NUM_VERTICIES, NUM_PARTITIONS, MAX_VERT_PER_GROUP
global MIN_VERT_PER_GROUP, MAX_OUTGOING_EDGES_PER_GROUP 
global MIN_OUTGOING_EDGES_PER_GROUP,MAX_INNER_EDGES_PER_GROUP
global MIN_INNER_EDGES_PER_GROUPm, Even, NUM_EDGES, randomGraph

"""
NUM_VERTICIES = 50

#Only if random graph
randomGraph = False
NUM_EDGES = 100

# if None random fill in the following:
Even = False # All partitions have the same number of Vert?
NUM_PARTITIONS =  10    # Only use if Even = True


MAX_VERT_PER_GROUP = 10
MIN_VERT_PER_GROUP = 10

MAX_OUTGOING_EDGES_PER_GROUP = 10
MIN_OUTGOING_EDGES_PER_GROUP = 10

MAX_INNER_EDGES_PER_GROUP = 10
MIN_INNER_EDGES_PER_GROUP = 10
"""

#Todo: add loop?
#####################
class partition():
	def __init__ (self):
		self.partID = -1
		self.vertices = [] 
		self.neighPart = [] #Neighbour partitions
		self.part = [] # Will store grouping info (optimal partition)
		self.edges = [] # all edges
		self.inner = [] # Internal edges 
		self.outter = [] # Out going edges
		self.numVert = -1
		self.numEdge = -1 

	#maks a chain of vertices
	def makeChain(self, size, start):
		for i in range(start, start+size):
			self.edges.append(str(i)+"\t"+str(i+1))
			self.vertices.append(i)
			self.part.append(str(i)+"\t"+str(self.partID))
		
		self.vertices.append(start+size)
		self.part.append(str(start+size)+"\t"+str(self.partID))

	#Will make connections b/w groups
	def InterGroup(self, size, graph):
		global MIN_OUTGOING_EDGES_PER_GROUP, MAX_OUTGOING_EDGES_PER_GROUP
		
		count = 0
		domain = rn.randint(MIN_OUTGOING_EDGES_PER_GROUP, MAX_OUTGOING_EDGES_PER_GROUP) 
		while count != domain:
			#print(count, "2", domain)
			partIndx = rn.randint(0,len(graph)-1) #Get target group index

			if graph[partIndx] != self:
				part = graph[partIndx] # target group
				
				nodeIndx = rn.randint(0, part.getSize()-1) #Get target node index
				nodeTo = part.vertices[nodeIndx]

				nodeIndx = rn.randint(0,self.getSize()-1) #Get self part's node index
				nodeFrom = self.vertices[nodeIndx]

				self.edges.append(str(nodeFrom)+"\t"+str(nodeTo))
				self.outter.append(str(nodeFrom)+"\t"+str(nodeTo))
				self.neighPart.append(part.partID)
				count += 1 
	
	#Will make connections within the partition:
	def IntraGroup(self):
		global MIN_INNER_EDGES_PER_GROUP, MAX_INNER_EDGES_PER_GROUP

		if self.getSize() == 1:
			self.edges.append(str(self.vertices[0])+"\t"+str(self.vertices[0]))
			self.inner.append(str(self.vertices[0])+"\t"+str(self.vertices[0]))
			return

		count = len(self.edges)  
		domain = rn.randint(MIN_INNER_EDGES_PER_GROUP, MAX_INNER_EDGES_PER_GROUP)
		while count != domain:
			#print(count,"3", domain)
			nodeIndx = rn.randint(0,self.getSize()-1) #Get target node index
			nodeTo = self.vertices[nodeIndx]
			nodeIndx = rn.randint(0,self.getSize()-1) #Get self part's node index
			nodeFrom = self.vertices[nodeIndx]
			#print(self.getSize(), nodeFrom, nodeTo)
			if nodeFrom == nodeTo:
				if (str(nodeFrom)+"\t"+str(nodeTo)) not in self.inner:
					self.edges.append(str(nodeFrom)+"\t"+str(nodeTo))
					self.inner.append(str(nodeFrom)+"\t"+str(nodeTo))
					count += 1
			else:
				self.edges.append(str(nodeFrom)+"\t"+str(nodeTo))
				self.inner.append(str(nodeFrom)+"\t"+str(nodeTo))
				count += 1

	def getSize(self):
		return len(self.vertices)

#XXXXXXXXXXXEnd ClassXXXXXXXXXXXXXXX
def error(x):
	print("\033[31mError\033[0m: "+ x)

def makeGraph():
	global NUM_VERTICIES, NUM_PARTITIONS, MAX_VERT_PER_GROUP
	global MIN_VERT_PER_GROUP, MAX_OUTGOING_EDGES_PER_GROUP 
	global MIN_OUTGOING_EDGES_PER_GROUP,MAX_INNER_EDGES_PER_GROUP
	global MIN_INNER_EDGES_PER_GROUPm, Even
	
	graph = []
	vertCount = 0 #Number of added verts
	start = 0 #First vert id
	
	if Even:
		if (NUM_VERTICIES%NUM_PARTITIONS) != 0:
			error("Can not distribute evenly!")
			exit(0)
		ratio = int(NUM_VERTICIES/NUM_PARTITIONS)-1
		for i in range(0,NUM_PARTITIONS):
			p = partition()
			p.partID = i
			p.makeChain(ratio, start)
			graph.append(p)
			start = start + ratio + 1
			
	else:
		i = 0
		while vertCount < NUM_VERTICIES:
			size = rn.randint(MIN_VERT_PER_GROUP, MAX_VERT_PER_GROUP)
			if start+size > NUM_VERTICIES:
				size = NUM_VERTICIES - start -1 #Fill last group till NUM_Vert
			p = partition()
			p.partID = i
			p.makeChain(size, start)
			graph.append(p)
			start = start + size + 1
			vertCount = vertCount + p.getSize()
			i += 1

	#"Make more connections"		
	for p in graph:
		size = MAX_VERT_PER_GROUP - p.getSize() - 1
		p.InterGroup(size, graph)
		p.IntraGroup()

	return graph


"""

def makeRandomGraph():
	global NUM_VERTICIES, NUM_EDGES
	edges = []
	graph = []
	for i in range(0, NUM_VERTICIES):
		for j in range(0, NUM_VERTICIES):
			if i != j:
				edges.append(str(i)+"\t"+str(j))
	i = 0 
	while(i < NUM_VERTICIES):
		j = rn.randint(0, len(edges))
		graph.append(edges.pop(j))
		i += 1
"""

def makeRandomGraphOld():
	global NUM_VERTICIES, NUM_EDGES
	vertices = []
	graph = []
	for i in range(0,NUM_VERTICIES):
		vertices.append(i)

	i = 0 
	while i != NUM_EDGES:
		nodeIndx = rn.randint(0,len(vertices)-1)
		nodeTo = vertices[nodeIndx]
		nodeIndx = rn.randint(0,len(vertices)-1)
		nodeFrom = vertices[nodeIndx]
		"""
		if nodeFrom == nodeTo:
			if [str(nodeFrom), str(nodeTo)] not in graph:
					graph.append([str(nodeFrom), str(nodeTo)])
					i += 1
		"""
		if nodeFrom != nodeTo and [str(nodeFrom), str(nodeTo)] not in graph:
			graph.append([str(nodeFrom), str(nodeTo)])
			i += 1

	with open("Random_Graph.dat","w") as f:
		for line in sorted(graph, key=lambda x: int(x[0])):
			f.write("\t".join(line)+"\n")

def toFile(graph, name = ""):
	f1 = "Edges%s.dat" %name
	f2 = "Optimal%s.dat" %name

	with open(f1,"w") as f:
		for p in graph: 
			if [] not in p.edges:
				for i in p.edges:
					f.write(i+"\n")
			else:
				error("toFile error")
	with open(f2,"w") as f:
		temp = []
		for p in graph:
			for i in p.part:
				temp.append(i.split("\t"))

		for i in sorted(temp,key=lambda x: int(x[0])):
			f.write("\t".join(i)+"\n")

def barabasi2(n, m0, m):
	#Load graph with initial nodes
	g = [] 
	total_edges = 0 
	for i in range(m0):
		temp = [int(x) for x in range(m0) if x != i] # Connect all nodes
		total_edges += len(temp)
		g.append(temp)

	for i in range(m0+1,n):
		sys.stdout.write("Progress: %d%%     \r" %(100*i/n))
		sys.stdout.flush()
		temp = [] #All connections will be saved here
		while(len(temp) < m):
			j = rn.randint(0,len(g)-1) #Pick random
			p = len(g[j])/total_edges
			R = rn.uniform(0,1)
			if p > R and j not in temp and j != i:
				temp.append(j)
				g[j].append(i)
				total_edges += 2
		g.append(temp)
	

	with open("Barabasi_Albert_Graph.dat","w") as f:
		for i, temp in enumerate(g):
			for j in temp:
				f.write(str(i)+"\t"+str(j)+"\n")
	#Test
	degree = {} 
	for i, temp in enumerate(g):
		#print(i, len(temp))
		if len(temp)in degree:
			degree[len(temp)] += 1
		else:
			degree[len(temp)] = 1

	with open("distribution_plot2.dat","w") as f:
		for key in degree:
			f.write(str(key)+"\t"+str(degree[key])+"\n")

	return g

def barabasi(n, m0, m):
	g = []
	deg = {}
	ti = time.time()
	for i in range(m0):
		for x in range(m0):
			if x != i:
				g.append([i, int(x)]) 
		deg[i] = m0 - 1

	for i in range(m0+1, n):
		sys.stdout.write("Progress: %d%%     \r" %(100*i/n))
		sys.stdout.flush()

		temp = 0
		arr = []
		deg[i] = 0
		while(temp < m):
			j = rn.randint(0,len(g)-1)
			j = g[j][1]
			if j != i and [i, j] not in arr:
				g.append([j, i])
				g.append([i, j])
				deg[i] += 1
				deg[j] += 1

				arr.append([i, j])
				temp += 1

	degree = {} 
	for key in deg:
		if deg[key] in degree:
			degree[deg[key]] += 1
		else:
			degree[deg[key]] = 1

	with open("distribution_plot%s.dat"%str(n),"w") as f:
		for key in degree:
			f.write(str(key)+"\t"+str(degree[key]/len(g))+"\n")
	print(time.time()-ti)


    
def options():
	global NUM_VERTICIES, NUM_PARTITIONS, MAX_VERT_PER_GROUP
	global MIN_VERT_PER_GROUP, MAX_OUTGOING_EDGES_PER_GROUP 
	global MIN_OUTGOING_EDGES_PER_GROUP,MAX_INNER_EDGES_PER_GROUP
	global MIN_INNER_EDGES_PER_GROUP, Even, NUM_EDGES, randomGraph
	print("Options")
	print("(0) RandomGraph | (1) OpimalRandomGraph |(2) EqualPartition")
	print("(3) Barabasi_Albert")

	op = eval(input("Input => "))
	randomGraph = False
	Even = False

	if 0 <= op <= 2:
		NUM_VERTICIES = eval(input("NUM_VERTICIES = "))
	if op == 0:
		randomGraph = True
		NUM_EDGES = eval(input("NUM_EDGES = "))
		if(NUM_VERTICIES**2 - NUM_VERTICIES) < NUM_EDGES:
			error("More edges requested than can be produced")
			print("(NUM_VERTICIES**2 - NUM_VERTICIES) < NUM_EDGES")
			exit(0)  
		if (NUM_EDGES < NUM_VERTICIES -1):
			print("Warning: (NUM_EDGES < NUM_VERTICIES -1)")
	elif op == 1 or op == 2:
		if op == 2:
			Even = True
			NUM_PARTITIONS = eval(input("NUM_PARTITIONS = "))
			if (NUM_VERTICIES%NUM_PARTITIONS) != 0:
				error("Can not distribute evenly!")
				exit(0)
			MAX_VERT_PER_GROUP = NUM_VERTICIES/NUM_PARTITIONS
			MIN_VERT_PER_GROUP = NUM_VERTICIES/NUM_PARTITIONS
		else:
			MAX_VERT_PER_GROUP = eval(input("MAX_VERT_PER_GROUP = "))
			MIN_VERT_PER_GROUP = eval(input("MIN_VERT_PER_GROUP = "))

		MAX_OUTGOING_EDGES_PER_GROUP = eval(input("MAX_OUTGOING_EDGES_PER_GROUP = "))
		MIN_OUTGOING_EDGES_PER_GROUP = eval(input("MIN_OUTGOING_EDGES_PER_GROUP = "))

		MAX_INNER_EDGES_PER_GROUP = eval(input("MAX_INNER_EDGES_PER_GROUP = "))
		MIN_INNER_EDGES_PER_GROUP = eval(input("MIN_INNER_EDGES_PER_GROUP = "))
	
	elif op == 3:
		n = eval(input("Total_num_nodes(n) = "))
		m0 = eval(input("Initial_num_nodes(m0) = "))
		m = eval(input("Degree_per_node(m < m0) = "))
		if m <= m0:
			barabasi(n, m0, m)
		else:
			error("m > m0 : Degree_per_node too large")
			op = input("Exit?(y | n) ")
			if op == "y" or op == "Y":
				exit(0)
			options()
		exit(0)
	else:
		error("Not valid choice!")
		op = input("Exit?(y | n) ")
		if op == "y" or op == "Y":
			exit(0)
		options()

print("Start")
print("GraphGen 1.0\n")

if len(argv) > 1:
	"""
	n = int(argv[1])
	m0 = int(argv[2])
	m = int(argv[3])
	if m <= m0:
		barabasi(n, m0, m)
	else:
		error("m > m0 : Degree_per_node too large")
		op = input("Exit?(y | n) ")
		if op == "y" or op == "Y":
			exit(0)
		options()
	exit(0)
	"""
	Even = True
	NUM_VERTICIES = int(argv[1])
	NUM_PARTITIONS = int(argv[2])
	MAX_VERT_PER_GROUP = NUM_VERTICIES/NUM_PARTITIONS
	MIN_VERT_PER_GROUP = NUM_VERTICIES/NUM_PARTITIONS
	MAX_OUTGOING_EDGES_PER_GROUP = int(argv[3])
	MIN_OUTGOING_EDGES_PER_GROUP = int(argv[3])
	MAX_INNER_EDGES_PER_GROUP = int(argv[4])
	MIN_INNER_EDGES_PER_GROUP = int(argv[4])
	graph = makeGraph()
	toFile(graph, argv[5])
	exit(0)


else:
	options()

if randomGraph == False:
	graph = makeGraph()
	toFile(graph)
	if Even:
		print("Number of Partitions: ",len(graph))
else:
	graph = makeRandomGraph()


