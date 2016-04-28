import random as rn 

### Parameters ###
global NUM_VERTICIES, NUM_PARTITIONS, MAX_VERT_PER_GROUP
global MIN_VERT_PER_GROUP, MAX_OUTGOING_EDGES_PER_GROUP 
global MIN_OUTGOING_EDGES_PER_GROUP,MAX_INNER_EDGES_PER_GROUP
global MIN_INNER_EDGES_PER_GROUPm, Even, NUM_EDGES


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

		count = 0 
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

def makeRandomGraph():
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
		if nodeFrom == nodeTo:
			if [str(nodeFrom), str(nodeTo)] not in graph:
					graph.append([str(nodeFrom), str(nodeTo)])
					i += 1
		else:
			graph.append([str(nodeFrom), str(nodeTo)])
			i += 1

	with open("Random_Graph.dat","w") as f:
		for line in sorted(graph, key=lambda x: int(x[0])):
			f.write("\t".join(line)+"\n")


def toFile(graph):
	f1 = "Edges.dat"
	f2 = "Optimal.dat"

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

    

print("Start")
if randomGraph == False:
	graph = makeGraph()
	toFile(graph)
	if Even:
		print("Number of Partitions: ",len(graph))
else:
	graph = makeRandomGraph()




