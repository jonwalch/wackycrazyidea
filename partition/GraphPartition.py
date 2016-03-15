import numpy as np
from sys import argv

class group(): #Contains group information
	def __init__(self):
		self.grp = -1
		self.inner = 0 
		self.outter = 0
		self.count = 0
		self.connceted_groups = []
		self.weight = []

	def connect(self, grp2):
		if grp2 == self.grp:
			self.inner += 1
		else:
			self.outter += 1
			if grp2 not in self.connceted_groups:
				self.connceted_groups.append(grp2)
				self.weight.append(1)
			else:
				indx = self.connceted_groups.index(grp2)
				self.weight[indx] += 1

def read_grps(filename):
	global vs, grps, edgs, groups
	
	with open(filename) as f:
		for line in f:
			line = line.split()
			if len(line) != 2:
				print("Error: ", line)
				continue
			vs.append(int(line[0]))
			grps.append(int(line[1]))
			edgs.append([])
			grp = int(line[1])
			if grp in groups:
				groups[grp].count += 1
			else:
				temp = group()
				temp.grp = grp
				temp.count = 1
				groups[grp] = temp

def read_edges(filename):
	global vs, grps, edgs, groups, total
	total = 0
	with open(filename) as f:
		for line in f:
			line = line.split()
			if len(line) != 2:
				print("Error: ", line)
				continue

			i1 = int(line[0])
			i2 = int(line[1])
			edgs[i1].append(vs[i2])
			groups[grps[i1]].connect(grps[i2]) 
			total += 1

def printGroupInfo():
	global vs, grps, edgs, groups, total
	c = 0
	c2 = 0
	print("======= Information =======")
	print("Graph:")
	print("Total number of vertices: ", len(vs))
	print("Total number of edges: ", total)
	print("Total number of groups: ", len(groups))
	
	print("")
	for grp in groups:
		print("Group: ", grp)
		group = groups[grp]
		print("Number of vertices: ", group.count)
		print("Number of internal edges: ", group.inner)
		print("Number of edges to other groups: ", group.outter)
		print("Connected groups: ", ",".join([str(i) for i in group.connceted_groups]))
		print("")
		c += group.outter
		c2 += group.inner
	print("Graph Summary:")
	print("Total number of vertices: ", len(vs))
	print("Total number of edges: ", total)
	print("Total number of groups: ", len(groups))
	print("Total number of edges b/w groups: ", c)
	print("% Edges that are b/w groups: ", 100*c/(c+c2))
	print("% Edges that are inside a group: ", 100*c2/(c+c2))

def writeGroupInfo():
	global vs, grps, edgs, groups, total
	f = open("report.dat","w")
	c = 0
	c2 = 0

	f.write("======= Information ======="+"\n")
	f.write("Graph:\n")
	f.write("Total number of vertices: "+str(len(vs))+"\n")
	f.write("Total number of edges: "+str(total)+"\n")
	f.write("Total number of groups: "+str(len(groups))+"\n")
	f.write("\n")

	for grp in groups:
		f.write("Group: "+str(grp)+"\n")
		group = groups[grp]
		f.write("Number of vertices: "+str(group.count)+"\n")
		f.write("Number of internal edges: "+str(group.inner)+"\n")
		f.write("Number of edges to other groups: "+str(group.outter)+"\n")
		f.write("Connected groups: "+ ",".join([str(i) for i in group.connceted_groups])+"\n")
		f.write("\n")
		c += group.outter
		c2 += group.inner

	f.write("Graph Summary:\n")
	f.write("Total number of vertices: "+str(len(vs))+"\n")
	f.write("Total number of edges: "+str(total)+"\n")
	f.write("Total number of groups: "+str(len(groups))+"\n")
	f.write("Total number of edges b/w groups: "+str(c)+"\n")
	f.write("% Edges that are b/w groups: "+str(100*c/(c+c2))+"\n")
	f.write("% Edges that are inside a group: "+str(100*c2/(c+c2))+"\n")
	f.close()

def mathematica():
	global vs, grps, edgs, groups, total

	filename = input("File name: ")
	with open(filename,"w") as f:	
		for i, v in enumerate(vs):
			for e in edgs[i]:
				f.write(str(v)+"->"+str(e)+" ")

	with open("color"+filename, "w") as f:
		for i, v in enumerate(vs):
			f.write(str(grps[i])+" ")
	
	with open("groups_"+filename, "w") as f:
		temp = groups
		for grp in temp:
			grp1 = temp[grp]
			for i in grp1.connceted_groups:
				f.write(str(grp1.grp)+"->"+str(i)+" ")

	with open("group_weight_"+filename, "w") as f:
		temp = groups
		for grp in temp:
			grp1 = temp[grp]
			for i in grp1.weight:
				f.write(str(i)+" ")

if __name__ == "__main__":
	global vs, grps, edgs, groups, total
	vs = []
	grps = []
	edgs = []
	groups = {}
	if len(argv) == 2:
		if argv[1] == "h":
				print("Input Format: GroupInfo EdgeInfo")
				exit(0)
	try:
		file1 = argv[1]
		file2 = argv[2]
	except:
		file1 = input("Group info: ")
		file2 = input("Edge info: ")

	print("Loading file1...")
	read_grps(file1)
	print("Loading file2...")
	read_edges(file2)

	printGroupInfo()
	writeGroupInfo()
	print("")
	#print("Making Mathematica Files")
	#mathematica()

	#input("\nENTERtoEXIT")