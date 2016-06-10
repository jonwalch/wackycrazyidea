from sys import argv

def loadfile (fn):
	data = []
	with open(fn) as f:
		for line in f:
			data.append([eval(i) for i in line.split()])
	return data 

def sort(data):
	grp = {}
	for l in data:
		if l[2] in grp:
			grp[l[2]].append(l)
		else:
			grp[l[2]] = [l]
	return grp

def findNumNode(grp):
	NodeCount = {}
	allnodes = []
	count = 0 
	for key in grp:
		arr = grp[key]
		temp = []
		for i in arr:
			temp.append(i[0])
			temp.append(i[1])
			allnodes.append(i[0])
			allnodes.append(i[1])
		NodeCount[key] = len(set(temp))
		count += len(set(temp))
	allnodes = set(allnodes)
	return [allnodes, NodeCount, count]



if len(argv) > 1:
	if argv[1] == "h":
		print("python kddOutPut")
	else:
		fn = argv[1]
		data = loadfile(fn)
		grp = sort(data)
		nodes, Ncount ,count = findNumNode(grp)
		print(count - len(nodes))

