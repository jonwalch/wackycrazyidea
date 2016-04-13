File = "Edges.dat"
graph = []
neb = {}
with open(File) as f:
	for line in f:
		graph.append(line)
		line = line.split()
		i1 = line[0]
		i2 = line[1]
		if i2 not in neb:
			neb[i2] = [i1]
		else:
			neb[i2].append(i1)

with open("NeighborInfoX.dat","w") as f:
    for key in sorted(neb, key = int):
        i = key +":"+" ".join(neb[key])
        f.write(i+"\n")