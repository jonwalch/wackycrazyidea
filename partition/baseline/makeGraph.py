import random as rn

NUM_VERTICIES = 100
MAX_VERT_PER_GROUP = 25
MAX_OUTGOING_EDGES_PER_GROUP = 5
MAX_INNER_EDGES_PER_GROUP = 30

graph = []
vlist = []
neb = {}
end = 0
c = 0
o = ""

while(len(vlist) < NUM_VERTICIES):
    temp = []
    x = rn.randint(1, MAX_VERT_PER_GROUP)
    if end+x > NUM_VERTICIES:
        x = NUM_VERTICIES - end -1
    for i in range(end, end+x):
        temp.append(str(i)+"\t"+str(i+1))
        vlist.append(str(i)+"\t"+str(c))
        neb[str(i)] = []
    vlist.append(str(end+x)+"\t"+str(c))
    neb[str(end+x)] = []
    graph.append(temp)
    end = end + x + 1
    c += 1
    print(str(c)," | ", len(vlist))
    #o = input(str(c))

print("")

InterGroup = []
for i, temp in enumerate(graph):
    j = 0
    while j < rn.randint(1, MAX_OUTGOING_EDGES_PER_GROUP):
        gID = rn.randint(0,len(graph)-1)
        nodeID = rn.randint(0,len(graph[gID])-1)
        if i != gID:
            nodeTO = graph[gID][nodeID].split("\t")[0]
            nodeFromID = rn.randint(0,len(graph[i])-1)
            nodeFrom = graph[i][nodeFromID].split()[0]
            InterGroup.append(str(nodeFrom)+"\t"+str(nodeTO))
            j += 1
    print(str(i),"/",str(len(graph)), j)
graph.append(InterGroup)

IntraGroup = []
for i, temp in enumerate(graph):
    j = 0
    while j < rn.randint(1, MAX_INNER_EDGES_PER_GROUP):
        gID = i#rn.randint(0,len(graph)-1)
        nodeID = rn.randint(0,len(graph[gID])-1)
    
        nodeTO = graph[gID][nodeID].split("\t")[0]
        nodeFromID = rn.randint(0,len(graph[i])-1)
        nodeFrom = graph[i][nodeFromID].split()[0]
        if nodeFrom == nodeTO:
         if (str(nodeFrom)+"\t"+str(nodeTO)) not in IntraGroup:
            IntraGroup.append(str(nodeFrom)+"\t"+str(nodeTO))
            j += 1
        else:
            IntraGroup.append(str(nodeFrom)+"\t"+str(nodeTO))
            j += 1

    print(str(i),"/",str(len(graph)), j)
graph.append(IntraGroup)

for i, temp in enumerate(graph):
    for line in temp:
        line = line.split()
        i1 = line[0]
        i2 = line[1]
        neb[i2].append(i1)


with open("NeighborInfo.dat","w") as f:
    for key in sorted(neb, key = int):
        i = key +":"+" ".join(neb[key])
        f.write(i+"\n")

with open("Edges.dat","w") as f:
    for g in graph:
        if [] not in g:
            for i in g:
                #print(i.replace("\t","->"))
                f.write(i+"\n")
            #print()

with open("Test.dat","w") as f:
    for i in vlist:
        f.write(i+"\n") 
