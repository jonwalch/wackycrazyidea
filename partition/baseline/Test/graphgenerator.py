import random as rn

NUM_VERTICIES = 1000000
MAX_OUTGOING_EDGES_PER_GROUP = 20000000
MAX_VERT_PER_GROUP = 100000

graph = []
vlist = []
end = 0
c = 0
o = ""

while(len(vlist) < NUM_VERTICIES):
    temp = []
    x = rn.randint(1, MAX_VERT_PER_GROUP)
    for i in range(end, end+x):
        temp.append(str(i)+"\t"+str(i+1))
        vlist.append(str(i)+"\t"+str(c))
    vlist.append(str(end+x)+"\t"+str(c))
        
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

with open("arg2.dat","w") as f:
    for g in graph:
        if [] not in g:
            for i in g:
                #print(i.replace("\t","->"))
                f.write(i+"\n")
            #print()

with open("arg1.dat","w") as f:
    for i in vlist:
        f.write(i+"\n") 
