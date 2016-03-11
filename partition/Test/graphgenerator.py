import random as rn

NUM_VERTICIES = 100000000
MAX_OUTGOING_EDGES_PER_GROUP = 5
MAX_VERT_PER_GROUP = 8000000
graph = []
vlist = []
end = 0
c = 0
o = ""

while(len(vlist)< NUM_VERTICIES):
    temp = []
    x = rn.randint(1,MAX_VERT_PER_GROUP)
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
for i, temp in enumerate(graph):
    j = 0
    while j < rn.randint(1,MAX_OUTGOING_EDGES_PER_GROUP):
        y = rn.randint(0,len(graph)-1)
        t = rn.randint(0,len(graph[y])-1)
        if i != y:
            y = graph[y][t].split("\t")[1]
            t = rn.randint(0,len(graph[i])-1)
            graph[i].append(str(y)+"\t"+str(graph[i][t].split()[0]))
            j += 1
    print(str(i),"/",str(len(graph)))

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
