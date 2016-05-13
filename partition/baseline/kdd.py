#usage python3.3 kdd.py partitioninputfile edgeinputfile
import sys
import collections
import math

def main():

  if len(sys.argv) < 3:
    print("wrong num cmd line args. usage is python3.3 kdd.py partitioninputfile edgeinputfile")
    return 1

  partitionInput = sys.argv[1]
  edgeList = sys.argv[2]

  vertexCluster = collections.defaultdict(list)
  edgeCluster = collections.defaultdict(list)
  cutEdges = collections.defaultdict(list)
  numVertex = 0

  with open(partitionInput) as f:
    for line in f:
      if line[0] != "#":
        split = line.split()
        vertex = int(split[0])
        partition = int(split[1])
        vertexCluster[vertex] = partition
        numVertex += 1
        #print("vertex " + str(vertex) + " is in parition " + str(partition))

  with open(edgeList) as f:
    for line in f:
      if line[0] != "#":
        split = line.split()
        u = int(split[0])
        v = int(split[1])

        partitionNum = vertexCluster[u]
        if vertexCluster[u] == vertexCluster[v] and (v,u) not in edgeCluster[partitionNum]:
          edgeCluster[partitionNum].append((u,v))

        elif (u,v) not in cutEdges[partitionNum] and (v,u) not in cutEdges[partitionNum]:
          cutEdges[vertexCluster[u]].append((u,v))
          cutEdges[vertexCluster[v]].append((u,v))

  pairList = []
  
  for i in range(numVertex - 1):
    for j in range(i + 1, numVertex - 1):
      pairList.append((i,j))

  for i,j in pairList:
    setA = []
    setB = []

    if len(cutEdges[i]) == 0 or len(cutEdges[j]) == 0:
      continue

    if len(cutEdges[i]) % 2 == 0:
      sizeSetA = len(cutEdges[i]) / 2
      sizeSetB = sizeSetA
    else:
      sizeSetA = floor(len(cutEdges[i]) / 2)
      sizeSetB = sizeSetA + 1

    #for each edge between clusters
      if len(setA) == setSizeA:
        setB.append(u,v)
      elif len(setB) == setsizeB:
        setA.append(u,v)
      else:
        rand = random.choice([True, False])
        if rand:
          setA.append(u,v)
        else:
          setB.append(u,v)

    edgeCluster[vertexCluster[u]].append(setA)
    edgeCluster[vertexCluster[v]].append(setB)

if __name__ == "__main__":
  main()
