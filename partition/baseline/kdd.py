#usage python3.3 kdd.py partitioninputfile edgeinputfile
import sys
import collections
import math
import random

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

  edgesDone = 0
  edgesToCut = 0

  with open(edgeList) as f:
    for line in f:
      if line[0] != "#":
        split = line.split()
        u = int(split[0])
        v = int(split[1])
        partitionNum = vertexCluster[u]

        if vertexCluster[u] == vertexCluster[v] and (v,u) not in edgeCluster[partitionNum]:
          edgeCluster[partitionNum].append((u,v))
          edgesDone += 1

        elif (u,v) not in cutEdges[partitionNum] and (v,u) not in cutEdges[partitionNum]:
          cutEdges[vertexCluster[u]].append((u,v))
          cutEdges[vertexCluster[v]].append((u,v))
          edgesToCut += 1

        #else: #duplicate or v,u already accounted for
          #print((u,v))

  #print(edgesDone)
  #print(edgesToCut)
  numEdges = edgesDone + edgesToCut

  pairList = []

  for i in range(len(cutEdges)):
    for j in range(i + 1, len(cutEdges)):
      pairList.append((i,j))

  for i,j in pairList:
    setA = []
    setB = []

    sharedEdges = set(cutEdges[i]).intersection(set(cutEdges[j])) #find common elements between the two
    sharedEdges = list(sharedEdges)

    if len(sharedEdges) == 0: #if no edges to partition go to next iteration
      continue

    if len(sharedEdges)% 2 == 0:
      sizeSetA = int(len(sharedEdges) / 2)
      sizeSetB = sizeSetA
    else:
      sizeSetA = int(math.floor(len(sharedEdges) / 2))
      sizeSetB = sizeSetA + 1

    for (u,v) in sharedEdges:
    #after each append here, must remove edge from opposite set's cutEdges
      if len(setA) == sizeSetA:
        setB.append((u,v))
        cutEdges[i].remove((u,v))
        cutEdges[j].remove((u,v))
      elif len(setB) == sizeSetB:
        setA.append((u,v))
        cutEdges[i].remove((u,v))
        cutEdges[j].remove((u,v))
      else:
        rand = random.choice([True, False])
        if rand:
          setA.append((u,v))
          cutEdges[i].remove((u,v))
          cutEdges[j].remove((u,v))
        else:
          setB.append((u,v))
          cutEdges[i].remove((u,v))
          cutEdges[j].remove((u,v))

    edgeCluster[i] += setA
    edgeCluster[j] += setB

  partitionedEdges = 0

if __name__ == "__main__":
  main()
