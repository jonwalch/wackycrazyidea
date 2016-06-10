#usage python3.3 kdd.py partitioninputfile edgeinputfile mode inputfiletype
#mode = 0 or 1, 0 is rand, 1 is optimal
#inputfiletype = 0 or 1, 0 is normal, 1 for .mtx
import sys
import collections
import math
import random
import numpy as np

def optimalPartition(pairList, cutEdges, edgeCluster):
  print("not done with optimal. on hold")
  maxs = []

  for i,j in pairList:
    setA = []
    setB = []

    sharedEdges = set(cutEdges[i]).intersection(set(cutEdges[j])) #find common elements between
    sharedEdges = list(sharedEdges)

    if len(sharedEdges) == 0: #if no edges to partition go to next iteration
      continue
    
    term1 = 0
    term2 = 0
    term3 = 0



    tonsOfGarbage = 0
    I = np.argmax(tonsOfGarbage)

  edgeCluster[i] += setA
  edgeCluster[j] += setB

def randPartition(pairList, cutEdges, edgeCluster):
  
  for i,j in pairList:
    setA = []
    setB = []

    sharedEdges = set(cutEdges[i]).intersection(set(cutEdges[j])) #find common elements between
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

def main():

  if len(sys.argv) < 4:
    print("wrong num cmd line args. usage is python3.3 kdd.py partitioninputfile edgeinputfile mode inputfiletype")
    return 1

  partitionInput = sys.argv[1]
  edgeList = sys.argv[2]
  mode = int(sys.argv[3])
  fileType = int(sys.argv[4])

  vertexCluster = collections.defaultdict(list)
  edgeCluster = collections.defaultdict(list)
  cutEdges = collections.defaultdict(list)
  numVertex = 0

  if fileType == 0:
    with open(partitionInput) as f:
      for line in f:
        if line[0] != "#" and line[0] != "%":
          split = line.split()
          vertex = int(split[0])
          partition = int(split[1])
          vertexCluster[vertex] = partition
          numVertex += 1
          #print("vertex " + str(vertex) + " is in parition " + str(partition))
  else:
    count = 1
    with open(partitionInput) as f:
      for line in f:
        if line[0] != "#" and line[0] != "%":
          split = line.split()
          vertex = count 
          partition = int(split[0])
          vertexCluster[vertex] = partition
          numVertex += 1
          count += 1

  edgesDone = 0
  edgesToCut = 0

  for i in vertexCluster:
    if type(i) != int:
      print(i)
      print("NOT INT")

  if fileType == 0:
    with open(edgeList) as f:
      for line in f:
        if line[0] != "#" and line[0] != "%":
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
  else:
    first = True
    offset = 0
    with open(edgeList) as f:
      for line in f:
        if line[0] != "#" and line[0] != "%":

          if first:
            split = line.split()
            offset = int(split[0])
            first = False
            continue

          split = line.split()
          u = int(split[0])
          v = int(split[1]) + offset
          partitionNum = vertexCluster[u]

          if type(vertexCluster[v]) != int:
            print(v)
            print(vertexCluster[v])
            return

          if vertexCluster[u] == vertexCluster[v] and (v,u) not in edgeCluster[partitionNum]:
            edgeCluster[partitionNum].append((u,v))
            edgesDone += 1

          elif (u,v) not in cutEdges[partitionNum] and (v,u) not in cutEdges[partitionNum]:
            cutEdges[vertexCluster[u]].append((u,v))
            cutEdges[vertexCluster[v]].append((u,v))
            edgesToCut += 1


  #print(edgesDone)
  #print(edgesToCut)
  #numEdges = edgesDone + edgesToCut

  pairList = []

  for i in range(len(cutEdges)):
    for j in range(i + 1, len(cutEdges)):
      pairList.append((i,j))

  if mode == 0:
    randPartition(pairList, cutEdges, edgeCluster)
  else:
    optimalPartition(pairList, cutEdges, edgeCluster)

  #partitionedEdges = 0
  output = open("kdd" + str(mode) + ".out", 'w+')

  for group in edgeCluster:
    for entry in edgeCluster[group]:
      output.write(str(entry[0]) + " " + str(entry[1]) + " " + str(group) + "\n")
      #format: edgeEndpoint1 edgeEndpoint2 group

if __name__ == "__main__":
  main()
