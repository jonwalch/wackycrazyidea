#usage python3.3 baseline.py file1 file2 num_partition
import sys

class Node:
  def __init__(self, num, neighbors):
    self.num = num
    self.neighbors = neighbors
    self.groupNum = -1

if len(sys.argv) != 4:
  print("wrong num cmd line args. usage is python3.3 baseline.py file1 file2 num_partition")
  return 1

def baseline1():
  edgeListFile = sys.argv[1]
  nodeNumNeigh = sys.argv[2]

  partition_num = sys.argv[3] 
  edge_num = sum(1 for line in open(edgeListFile))

  nodes = [] #each node has .neighbor

  with open(nodeNumNeigh) as f:
    for line in f:
      split = line.split(":")
      tempNodeNum = split[0]
      neighbors = []
      for i in range (1, len(split)):
        neighbors.append(split[i])
      newnode = Node(tempNodeNum, neighbors)

  curGroupEdges = 0
  group_num = 0

  for i,node in enumerate(nodes):
    while curGroupEdges < edge_num / partition_num:
      nodes[i].groupNum = group_num
      curGroupEdges += len(node.neighbors)
    group_num += 1
  curGroupEdges = 0

if __name__ == "__main__":
  baseline1()
