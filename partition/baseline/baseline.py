#usage  with input file format vertex1 vertex2 python3.3 baseline.py neighborinfo num_partition
import sys

class Node:
  def __init__(self, num, neighbors):
    self.num = num
    self.neighbors = neighbors
    self.groupNum = -1

def baseline1():

  if len(sys.argv) < 3:
    print("wrong num cmd line args. usage is python3.3 baseline.py file1 num_partition")
    return 1

  nodeNumNeigh = sys.argv[1]
  partition_num = int(sys.argv[2])
  edge_num = 0
  nodes = {}

  with open(nodeNumNeigh) as f:
    for line in f:
      split = line.split()
      curNodeNum = split[0]
      destNode = split[1]
      edge_num += 1

      if curNodeNum in nodes:
        nodes[curNodeNum].neighbors.append(destNode)
      else:
        newnode = Node(curNodeNum, [destNode])
        nodes[curNodeNum] = newnode

  curGroupEdges = 0
  group_num = 0
  threshold = float(edge_num) / float(partition_num)
  output = open("baseline1.out", 'w+')

  for i,node in nodes.items():
    #print("node neighbors # = " + str(len(node.neighbors)))
    if curGroupEdges + len(node.neighbors) < threshold: 
      nodes[i].groupNum = group_num
      curGroupEdges += len(node.neighbors)
      output.write(str(i) + " " + str(group_num) + "\n")
    else:
      group_num += 1
      curGroupEdges = 0
      nodes[i].groupNum = group_num
      curGroupEdges += len(node.neighbors)
      output.write(str(i) + " " + str(group_num) + "\n")

if __name__ == "__main__":
  baseline1()
