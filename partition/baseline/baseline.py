#usage python3.3 baseline.py edgelist neighborinfo num_partition
import sys
import pprint

class Node:
  def __init__(self, num, neighbors):
    self.num = num
    self.neighbors = neighbors
    self.groupNum = -1

def baseline1():
  pp = pprint.PrettyPrinter(indent=4)

  if len(sys.argv) != 4:
    print("wrong num cmd line args. usage is python3.3 baseline.py file1 file2 num_partition")
    return 1

  edgeListFile = sys.argv[1]
  nodeNumNeigh = sys.argv[2]

  partition_num = int(sys.argv[3])
  edge_num = sum(1 for line in open(edgeListFile))
  threshold = edge_num / partition_num

  nodes = [] #each node has .neighbor

  with open(nodeNumNeigh) as f:
    for line in f:
      split = line.split(":")
      curNodeNum = split[0]
      neighbors = []
      n = ""
      for i in range (1, len(split)):
        n += str(split[i]) + " "
        neighbors.append(split[i])
      print(str(curNodeNum) + ": " + n)
      n = ""
      newnode = Node(curNodeNum, neighbors)
      nodes.append(newnode)

  curGroupEdges = 0
  group_num = 0

  output = open("baseline1.out", 'w+')

  for i,node in enumerate(nodes):
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
