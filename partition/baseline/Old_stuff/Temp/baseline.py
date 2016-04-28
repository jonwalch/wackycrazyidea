#usage python3.3 baseline.py neighborinfo num_partition
import sys
#T
class Node:
  def __init__(self, num, neighbors):
    self.num = num
    self.neighbors = neighbors
    self.groupNum = -1

def baseline1():
  """
  if len(sys.argv) != 3:
    print("wrong num cmd line args. usage is python3.3 baseline.py file1 num_partition")
    return 1
  """
  nodeNumNeigh = "NeighborInfo.dat"#sys.argv[1]
  partition_num = 9#int(sys.argv[2])
  edge_num = 0
  nodes = [] #each node has .neighbor

  with open(nodeNumNeigh) as f:
    for line in f:
      split = line.split(":")

      curNodeNum = split[0]
      neighbors = []

      for i in split[1].split():
        edge_num += 1
        neighbors.append(i)

      newnode = Node(curNodeNum, neighbors)
      nodes.append(newnode)

  curGroupEdges = 0
  group_num = 0
  threshold = float(edge_num) / float(partition_num)

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
  input('ok')
