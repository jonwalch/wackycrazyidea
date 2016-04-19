#usage  with input file format vertex1 vertex2 python3.3 baseline.py neighborinfo num_partition baselineNum
import sys
from random import shuffle
import re
import collections

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

class Node:
  def __init__(self, num, neighbors):
    self.num = num
    self.neighbors = neighbors
    self.groupNum = -1

def main():

  if len(sys.argv) < 4:
    print("wrong num cmd line args. usage is python3.3 baseline.py file1 num_partition")
    return 1

  nodeNumNeigh = sys.argv[1]
  partition_num = int(sys.argv[2])
  edge_num = 0
  mode = sys.argv[3]
  nodes = collections.OrderedDict()

  if mode == "1":
    edge_num = baseline1(nodes, nodeNumNeigh)

  elif mode == "2":
    edge_num = baseline2(nodes, nodeNumNeigh)

  curGroupEdges = 0
  group_num = 0
  threshold = float(edge_num) / float(partition_num)
  result = []
  output = open("baseline" + mode + ".out", 'w+')

  for i,node in nodes.items():
    if curGroupEdges + len(node.neighbors) < threshold: 
      nodes[i].groupNum = group_num
      curGroupEdges += len(node.neighbors)
      result.append(str(i) + " " + str(group_num) + "\n")
    else:
      group_num += 1
      curGroupEdges = 0
      nodes[i].groupNum = group_num
      curGroupEdges += len(node.neighbors)
      result.append(str(i) + " " + str(group_num) + "\n")

  result = natural_sort(result)

  for i in result:
    output.write(i)

def baseline1(nodes, nodeNumNeigh): #default ordering
  edges = []

  with open(nodeNumNeigh) as f:
    for line in f:
      edges.append(line)

  #edges = list(set(edges)) #remove duplicate edges
  edges = natural_sort(edges)

  for i in range(len(edges)):
    split = edges[i].split()
    curNodeNum = split[0]
    destNode = split[1]

    if curNodeNum in nodes:
      nodes[curNodeNum].neighbors.append(destNode)
    else:
      newnode = Node(curNodeNum, [destNode])
      nodes[curNodeNum] = newnode

    if destNode not in nodes:
      newnode = Node(destNode, [])
      nodes[destNode] = newnode

  return len(edges)

def baseline2(nodes, nodeNumNeigh): #random ordering
  edges = []

  with open(nodeNumNeigh) as f:
    for line in f:
      edges.append(line)

  shuffle(edges)

  for i in range(len(edges)):
    split = edges[i].split()
    curNodeNum = split[0]
    destNode = split[1]

    if curNodeNum in nodes:
      nodes[curNodeNum].neighbors.append(destNode)
    else:
      newnode = Node(curNodeNum, [destNode])
      nodes[curNodeNum] = newnode

    if destNode not in nodes:
      newnode = Node(destNode, [])
      nodes[destNode] = newnode

  return len(edges)

if __name__ == "__main__":
  main()
