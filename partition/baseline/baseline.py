
partition_num = 100 #input parameter
edge_num = 100 #count num line sin input file

curGroupEdges = 0
group_num = 0

nodes = [] #each node has .neighbor

for i,node in enumerate(nodes):
  while curGroupEdges < edge_num / partition_num:
    nodes[i] = group_num
    curGroupEdges += len(node.neighbors)
  group_num += 1
  curGroupEdges = 0
