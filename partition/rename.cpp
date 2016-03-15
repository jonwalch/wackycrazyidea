#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;

//This program renames vertices based on their partitioning,
// then outputs a sorted list of edges.
int main(int argc, char *argv[]) {
	char garbage[2048];
	
	//NOTE: arguments are: <edge file> <partition file> <# edges in edge file>
	
	//Open partition file:
	ifstream partFile(argv[2]);
	if(partFile == NULL) {
		cerr << "Cannot open vertex partition file." << endl;
		return 0;
	}
	
	//Get number of vertices:
	int numVertices;
	partFile >> numVertices;
	
	//Read in vertex map:
	int * map = (int *) malloc(sizeof(int) * numVertices);
	for(int x = 0; x < numVertices; x++) {
		int vertex, garbage;
		partFile >> vertex;
		partFile >> garbage;
		map[vertex] = x;
	}
	
	//Open edge file and skip past comment lines:
	ifstream edgeFile(argv[1]);
	if(edgeFile == NULL) {
		cerr << "Cannot open original edge list file." << endl;
		return 0;
	}
	while(edgeFile.peek() == '#') {
		edgeFile.getline(garbage, 2048);
	}
	
	//Read in the list of edges, remapping their vertices as we go:
	int numEdges = 0;
	pair<int, int> * edges = (pair<int, int> *) malloc(sizeof(pair<int, int>) * atoi(argv[3]));
	while(!edgeFile.eof()) {
		int v1, v2;
		edgeFile >> v1;
		if(!edgeFile.fail()) {
			edgeFile >> v2;
			edges[numEdges].first = map[v1];
			edges[numEdges].second = map[v2];
			numEdges++;
		}
	}
	
	//Sanity check:
	if(numEdges != atoi(argv[3])) {
		cerr << "SANITY CHECK FAILED ~53: mismatched edge count (" << numEdges << " vs " << atoi(argv[3]) << ").\n";
	}
	
	//Sort edges:
	sort(edges, edges + numEdges);
	
	//Output:
	for(int x = 0; x < numEdges; x++) {
		cout << edges[x].first << "\t" << edges[x].second << "\n";
	}
	
	//Cleanup and exit:
	free(edges);
	free(map);
	return 0;
}
