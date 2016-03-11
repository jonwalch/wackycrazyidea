#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>
using namespace std;

typedef struct {
	vector<int> incomingEdges;
	vector<int> outgoingEdges;
	
	int totalEdges;
	int newStart;
} vertex;

//Input: File with list of edges in directed graph.
//Output: Weighted graph with vertices split into multiple vertices, and file with vertex mapping.
int main(int argc, char *argv[]) {
	//NOTE: arguments are: <edge file> <# vertices>
	
	//Notes to self regarding conventions: put original edge first for each vertex in new graph.
	//			map vertices so that ones with incoming edges have consecutive IDs.
	
	//Get number of vertices
	int numVertices = atoi(argv[2]);
	
	//Open partition file:
	ifstream edgeFile(argv[1]);
	if(edgeFile == NULL) {
		cerr << "ERROR: Cannot open edge file." << endl;
		return 0;
	}
	
	//Skip comment lines:
	char garbage[2048];
	while(edgeFile.peek() == '#') {
		edgeFile.getline(garbage, 2048);
	}
	
	//Malloc original graph:
	vertex * G = (vertex *) malloc(sizeof(vertex) * numVertices);
	
	//Read in the list of edges:
	int orgEdges = 0;
	while(!edgeFile.eof()) {
		int v1, v2;
		edgeFile >> v1;
		if(!edgeFile.fail()) {
			edgeFile >> v2;
			G[v1].outgoingEdges.push_back(v2);
			orgEdges++;
		}
	}
	
	//Figure out number of edges (in either direction) for each vertex:
	for(int x = 0; x < numVertices; x++) {		
		//Add number of bi-directional and outgoing edges:
		G[x].totalEdges += G[x].outgoingEdges.size();
		
		//For target vertices, deal with incoming edge from here:
		for(int y = 0; y < G[x].outgoingEdges.size(); y++) {
			int z = G[x].outgoingEdges[y];
			//G[z].incomingEdges.push_back(x);
			G[z].totalEdges++;
		}
	}
	
	//Get starting index for each vertex's new vertex set:
	int newStart = 0;
	int totalEdgesInGraph = 0;
	for(int x = 0; x < numVertices; x++) {
		G[x].newStart = newStart;
		
		newStart += G[x].totalEdges;
		totalEdgesInGraph += (G[x].totalEdges * 2) - 1;
	}
	
	//Figure out vertex mapping for edge lists:
	for(int x = 0; x < numVertices; x++) {
		int offset = G[x].newStart + G[x].incomingEdges.size();
		for(int y = 0; y < G[x].outgoingEdges.size(); y++) {
			int z = G[x].outgoingEdges[y];
			int offset2 = G[z].newStart;
			G[x].outgoingEdges[y] = (G[z].incomingEdges.size() + offset2);
			G[z].incomingEdges.push_back(offset++);
		}
	}
	
	//Output vertex mapping
	stringstream ss1;
	ss1 << argv[1] << ".mapping";
	ofstream mapFile(ss1.str().c_str());
	mapFile << numVertices << "\n";
	for(int x = 0; x < numVertices; x++) {
		mapFile << G[x].newStart << "\n";
	}
	mapFile.close();
	
	//Output new graph
	stringstream ss2;
	ss2 << argv[1] << ".transgraph";
	ofstream graphFile(ss2.str().c_str());
	graphFile << newStart << " " << totalEdgesInGraph << " 001\n";
	for(int x = 0; x < numVertices; x++) {
		int v = G[x].newStart;
		
		for(int y = 0; y < G[x].incomingEdges.size(); y++) {
			graphFile << G[x].incomingEdges[y] << " 99999 ";
			if(y > 0) {
				graphFile << v - 1 << " 99999 ";
			}
			if(y < G[x].incomingEdges.size() - 1) {
				graphFile << v + 1 << " 99999 ";
			} else if(G[x].outgoingEdges.size() > 0) {
				graphFile << v + 1 << " 1 ";
			}
			graphFile << "\n";
			v++;
		}
		for(int y = 0; y < G[x].outgoingEdges.size(); y++) {
			graphFile << G[x].outgoingEdges[y] << " 99999 ";
			if(y > 0 || G[x].incomingEdges.size() > 0) {
				graphFile << v - 1 << " 1 ";
			}
			if(y < G[x].outgoingEdges.size() - 1) {
				graphFile << v + 1 << " 1 ";
			}
			graphFile << "\n";
			v++;
		}
	}
	graphFile.close();
	
	//Cleanup and exit:
	free(G);
	return 0;
}
