#pragma once

#include <lemon/smart_graph.h>
#include <lemon/preflow.h>
#include <lemon/edmonds_karp.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/elevator.h>
#include <lemon/list_graph.h>

#define LARGE_WEIGHT 10000

class _tree
{
	public:

	void *instance; // this is going to be the tree

	// num vertices and edges
	int num_vertices;
	int num_edges; 

	uint32_t* bin_vertices; // binary number for vertices in the tree
	
	int* edges; // contains the indices of the edges in the tree
	int* vertices; // contains the indices of the vertices in the tree

	int* degree; // contains the degree of the vertices in the tree

	// decomposition
	int* decomp; // contains the decomp index of the vertices in the tree

	int weight; // the weight of the tree

	double reduced_cost; // the reduced cost of the tree

	bool isSpanningTree; // if the tree is a spanning tree over its vertices

	int num_forbidden_edges_in_mst; // the number of forbidden edges in the minumum spanning tree


	// constructor and destructor
	_tree(void *instance);
	~_tree();

	// print methods
	void PrintTree();

	// compute 
	
	// tree init methods 
	void SingltonTree(int vertex);

	void CopyVertices(uint32_t* bin_vertices);
	void CopyVertices(int num_vertices, int* vertices, uint32_t* bin_vertices);
	void CopyVertices(int num_vertices, int* vertices);
	
	void CopyEdges( int num_edges, int* edges);
	void InitTreeByEdges(int num_edges, int* edges);

	// Modify the tree
	void AddEdge(int edge);
	void RemoveEdge(int edge);
	void AddVertex(int vertex);
	void RemoveVertex(int vertex);


	// compute 
	int ComputeWeight();
	void RecomputeDegree();
	void ComputeMST(bool* forbidden_edges = NULL);

	// determine 
	bool IsSpanningTree(int edge);
	bool IsEdgeInTree(int edge);
	bool GetAnyIncidentVertex(int edge); // returns a vertex of the edge that is connected to the tree; if both vertices are connected, returns any of the two
	int IncidentVertex(int edge);  // returns a vertext  of the edge that is connected to the tree if the edge is connected; returns -1 if both vertices are connected
	bool IsConnected(); // is tree connected

	// split the tree into k trees
	_tree** SplitIntoKTrees(int k, bool* forbidden_edges);
	int TraverseToDecompose(int vertex, int prev_vertex, int i, int largest_i, bool* edge_to_be_removed);

	// 
	int ComputeSpanningKForest(int k, void* instance, int** vertices, int** edges, bool* edge_to_be_removed, int* num_vertices, int* num_edges, int* tree_weights, uint32_t** bin_vertices, bool *forbidden_edges);


};

