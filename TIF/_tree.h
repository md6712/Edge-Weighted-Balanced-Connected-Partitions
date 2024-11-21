#pragma once

#include <lemon/smart_graph.h>
#include <lemon/preflow.h>
#include <lemon/edmonds_karp.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/elevator.h>
#include <lemon/list_graph.h>


class _tree
{
	public:

	void *instance; // this is going to be the tree
	int num_vertices;
	uint32_t* bin_vertices; // binary number for vertices in the tree		
	
	int* edges; // contains the indices of the edges in the tree
	int* vertices; // contains the indices of the vertices in the tree

	// decomposition
	int* decomp; // contains the indices of the vertices in the tree

	int weight; // the weight of the tree

	bool isSpanningTree; // if the tree is a spanning tree over its vertices

	_tree(int *vertices, int num_vertices, uint32_t* bin_vertices, void *instance);
	~_tree();

	void computeMST();
	void printTree();
	void singltonTree(int vertex);
	void fixEdges(int* edges, int num_edges);

	_tree** splitIntoKTrees(int k);
	void traverseToDecompose(int vertex, int prev_vertex, int i, bool* edge_to_be_removed);

	int computeSpanningKForest(int k, void* instance, int** vertices, int** edges, bool* edge_to_be_removed, int* num_vertices, int* num_edges, int* tree_weights, uint32_t** bin_vertices);

};

