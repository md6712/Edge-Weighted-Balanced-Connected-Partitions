#pragma once

#include "_g.h"
#include "bmemory.h"

#define SIZE_OF_EDGES_BINARY 4

struct _node_t {
	int UB; 
	int LB;
	uint32_t forbidden_edges[SIZE_OF_EDGES_BINARY]; // size 1 for each 32 edges; 2 for each 64 edges; 3 for each 96 edges; 4 for each 128 edges
	uint32_t mst_edges[SIZE_OF_EDGES_BINARY];
	_node_t* next;
};


class bb
{
	_g* instance;

	// linked list to store all nodes
	_node_t* root;

	// benv	
	_benv_t* benv;

	// UB components 
	int** matrix;
	uint32_t edges_to_be_removed[SIZE_OF_EDGES_BINARY];
	uint32_t edges_selected[SIZE_OF_EDGES_BINARY];
	int* tree_weights; 
	int* num_vertices;
	int* num_edges;
	int** vertices; 
	int** edges_ub; 
	uint32_t** bin_vertices;
	int *decomp;
		
public:

	// MST components 
	ListGraph* graph;
	ListGraph::Node* nodes;
	ListGraph::Edge* edges;
	ListGraph::EdgeMap<int>* edge_weights;
	ListGraph::EdgeMap<bool>* mst;

	bb(_g* instance);
	~bb();

	// run the branch and bound
	bb* Run();

	// branch the node
	void branch(_node_t* node);

	// create root node; 
	_node_t* create_root_node();

	// compute mst
	void compute_mst(_node_t* node);

	// compute upper bound
	void compute_upper_bound(_node_t* node);

	// compute spanning k forest
	int compute_spanning_k_forest(_node_t* node);

	// traverse to decompose
	int traverse_to_decompose(_node_t* node, int vertex, int prev_vertex, int i, int largest_i, uint32_t* edge_to_be_removed);

	// print
	void print_node(_node_t* node);
	

};

