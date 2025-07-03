#pragma once

#include "_g.h"
#include "bmemory.h"
#include "sbbt.h"

#define SIZE_OF_EDGES_BINARY 20

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

	// benv
	_benv_t* benv_trees;

	// ssbt
	_multi_hashed_sbbt* sbbt;

	// ssbt for trees
	_multi_hashed_sbbt* sbbt_trees;

	// counters 
	int num_nodes;
	int num_kruskal_calls; 
	int num_ub_calls;

	// clock ticks
	clock_t kruskal_total_time;
	clock_t ub_total_time;



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

	// sbbt components
	_small_tree* temp_tree;
	_sbbt_node* curr_node, * new_node, * mother_node, * tmp_node = NULL;
	int l;
		
public:

	// MST components 
	ListGraph* graph;
	ListGraph::Node* nodes;
	ListGraph::Edge* edges;
	ListGraph::EdgeMap<int>* edge_weights;
	ListGraph::EdgeMap<bool>* mst;

	long long int num_trees_generated;

	bb(_g* instance);
	~bb();

	// run the branch and bound
	bb* Run(int time_limit);

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
	
	// hash tree
	__inline uint16_t hash_tree(_small_tree* tree);

	// populate all trees
	void populate_all_trees();

	// populate all trees investigate
	void populate_all_trees_investigate(_sbbt_node* node);

	// add singleton trees
	void add_singlton_trees();

	// add single-edge trees
	void add_single_edge_trees();

};

