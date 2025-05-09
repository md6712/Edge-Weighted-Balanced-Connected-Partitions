#pragma once
#define PCST_LARGE 100000

class _pcst
{
public:

	bool log = false; 

	void* instance; // this is going to be the graph

	double UB; // upper bound for the weight of the tree

	int num_vertices;    // num_vertices of the original graph
	int num_edges;       // num_edges of the original graph

	int num_vertices_reduced; // num_vertices of the reduced graph
	int num_edges_reduced;    // num_edges of the reduced graph

	int(*edges)[2];				// 2D array to store the edges
	
	double* vertex_prize;		// array to store the prize of each vertex
	int* vertex_degree; 		// array to store the degree of each vertex
	double* edge_cost;			// array to store the cost of each edge
	double* edge_weight;		// array to store the weight of each edge
	
	bool* edge_active;			// array to store the active status of each edge
	bool* vertex_active;		// array to store the active status of each vertex	

	double** all_pairs_shortest_path; // 3D array to store the all pairs shortest path
	double*** all_pairs_shortest_path_with_capacity; // 3D array to store the all pairs shortest path

	bool* roots;					// array to store the roots of the tree


	// Aborecence instance
	// // // //  // // // // // // 
	int num_vertices_aborescence;    // num_vertices of the aborescence
	int num_arcs_aborescence;        // num_edges of the aborescence

	int num_vertices_aborescence_reduced; // num_vertices of the reduced aborescence
	int num_arcs_aborescence_reduced;    // num_edges of the reduced aborescence

	int (*arcs)[2];				// 2D array to store the arcs
	double* arc_cost;			// array to store the cost of each arc
	double* arc_weight;			// array to store the weight of each arc

	bool* arc_active;						// array to store the active status of each arc
	bool* vertex_aborescence_active;		// array to store the active status of each vertex


	// // Methods
	_pcst(void* g);
	~_pcst(void);

	// print instane 
	_pcst* print_instance();

	// set upper bound on the weight of the tree
	_pcst* set_upper_bound(double);

	// reset the number of vertices and edges
	_pcst* reset_num_vertices_edges();

	// reset the edge and vertex active status
	_pcst* reset_active_status();

	// set the prizes for the vertices
	_pcst* set_vertex_prizes(double*);

	// set the costs for the edges
	_pcst* set_edge_costs(double*);

	// set roots	
	_pcst* set_roots(int ,int*);

	// compute the all pairs shortest path
	_pcst* floyd_warshal();

	// check select budgeted shortest paths 
	_pcst* check_select_budgeted_shortest_paths();

	// computed budgeted shortest path between two vertices
	double budgeted_shortest_path(int u, int v, int w);

	// reduce the graph
	_pcst* reduce_graph();

	// compute degrees
	_pcst* compute_degree();

	// degree one prune 
	_pcst* degree_one_prune();

	// degree two prune
	_pcst* degree_two_prune();

	// solve
	_pcst* solve();	

	// print the graph
	_pcst* print();

	// print all pairs shortest path
	_pcst* print_all_pairs_shortest_path();

	// make the aborescence instance
	_pcst* make_aborescence_instance();

	// print the aborescence instance
	_pcst* print_aborescence_instance();



	// set log	
	_pcst* set_log(bool log) {
		this->log = log;
		return this;
	}


};

