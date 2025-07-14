#include "_pcst.h"


#include <lemon/kruskal.h>
#include <lemon/list_graph.h>
#include <algorithm>
#include "_g.h"


// constructor
_pcst::_pcst(void* g) {
	// read graph 
	this->instance = g;
	_g* instance = (_g*)g;

	// set the number of vertices and edges
	this->num_vertices = instance->num_vertices;
	this->num_edges = instance->num_edges;

	// set the upper bound
	this->UB = instance->UB;

	// allocate memory for the edges and vertices
	this->edges = new int [this->num_edges][2];
	this->vertex_prize = new double[this->num_vertices];
	this->vertex_degree = new int[this->num_vertices];
	this->edge_cost = new double[this->num_edges];
	this->edge_weight = new double[this->num_edges];
	this->edge_active = new bool[this->num_edges];
	this->vertex_active = new bool[this->num_vertices];

	// allocate memory for the roots
	this->roots = new bool[this->num_vertices];

	// copy edges and their weights
	for (int e = 0; e < this->num_edges; e++) {
		this->edges[e][0] = instance->edges[e][0];
		this->edges[e][1] = instance->edges[e][1];
		this->edge_weight[e] = instance->edges[e][2];
	}
	
	// allocate memory for the all pairs shortest path
	this->all_pairs_shortest_path = new double* [this->num_vertices];
	this->all_pairs_shortest_path_with_capacity = new double** [this->num_vertices];
	for (int i = 0; i < this->num_vertices; i++) {
		this->all_pairs_shortest_path[i] = new double [this->num_vertices];
		this->all_pairs_shortest_path_with_capacity[i] = new double* [this->num_vertices];
		for (int j = 0; j < this->num_vertices; j++) {
			this->all_pairs_shortest_path_with_capacity[i][j] = new double[this->UB+1];
		}
	}

	// set the number of vertices and arcs for aborescence instance
	this->num_vertices_aborescence = this->num_vertices + 1 ;
	this->num_arcs_aborescence = this->num_edges *2 + this->num_vertices;

	// allocate memory for the arcs and vertices
	this->arcs = new int[this->num_arcs_aborescence][2];
	this->arc_cost = new double[this->num_arcs_aborescence];
	this->arc_weight = new double[this->num_arcs_aborescence];

	// allocate memory for the active status of the arcs and vertices
	this->arc_active = new bool[this->num_arcs_aborescence];
	this->vertex_aborescence_active = new bool[this->num_vertices_aborescence];

	// copy the arcs and their weights
	for (int e = 0; e < this->num_edges; e++) {
		// first arc in position e
		this->arcs[e][0] = instance->edges[e][0];
		this->arcs[e][1] = instance->edges[e][1];
		this->arc_weight[e] = instance->edges[e][2];
		
		// second arc in position e+ num_edges
		this->arcs[e + this->num_edges][0] = instance->edges[e][1];
		this->arcs[e + this->num_edges][1] = instance->edges[e][0];
		this->arc_weight[e + this->num_edges] = instance->edges[e][2];
	}

	// copy the arcs from root to vertices
	for (int v = 0; v < this->num_vertices; v++) {
		this->arcs[v + 2 * this->num_edges][0] = this->num_vertices;
		this->arcs[v + 2 * this->num_edges][1] = v;
		this->arc_weight[v + 2 * this->num_edges] = 0;
	}



	// define graph
	this->graph_mst = new ListGraph();

	// mst components
	// define nodes 
	nodes_mst = new ListGraph::Node[num_vertices];

	// include the nodes
	for (int i = 0; i < num_vertices; i++) {
		nodes_mst[i] = this->graph_mst->addNode();
	}

	// define edges
	edges_mst = new ListGraph::Edge[instance->num_vertices * (instance->num_vertices - 1) / 2];

	// define edge_weight
	edge_weights_mst = new ListGraph::EdgeMap<double>(*graph_mst);

	int e = 0; // edge counter

	// for each pair of vertices there is an edge 
	for (int u = 0; u < num_vertices; u++) {
		for (int v = u + 1; v < num_vertices; v++) {
			edges_mst[e++] = this->graph_mst->addEdge(nodes_mst[u], nodes_mst[v]);			
		}
	}

	// define mst 
	mst = new ListGraph::EdgeMap<bool>(*graph_mst, false);

	// shortest path vertex choice
	this->shortest_path_vertex_choice = new int[this->num_vertices];

}

// destructor
_pcst::~_pcst() {
	delete[] this->edges;
	delete[] this->vertex_prize;
	delete[] this->vertex_degree;
	delete[] this->edge_cost;
	delete[] this->edge_weight;
	delete[] this->edge_active;
	delete[] this->vertex_active;

	delete[] this->roots;
	for (int i = 0; i < this->num_vertices; i++) {
		for (int j = 0; j < this->num_vertices; j++) {
			delete[] this->all_pairs_shortest_path_with_capacity[i][j];
		}
		delete[] this->all_pairs_shortest_path[i];
		delete[] this->all_pairs_shortest_path_with_capacity[i];
	}
	delete[] this->all_pairs_shortest_path;
	delete[] this->all_pairs_shortest_path_with_capacity;

	// delete aborescence instance
	delete[] this->arcs;
	delete[] this->arc_cost;
	delete[] this->arc_weight;
	delete[] this->arc_active;
	delete[] this->vertex_aborescence_active;


	// delete the graph_mst components
	delete[] this->nodes_mst;
	delete[] this->edges_mst;
	delete this->edge_weights_mst;
	delete this->mst;
	delete this->graph_mst;

	// delete shortest path vertex choice
	delete[] this->shortest_path_vertex_choice;

}

// set upper bound on the weight of the tree
_pcst* _pcst::set_upper_bound(double UB) {
	this->UB = UB;
	return this;
}

// reset the edge and vertex active status
_pcst* _pcst::reset_active_status() {
	memset(this->edge_active, 1, sizeof(bool) * this->num_edges);
	memset(this->vertex_active, 1, sizeof(bool) * this->num_vertices);
	this->num_vertices_reduced = this->num_vertices;
	this->num_edges_reduced = this->num_edges;
	return this;
}

// set vertex prizes
_pcst* _pcst::set_vertex_prizes(double* vertex_prizes) {
	memcpy(this->vertex_prize, vertex_prizes, sizeof(double) * this->num_vertices);
	return this;
}

// set edge costs
_pcst* _pcst::set_edge_costs(double* edge_costs) {
	memcpy(this->edge_cost, edge_costs, sizeof(double) * this->num_edges);
	return this;
}

// set roots
_pcst* _pcst::set_roots(int num_roots, int* roots) {
	memset(this->roots, 0, sizeof(bool) * this->num_vertices);
	// for each root set the active status to true
	for (int r = 0; r < num_roots; r++) {
		this->roots[roots[r]] = true;
	}

	return this;
}

// reset the number of vertices and edges
_pcst* _pcst::reset_num_vertices_edges() {
	this->num_vertices_reduced = this->num_vertices;
	this->num_edges_reduced = this->num_edges;
	return this;
}

// print 
_pcst* _pcst::print() {

	// print tree weight upper bound
	printf("Upper bound: %4.2f\n", this->UB);

	// show the number of vertices and edges
	printf("Number of vertices: %d\n", this->num_vertices);
	printf("Number of edges: %d\n", this->num_edges);

	// roots of the tree
	printf("Roots: ");
	for (int v = 0; v < this->num_vertices; v++) {
		if (this->roots[v]) {
			printf("%d ", v);
		}		
	}
	printf("\n");

	// print number of vertices and edges in the reduced graph
	printf("Number of vertices in the reduced graph: %d\n", this->num_vertices_reduced);
	printf("Number of edges in the reduced graph: %d\n", this->num_edges_reduced);
	
	// print vertices 
	printf("Vertices:\n");
	for (int v = 0; v < this->num_vertices; v++) {
		if (this->vertex_active[v]) {
			printf("(%d %4.2f)\n", v, this->vertex_prize[v]);
		}		
	}

	// print the edges
	printf("Edges:\n");
	for (int e = 0; e < this->num_edges; e++) {
		if (this->edge_active[e]) {
			printf("(%d %d %4.2f %4.2f)\n", this->edges[e][0], this->edges[e][1], this->edge_weight[e], this->edge_cost[e]);
		}
	}

	return this;
}

// reduce the graph
_pcst* _pcst::reduce_graph() {
	
	// iterator
	int itr = 0;

	// track the number of vertices and edges before reduction
	int n_vertices_before = this->num_vertices;
	int n_edges_before = this->num_edges;

	// n_vertices and edges after reduction
	int n_vertices_after = this->num_vertices;
	int n_edges_after = this->num_edges;

	// delete all edges that are heavier than UB
	for (int e = 0; e < this->num_edges; e++) {
		if (this->edge_weight[e] > this->UB) {
			this->edge_active[e] = false;
			this->num_edges_reduced--;
			// print edge removed
			if (log) {
				printf("Delete edge %d %d with weight %4.2f > UB %4.2f\n", this->edges[e][0], this->edges[e][1], this->edge_weight[e], this->UB);
			}
		}
	}

	// delete all edges that joint with any of its adjacent edges is heavier than UB
	for (int e = 0; e < this->num_edges; e++) {
		if (this->edge_active[e]) {
			bool candidate_for_reduction = true;
			int u = this->edges[e][0];
			int v = this->edges[e][1];
			for (int ep = 0; ep < this->num_edges; ep++) {
				if (this->edge_active[ep]) {
					if (this->edges[ep][0] == u || this->edges[ep][1] == u || this->edges[ep][0] == v || this->edges[ep][1] == v) {
						if (this->edge_weight[e] + this->edge_weight[ep] <= this->UB) {
							candidate_for_reduction = false;
						}
					}
				}
			}

			if (candidate_for_reduction) {
				this->edge_active[e] = false;
				this->num_edges_reduced--;
				// print edge removed
				if (log) {
					printf("Delete edge (%d %d) with weight %4.2f > UB %4.2f\n", this->edges[e][0], this->edges[e][1], this->edge_weight[e], this->UB);
				}
			}
		}
	}
	
	while (itr++ < 100) {

		// remember the number of vertices and edges before reduction
		n_vertices_before = this->num_vertices_reduced;
		n_edges_before = this->num_edges_reduced;

		// compute the all pairs shortest path
		floyd_warshal();

		// compute budgeted shortest path for the heavy pairs (here we select those pairs with shortest path larger than their cost)
		check_select_budgeted_shortest_paths();

		// for each edge compute path with highest reward 


		// print the graph
		if (log) this->print();

		// compute degrees
		this->compute_degree();

		// prune degree one
		this->degree_one_prune();

		// print the graph
		if (log)  this->print();

		// prune degree two	
		this->degree_two_prune();

		// print the graph
		if (log)  this->print();

		// set number of vertices and edges after reduction
		n_vertices_after = this->num_vertices_reduced;
		n_edges_after = this->num_edges_reduced;

		// check if the number of vertices and edges has changed
		if (n_vertices_before == n_vertices_after && n_edges_before == n_edges_after) {
			break;
		}
	}

	// print the number of vertices and edges after reduction
	if (true) {
		printf("Number of vertices after reduction: %d\n", this->num_vertices_reduced);
		printf("Number of edges after reduction: %d\n", this->num_edges_reduced);
	}

	return this;
}

// compute the all pairs shortest path
_pcst* _pcst::floyd_warshal() {

	// memset the all pairs shortest path to a very large number
	for (int i = 0; i < this->num_vertices; i++) {
		for (int j = 0; j < this->num_vertices; j++) {
			this->all_pairs_shortest_path[i][j] = PCST_LARGE;
			for (int w = 0; w < this->UB + 1; w++) {				
				this->all_pairs_shortest_path_with_capacity[i][j][w] = PCST_LARGE;
			}
		}
	}

	// loop over edges and update the shortest paths
	for (int e = 0; e < this->num_edges; e++) {		
		if (this->edge_active[e]) {
			if (this->edge_cost[e] < this->all_pairs_shortest_path[this->edges[e][0]][this->edges[e][1]]) {
				this->all_pairs_shortest_path[this->edges[e][0]][this->edges[e][1]] = this->edge_cost[e];
				this->all_pairs_shortest_path[this->edges[e][1]][this->edges[e][0]] = this->edge_cost[e];
			}
		}		
	}

	if (log) print_all_pairs_shortest_path();

	// loop over u, v , k and update the shortest path
	for (int k = 0; k < this->num_vertices; k++) {
		for (int u = 0; u < this->num_vertices; u++) {
			for (int v = 0; v < this->num_vertices; v++) {
				if (this->all_pairs_shortest_path[u][k] + this->all_pairs_shortest_path[k][v] < this->all_pairs_shortest_path[u][v]
					) {
					this->all_pairs_shortest_path[u][v] = this->all_pairs_shortest_path[u][k] + this->all_pairs_shortest_path[k][v];
				}
			}
		}
	}

	if (log) print_all_pairs_shortest_path();

	return this;
}

// check select budgeted shortest paths
_pcst* _pcst::check_select_budgeted_shortest_paths() {

	// for each edge, check if we want to compute budgeted shortest path. 
	int u, v;

	for (int e = 0; e < this->num_edges; e++) {
		if (this->edge_active[e]) {
			// define u, v of the edge
			u = this->edges[e][0];
			v = this->edges[e][1];
			// check if cost of the edge is larger than the shortest path computed 

			if (this->edge_cost[e] > this->all_pairs_shortest_path[u][v]) {

				// compute the budgeted shortest path
				double cost = budgeted_shortest_path(u, v, this->edge_weight[e], false);

				if (cost < this->edge_cost[e]) {
					// if the cost is less than the edge cost, delete the edge
					this->edge_active[e] = false;
					this->num_edges_reduced--;
					// print edge removed
					if (log)
						printf("Delete edge %d %d with cost %4.2f > %4.2f\n", u, v, this->edge_cost[e], cost);
				}				
			}
		}
	}

	return this;
}

// compute budgeted shortest path between two vertices
double _pcst::budgeted_shortest_path(int u, int v, int w, bool positive_prizes) {
	
	// start from u check all edges adjacent to u
	// check if the edge is active
	if (this->all_pairs_shortest_path_with_capacity[u][v][w] < PCST_LARGE) {
		return this->all_pairs_shortest_path_with_capacity[u][v][w];
	}

	double best_cost = PCST_LARGE;

	// the edge choice for vertex e is not set
	shortest_path_vertex_choice[u] = -1;

	for (int e = 0; e < this->num_edges; e++) {
		if (this->edges[e][0] == u || this->edges[e][1] == u) {
			// if the edge is active, check the edge
			if (this->edge_active[e]) {
				if (this->edge_weight[e] > w) {					
					continue; 
				}

				// compute the other side of edge
				int vp = (this->edges[e][0] == u) ? this->edges[e][1] : this->edges[e][0];

				if (v == vp) {
					// if the edge is adjacent to v, return the edge cost					
					if (this->edge_cost[e] < best_cost) {
						best_cost = this->edge_cost[e];
						shortest_path_vertex_choice[u] = e;
					}
				}
				else {
					// first compute the cost to go
					double cost_to_go = budgeted_shortest_path(vp, v, w - this->edge_weight[e],false);					

					// add the cost of the edge and shared prize of the edge
					if (cost_to_go < PCST_LARGE) {
						cost_to_go += this->edge_cost[e];
						if (this->vertex_prize[vp] < 0) // we cannot add positive prizes as it would not be dominant
							cost_to_go -= this->vertex_prize[vp];
						else {
							if (positive_prizes) { // when positive prizes are allowed, we can add them. this can be true in some procedures.
								cost_to_go -= this->vertex_prize[vp];
							}							
						}
					}

					

					// update the best cost
					if (cost_to_go < best_cost) {
						best_cost = cost_to_go;
						
						// update the edge choice						
						shortest_path_vertex_choice[u] = e;
					}					
				}
			}
		}
	}

	if (best_cost == PCST_LARGE) {
		best_cost -= 0.00001;
	}

	// update the best cost
	if (best_cost < PCST_LARGE) {
		this->all_pairs_shortest_path_with_capacity[u][v][w] = best_cost;		
		this->all_pairs_shortest_path_with_capacity[v][u][w] = best_cost;
	}

	//printf("Best cost from %d to %d with capacity %d is %4.2f\n", u, v, w, best_cost);

	return best_cost;
}

// print all pairs shortest path
_pcst* _pcst::print_all_pairs_shortest_path() {
	printf("All pairs shortest path:\n");
	printf("-> \t");
	for (int v = 0; v < this->num_vertices; v++) {		
		if (vertex_active[v]) {
			printf("%d: \t", v);
		}
	}
	printf("\n");

	for (int u = 0; u < this->num_vertices; u++) {
		if (vertex_active[u]) {
			printf("%d: \t", u);
			for (int v = 0; v < this->num_vertices; v++) {
				if (vertex_active[v]) {
					if (this->all_pairs_shortest_path[u][v] < 1000)
						printf("%4.2f \t", this->all_pairs_shortest_path[u][v]);
					else {
						printf("inf \t");
					}
				}
			}
			printf("\n");
		}
	}
	return this;
}

// compute vertex degree
_pcst* _pcst::compute_degree() {

	// memset the vertex degree to 0
	memset(this->vertex_degree, 0, sizeof(int) * this->num_vertices);	

	// loop over the edges
	for (int e = 0; e < this->num_edges; e++) {
		if (this->edge_active[e]) {
			this->vertex_degree[this->edges[e][0]]++;
			this->vertex_degree[this->edges[e][1]]++;
		}
	}

	//// print degrees 
	//printf("Degrees:\n");
	//for (int v = 0; v < this->num_vertices; v++) {
	//	if (this->vertex_active[v]) {
	//		printf("%d (%d)\n", v, this->vertex_degree[v]);
	//	}
	//}

	return this;
	
}

// degree 1 prune 
_pcst* _pcst::degree_one_prune() {
	// loop over the vertices
	for (int v = 0; v < this->num_vertices; v++) {

		// check if v is root
		bool is_root = false;
		
		// if v is root, continue
		if (this->roots[v]) {
			continue;
		}

		if (this->vertex_active[v] && this->vertex_degree[v] == 1) {
			// loop over the edges to find adjacent edge
			for (int e = 0; e < this->num_edges; e++) {
				if (this->edges[e][0] == v || this->edges[e][1] == v) {
					// if the adjacent edge is active, check the edge
					if (this->edge_active[e]) {
						// if the prize of the vertex is less than the cost of the edge, deactivate the edge
						if (this->vertex_prize[v] < this->edge_cost[e]) {
							this->edge_active[e] = false;
							this->num_edges_reduced--;
							this->vertex_active[v] = false;
							this->num_vertices_reduced--;	

							// print the delete
							if (log) 
								printf("Delete edge %d %d with cost %4.2f > %4.2f\n", this->edges[e][0], this->edges[e][1], this->edge_cost[e], this->vertex_prize[v]);
						}
					}					
					break;
				}				
			}
		}
	}
	return this;
}

// degree 2 prune 
// TODO : This is not correct
_pcst* _pcst::degree_two_prune() {
	// loop over the vertices
	for (int v = 0; v < this->num_vertices; v++) {
		
		// if v is root, continue
		if (this->roots[v]) {
			continue;
		}


		if (this->vertex_active[v] && this->vertex_degree[v] == 2) {
			int e1 = -1, e2 = -1; 
			// loop over the edges to find adjacent edge
			for (int e = 0; e < this->num_edges; e++) {
				if (this->edges[e][0] == v || this->edges[e][1] == v) {
					// if the adjacent edge is active, check the edge
					if (this->edge_active[e]) {
						if (e1 == -1) {
							e1 = e;
						}
						else {
							e2 = e;
						}
					}
				}
			}

			// check if e1 + v is profitable
			if (edge_cost[e1] - vertex_prize[v] < 0) {
				continue;
			}

			// check if e2 + v is profitable
			if (edge_cost[e2] - vertex_prize[v] < 0) {
				continue;
			}

			// check if the edges are found
			if (e1 != -1 && e2 != -1) {
				// find the vertices of u and w
				int u = this->edges[e1][0] == v ? this->edges[e1][1] : this->edges[e1][0];
				int w = this->edges[e2][0] == v ? this->edges[e2][1] : this->edges[e2][0];
				
				if (this->edge_cost[e1] + this->edge_cost[e2] - this->vertex_prize[v] > this->all_pairs_shortest_path[u][w]) {					
					//printf("Potentially Delete %d %d %d\n", u, v, w);
					// compute cost of budgeted shortest path
					double cost = budgeted_shortest_path(u, w, this->edge_weight[e1] + this->edge_weight[e2], false);

					// check if the cost of budgeted shortest path is less than the cost of the edge
					if (cost < this->edge_cost[e1] + this->edge_cost[e2]) {
						// delete the edges
						this->edge_active[e1] = false;
						this->edge_active[e2] = false;
						this->num_edges_reduced -= 2;
						this->vertex_active[v] = false;
						this->num_vertices_reduced--;
						
						// print the delete
						if (log)
							printf("Delete edges %d %d and %d %d with cost %4.2f + %4.2f > %4.2f\n", this->edges[e1][0], this->edges[e1][1], this->edges[e2][0], this->edges[e2][1], this->edge_cost[e1], this->edge_cost[e2], cost);						
					}					
				}

			}
		}
	}
	return this;
}

// print instance 
_pcst* _pcst::print_instance() {
	// print the instance
	printf("PCST instance:\n");
	printf("Number of vertices: %d\n", this->num_vertices);
	printf("Number of edges: %d\n", this->num_edges);
	printf("Upper bound: %4.2f\n", this->UB);

	// print the vertices
	for (int v = 0; v < this->num_vertices; v++) {
		printf("Vertex %d: %4.2f\n", v, this->vertex_prize[v]);
	}
	// print the edges
	for (int e = 0; e < this->num_edges; e++) {
		printf("Edge %d: (%d %d) %4.2f\n", e, this->edges[e][0], this->edges[e][1], this->edge_cost[e]);
	}
	return this;
}

// make the aborescence instance
_pcst* _pcst::make_aborescence_instance() {
	// loop over the edges

	// reduced number of vertices and arcs are zero
	this->num_vertices_aborescence_reduced = 0;
	this->num_arcs_aborescence_reduced = 0;

	// for each arc, if the associated edge is active, add make the arc active and compute its cost 
	for (int e = 0; e < this->num_edges; e++) {
		if (this->edge_active[e]) {
			int a = e; 						
			this->arc_active[a] = true;
			this->arc_cost[a] = this->edge_cost[e] - this->vertex_prize[arcs[a][1]];
			this->num_arcs_aborescence_reduced++;
			

			a = e + this->num_edges;			
			this->arc_active[a] = true;
			this->arc_cost[a] = this->edge_cost[e] - this->vertex_prize[arcs[a][1]];
			this->num_arcs_aborescence_reduced++;				
		}
		else {
			// are not active
			this->arc_active[e] = false; 
			this->arc_active[e + this->num_edges] = false;
		}
	}

	// for each vertex, if the vertex is active, add make the vertex active
	for (int v = 0; v < this->num_vertices; v++) {
		if (this->vertex_active[v]) {
			this->vertex_aborescence_active[v] = true;
			this->num_vertices_aborescence_reduced++;
		}
		else {
			this->vertex_aborescence_active[v] = false;
		}
	}

	// create aborescence root
	this->vertex_aborescence_active[this->num_vertices] = true;
	this->num_vertices_aborescence_reduced++;

	// create arcs from the root to the vertices
	for (int v = 0; v < this->num_vertices; v++) {
		if (this->vertex_aborescence_active[v] /*&& this->vertex_prize[v] > 0*/) {
			int a = this->num_edges * 2 + v;
			this->arc_active[a] = true;
			this->arc_cost[a] = -this->vertex_prize[v];
			this->num_arcs_aborescence_reduced++;
		}
		else {
			this->arc_active[this->num_edges * 2 + v] = false;
		}
	}

	return this;
}

// print the aborescence instance
_pcst* _pcst::print_aborescence_instance() {
	printf("Aborescence instance:\n");
	printf("Number of vertices: %d\n", this->num_vertices_aborescence_reduced);
	printf("Number of arcs: %d\n", this->num_arcs_aborescence_reduced);
	printf("Vertices: ");
	for (int v = 0; v < this->num_vertices_aborescence; v++) {
		if (this->vertex_aborescence_active[v]) {
			printf("%d ", v);
		}		
	}
	printf("\n");
	printf("Arcs:\n");
	for (int a = 0; a < this->num_arcs_aborescence; a++) {
		if (this->arc_active[a]) {
			printf("(%d %d %4.2f)\n", this->arcs[a][0], this->arcs[a][1], this->arc_cost[a]);
		}		
	}
	return this;
}


// heuristic
_pcst* _pcst::heuristic(void* g) {	

	// instance alias
	_g* instance = (_g*)g;

	// using the graph_mst; first we update the edge weights
	int e = 0; 
	for (int u = 0; u < this->num_vertices; u++) {
		for (int v = u + 1; v < this->num_vertices; v++) {			
			if (this->vertex_prize[u] > 0 && this->vertex_prize[v] > 0) {
				(*edge_weights_mst)[edges_mst[e]] = this->all_pairs_shortest_path[u][v];
			}
			else {
				(*edge_weights_mst)[edges_mst[e]] = 1000000;
			}
			e++;
		}
	}

	// make mst all false; 
	for (int i = 0; i < this->num_vertices; i++) {
		(*mst)[edges_mst[i]] = false;
	}

	// compute mst 
	kruskal(*graph_mst, *edge_weights_mst, *mst);
	
	struct _pair_ratio {
		int u;
		int v;
		double mst; 
		double cost; 
		double ratio;
		int shortest_path_weight;
	};

	_pair_ratio* pair_ratio = new _pair_ratio[this->num_vertices - 1];

	e = 0;
	int pr = 0;
	for (int u = 0; u < this->num_vertices; u++) {
		for (int v = u + 1; v < this->num_vertices; v++) {
			if ((*mst)[edges_mst[e]]) {				
				// copy the vertices
				pair_ratio[pr].u = u;
				pair_ratio[pr].v = v;
				pair_ratio[pr].mst = (*edge_weights_mst)[edges_mst[e]];

				// compute the cost 
				double cost = this->vertex_prize[u] + this->vertex_prize[v] - (*edge_weights_mst)[edges_mst[e]];				

				pair_ratio[pr].cost = cost;			
				pair_ratio[pr].shortest_path_weight = instance->shortest_path_weight[u][v];
				pair_ratio[pr++].ratio = cost / instance->shortest_path_weight[u][v];
			}
			e++;
		}
	}

	// sort the pair ratios based on ratio
	std::sort(pair_ratio, pair_ratio + pr, [](const _pair_ratio& a, const _pair_ratio& b) {
		return a.ratio > b.ratio;
		});

	// print all pair_ratios
	printf("Pair ratios:\n");
	for (int i = 0; i < pr; i++) {
		if (pair_ratio[i].cost < 0) {
			continue;
		}
		printf("(%3d %3d) %10.2f %10.2f %5d %10.2f\n", pair_ratio[i].u, pair_ratio[i].v, pair_ratio[i].mst, pair_ratio[i].cost, pair_ratio[i].shortest_path_weight, pair_ratio[i].ratio);
	}		

	// keep track of vertices in the small tree
	int* vertices_in_small_tree = new int[this->num_vertices];
	int num_vertices_in_small_tree = 0;

	// total net reward 
	double total_net_reward = 0;

	// create a small tree
	_small_tree* small_tree = new _small_tree();

	// set the weight of tree
	small_tree->weight = 0;	

	// add the vertices to the small tree
	int u = pair_ratio[0].u;
	
	// add u to the small tree
	addbin(small_tree->bin_vertices, u);
	vertices_in_small_tree[num_vertices_in_small_tree++] = u;

	// print the small tree	
	small_tree->print_vertices(instance);

	// update the net reward 
	total_net_reward += this->vertex_prize[u];

	// find the path with the lowest weight between u and v	
	// here it is better to already compute the shortest path cost and reward between u and v
	// and save the best option in each iteration, that would help with reproducing the path	
	double bsp = budgeted_shortest_path(u, pair_ratio[0].v, instance->UB, true);

	// start from u, each time go to the edge choice of the vertex, move and compute cost
	// add the vertex to the small tree
	int v = pair_ratio[0].v;

	// get best edge out of the vertex
	int e1 = shortest_path_vertex_choice[u];

	// loop until we reach the vertex v
	while (e1 != -1) {		
		// consider the weight of the edge in the tree
		small_tree->weight += this->edge_weight[e1];

		// update net reward
		total_net_reward -= this->edge_cost[e1];

		// choose the other vertex 
		u = this->edges[e1][0] == u ? this->edges[e1][1] : this->edges[e1][0];

		// add the vertex to the tree
		addbin(small_tree->bin_vertices, u);
		vertices_in_small_tree[num_vertices_in_small_tree++] = u;

		// update the net reward
		total_net_reward += this->vertex_prize[u];		

		// print the small tree	
		small_tree->print_vertices(instance);

		// print net reward
		printf("Net reward: %4.2f\n", total_net_reward);

		// get the best edge out of the vertex
		e1 = shortest_path_vertex_choice[u];		

		if (u == v) {
			break;
		}
	}

	

	// delete the pair ratio
	delete[] pair_ratio;
	delete[] vertices_in_small_tree;

	return this;
}