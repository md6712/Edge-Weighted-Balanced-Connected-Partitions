#define _CRT_SECURE_NO_WARNINGS

#include "_g.h"
#include "_tree.h"
#include <iostream>
#include <Windows.h>
#include<direct.h>
#include "binary.h"
#include "AVLTree.h"

#include <lemon/edmonds_karp.h>

_g::_g(int num_vertices, int num_edges, int num_trees, int duplication_num)
{
	// initialize _opt and _lp_bound
	_opt = 0;
	_lp_bound = 0;
	_lb_root = 0;
	_gap = 0;


	// initialize the render object
	render = new _render();

	// allocate memory for the edges
	edges = new int[num_edges+ num_trees-1][3]; // we allow more edges in case if vertices have degree 0

	// set the number of vertices and edges
	this->num_vertices = num_vertices;
	this->num_edges = num_edges;	
	this->num_trees = num_trees;
	this->duplication_num = duplication_num;

	// allocate memory for the coordinates
	coords = new _coord[num_vertices];

	subsetsComputed = false;
	
	// allocate memory for the cycles
	cycles = new bool* [100];
	for (int k = 0; k<100; k++){
		cycles[k] = new bool[num_edges+10];
	}

	// allocate memory for the visited array
	visited = new bool[num_vertices];

	// allocate memory for the edge visited array
	edge_visited = new bool[num_edges];

	// set the filename
	filename = new char[10000];
	this->setFilename();


	// allocate memory for the start edge tree
	start_edge_tree = new int[num_trees];

	
	// allocate memory for the opt sol
	opt_edges = new int [num_edges];	


	// alloacte memory for the opt vertices
	opt_vertices = new int* [num_trees];
	num_opt_vertices = new int[num_trees];
	for (int i = 0; i < num_trees; i++) {
		opt_vertices[i] = new int[num_vertices];
	}

	// allocate memory for the parent array
	parent_mst = new int[num_vertices];

	// allocate memory for the bin vertices
	bin_vertices = new uint32_t[num_vertices];

	// allocate memory for the sum adjacent weight
	sum_adjacent_weight = new int[num_vertices];

	// allocate memory for the seperator
	seperator = new int[num_vertices];	

	// allocate memory for the colors
	v_colors = new color[num_vertices];

	// set up number of cuts 
	n_user_cuts = 0;
	n_lazy_cuts = 0;

	// trees for ub
	trees_ub = std::vector<_tree*>();

}

_g::~_g()
{
	// deallocate memory for the render object
	delete render;

	// deallocate memory for the edges	
	delete[] edges;	

	// deallocate memory for the arcs
	delete[] arcs;

	// deallocate memory for the opt edges
	delete[] opt_edges;

	// deallocate memory for the coordinates
	delete[] coords;

	// deallocate memory for the filename
	delete[] filename;

	// deallocate memory for the cycles
	for (int k = 0; k < 100; k++) {
		delete[] cycles[k];
	}
	delete[] cycles;

	// deallocate memory for the visited array
	delete[] visited;
	delete [] edge_visited;

	
	// deallocate memory for the start edge tree
	delete[] start_edge_tree;


	// deallocate memory for the subsets
	if (subsetsComputed) {
		delete[] nSubsets;

		for (int i = 0; i < pow(2, num_vertices); i++) {
			delete[] subsets[i];
		}
		delete[] subsets;

	}

	// deallocate memory for the opt vertices
	for (int i = 0; i < num_trees; i++) {
		delete[] opt_vertices[i];	
	}
	delete[] opt_vertices;
	delete[] num_opt_vertices;

	// deallocate memory for the parent array
	delete[] parent_mst;

	// deallocate memory for the bin vertices
	delete[] bin_vertices;

	// deallocate memory for the sum adjacent weight
	delete[] sum_adjacent_weight;

	// deallocate memory for the seperator
	delete[] seperator;

	// deallocate memory for the colors
	delete[] v_colors;

	// deallocate memory for the trees
	for (int i = 0; i < trees.size(); i++) {
		delete trees[i];
	}
	trees.clear();

	// for each tree in the trees_ub vector
	for (int i = 0; i < trees_ub.size(); i++) {
		delete trees_ub[i];
	}
	trees_ub.clear();

	// deallocate memory for the digraph
	dg.clear();

	// delete the select trees for CG if size is greater than 0
	for (int i = 0; i < select_trees_for_CG.size(); i++) {
		delete select_trees_for_CG[i];
	}

	// deallocate memory for the select trees for CG	
	select_trees_for_CG.clear();

	// delete mwcs 
	if (mwcs != nullptr) {
		delete mwcs;
	}

	// delete pcst
	if (pcst != nullptr) {
		delete pcst;
	}
}

void _g::readGraph()
{
	
	// Open a file and write the string to the file
	FILE* f = fopen(filename, "r");
	if (f == NULL) {
		printf("Error opening file!\n");
	}
	else {
		// read the first line of the file with scanf
		fscanf(f, "%d %d %d", &num_vertices, &num_edges, &num_trees);
		
		// read edges
		for (int i = 0; i < num_edges; i++) {
			fscanf(f, "%d %d %d", &edges[i][0], &edges[i][1], &edges[i][2]);

			//without loss of generality, we assume that $u < v$ for each edge $(u,v) \in E$.
			if (edges[i][0] > edges[i][1]) {
				int temp = edges[i][0];
				edges[i][0] = edges[i][1];
				edges[i][1] = temp;
			}
		}
	}


	// if there is a vertex with degree 0, add an edge to the last vertex	
	for (int i = 0; i < num_vertices; i++) {
		int degree = 0;
		for (int e = 0; e < num_edges; e++) {
			if (edges[e][0] == i || edges[e][1] == i) {
				degree++;
			}
		}
		if (degree == 0) {
			edges[num_edges][0] = i;
			edges[num_edges][1] = num_vertices - 1;
			edges[num_edges][2] = LARGE_WEIGHT;
			num_edges++;
		}
	}


	// compute the number of arcs
	this->num_arcs = num_edges * 2 + num_vertices;

	// set arcs
	arcs = new int[num_arcs][3];

	// add arcs

	for (int i = 0; i < num_edges; i++) {
		// add the arcs
		arcs[i][0] = edges[i][0];
		arcs[i][1] = edges[i][1];
		arcs[i][2] = edges[i][2];

		arcs[i + num_edges][0] = edges[i][1];
		arcs[i + num_edges][1] = edges[i][0];
		arcs[i + num_edges][2] = edges[i][2];
	}

	// add arcs from s to all vertices in the first tree
	for (int i = 0; i < num_vertices; i++) {
		arcs[i + 2 * num_edges][0] = num_vertices;
		arcs[i + 2 * num_edges][1] = i;
		arcs[i + 2 * num_edges][2] = 0;
	}

	// close the file
	fclose(f);
	
}

void _g::createMWCS() {
	this->mwcs = new _mwcs(this);
}

void _g::createPCST() {
	this->pcst = new _pcst(this);
}


void _g::printGraph()
{
	printf("Number of vertices: %d\n", num_vertices);
	printf("Number of edges: %d\n", num_edges);
	printf("Number of trees: %d\n", num_trees);
	printf("Edges:\n");
	for (int e = 0; e < num_edges; e++) {
		printf("(%3d) [%3d %3d] : w %3d\n", e, edges[e][0], edges[e][1], edges[e][2]);
	}	

	/*printf("Arcs:\n");
	for (int a = 0; a < num_arcs; a++) {
		printf("[%2d] [%2d -> %2d]--%3d\n", a, arcs[a][0], arcs[a][1], arcs[a][2]);
	}*/
}

void _g::setFilename()
{
	sprintf_s(this->filename, sizeof(char)*200, "instances//test_%d_%d_%d_%d.txt", num_vertices, num_edges, num_trees, duplication_num);
}

char* _g::getFilename()
{
	return filename;
}

void _g::computeSubsets()
{


	// allocate memory for the subsets
	subsets = new uint32_t * [pow(2, num_vertices)];
	for (int i = 0; i < pow(2, num_vertices); i++) {
		subsets[i] = new uint32_t[binaryArrlength(num_vertices)];
		memset(subsets[i], 0, sizeof(uint32_t) * binaryArrlength(num_vertices));
	}

	// allocate memory for the number of subsets
	nSubsets = new int[pow(2, num_vertices)];

	// enumerate all subsets of a set from 1 to num_vertices
	for (int i = 1; i < pow(2, num_vertices); i++) {
		// convert i to binary
		int* binary = new int[num_vertices];
		int k = i;
		for (int j = 0; j < num_vertices; j++) {
			binary[j] = k % 2;
			k = k / 2;
		}

		// add the subset to the list of subsets
		int count = 0;
		for (int j = 0; j < num_vertices; j++) {
			if (binary[j] == 1) {
				addbin(subsets[i], j);
				count++;
			}
		}

		nSubsets[i] = count;
	}	

	subsetsComputed = true;
}

void _g::printSubsets()
{
	printf("\nSubsets\n");
	for (int i = 0; i < pow(2, num_vertices); i++) {
		// print
		printf("%d : %d \n", subsets[i][0], nSubsets[i]);		
	}
}

void _g::outputOPTEdges() {
	FILE* f = fopen("opt_edges.txt", "w");
	if (f == NULL) {
		printf("Error opening file!\n");
	}
	else {
		for (int i = 0; i < num_trees; i++) {
			for (int e = 0; e < num_vertices - num_trees; e++) {
				fprintf(f, "%d %d\n", edges[opt_edges[e]][0], edges[opt_edges[e]][1]);
			}
		}
	}
	fclose(f);	
}

// check if there are cycles in opt_edges
bool _g::CheckCyclesInOptEdges() {

	// initialize cycle to false
	bool cycle = false;

	// initialize cycles
	for (int i = 0; i < num_trees; i++) {
		memset(cycles[i], 0, sizeof(bool) * num_edges);		
	}

	// initialize visited array
	memset(visited, 0, sizeof(bool) * num_vertices);	

	// initialize edge visited array
	memset(edge_visited, 0, sizeof(bool) * num_edges);

	// init number of cycles
	num_cycles = 0;	

	// loop over the trees
	for (int i = 0; i < num_trees; i++) {
		// check if there are no edges in the tree
		if (i < num_trees - 1) {
			if (start_edge_tree[i] == start_edge_tree[i + 1]) {
				continue;
			}
		}
		
		// initialize in_cycle to false;  we are not yet in a cycle
		in_cycle = false;

		int ep = start_edge_tree[i];
	
		while (true) {
			// first edge of each tree			 
			while (edge_visited[ep]){
				ep++;
				if ((i < num_trees - 1 && ep == start_edge_tree[i+1]) || ep >= num_opt_edges - 1) {					
					break;
				}				
			}		

			// leave the upper loop if we are at the end of the edges
			if ((i < num_trees - 1 && ep == start_edge_tree[i + 1]) || ep >= num_opt_edges - 1) {
				break;
			}

			// start from an edge in the tree and check recursively for cycles		
			if (CheckCyclesInOptEdgesRec(i, ep, -1)) {
				num_cycles++; // increment the number of cycles
				cycle = true; // we have found a cycle
				in_cycle = false; // we are no longer in a cycle // setting this to false is important as we want to find a new cycle
				break;
			}
		}

		// check if there is singleton tree
		if (ep >= num_opt_edges - 1) {
			break;  // there is a singleton tree
		}

	}	

	// return true if a cycle is found and false otherwise
	return cycle;
}

// recursive function to check for cycles in tree i 
bool _g::CheckCyclesInOptEdgesRec(int i, int ep, int v) {
	// compute e from ep
	int e = opt_edges[ep];

	// mark edge e as visited
	edge_visited[ep] = true;

	// compute the adjacent vertex
	int v2 = (edges[e][0] == v) ? edges[e][1] : edges[e][0];

	// return true if edge e is already visited
	if (visited[v2]) {
		// add e to visited edges
		cycles[num_cycles][e] = true;
		// define the start edge of the cycle
		cycle_start_vertex = v2;
		// state that we are in a cycle
		in_cycle = true;
		// return true
		return true;
	}

	// initialize cycle
	bool cycle = false;

	// mark vertices v and v2 as visited	
	if (v != -1) visited[v] = true;
	if (v2 != -1) visited[v2] = true;

	for (int ep2 = 0; ep2 < num_opt_edges; ep2++) {
		if (ep2 != ep) {
			// compute e2 from ep2
			int e2 = opt_edges[ep2];

			// check if e and e2 are adjacent
			if (edges[e][0] == edges[e2][0] || edges[e][0] == edges[e2][1] || edges[e][1] == edges[e2][0] || edges[e][1] == edges[e2][1]) {				
  
				// compute the adjacent vertex
				int nv = 0 ;
				if (edges[e][0] == edges[e2][0] || edges[e][0] == edges[e2][1]) nv = edges[e][0];
				if (edges[e][1] == edges[e2][0] || edges[e][1] == edges[e2][1]) nv = edges[e][1];
				
				// check if the vertex is not the previous vertex
				if (nv != v) {
					cycle = CheckCyclesInOptEdgesRec(i, ep2, nv);

					// if a cycle is found,  return true
					if (cycle) {
						// if we are in cycle, store the edge
						if (in_cycle) {
							// add e to visited edges
							cycles[num_cycles][e] = true;
							// if e is the start edge of the cycle, we are no longer in a cycle
							if (edges[e][0] == cycle_start_vertex || edges[e][1] == cycle_start_vertex) {
								in_cycle = false;
							}
						}
						return true;
					}
				}
			}										
		}
	}
	
	// return false if no cycle is found
	return false;
}

// print the cycles
void _g::PrintCycles() {
	for (int i = 0; i < num_cycles; i++) {
		printf("Cycle %d: ", i);
		for (int e = 0; e < num_edges; e++) {
			if (cycles[i][e]) {
				printf("(%d %d)", edges[e][0], edges[e][1]);
			}
		}		
		printf("\n");
	}
	
}

// print the opt edges
void _g::PrintOptEdges(bool printWeight) {
	// choose the first tree
	int i = 0;
	
	// skip trees with no edges
	while (i < num_trees - 1 && start_edge_tree[i] == start_edge_tree[i + 1]) { i++;}
	
	// initialize printed
	bool *printed = new bool [num_edges];
	int w_edges = 0;
	for (int ep = 0; ep < num_opt_edges; ep++) {
		
		if (ep == start_edge_tree[i]) {
			if (i > 0) {
				// print w_edges
				if (printWeight)
				printf("\t Weight: %d", w_edges);					
			}
			printf("\nTree %d: ", i);
			i++;

			// skip trees with no edges
			while (i < num_trees - 1 && start_edge_tree[i] == start_edge_tree[i + 1]) { i++; }

			w_edges = 0;
		}
		int e = opt_edges[ep];
		w_edges += edges[e][2];
		printf("(%d %d)", edges[e][0], edges[e][1]);
		printed[e] = true;		
	}
	if (printWeight)
	printf("\t Weight: %d", w_edges);
	printf("\n");	

	delete [] printed;
	// print the remaining edges
	//// print the remaining edges
	//printf("\n Tree -1: ");
	//for (int e = 0; e < num_edges; e++) {
	//	if (!printed[e]) {
	//		printf("(%d %d)", edges[e][0], edges[e][1]);
	//	}
	//}		
}

_g* _g::SortEdges() {
	
	// sort the edges using qsort
	qsort(edges, num_edges, sizeof(edges[0]), [](const void* a, const void* b) {
		int* x = (int*)a;
		int* y = (int*)b;
		return x[2] - y[2];
	});

	// recompute arcs
	for (int i = 0; i < num_edges; i++) {
		arcs[i][0] = edges[i][0];
		arcs[i][1] = edges[i][1];
		arcs[i][2] = edges[i][2];
		arcs[i + num_edges][0] = edges[i][1];
		arcs[i + num_edges][1] = edges[i][0];
		arcs[i + num_edges][2] = edges[i][2];
	}
	
	return this;
}

// compute the minimum spanning tree using Kruskal's algorithm
int _g::MST(int tree, int* mst_vertices, int n) {
	// set the start_edge_tree 
	start_edge_tree[tree] = num_opt_edges;

	if (n <= 1) {
		return 0;
	}	
	// initialize the total weight of the MST
	int total_weight = 0;

	// initialize number of edges in the MST
	int num_edges_mst = 0;

	// initialize the parent array	
	for (int i = 0; i < num_vertices; i++) {
		parent_mst[i] = i;
	}

	// initiate bin vertices
	memset(bin_vertices, 0, sizeof(uint32_t) * binaryArrlength(num_vertices));
	
	// add the vertices to the bin vertices
	for (int i = 0; i < n; i++) {
		addbin(bin_vertices, mst_vertices[i]);
	}	

	// memset the sum adjecent weight for vertices 
	memset(sum_adjacent_weight, 0, sizeof(int) * num_vertices);
	
	
	//printf("\nTree: ");

	// loop over the edges
	for (int e = 0; e < num_edges; e++) {
		// check if the edge is in the selected vertices
		if (checkbin(bin_vertices, edges[e][0]) && checkbin(bin_vertices, edges[e][1])) {

			// get the vertices
			int u = edges[e][0];
			int v = edges[e][1];

			// sum the weight of the adjacent edges
			// this is used as a coeficient for logic based cuts
			sum_adjacent_weight[u] += edges[e][2];
			sum_adjacent_weight[v] += edges[e][2];

			// get the weight
			int w = edges[e][2];

			// find the parent of u and v
			int pu = u;
			while (parent_mst[pu] != pu) {
				pu = parent_mst[pu];
			}

			int pv = v;
			while (parent_mst[pv] != pv) {
				pv = parent_mst[pv];
			}

			// if the parent of u and v are different, add the edge to the MST
			if (pu != pv) {
				//print edge
				//printf("(%d %d)", u, v);
				
				// add the edge to the optimal edges array for tree 
				opt_edges[num_opt_edges++] = e;

				total_weight += w;
				num_edges_mst++; 
				parent_mst[pu] = pv;
			}
		}
	}

	if (num_edges_mst != n - 1) {
		// return -1 if the number of edges in the MST is not equal to the number of vertices - 1 
		// not a tree
		return -1;
	}
	else {
		// return the total weight of the MST
		return total_weight;
	}
}

// compute the u-v separator for vertices u and v
bool _g::UVSeperator(int u, int v) {
	
	// print u and v	
	for (int e = 0; e < num_edges; e++) {
		if ((edges[e][0] == u && edges[e][1] == v) || (edges[e][1] == u && edges[e][0] == v)) {			
			return false;
		}
	}

	// initialize num_seperators
	num_seperators = 0;

	// initialize the seperator
	memset(seperator, 0, sizeof(int) * num_vertices);

	// initialize the colors
	memset(v_colors, white, sizeof(color) * num_vertices);

	// travers neighbors of u and color all neighbors blue
	ColorNeighbors(u, blue, false);	

	// color all reachable vertices from v red; include also transitive neighbors	
	ColorNeighbors(v, red, true);

	//printf("\nV_Colors: ");
	//for (int i = 0; i < num_vertices; i++) {
	//	printf("%d ", v_colors[i]);
	//}

	// check if seperated 
	memset(visited, false, sizeof(bool) * num_vertices);
	bool seperated = Seperated(u, v);

	for (int i = 0; i < num_vertices; i++) {
		if (v_colors[i] == yellow) {
			v_colors[i] = white;
			memset(visited, false, sizeof(bool) * num_vertices);
			if (!Seperated(u, v)) {
				v_colors[i] = yellow;
			}
		}
	}

	return true;
}

// color neigbors of v color c 
void _g::ColorNeighbors(int v, color c, bool transitive) {
	for (int e = 0; e < num_edges; e++) {
		int v2 = -1;

		// select the adjecent vertex
		if (edges[e][0] == v) v2 = edges[e][1];
		if (edges[e][1] == v) v2 = edges[e][0];

		if (v2 >= 0)
		{
			// check if the vertex is not colored
			if (v_colors[v2] == white) {
				v_colors[v2] = c;
				if (transitive)
					ColorNeighbors(v2, c, transitive);
			}
		// if colored, use yellow
			else {
				v_colors[v2] = yellow;
			}
		}
	}
}

// check u-v seperated with S (yellow vertices)
bool _g::Seperated(int u, int v) {
	// mark u and v as visited
	visited[u] = visited[v] = true;

	// check if u and v are directly connected
	for (int e = 0; e < num_edges; e++) {
		if (edges[e][0] == u && edges[e][1] == v) {
			return false;
		}
		if (edges[e][0] == v && edges[e][1] == u) {
			return false;
		}
	}
	
	// check if neigbors of u are seperated from v
	bool seperated = true;
	for (int e = 0; e < num_edges; e++) {
		int w = -1;

		// select the adjecent vertex
		if (edges[e][0] == u) w = edges[e][1];
		if (edges[e][1] == u) w = edges[e][0];

		if (w >= 0) {
			if (v_colors[w] == yellow) {
				continue;
			}
			else {
				if (!visited[w]) {
					seperated = Seperated(w, v);
				}
			}

			if (!seperated) {
				break;
			}
		}
	}

	return seperated;
}

// print all pairs min seperator
void _g::PrintMinSeperators() {
	for (int u = 0; u < num_vertices; u++) {
		for (int v = u + 1; v < num_vertices; v++) {
			if (UVSeperator(u, v)) {
				printf("\nSeperator (%d %d): ", u, v);
				for (int i = 0; i < num_vertices; i++) {
					if (v_colors[i] == yellow) {
						printf("%d ", i);
					}
				}				
			}
		}
	}
}

void _g::setDiGraph() {

	// for each node, create a node 
	for (int i = 0; i < num_vertices; i++) {		
		dg.addNode();
	}	

	dg.addNode(); // add node s;  // this is node num_vertices
	dg.addNode(); // add node t;  // this is node num_vertices + 1

	// add arcs
	for (int e = 0; e < num_edges; e++) {
		dg.addArc(dg.nodeFromId(edges[e][0]), dg.nodeFromId(edges[e][1]));
	}

	for (int e = 0; e < num_edges; e++) {
		dg.addArc(dg.nodeFromId(edges[e][1]), dg.nodeFromId(edges[e][0]));
	}

	// add arcs from s to all vertices in the first tree
	for (int i = 0; i < num_vertices; i++) {
		dg.addArc(dg.nodeFromId(num_vertices), dg.nodeFromId(i));
	}

	// add arcs from all vertices in the last tree to t
	for (int i = 0; i < num_vertices; i++) {
		dg.addArc(dg.nodeFromId(i), dg.nodeFromId(num_vertices + 1));		
	}
	
	//// print the graph
	//for (ListDigraph::ArcIt a(dg); a != INVALID; ++a) {
	//	printf("Arc %d -> %d\n", dg.id(dg.source(a)), dg.id(dg.target(a)));
	//}
}


// generate trees
// the method only works for num_vertices < 32
void _g::generateTrees() {

	// check if num_vertices is less than 32
	if (num_vertices > 32) {
		printf("Number of vertices is greater than 32\n");
		throw ("Number of vertices is greater than 32\n");
		return;
	}


	// allocate memory for the vertices
	int* vertices = new int[num_vertices];

	// binary number for the vertices
	uint32_t* bin_vertices = new uint32_t[binaryArrlength(this->num_vertices)];
	memset(bin_vertices, 0, sizeof(uint32_t) * binaryArrlength(this->num_vertices));
	
	// compute the maximum number of subsets
	long maxBinNum = pow(2, num_vertices) - 1;
	
	// tree counter
	int treeCounter = 0;


	// print trees 
	//std::cout << "\n Trees: \n";

	// loop over all subsets of vertices and generate vertices vector
	for (uint32_t i = 1; i <= maxBinNum; i++) {
		memset(bin_vertices, 0, sizeof(uint32_t) * binaryArrlength(this->num_vertices));
		int count = 0;

		// check what vertices are in the subset
		for (int j = 0; j < num_vertices; j++) {
			if (checkbinSingle(i, j)) {
				vertices[count++] = j;
				addbin(bin_vertices, j);
			}
		}

		// create a minimum spanning tree 
		_tree* tree = new _tree( this);
		tree->CopyVertices(count, vertices, bin_vertices);
		if (tree->num_vertices >= 2) {
			tree->ComputeMST();
		}
		else
		{
			tree->SingltonTree(vertices[0]);
		}

		if (tree->isSpanningTree 
				&& tree->weight <= UB
			) {
			trees.push_back(tree);			
			tree->RecomputeDegree();
			treeCounter++;
			/*std::cout << treeCounter++ << ": ";
			tree->PrintTree();
			*/
			if (treeCounter % 10000 == 0) {
				std::cout << treeCounter << std::endl;
			}
		}
		else {
			delete tree;
		}
	}

	// delete local arrays
	delete[] vertices;
	delete[] bin_vertices;

	// print the number of trees
	std::cout << "Number of trees: " << trees.size() << std::endl;	
}

// generate select trees
void _g::generateSelectTrees() {

	// check if num_vertices is less than 32
	/*if (num_vertices > 32) {
		printf("Number of vertices is greater than 32\n");
		throw ("Number of vertices is greater than 32\n");
		return;
	}*/

	// allocate memory for the vertices
	int* vertices = new int[num_vertices];

	// binary number for the vertices
	uint32_t* bin_vertices = new uint32_t[binaryArrlength(this->num_vertices)];

	// compute the maximum number of subsets
	long maxBinNum = pow(2, num_vertices) - 1;

	// tree counter
	int treeCounter = 0;


	//create all singlton trees 
	for (int i = 0; i < num_vertices; i++) {
		memset(bin_vertices, 0, sizeof(uint32_t) * binaryArrlength(this->num_vertices));
		addbin(bin_vertices, i);
		// &i treat the single number as an array of size 1
		// 1 << i is the binary number for the vertex i		

		_tree* tree = new _tree(this);
		tree->CopyVertices(1, &i, bin_vertices);
		//tree->singltonTree(i);
		trees.push_back(tree);
		treeCounter++;
	}


	// reset bin_vertices
	memset(bin_vertices, 0, sizeof(uint32_t) * binaryArrlength(this->num_vertices));

	// create a tree including all vertices
	for (int v = 0; v < num_vertices; v++) {
		vertices[v] = v;
		addbin(bin_vertices, v);
	}
	_tree* tree = new _tree(this);
	tree->CopyVertices(num_vertices, vertices, bin_vertices);

	// set a bool vector of forbidden edges
	bool* forbidden_edges = new bool[num_edges];
	memset(forbidden_edges, 0, sizeof(bool) * num_edges);

	AVLTree* avl = new AVLTree();
	
	while (true) {
		tree->ComputeMST(forbidden_edges);

		_tree** trees_to_be_added = NULL;

		bool new_tree_added = false; // check if we add a tree in this iteration

		if (tree->num_forbidden_edges_in_mst <= num_trees - 1) {
			trees_to_be_added = tree->SplitIntoKTrees(this->num_trees, forbidden_edges );
			for (int i = 0; i < num_trees; i++) {
				if (trees_to_be_added[i]->weight > UB) continue;
				if (trees_to_be_added[i]->weight < 1) continue;
				if (avl->search(trees_to_be_added[i]->bin_vertices[0]) == nullptr) {					
					avl->insert(trees_to_be_added[i]->bin_vertices[0], trees_to_be_added[i]);
					trees.push_back(trees_to_be_added[i]);
					new_tree_added = true;

					// print tree 

					trees_to_be_added[i]->PrintTree(); 

					// update counter
					treeCounter += 1;
				}
				else {
					delete trees_to_be_added[i];
				}
			}

			// delete local arrays
			delete[] trees_to_be_added;

			if (!new_tree_added) {
				break;
			}
		}
		else {
			break;
		}

		// take the first edge of the tree that is not forbidden already			
		int k =0; 
		int e = tree->edges[k++];
		while (forbidden_edges[e]) {
			e = tree->edges[k++];
		}

		// print selected edge
		//printf("Selected edge: [%2d %2d] -> %2d\n", edges[e][0], edges[e][1], edges[e][2]);

		forbidden_edges[e] = true;
	}
	


	delete tree;
	delete avl;

	// delete local arrays
	delete[] vertices;
	delete[] bin_vertices;
	delete[] forbidden_edges;

	// print the number of trees
	std::cout << "Number of trees: " << trees.size() << std::endl;
}

// recompute LB 
void _g::recomputeLB() {
	// set the lower bound to 0
	LB = 0;

	LB = MST_weight - (num_trees - 1) * UB;

	for (int i = 0; i < num_trees-1; i++) {
		LB -= edges[num_edges - i - 1][2];
	}	

	if (LB < 0) {
		LB = 0;
	}
}

// compute UB 
void _g::computeUB() {
	
	// allocate memory for the vertices
	int* vertices = new int[num_vertices];

	// binary number for the vertices
	uint32_t* bin_vertices = new uint32_t[binaryArrlength(this->num_vertices)];

	for (int v = 0; v < num_vertices; v++) {
		vertices[v] = v;
		addbin(bin_vertices, v);
	}
	_tree* tree = new _tree(this);
	tree->CopyVertices(num_vertices, vertices, bin_vertices);

	// set a bool vector of forbidden edges
	bool* forbidden_edges = new bool[num_edges];
	memset(forbidden_edges, 0, sizeof(bool) * num_edges);
	
	// compute MST 
	tree->ComputeMST(forbidden_edges);

	// print the MST
	//tree->PrintTree();

	// recompute degrees
	tree->RecomputeDegree();

	// save the weight of the MST
	MST_weight = tree->weight;

	// split the trees into k trees	
	_tree** trees = tree->SplitIntoKTrees(this->num_trees, forbidden_edges);;

	// set the upper bound
	UB = 0;
	
	for (int i = 0; i < num_trees; i++) {
		//trees[i]->PrintTree();

		if (trees[i]->weight > UB) 
			UB = trees[i]->weight;
	}

	UB_naive = UB;

	 //print UB
	std::cout << "UB: " << UB << std::endl;

	// recompute LB
	recomputeLB();

	// print LB
	//std::cout << "LB: " << LB << std::endl;

	// delete local arrays
	delete[] vertices;
	delete[] bin_vertices;
	delete[] forbidden_edges;
	delete tree;

	// delete trees
	for (int i = 0; i < num_trees; i++) {
		delete trees[i];
	}
	delete[] trees;	
}


// forced directed layout 
void  _g::ForcedDirectedLayout() {

	// width and height of the window
	int width = 1200;
	int height = 800;

	// srand for random number generation
	//srand(time(NULL));
	srand(0);

	// give random coordinations to the vertices
	for (int i = 0; i < num_vertices; i++) {
		coords[i].x = (rand() % width);
		coords[i].y = (rand() % height);
	}

	// simulation constants 
	double k = sqrt(width * height) / num_vertices; // optimal distance between vertices
	double cooling_factor = 0.99; // cooling factor

	//// print the coordinates
	//for (int i = 0; i < num_vertices; i++) {
	//	printf("Vertex %d: (%f, %f)\n", i, coords[i].x, coords[i].y);
	//}

	// for each iteration
	for (int itr = 0; itr < 1000; itr++) {
		
		// reset dx and dy
		for (int i = 0; i < num_vertices; i++) {
			coords[i].dx = 0;
			coords[i].dy = 0;
		}

		// compute the repulsive forces
		for (int i = 0; i < num_vertices; i++) {
			for (int j = i + 1; j < num_vertices; j++) {
				double dx = coords[i].x - coords[j].x;
				double dy = coords[i].y - coords[j].y;
				double d = sqrt(dx * dx + dy * dy);
				if (d == 0) d = 0.0001;
				double force = k * k / d;
				coords[i].dx += force * dx / d;
				coords[i].dy += force * dy / d;
				coords[j].dx -= force * dx / d;
				coords[j].dy -= force * dy / d;
			}
		}

		// compute the attractive forces
		for (int e = 0; e < num_edges; e++) {
			int i = edges[e][0];
			int j = edges[e][1];
			double dx = coords[i].x - coords[j].x;
			double dy = coords[i].y - coords[j].y;
			double d = sqrt(dx * dx + dy * dy);
			if (d == 0) d = 0.0001;
			double force = d * d / k;
			coords[i].dx -= force * dx / d;
			coords[i].dy -= force * dy / d;
			coords[j].dx += force * dx / d;
			coords[j].dy += force * dy / d;
		}


		// update the coordinates
		for (int i = 0; i < num_vertices; i++) {
			coords[i].x += coords[i].dx * cooling_factor;
			coords[i].y += coords[i].dy * cooling_factor;			
		}

		// introduce min_x, min_y, max_x, max_y
		double min_x = 1200, min_y = 0, max_x = width, max_y = height;

		// map the coordinates to the window
		for (int i = 0; i < num_vertices; i++) {
			// compute lowest x and y coordinates
			min_x = min(coords[i].x, min_x);
			min_y = min(coords[i].y, min_y);

			// compute highest x and y coordinates
			max_x = max(coords[i].x, max_x);
			max_y = max(coords[i].y, max_y);			
		}

		// current width and height
		double cwidth = max_x - min_x;
		double cheight = max_y - min_y;

		// translate each coordinate to a range of 0 to width and 0 to height
		for (int j = 0; j < num_vertices; j++) {
			coords[j].x = (coords[j].x - min_x) / cwidth * width;
			coords[j].y = (coords[j].y - min_y) / cheight * height;
		}

		// reduce the temperature
		cooling_factor *= 0.99;				
	}

	// print the coordinates
	/*for (int i = 0; i < num_vertices; i++) {
		printf("Vertex %d: (%f, %f)\n", i, coords[i].x, coords[i].y);
	}*/

	// compute minimum x and y coordinates
	int min_x = width, min_y = height;
	int max_x = 0, max_y = 0;	
	for (int v = 0; v < num_vertices; v++) {
		min_x = (coords[v].x < min_x) ? coords[v].x : min_x;
		min_y = (coords[v].y < min_y) ? coords[v].y : min_y;
		max_x = (coords[v].x > max_x) ? coords[v].x : max_x;
		max_y = (coords[v].y > max_y) ? coords[v].y : max_y;
	}

	// center the coordinates
	int diff_x = width / 2 - (max_x + min_x) / 2; 
	int diff_y = height / 2 - (max_y + min_y) / 2;

	

	// shift the coordinates
	for (int v = 0; v < num_vertices; v++) {
		coords[v].x += diff_x;
		coords[v].y += diff_y;
	}
		
	// print the coordinates
	/*for (int i = 0; i < num_vertices; i++) {
		printf("Vertex %d: (%f, %f)\n", i, coords[i].x, coords[i].y);
	}*/

}

// draw the graph
void _g::DrawGraph(int* highlighted_edges, int num_highlighted_edges) {

	for (int i = 0; i < num_edges; i++) {
		int u = edges[i][0];
		int v = edges[i][1];
		double x1 = coords[u].x;
		double y1 = coords[u].y;
		double x2 = coords[v].x;
		double y2 = coords[v].y;
		render->addEdge(x1, y1, x2, y2, sfLinesStrip);
	}

	if (num_highlighted_edges > 0) {
		for (int i = 0; i < num_highlighted_edges; i++) {
			int u = edges[highlighted_edges[i]][0];
			int v = edges[highlighted_edges[i]][1];
			double x1 = coords[u].x;
			double y1 = coords[u].y;
			double x2 = coords[v].x;
			double y2 = coords[v].y;
			render->highlightEdge(x1, y1, x2, y2);
		}
		
	}

	for (int i = 0; i < num_vertices; i++) {
		double x = coords[i].x;
		double y = coords[i].y;
		render->addNode(i, x, y);
	}

	render->draw();

}



void _g::populate_trees_ub_from_select_trees() {

	// search all trees in the selected trees, and if tree is part of optimal, 
	// add it to the trees_ub

	for (int i = 0; i < select_trees_for_CG.size(); i++) {
		_small_tree* st = select_trees_for_CG[i];
		

		// check if the tree is in the best ub trees
		if (st->part_of_optimal == num_upper_bound_updates) {						

			_tree* tree = new _tree(this);

			tree->CopyVertices(st->bin_vertices);	

			if (tree->num_vertices >= 2) {
				tree->ComputeMST();
			}
			else
			{
				tree->SingltonTree(tree->vertices[0]);
			}

			trees_ub.push_back(tree);				
		}		
	}

	// sort trees based on their weight
	std::sort(trees_ub.begin(), trees_ub.end(), [](const _tree* a, const _tree* b) {
		return a->weight > b->weight;
		});


	//// print trees	
	//for (int i = 0; i < trees_ub.size(); i++) {		
	//	trees_ub[i]->PrintTree();	
	//}

}