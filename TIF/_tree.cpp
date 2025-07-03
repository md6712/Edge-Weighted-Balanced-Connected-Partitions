#include "_tree.h"
#include <lemon/kruskal.h>
#include <lemon/list_graph.h>
#include <iostream>
#include <algorithm>
#include "binary.h"

#include "_g.h"

void _small_tree::print_vertices(void *g)
{
	_g* graph = (_g*)g;
	printf("Tree: ");
	printf("\t Weight: %d \t Opt: %d \t", weight, part_of_optimal);
	for (int v = 0; v < graph->num_vertices; v++)
	{
		if (checkbin(bin_vertices, v))
		{
			printf("%d ", v);
		}		
	}
	// print weight
	printf("\n");	
}

// constructor for the tree
_small_tree::_small_tree() {
	part_of_optimal = 0;
	memset(bin_vertices, 0, sizeof(uint32_t) * SIZE_OF_VERTICES_BINARY);
}

// destructor for the tree
_small_tree::~_small_tree() {
	// nothing to do here, since we don't allocate any memory in the constructor
}

// compute MST 
void _small_tree::computeMST(void* g)
{
	// convert to a _tree
	_tree* tree = new _tree(g);

	// copy the vertices from the binary representation
	tree->CopyVertices(this->bin_vertices);	

	// compute the minimum spanning tree		
	tree->ComputeMST();		

	// if weight is zero, then the tree is empty
	if (tree->isSpanningTree) {
		if (tree->weight == 0) {
			this->weight = MAXINT;
		}
		else {
			// set the weight of the tree
			this->weight = tree->weight;
		}
	}
	else {
		this->weight = MAXINT; // if the tree is not spanning, set the weight to MAXINT
	}

	// delete tree
	delete tree;
}

_tree::_tree(void *instance)
{
	// get num_vertices in the instance graph
	int num_vertices = ((_g*)instance)->num_vertices;

	// tree is empty 
	this->num_vertices = 0;
	this->weight = 0;

	// the tree is not a spanning tree that means not all vertices in the graph is covered. an empty tree is also not spanning tree. 
	this->isSpanningTree = false;

	this->edges = new int[num_vertices - 1];
	memset(this->edges, -1, sizeof(int) * (num_vertices - 1)); // initialize the edges array

	this->vertices = new int[num_vertices];
	this->decomp = new int[num_vertices];
	this->degree = new int[num_vertices];	
	this->bin_vertices = new uint32_t[binaryArrlength(num_vertices)];
	memset(this->bin_vertices, 0, sizeof(uint32_t)*binaryArrlength(num_vertices));

	// initialize the visited and stack arrays	
	this->visited = new bool[num_vertices];
	this->stack = new int[num_vertices];
		
	this->instance = instance;
	
	// this is used for mst
	this->num_forbidden_edges_in_mst = 0;
}

//memcpy(this->vertices, vertices, sizeof(int)* num_vertices);
//this->bin_vertices = new uint32_t[binaryArrlength(((_g*)instance)->num_vertices)];
//memcpy(this->bin_vertices, bin_vertices, sizeof(uint32_t)* binaryArrlength(((_g*)instance)->num_vertices));

_tree::~_tree()
{
	delete[] this->edges;
	delete[] this->vertices;
	delete[] this->decomp;
	delete[] this->bin_vertices;
	delete[] this->degree;

	// delete visited and stack arrays
	delete[] this->visited;
	delete[] this->stack;
}

// singleton tree
void _tree::SingltonTree(int vertex)
{
	this->num_vertices = 1;
	this->vertices[0] = vertex;
	this->weight = 0;
	this->isSpanningTree = true;
}

// copy vertices
void _tree::CopyVertices(uint32_t* bin_vertices)
{	
	memcpy(this->bin_vertices, bin_vertices, sizeof(uint32_t)*binaryArrlength(((_g*)instance)->num_vertices));

	// for all vertices in the graph
	for (int i = 0; i < ((_g*)instance)->num_vertices; i++)
	{
		// check if the vertex is in the tree
		if (checkbin(bin_vertices, i))
		{
			this->vertices[this->num_vertices++] = i;
		}
	}
}
void _tree::CopyVertices(int num_vertices, int* vertices, uint32_t* bin_vertices)
{
	memcpy(this->bin_vertices, bin_vertices, sizeof(uint32_t)*binaryArrlength(((_g*)instance)->num_vertices));
	memcpy(this->vertices, vertices, sizeof(int)*num_vertices);
	this->num_vertices = num_vertices;
}
void _tree::CopyVertices(int num_vertices, int* vertices)
{
	// for all vertices in the graph
	for (int i = 0; i < ((_g*)instance)->num_vertices; i++)
	{
		// check if the vertex is in the tree
		if (checkbin(bin_vertices, i))
		{
			this->vertices[this->num_vertices++] = i;
		}
	}
}

// build graph by edges
void _tree::InitTreeByEdges(int num_edges, int* edges) {
		// get the instance graph
	_g* g = (_g*)this->instance;

	// for all edges, add them to tree
	for (int e = 0; e < num_edges; e++) {
		AddEdge(edges[e]);
	}
}

// compute the minimum spanning tree
void _tree::ComputeMST(bool *forbidden_edges)
{
	// set the number of forbidden sets in MST
	num_forbidden_edges_in_mst = 0;

	// get the original graph
	_g* g = (_g*)this->instance;

	// create a listGraph
	ListGraph graph;	

	// node one to one 
	int *nodeOneToOneMap = new int[g->num_vertices];

	// edge one to one
	int *edgeOneToOneMap = new int[g->num_edges];
		
	// add the vertices
	ListGraph::Node* nodes = new ListGraph::Node[this->num_vertices];
	for (int i = 0; i < this->num_vertices; i++)
	{
		nodes[i] = graph.addNode();		
		nodeOneToOneMap[this->vertices[i]] = i;
	}

	//ListGraph::NodeMap<int> nodeMap(graph);
	
	// add the edges
	// for all edges in g, if the edge vertices is one of the vertices in this tree, add the edge to the graph
	ListGraph::Edge* edges = new ListGraph::Edge[g->num_edges];
	int edgeCount = 0;
	for (int i = 0; i < g->num_edges; i++)
	{			

		// get the vertices of the edge
		int u = g->edges[i][0];
		int v = g->edges[i][1];

		// get the node id
		int nodeU = nodeOneToOneMap[u];
		int nodeV = nodeOneToOneMap[v];

		// check if both vertices of the edges is in this tree
		bool VisInTree = false;
		bool UisInTree = false;
		for (int s = 0; s < this->num_vertices; s++)
		{
			if (VisInTree && UisInTree)
				break;

			if (v == this->vertices[s])
			{
				VisInTree = true;
				continue;
			}

			else if (u == this->vertices[s])
			{
				UisInTree = true;
				continue;
			}
		}

		if (VisInTree && UisInTree)
		{
			edges[edgeCount] = graph.addEdge(nodes[nodeU], nodes[nodeV]);
			edgeOneToOneMap[edgeCount++] = i;
		}
	}

	if (edgeCount >= this->num_vertices - 1) {
		// add a weight map
		ListGraph::EdgeMap<int> edgeWeight(graph);

		// set the weight of the edges
		for (int i = 0; i < edgeCount; i++)
		{
			int edgeId = edgeOneToOneMap[i];
			if (forbidden_edges != nullptr && forbidden_edges[edgeId]) {
				edgeWeight[edges[i]] = LARGE_WEIGHT;
			}
			else {
				edgeWeight[edges[i]] = g->edges[edgeId][2];
			}
			
		}

		// create a mst map
		ListGraph::EdgeMap<bool> mst(graph, false);

		// compute the minimum spanning tree
		kruskal(graph, edgeWeight, mst);


		// counter for the edges
		num_edges = 0;

		// weight of the tree
		this->weight = 0;

		// loop over the edges in mst and print the edges
		for (ListGraph::EdgeIt e(graph); e != INVALID; ++e)
		{
			if (mst[e])
			{
				// get the edge id
				int edgeId = edgeOneToOneMap[graph.id(e)];								

				// add the edge to the tree
				this->edges[num_edges] = edgeId;

				// check if e is forbidden 
				if (forbidden_edges!= nullptr && forbidden_edges[edgeId]) {
					num_forbidden_edges_in_mst++;
				}
				
				// add the weight of the edge to the tree
				this->weight += g->edges[edgeId][2];

				// increment the edge counter
				num_edges++;
			}
		}

		// check if the number of edges in the tree is equal to the number of vertices - 1
		if (num_edges != this->num_vertices - 1)
		{
			this->isSpanningTree = false;
		}
		else
		{
			this->isSpanningTree = true;		
		}
	}
	else
	{
		this->isSpanningTree = false;
	}

	// delete all arrays in this function
	delete[] nodes;
	delete[] edges;
	delete[] nodeOneToOneMap;
	delete[] edgeOneToOneMap;
}

void _tree::PrintTree()
{

	// get the instance graph
	_g* g = (_g*)this->instance;

	std::cout << "Vertices: ";
	for (int i = 0; i < this->num_vertices; i++)
	{
		std::cout << this->vertices[i] << " ";
	}	

	std::cout << "\t Edges: ";
	for (int i = 0; i < this->num_edges; i++)
	{
		// get the edge
		int edgeId = this->edges[i];

		// vertices of the edge
		int u = g->edges[edgeId][0];
		int v = g->edges[edgeId][1];

		std::cout << "(" << u << "," << v << ")";
	}

	std::cout << "\t Degrees:"; 
	for (int i = 0; i < this->num_vertices; i++)
	{
		std::cout << this->degree[i] << " ";
	}

	std::cout << "\t Weight: " << this->weight << std::endl;
}

// print vertices and weights
void _tree::PrintVerticesWeight() {
	// get the instance graph
	_g* g = (_g*)this->instance;
	std::cout << "Vertices: ";
	for (int i = 0; i < this->num_vertices; i++)
	{
		std::cout << this->vertices[i] << " ";
	}
	std::cout << "\t Weight: " << this->weight << std::endl;
}

// split into k trees; 
_tree** _tree::SplitIntoKTrees(int k, bool *forbidden_edges) {

	// get the instance graph
	_g* g = (_g*)this->instance;

	// print edges 
	for (int i = 0; i < this->num_vertices - 1; i++) {
		int e = this->edges[i];
		int u = g->edges[e][0];
		int v = g->edges[e][1];
		int w = g->edges[e][2];
		bool f = forbidden_edges[e];
		//std::cout << "Edge: " << u << " " << v << " " << w;
		//if (f) std::cout << " forbidden";	
		//std::cout << std::endl;		
	}

	// sort the edges based on their weights
	std::sort(this->edges, this->edges + (this->num_vertices - 1), [g, forbidden_edges](int a, int b) {
		if (forbidden_edges != NULL)
			return (g->edges[a][2] > g->edges[b][2] || forbidden_edges[a]);
		else
			return g->edges[a][2] > g->edges[b][2];
		});				

	// print edges 
	/*for (int i = 0; i < this->num_vertices - 1; i++) {
		int u = g->edges[this->edges[i]][0];
		int v = g->edges[this->edges[i]][1];
		int w = g->edges[this->edges[i]][2];
		bool f = forbidden_edges[this->edges[i]];
		std::cout << "Edge: " << u << " " << v << " " << w;
		if (f) std::cout << " forbidden";	
		std::cout << std::endl;		
	}*/

	// create a Matrix k-1 to n-k
	int ** matrix = new int* [k - 1];
	for (int i = 0; i < k - 1; i++) {
		matrix[i] = new int[this->num_vertices - k];
	}

	// create an array to store the edge removals ; -1 means the edge is removed
	bool* edge_to_be_removed = new bool[num_vertices - 1];	
	memset(edge_to_be_removed, 0, sizeof(bool)*num_vertices - 1);

	// records if the edge has been selected to be removed
	bool* edge_selected = new bool[this->num_vertices - 1];

	// weights of trees 
	int* tree_weights = new int[k];
	memset(tree_weights, 0, sizeof(int)*k);

	// max forest weight
	int maxForestWeight = 0;

	// number of vertices in each tree
	int* num_vertices = new int[k];

	// vertices in each tree
	int** vertices = new int* [k];
	uint32_t** bin_vertices = new uint32_t * [k];
	for (int i = 0; i < k; i++) {
		vertices[i] = new int[this->num_vertices];
		bin_vertices[i] = new uint32_t[binaryArrlength(this->num_vertices)];
	}

	// edges for trees
	int** edges = new int* [k];
	int* num_edges = new int[k];
	for (int i = 0; i < k; i++) {
		edges[i] = new int[this->num_vertices - 1];
	}


	// select the first k-1 edges to be removed
	for (int e = 0; e < k - 1; e++) {
		edge_to_be_removed[e] = true;
		edge_selected[e] = true;
	}
	
	//// print edges to be removed
	//for (int i = 0; i < this->num_vertices - 1; i++) {
	//	if (edge_to_be_removed[i]) {
	//		int u = g->edges[this->edges[i]][0];
	//		int v = g->edges[this->edges[i]][1];
	//		int w = g->edges[this->edges[i]][2];
	//		std::cout << "R Edge: " << u << " " << v << " " << w << std::endl;
	//	}
	//}

	// compute the max forest weight
	maxForestWeight = ComputeSpanningKForest(k, g, vertices, edges, edge_to_be_removed, num_vertices, num_edges, tree_weights, bin_vertices, forbidden_edges);

	int minMaxForestWeight = maxForestWeight;	

	while (true) {
		// min max forest weight
		int a_min = -1;
		int d_min = -1;

		int aa = 0;
		for (int a = 0; a < this->num_vertices - 1; a++) {
			if (!edge_to_be_removed[a]) continue;
			int dd = 0;
			for (int d = 0; d < this->num_vertices - 1; d++) {
				if (edge_to_be_removed[d]) continue;				

				edge_to_be_removed[a] = false;
				edge_to_be_removed[d] = true;

				//// print a and d
				//std::cout << "A: " << a << " D: " << d << std::endl;

				//// print edges to be removed
				//for (int i = 0; i < this->num_vertices - 1; i++) {
				//	if (edge_to_be_removed[i]) {
				//		int u = g->edges[this->edges[i]][0];
				//		int v = g->edges[this->edges[i]][1];
				//		int w = g->edges[this->edges[i]][2];
				//		std::cout << "R Edge: " << u << " " << v << " " << w << std::endl;
				//	}
				//}

				// remove the edge
				matrix[aa][dd] = ComputeSpanningKForest(k, g, vertices, edges, edge_to_be_removed, num_vertices, num_edges, tree_weights, bin_vertices, forbidden_edges);				

				

				//std::cout << "Max Forest Weight: " << matrix[aa][dd] << std::endl;

				// check if the weight is less than the min max forest weight
				if (matrix[aa][dd] < minMaxForestWeight) {
					minMaxForestWeight = matrix[aa][dd];
					a_min = a;
					d_min = d;
				}

				else if (matrix[aa][dd] == minMaxForestWeight) {
					if (!edge_selected[d]) {
						a_min = a;
						d_min = d;
					}
				}

				
				dd++;

				//reverse
				edge_to_be_removed[a] = true;
				edge_to_be_removed[d] = false;
			}
			aa++;
		}


		//// print the matrix
		//for (int a = 0; a < k - 1; a++) {
		//	for (int d = 0; d < this->num_vertices - k; d++) {
		//		std::cout << matrix[a][d] << " ";
		//	}
		//	std::cout << std::endl;
		//}

		if (a_min == -1) {			
			break;
		}
		else {
			edge_to_be_removed[a_min] = false;
			edge_to_be_removed[d_min] = true;
			edge_selected[d_min] = true;			
		}
			
	}
	maxForestWeight = ComputeSpanningKForest(k, g, vertices, edges, edge_to_be_removed, num_vertices, num_edges, tree_weights, bin_vertices, forbidden_edges);

	// update graph upper bound: Min (UB, maxForestWeight)
	if (maxForestWeight < g->UB) {
		g->UB = maxForestWeight;
		g->recomputeLB();
	}
	
	// Now we create the trees and output them
	// create an array of trees	

	_tree** trees = new _tree * [k];

	
	// create the trees
	for (int i = 0; i < k; i++) {
		// print bin vertices

		//if (num_edges[i] > 0) {
			trees[i] = new _tree(this->instance);		
			trees[i]->InitTreeByEdges(num_edges[i], edges[i]);
			//trees[i]->PrintTree();
		//}
	}

	// delete arrays
	delete[] num_vertices;
	for (int i = 0; i < k; i++) {
		delete[] vertices[i];
		delete[] bin_vertices[i];
		delete[] edges[i];
	}

	for (int i = 0; i < k - 1; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;

	delete[] vertices;
	delete[] bin_vertices;

	delete[] edges;
	delete[] num_edges;

	delete[] edge_to_be_removed;
	delete[] edge_selected;
	delete[] tree_weights;

	// return the trees
	return trees;
}

int _tree::ComputeSpanningKForest(int k, void* instance, int** vertices, int** edges, bool* edge_to_be_removed, int* num_vertices, int* num_edges, int* tree_weights, uint32_t ** bin_vertices, bool *forbidden_edges) {	

	// get the instance graph
	_g* g = (_g*)instance;

	// traverse the tree to decompose the tree
	TraverseToDecompose(0, -1, 0, 0, edge_to_be_removed);

	memset(num_vertices, 0, sizeof(int) * k);
	memset(num_edges, 0, sizeof(int) * k);
	memset(tree_weights, 0, sizeof(int) * k);

	// reset bin_vertices
	for (int i = 0; i < k; i++) {
		memset(bin_vertices[i], 0, sizeof(uint32_t) * binaryArrlength(this->num_vertices));
	}

	for (int v = 0; v < this->num_vertices; v++) {
		int i = this->decomp[v];
		int vertex = this->vertices[v];
		vertices[i][num_vertices[i]++] = vertex;
		addbin(bin_vertices[i], vertex);
	}

	// compute the edges for each tree
	for (int e = 0; e < this->num_vertices - 1; e++) {
		// get the edge
		int edgeId = this->edges[e];

		// get the vertices of the edge
		int u = g->edges[edgeId][0];
		int v = g->edges[edgeId][1];	

		// get the tree of the vertices
		int treeU = this->decomp[u];
		int treeV = this->decomp[v];

		// if the vertices are in the same tree
		if (treeU == treeV) {
			// add the edge to the tree
			edges[treeU][num_edges[treeU]++] = edgeId;
			// weight of tree
			if (!forbidden_edges[edgeId])
				tree_weights[treeU] += g->edges[edgeId][2];
			else
				tree_weights[treeU] += LARGE_WEIGHT;
		}		
	}

	int maxForestWeight = 0;
	// compute maximum Forest weight
	for (int i = 0; i < k; i++) {

		if (tree_weights[i] > maxForestWeight) {
			maxForestWeight = tree_weights[i];
		}

		//// print the tree
		//std::cout << "Tree " << i << " Weight: " << tree_weights[i] << " ";

		//// print edges
		//for (int j = 0; j < num_edges[i]; j++) {
		//	int edgeId = edges[i][j];
		//	int u = g->edges[edgeId][0];
		//	int v = g->edges[edgeId][1];
		//	std::cout << "(" << u << "," << v << ")" << " "; 
		//}

		//std::cout << std::endl;	
	}

	//// print the maximum forest weight
	//std::cout << "Max Forest Weight: " << maxForestWeight << std::endl;

	return maxForestWeight;
}

// traverse to decompose
int _tree::TraverseToDecompose(int vertex, int prev_vertex, int i, int largest_i,  bool* edge_to_be_removed) {
	// get the instance graph
	_g* g = (_g*)this->instance;

	// choose the vertex to tree i
	this->decomp[vertex] = i;

	// print decompsition	
	//std::cout << "Vertex: " << vertex << " Tree: " << i << std::endl;

	int connector; // the vertex to connect to

	int new_tree_vertex[4]; // the vertex to create a new tree
	int new_tree_vertex_count = 0;		

	// check all edges
	for (int e = 0; e < num_vertices - 1; e++) {

		// get the edge
		int edgeId = this->edges[e];

		// get the vertices of the edge
		int u = g->edges[edgeId][0];
		int v = g->edges[edgeId][1];


		connector = -1;

		// if edge is adjacent to the vertex
		if (vertex == u) {
			connector = v;
		}
		else if (vertex == v) {
			connector = u;
		}

		// check if the vertex is connected
		if (connector != -1 && connector != prev_vertex) {

			// if edge is to be excluded
			if (edge_to_be_removed[e]) {
				// create a new tree
				new_tree_vertex[new_tree_vertex_count++] = connector;
			}
			else {
				largest_i = TraverseToDecompose(connector, vertex, i, largest_i, edge_to_be_removed);

				// print connector vertex i largest i 
				//std::cout << "Connector: " << connector << " Vertex: " << vertex << " i: " << i << " Largest i: " << largest_i << std::endl;

			}

		}

	}

	// if a new tree is to be created
	for (int s = 0; s < new_tree_vertex_count; s++) {
		largest_i = TraverseToDecompose(new_tree_vertex[s], vertex, largest_i +1, largest_i+1, edge_to_be_removed);
	}

	return largest_i;
}

// add vertex
void _tree::AddVertex(int vertex)
{
	// if vertex already is in the tree, increase the degree
	if (checkbin(this->bin_vertices, vertex)) {
		// find the vertex in the tree
		for (int i = 0; i < this->num_vertices; i++) {
			if (this->vertices[i] == vertex) {
				degree[i]++;
				break;
			}
		}
		return;
	}

	// if vertex is new, add it to the tree
	degree[this->num_vertices] = 1;
	this->vertices[this->num_vertices++] = vertex;	
	addbin(this->bin_vertices, vertex);	
}

// remove vertex
void _tree::RemoveVertex(int vertex)
{
	// get the instance graph
	_g* g = (_g*)this->instance;

	// get the vertex index in the tree
	int vertexInTree = -1;
	for (int i = 0; i < num_vertices; i++) {
		if (vertices[i] == vertex) {
			vertexInTree = i;
			break;
		}
	}

	// if the vertex is not in the tree
	if (vertexInTree == -1) {
		return;
	}

	// if vertex will be still connected to the tree, just decrease the degree
	if (degree[vertexInTree] > 1) {
		degree[vertexInTree]--;
		return;
	}

	// if degree becomes zero remove the vertex
	for (int i = vertexInTree; i < num_vertices - 1; i++) {
		vertices[i] = vertices[i + 1];
		degree[i] = degree[i + 1];
	}

	// remove the vertex from the binary array
	removebin(this->bin_vertices, vertex);

	// decrement the number of vertices
	num_vertices--;
}

// add edge
void _tree::AddEdge(int edge)
{	
	if (IsEdgeInTree(edge)) {
		return;
	}

	this->edges[num_edges++] = edge;
	int u = ((_g*)this->instance)->edges[edge][0];
	int v = ((_g*)this->instance)->edges[edge][1];

	AddVertex(u);	
	AddVertex(v);	

	this->weight += ((_g*)this->instance)->edges[edge][2];

	//// print 
	//std::cout << "Edge Added: " << u << " " << v << std::endl;
}

// remove edge
void _tree::RemoveEdge(int edge)
{
	// get the instance graph
	_g* g = (_g*)this->instance;

	// get the vertices of the edge
	int u = g->edges[edge][0];
	int v = g->edges[edge][1];	

	RemoveVertex(u);
	RemoveVertex(v);

	// remove the edge
	for (int i = 0; i < num_edges; i++) {
		if (edges[i] == edge) {
			while (i < num_edges - 1) {
				edges[i] = edges[i + 1];
				i++;
			}		
			break;
		}
	}
	num_edges--;
}

// check if the edge is in the tree
bool _tree::IsEdgeInTree(int edge)
{
	for (int i = 0; i < num_vertices - 1; i++) {
		if (edges[i] == edge) {
			return true;
		}
	}
	return false;
}

// incident vertex
// this function returns that is connected to the tree. If both vertices are connected, it will return -1 
int _tree::IncidentVertex(int edge) 
{
	int u = ((_g*)this->instance)->edges[edge][0];
	int v = ((_g*)this->instance)->edges[edge][1];

	bool uInTree = false;
	bool vInTree = false;

	if (checkbin(this->bin_vertices, u)) {
		uInTree = true;
	}
	if (checkbin(this->bin_vertices, v)) {
		vInTree = true;
	}
	
	if (uInTree && !vInTree) {
		return u;
	}
	else if (!uInTree && vInTree) {
		return v;
	}	
	else {
		return -1;
	}

}

// any incident vertex 
bool _tree::GetAnyIncidentVertex(int edge) {  // returns a vertex of the edge that is connected to the tree; if both vertices are connected, returns any of the two
	int u = ((_g*)this->instance)->edges[edge][0];
	int v = ((_g*)this->instance)->edges[edge][1];

	bool uInTree = false;
	bool vInTree = false;

	if (checkbin(this->bin_vertices, u)) {
		uInTree = true;
	}
	if (checkbin(this->bin_vertices, v)) {
		vInTree = true;
	}

	if (uInTree || vInTree) {
		return true;
	}
	else {
		return false;
	}

}

// compute the degree of the vertices in the tree
void _tree::RecomputeDegree()
{
	// get the instance graph
	_g* g = (_g*)this->instance;

	// initialize the degree
	memset(degree, 0, sizeof(int)*num_vertices);

	// loop over the edges
	for (int i = 0; i < num_edges; i++) {
		// get the edge
		int edgeId = edges[i];

		// get the vertices of the edge
		int u = g->edges[edgeId][0];
		int v = g->edges[edgeId][1];

		// get the vertices in the tree
		int uInTree = -1;
		int vInTree = -1;

		for (int j = 0; j < num_vertices; j++) {
			if (vertices[j] == u) {
				uInTree = j;
			}
			if (vertices[j] == v) {
				vInTree = j;
			}
		}

		// increment the degree
		degree[uInTree]++;
		degree[vInTree]++;
	}
}

// check if the tree is a spanning tree
bool _tree::IsSpanningTree()
{
	// get the instance graph
	_g* g = (_g*)this->instance;
	// check if the number of edges is equal to the number of vertices - 1
	if (num_edges == num_vertices - 1) {
		// check if the tree is connected
		memset(visited, 0, sizeof(bool) * num_vertices);
		// create a stack		
		int stackCount = 0;
		// push the first vertex to the stack
		stack[stackCount++] = vertices[0];
		visited[0] = true;
		// loop over the stack
		while (stackCount > 0) {
			// pop the vertex
			int vertex = stack[--stackCount];
			// loop over the edges
			for (int i = 0; i < num_edges; i++) {
				// get the edge
				int edgeId = edges[i];
				// get the vertices of the edge
				int u = g->edges[edgeId][0];
				int v = g->edges[edgeId][1];
				// check if the edge is incident to the vertex
				if (u == vertex || v == vertex) {
					// get the other vertex
					int otherVertex = u == vertex ? v : u;
					// get the index of the other vertex
					int otherVertexIndex = -1;
					for (int j = 0; j < num_vertices; j++) {
						if (vertices[j] == otherVertex) {
							otherVertexIndex = j;
							break;
						}
					}
					// check if the other vertex is not visited
					if (!visited[otherVertexIndex]) {
						visited[otherVertexIndex] = true;
						stack[stackCount++] = otherVertex;
					}
				}
			}
		}
		// check if all vertices are visited
		for (int i = 0; i < num_vertices; i++) {
			if (!visited[i]) {				
				return false;
			}
		}		
		return true;
	}
	else {
		return false;
	}
}

