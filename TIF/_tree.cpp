#include "_tree.h"
#include <lemon/kruskal.h>
#include <lemon/list_graph.h>
#include <iostream>
#include <algorithm>
#include "binary.h"

#include "_g.h"


_tree::_tree(int* vertices, int num_vertices, uint32_t* bin_vertices, void *instance)
{
	this->num_vertices = num_vertices;
	this->edges = new int[num_vertices - 1];
	this->vertices = new int[num_vertices];
	this->decomp = new int[num_vertices];
	memcpy(this->vertices, vertices, sizeof(int)*num_vertices);	
	this->bin_vertices = bin_vertices;
	this->weight = 0;
	this->instance = instance;
	this->isSpanningTree = false;
}

_tree::~_tree()
{
	delete[] this->edges;
	delete[] this->vertices;
}

void _tree::singltonTree(int vertex)
{
	this->num_vertices = 1;
	this->vertices = new int[1];
	this->vertices[0] = vertex;	
	this->weight = 0;
	this->isSpanningTree = true;
}

void _tree::fixEdges(int* edges, int num_edges)
{
	this->edges = new int[num_edges];
	memcpy(this->edges, edges, sizeof(int)*num_edges);
	
	// compute weight 
	_g* g = (_g*)this->instance;
	for (int i = 0; i < num_edges; i++)
	{
		this->weight += g->edges[edges[i]][2];
	}
}

// compute the minimum spanning tree

void _tree::computeMST()
{

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
		int u = g->edges[i][0];
		int v = g->edges[i][1];

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
			edgeWeight[edges[i]] = g->edges[edgeId][2];
		}

		// create a mst map
		ListGraph::EdgeMap<bool> mst(graph, false);

		// compute the minimum spanning tree
		kruskal(graph, edgeWeight, mst);


		// counter for the edges
		int edgeCounter = 0;

		// loop over the edges in mst and print the edges
		for (ListGraph::EdgeIt e(graph); e != INVALID; ++e)
		{
			if (mst[e])
			{
				// get the edge id
				int edgeId = edgeOneToOneMap[graph.id(e)];								

				// add the edge to the tree
				this->edges[edgeCounter] = edgeId;
				
				// add the weight of the edge to the tree
				this->weight += g->edges[edgeId][2];

				// increment the edge counter
				edgeCounter++;
			}
		}

		// check if the number of edges in the tree is equal to the number of vertices - 1
		if (edgeCounter != this->num_vertices - 1)
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


void _tree::printTree()
{

	// get the instance graph
	_g* g = (_g*)this->instance;

	std::cout << "Vertices: ";
	for (int i = 0; i < this->num_vertices; i++)
	{
		std::cout << this->vertices[i] << " ";
	}	

	std::cout << "\t Edges: ";
	for (int i = 0; i < this->num_vertices - 1; i++)
	{
		// get the edge
		int edgeId = this->edges[i];

		// vertices of the edge
		int u = g->edges[edgeId][0];
		int v = g->edges[edgeId][1];

		std::cout << "(" << u << "," << v << ")";
	}

	std::cout << "\t Weight: " << this->weight << std::endl;
}


// split into k trees; 
_tree** _tree::splitIntoKTrees(int k) {

	// get the instance graph
	_g* g = (_g*)this->instance;

	// sort the edges based on their weights
	std::sort(this->edges, this->edges + (this->num_vertices - 1), [g](int a, int b) {
		return g->edges[a][2] > g->edges[b][2];
		});

	// print edges 
	for (int i = 0; i < this->num_vertices - 1; i++) {
		int u = g->edges[this->edges[i]][0];
		int v = g->edges[this->edges[i]][1];
		int w = g->edges[this->edges[i]][2];
		std::cout << "Edge: " << u << " " << v << " " << w << std::endl;		
	}

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
	
	// compute the max forest weight
	maxForestWeight = computeSpanningKForest(k, g, vertices, edges, edge_to_be_removed, num_vertices, num_edges, tree_weights, bin_vertices);

	while (true) {
		// min max forest weight
		int minMaxForestWeight = maxForestWeight;
		int a_min = -1;
		int d_min = -1;

		for (int a = 0; a < k - 1; a++) {
			if (!edge_to_be_removed[a]) continue;
			int b = 0;
			for (int d = 0; d < this->num_vertices - 1; d++) {
				if (edge_to_be_removed[d]) continue;				

				edge_to_be_removed[a] = false;
				edge_to_be_removed[d] = true;

				// remove the edge
				matrix[a][b] = computeSpanningKForest(k, g, vertices, edges, edge_to_be_removed, num_vertices, num_edges, tree_weights, bin_vertices);

				// check if the weight is less than the min max forest weight
				if (matrix[a][b] < minMaxForestWeight) {
					minMaxForestWeight = matrix[a][b++];
					a_min = a;
					d_min = d;
				}

				else if (matrix[a][b] == minMaxForestWeight) {
					if (!edge_selected[d]) {
						a_min = a;
						d_min = d;
					}
				}

				b++;

				//reverse
				edge_to_be_removed[a] = true;
				edge_to_be_removed[d] = false;
			}
		}

		if (a_min == -1) {
			break;
		}
		else {
			edge_to_be_removed[a_min] = false;
			edge_to_be_removed[d_min] = true;

			maxForestWeight = computeSpanningKForest(k, g, vertices, edges, edge_to_be_removed, num_vertices, num_edges, tree_weights, bin_vertices);
		}
	}

	// print the matrix
	for (int a = 0; a < k - 1; a++) {
		for (int d = 0; d < this->num_vertices - k; d++) {
			std::cout << matrix[a][d] << " ";
		}
		std::cout << std::endl;
	}

	// Now we create the trees and output them
	// create an array of trees
	_tree** trees = new _tree*[k];
	
	// create the trees
	for (int i = 0; i < k; i++) {
		// create a new tree
		trees[i] = new _tree(vertices[i], num_vertices[i], bin_vertices[i], this->instance);
		trees[i]->fixEdges(edges[i], num_edges[i]);
		trees[i]->printTree();
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

	delete[] edge_to_be_removed;
	delete[] edge_selected;
	delete[] tree_weights;



	// return the trees
	return trees;
}

int _tree::computeSpanningKForest(int k, void* instance, int** vertices, int** edges, bool* edge_to_be_removed, int* num_vertices, int* num_edges, int* tree_weights, uint32_t ** bin_vertices) {

	// get the instance graph
	_g* g = (_g*)instance;

	// traverse the tree to decompose the tree
	traverseToDecompose(0, -1, 0, edge_to_be_removed);

	memset(num_vertices, 0, sizeof(int) * k);
	memset(num_edges, 0, sizeof(int) * k);
	memset(tree_weights, 0, sizeof(int) * k);

	for (int v = 0; v < this->num_vertices; v++) {
		int i = this->decomp[v];
		int vertex = this->vertices[v];
		vertices[i][num_vertices[i]++] = vertex;
		addbin(bin_vertices[i], vertex);
	}

	// compute the edges for each tree
	for (int e = 0; e <= this->num_vertices - 2; e++) {
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
			tree_weights[treeU] += g->edges[edgeId][2];
		}
	}

	int maxForestWeight = 0;
	// compute maximum Forest weight
	for (int i = 0; i < k; i++) {

		if (tree_weights[i] > maxForestWeight) {
			maxForestWeight = tree_weights[i];
		}

		std::cout << "Tree " << i << " Weight: " << tree_weights[i] << std::endl;
	}

	// print the maximum forest weight
	std::cout << "Max Forest Weight: " << maxForestWeight << std::endl;

	return maxForestWeight;
}

// traverse to decompose
void _tree::traverseToDecompose(int vertex, int prev_vertex, int i, bool* edge_to_be_removed) {
	// get the instance graph
	_g* g = (_g*)this->instance;

	// choose the vertex to tree i
	this->decomp[vertex] = i;

	//// print decompsition	
	//std::cout << "Vertex: " << vertex << " Tree: " << i << std::endl;

	int connector; // the vertex to connect to

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
				traverseToDecompose(connector, vertex, i + 1, edge_to_be_removed);
			}
			else {
				traverseToDecompose(connector, vertex, i, edge_to_be_removed);
			}

		}

	}
}


