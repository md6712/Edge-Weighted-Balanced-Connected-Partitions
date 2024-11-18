#include "_tree.h"
#include <lemon/kruskal.h>
#include <lemon/list_graph.h>
#include <iostream>
#include <algorithm>

#include "_g.h"


_tree::_tree(int* vertices, int num_vertices, uint32_t bin_vertices, void *instance)
{
	this->num_vertices = num_vertices;
	this->edges = new int[num_vertices - 1];
	this->vertices = new int[num_vertices];
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