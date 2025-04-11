#pragma once
#include <cstdint>


#include <lemon/smart_graph.h>
#include <lemon/preflow.h>
#include <lemon/edmonds_karp.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/elevator.h>
#include <lemon/list_graph.h>

#include "_tree.h"
#include "_mwcs.h"
#include "_pcst.h"

#include "_render.h"

using namespace lemon;

using Graph = ListDigraph;

typedef Graph::Node Node;	
typedef Graph::Arc DEdge;

enum color {
	white,
	blue, 
	red, 
	green, 
	yellow
};

struct _coord {
	double x;
	double y;
	double dx;
	double dy;
};



class _g
{
	
	char* filename; 

	public:	

		Graph dg;
		_render* render;

		int num_vertices;
		int num_edges;
		int num_arcs; 
		int num_trees; // number of trees in the solution
		int(* edges)[3]; // 2D array to store the edges
		int(* arcs)[3]; // 2D array to store the arcs
		int *opt_edges;
		int num_opt_edges;
		int *start_edge_tree;

		_coord* coords; // array to store the coordinates of the vertices

		int** opt_vertices; // 2D array to store the vertices of the optimal solution - for each tree
		int* num_opt_vertices;

		bool subsetsComputed;

		bool** cycles;
		bool* visited;
		bool* edge_visited;
		int num_cycles;
		int cycle_start_vertex;
		bool in_cycle;

		_mwcs* mwcs = nullptr; // the associated minimum weight connected subgraph problem

		_pcst* pcst = nullptr; // the associated prize collecting steiner tree problem		

		std::vector <_tree*> trees; // vector to store the trees in the forest

		std::vector <_small_tree*> select_trees_for_CG; // vector to store the selected trees in the forest

		uint32_t** subsets; // 2D array to store the subsets of vertices
		int* nSubsets; // array to store the number of subsets for each tree

		// MST 
		int* parent_mst;  // array to store the parent of each vertex in the MST
		uint32_t* bin_vertices; // array to store the vertices in binary format
		int* sum_adjacent_weight; // array to store the sum weight of the adjacent edges for each vertex

		// UV Seperator
		int* seperator;
		int num_seperators;	
		color* v_colors;


		int UB = INT32_MAX;  // upper bound on the size of maximum tree in the forest
		int LB = 0; // lower bound on the size of smallest tree in the forest
		int MST_weight = 0; // weight of the minimum spanning tree
		

		double _gap;
		double _opt;
		_g(int num_vertices, int num_edges, int num_trees);
		~_g();

		void readGraph();
		void createMWCS();
		void createPCST();
		void printGraph();

		void setFilename();

		void computeSubsets();
		void printSubsets();
		void outputOPTEdges();

		char* getFilename();
		
		bool CheckCyclesInOptEdges();
		bool CheckCyclesInOptEdgesRec(int i, int ep, int v);

		void PrintCycles();
		void PrintOptEdges(bool printWeight = false);

		_g* SortEdges();

		int  MST(int tree, int* mst_vertices, int n);
		bool UVSeperator(int u, int v);
		void ColorNeighbors(int v, color c, bool transitive);
		bool Seperated(int u, int v);
		void PrintMinSeperators();

		void setDiGraph();	

		// generate trees
		void generateTrees();
		void generateSelectTrees();

		void recomputeLB();
		void computeUB();

		void ForcedDirectedLayout(); // force directed layout algorithm to define the coordinates of the vertices in the graph
		void DrawGraph(int* highlighted_edges = nullptr, int num_highlighted_edges = 0); // draw the graph

};

