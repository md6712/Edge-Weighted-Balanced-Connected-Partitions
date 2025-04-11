#pragma once
#include "_tree.h"
#include <windows.h>
#include <iostream>
#include <sstream>
#include <string>

class _mwcs
{
public:

	void* instance; // this is going to be the graph

	int num_vertices;   // this is equal to num_vertices + num_edges of original graph
	int num_edges;  // this is equal double the size of the original graph

	double* vertex_weight; // array to store the weight of each vertex
	int(* edges)[2]; // 2D array to store the edges

	int* vertex_map;  // maps the associated edge in the original graph; -1 for vertices
	char* stringfile; // string to store file data	

	_small_tree* tree; // the tree associated to the solution

	int* vertices_in_tree; // vertices in the tree

	double theta; // dual value of the knapsack constraint

	double objective_value; // objective value of the solution



	// I didnt finish this PIPE IDEA!
	// pipes
	HANDLE hNamedPipe_input; // handle to the named pipe for input
	HANDLE hNamedPipe_output; // handle to the named pipe for output

	HANDLE hStdOutRead, hStdOutWrite; // handles to the standard output
	HANDLE hStdInRead, hStdInWrite; // handles to the standard input

	SECURITY_ATTRIBUTES saAttr;

	const char pipeName_input[50] = "\\\\.\\pipe\\mwcs_input"; // name of the named pipe for input
	const char pipeName_output[50] = "\\\\.\\pipe\\mwcs_output"; // name of the named pipe for output

	
	
	/// METHODS
	
	// constructor 
	_mwcs(void *g);

	// destructor 
	~_mwcs();

	// set the weights for the vertices 
	_mwcs* set_vertex_weights(double *, double);

	_mwcs* solve();

	_small_tree* get_tree();

	_mwcs* print();

	_mwcs* print_in_file();

	_mwcs* read_solution_from_file();


	// I didn't finish this PIPE IDEA! Nust be implemented in the future if needed.
	// Instead we try to use PCST implementation
	_mwcs* pipes_ready();

	_mwcs* pipes_close();

	_mwcs* pipes_write();

	_mwcs* pipes_solve();

	_mwcs* pipes_read();
};

