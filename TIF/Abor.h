#pragma once
#include "Cplex.h"
class Abor :
    public Cplex
{
public:

	// variables
	IloNumVarArray x; // x variables
	IloNumVarArray y; // y variables

	// costs 
	IloNumArray costs; // cost of each arc
	IloNumArray vertex_prizes; // prize of each vertex

	// upper bound on weight
	IloNum upper_bound_weight; // upper bound on the weight of each vertex


	// objective
	IloObjective objective; // objective function

	// upper_bound constraint
	IloRange upper_bound; // upper bound constraint


	// arrays used in callback
	double* x_value; // value of x variables
	double* y_value; // value of y variables

	int* S_r; 
	int* S_v;

	// constructor
	Abor(_g* g, bool redirect, bool linear);
	~Abor();
	// run the algorithm
	Abor* Run();
	// define variables
	void DefVars();
	void DefVarX();
	void DefVarY();
	
	// add objective function
	void AddObj();
	
	// add constraints
	void AddCons();	

	// constraints that makes sure that one incoming arc of each vertex is selected
	void AddConsOneIncomingArc();

	// constraints for imposing upper bounds of the tree 
	void AddConsUpperBound();

	// constraints for imposing exactly one arc from root
	void AddConsExactlyOneArcFromRoot();

	// constraints for imposing that if a not terminal vertex is selected, then at least one of its outgoing arcs is selected
	void AddConsAtLeastOneOutgoingArcTerminal();

	// we impose loops of size two
	void AddConsLoopsOfSizeTwo();


	Abor* Init(_pcst* pcst);

	Abor* PrintSol();
	Abor* PrintModel();
};



