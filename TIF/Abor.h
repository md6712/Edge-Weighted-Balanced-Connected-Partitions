#pragma once
#include "Cplex.h"
#include "CG_branch.h"

struct AborSol {	
	double value; // optimal value
	_small_tree* tree; // tree associated to the solution
};

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

	// objective constraint > 0
	IloRange objective_constraint; // objective constraint

	// upper_bound constraint
	IloRange upper_bound; // upper bound constraint

	// range array for branch constraints
	IloRangeArray branch_constraints; // branch constraints

	// arrays used in callback
	double* x_value; // value of x variables
	double* y_value; // value of y variables

	// arrays used in the model : in the callback
	int* S_r; 
	int* S_v;

	// solution 
	AborSol *opt; 
	bool solution_found; // solution found
	

	// constructor
	Abor(_g* g, bool redirect, bool linear);
	~Abor();
	// run the algorithm
	AborSol* Run(bool printCuts);
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


	// void save solution
	void SaveOpt();

	// void reset solution
	void ResetOpt();


	Abor* Init(_pcst* pcst, double fixed_cost);	

	Abor* PrintSol();
	Abor* PrintModel();


	Abor* AddConstraintsBranching(BP_node* node);

};



