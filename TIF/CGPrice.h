#pragma once
#include "Cplex.h"
#include "Abor.h"

class CGPrice :
    public Cplex
{
    public:
	double* y_value;
    double* x_value;
	bool halted = false;
    double ub_positive = 0;

    // branch and bound 
    int num_vertices_pos_zeta; 
    int* vertcies_pos_zeta; 
    double* vertex_weights; 
    double* vertex_prizes; 
	double* edge_costs;

    int* roots_pcst;

	double best_mwcs = -INFINITY;
	double best_pcst = -INFINITY;   

	// Aborescence instance
	Abor* aborescence = NULL;

    IloNum theta;
    IloNumArray eta; 
    IloNumArray zeta;

    IloNumVarArray x;
    IloNumVarArray y;

    // martix variables
    NumVarMatrix phi;

    // objective function
    IloObjective objective;        

    CGPrice(_g*, bool);
    ~CGPrice();

    CGPrice* Run();

    void DefVar();
    void DefVarX();
    void DefVarY();
    void DefVarPhi();

    void AddObj();
    void AddObj2();

    void AddCons();
    void AddConsNumEdgesVerticesMinusOne();
    void XYRelation();
    void PhiRelation();
    void AtleastOneEdge();

    // cost of tree extra constraint
    void AddConsCostOfTree();

    // print the model
    CGPrice* PrintModel();

    // print the solution
    CGPrice* PrintSol();

    // set dual values
    CGPrice* setTheta(IloNum);
    CGPrice* setEta(IloNumArray);
    CGPrice* setZeta(IloNumArray);
    CGPrice* UpdateObjectiveCoefficients();

    // get the tree associated to the optimal solution
    _small_tree* GetTree();

    // set print cuts and cycles
    CGPrice* SetPrintCuts(bool);
    CGPrice* SetPrintCycles(bool);

	// fix the solution
	CGPrice* FixSol();

    // Heuritic 
    _tree* heuristic();

	double ComputeReward(_tree * tree);

   
    // solve using maximum weight connected subgraph
    CGPrice* solve_mwcs();

	// solve using prize collecting steiner tree
	CGPrice* solve_pcst();
};

