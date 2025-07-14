#pragma once
#include "Cplex.h"
#include "CGPrice.h"
#include "CG_branch.h"
#include "timer.h"



class CG :
    public Cplex
{

    public:
    //variables
    IloNumVarArray x;
    IloNumVar z;        

    // constraints 
    IloRange Knapsack;
    IloRangeArray Cover;
    IloRangeArray Z;

    // dual values for constraints
    IloNum dualKnapsack;   // theta in the model
    IloNumArray dualCover;      // eta in the model
    IloNumArray dualZ;   	    // zeta in the model   

    IloObjective obj;

    // pricing problem
    CGPrice* pricing_problem;

    // run the model
    CG* Run_CG(BP_node* node, Timer timer_CG);

	CG* Run_BP();

	CG* Impose_Braching(BP_node* node);


    // matrix total assignment on pairs of vertices 
	double** matrix_total_assignment; // matrix of total assignment
	int best_pair[2]; // best pair of vertices

	// sorted linked list to keep the BP nodes
	//BP_node* root = nullptr;

	// vector of active nodes   
	std::vector<BP_node*> active_nodes;
    
    // vector or premature nodes
	std::vector<BP_node*> premature_nodes;

	// vector of nodes to be removed
	std::vector<BP_node*> nodes_to_remove;




    // sorted indices for + reduced cost reduction
	int* sorted_indices;

    // constructor and destructor   
    CG(_g*, bool, bool);
    ~CG();

    // print the solution
    CG* PrintSol();

    // print the model
    CG* PrintModel();

    // define constraints
    void AddCons();

    // add variables 
    void AddVars(bool integer);
	
    void AddVar(int i, bool integer);

	// add x variables
	void AddXVars(bool integer);

	// remove variables
	void RemoveXVars();

    // update UB
    void UpdateUB();

    // remove variables with large positive reduced cost 
	void RemoveLargePositiveReducedCost(IloNumArray reducedCosts);

    // update recursive prematurity
	void UpdatePrematureParent(BP_node* node);

	// update confirmed lower bound
	void UpdateConfirmedLB(BP_node* node);

    // add objective function
    void AddObj();

	// void compute total assignment
	void ComputeTotalAssignment(BP_node* node);

	// void print total assignment
	void PrintTotalAssignment();

    // add node to linked list
	void AddNodeSorted(BP_node* node);

	// add artificial columns
	void AddArtificialColumns(BP_node* node);

    // print NodeList   
	void PrintNodeList();

    // print node
	void PrintNode(BP_node* node);
    

};

