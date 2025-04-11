#pragma once
#include "Cplex.h"
#include "CGPrice.h"
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
    CG* Run();

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
    void AddVars();
    void AddVar(int i);

    // add objective function
    void AddObj();

};

