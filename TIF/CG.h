#pragma once
#include "Cplex.h"
class CG :
    public Cplex
{

    public:
    //variables
    IloNumVarArray x;
    IloNumVar z;

    // constraints 
    IloRangeArray Knapsack;
    IloRangeArray Cover;
    IloRangeArray Z;

    IloObjective obj;


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
    void AddVar();

    // add objective function
    void AddObj();

};

