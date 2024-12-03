#pragma once
#include "Cplex.h"
class CGPrice :
    public Cplex
{
    public:
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
    _tree* GetTree(int *vertices, uint32_t *bin_vertices);

    // set print cuts and cycles
    CGPrice* SetPrintCuts(bool);
    CGPrice* SetPrintCycles(bool);

    // Heuritic 
    _tree* heuristic();

};

