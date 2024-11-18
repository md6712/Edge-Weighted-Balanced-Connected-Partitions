#pragma once
#include "Cplex.h"
class ExF :
    public Cplex
{

    public:
	NumVarMatrix x; // variable is used to assign edges to trees
	NumVarMatrix xbar;  // variable is used to assign arcs to trees
	NumVarMatrix y; // variable is used to define the root vertex of each tree
	NumVarMatrix z; // variable is used to assign vertices to trees

	// define flow variable cubic for each arc, each vertex, and each tree
	NumVarCube f; // variable is used to define the flow of each arc

	// define variable for the objective function
	IloNumVar eta; // variable is used to define the objective function

	

	// define all variables
	void DefVar(); 
	void DefVarX();  
	void DefVarXbar();  
	void DefVarY();  
	void DefVarZ();  
	void DefVarEta();  
	void DefVarF();

	// define the objective function
	void AddObj();

	// define the constraints
	void AddCons();
	void AssignVerticesToTrees();
	void DecideRootOfTrees();
	void SetYZRelation();
	void SetXZRelation();
	void SumOfFlowFromSToVEqualsZui();
	void FlowConservation(); 
	void FlowConservationUV();
	void FlowXBarRelation();
	void Flow0vYRelation();
	void XbarZRelation();
	void XXbarRelation();






	// constructor and destructor	
	ExF(_g*, bool, bool);
	~ExF();

	// run the model
	ExF* Run();

	// print the solution
	ExF* PrintSol();
	ExF* PrintModel();

	// set the initial solution
	ExF* SetInitSol();

	// force the solution
	ExF* ForceSol();

	// save the optimal solution
	void SaveOpt();

};

