#pragma once
#include "Cplex.h"
class SetCoverF :
    public Cplex
{
public:
	IloNumVarArray x;  // variable which is one if tree t is selected
	
	
	IloNumVar z; // computes the heaviest tree


	// define all variables
	void DefVar();
	void DefVarX();	
	void DefVarZ();	

	// define the objective function
	void AddObj();

	// define the constraints
	void AddCons();
	void SetCover();
	void SetZ();
	void KnapSack();



	// constructor and destructor
	SetCoverF(_g*, bool, bool);
	~SetCoverF();


	
	// run the model
	SetCoverF* Run();

	// print the solution
	SetCoverF* PrintSol();
	SetCoverF* PrintModel();

	// set the initial solution
	SetCoverF* SetInitSol();

	// force the solution
	SetCoverF* ForceSol();

	// save the optimal solution
	void SaveOpt();

	// set integer or linear
	SetCoverF* SetInteger();
	SetCoverF* SetLinear();
};

