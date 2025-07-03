#pragma once
#include "Cplex.h"
#include "_g.h"

class FlowF :
    public Cplex
{

public:

	// counts for the user and lazy callback
	int usercallback_count = 0;
	int lazycallback_count = 0;

	bool user_cuts_active = false;
	bool mutual_exclusion_cycles = false;



	// variables for cplex
	IloNumVarArray x; 
	IloNumVarArray y;
	IloNumVarArray f;
	IloNumVarArray f0;

	// theta
	NumVarMatrix theta;

	IloNumVar eta; // the objective function

	// constructor and destructor
	FlowF(_g*, bool, bool, bool);
	~FlowF();

	// define variables
	void DefVar();
	void DefVarX();
	void DefVarY();
	void DefVarF();
	void DefVarF0();
	void DefVarEta();
	void DefVarTheta();


	// add objective function
	void AddObj();

	// add constraints
	void AddCons();
	void AddConsComputeEta();
	void AddConsSelectKTreeRepresentatives();
	void AddConsSelectArcForRepresentatives();
	void AddConsFlow();
	void AddConsBoundFByX();
	void AddConsBoundF0ByY();
	void AddConsSymmetryRootVertex();

	void AddHandeSuggestion();
	
	void AddConsCycleXY();

	// constraints for theta 
	// constraints for theta and f0
	void AddConsThetaF0();
	// constraints for theta and eta
	void AddConsThetaEta();


	// branching
	void AddPriorityY();

	void AddConsBoundFByXandY();

	// run the model
	FlowF* Run();

	// set initial solution
	FlowF* SetInitSol();
	int TraveseInitSol(int* sol_x, int* sol_f, _tree* tree, int v, int vp);


	// print 
	FlowF* PrintSol();
	FlowF* PrintModel();


	void SaveOpt();


	// set user cuts active
	void SetUserCutsActive(bool active) {
		user_cuts_active = active;
	}


	// print the solution



};

