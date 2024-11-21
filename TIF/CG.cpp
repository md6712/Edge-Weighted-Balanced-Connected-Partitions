#include "CG.h"
#include <ilcplex/ilocplex.h>
#include "binary.h"

// constructor
CG::CG(_g* g, bool redirect, bool linear) :Cplex(g) {
	printCuts = false;
	printCycles = false;
	if (linear) {
		SetLinear();
	}
	else {
		SetInteger();
	}

	if (!redirect) {

		// Add Objective function
		AddObj();

		// Add Constraints
		AddCons();

		// Define Variables 
		AddVar();				
	}
}

// destructor
CG::~CG() {

}


// run


CG* CG::Run() {

	Cplex::Run();

	//SaveOpt();
	//instance->PrintOptEdges();

	return this;
}


CG* CG::PrintModel() {
	cplex.exportModel("ModelCG.lp");
	return this;
}


CG* CG::PrintSol() {
	//printf("\n
	return this;
}

// add objective function
void CG::AddObj() {
	// set the objective function
	obj = IloObjective(env, 0, IloObjective::Minimize, "obj");
	model.add(obj);
}


// add constraints
void CG::AddCons() {
	// Create arrays
	Knapsack = IloRangeArray(env, 1);
	Cover = IloRangeArray(env, this->instance->num_vertices);
	Z = IloRangeArray(env, this->instance->num_vertices);

	// set bounds 
	Knapsack[0] = IloRange(env, -this->instance->num_trees, IloInfinity);
	for (int v = 0; v < this->instance->num_vertices; v++) {
		Cover[v] = IloRange(env, 1, 1);
		Z[v] = IloRange(env, 0, IloInfinity);
	}
	
	// add to model
	model.add(Knapsack);
	model.add(Cover);
	model.add(Z);
}

// add variables	
void CG::AddVar() {
	// instance g
	_g* g = instance;	

	// define the variable Z
	sprintf(name, "z");
	z = IloNumVar(env, 0, IloInfinity, ILOFLOAT, name);
	
	obj.setLinearCoef(z, 1);

	for (int v = 0; v < g->num_vertices; v++) {
					
		Z[v].setLinearCoef(z, 1);						
	}

	// define the variable X
	x = IloNumVarArray(env);

	for (int i = 0; i < g->trees.size(); i++) {
		// name 
		sprintf(name, "x_%d", i);
		IloNumVar Xvar(env, 0, IloInfinity, ILOFLOAT, name);
		
		obj.setLinearCoef(Xvar, 0);
		Knapsack[0].setLinearCoef(Xvar, -1);
		
		for (int v = 0; v < g->num_vertices; v++) {
			if (checkbin(g->trees[i]->bin_vertices,v)) {				
				Cover[v].setLinearCoef(Xvar, 1);				
			}
			else {				
				Cover[v].setLinearCoef(Xvar, 0);				
			}

			if (g->trees[i]->vertices[0] == v) {
				Z[v].setLinearCoef(Xvar, -g->trees[i]->weight);
			}
			else {
				Z[v].setLinearCoef(Xvar, 0);
			}
		}		

		x.add(Xvar);
	}

	// add variables to the model
	model.add(z);
	model.add(x);	
}


