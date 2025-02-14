#include "SetCoverF.h"

// constructor
SetCoverF::SetCoverF(_g* g, bool redirect, bool linear) :Cplex(g) {
	printCuts = false;
	printCycles = false;
	if (linear) {
		SetLinear();
	}
	else {
		SetInteger();
	}

	if (!redirect) {
		// Define Variables 
		DefVar();

		// Add Objective function 
		AddObj();

		// Add Constraints
		AddCons();
	}
}

// destructor
SetCoverF::~SetCoverF() {

}

// run
SetCoverF* SetCoverF::Run() {

	Cplex::Run();

	instance->_opt = opt;
	instance->_gap = gap;

	//SaveOpt();
	//instance->PrintOptEdges();

	return this;
}

SetCoverF* SetCoverF::PrintModel() {
	cplex.exportModel("ModelSetCoverF.lp");
	return this;
}

void SetCoverF::SaveOpt() {
	//this->instance->num_opt_edges = 0;
	//for (int i = 0; i < this->instance->num_trees; i++) {
	//	this->instance->start_edge_tree[i] = this->instance->num_opt_edges;
	//	for (int e = 0; e < this->instance->num_edges; e++) {
	//		if (this->cplex.isExtracted(this->x[e][i])) {
	//			double value = this->cplex.getValue(this->x[e][i]);
	//			if (value > 0.5) {
	//				this->instance->opt_edges[this->instance->num_opt_edges++] = e;
	//			}
	//		}
	//	}
	//}
}

// define all variables
void SetCoverF::DefVar() {
	DefVarX();
	DefVarZ();
}

// define variable x
void SetCoverF::DefVarX() {

	// read instance
	_g* g = instance;

	// number of generated trees
	int num_trees = g->trees.size();

	// define the variable x
	x = IloNumVarArray(env, num_trees, 0, 1, integer ? ILOINT : ILOFLOAT);
	for (int T = 0; T < num_trees; T++) {
		// name 
		sprintf(name, "x_%d", T);
		x[T].setName(name);
	}
}

// define variable z
void SetCoverF::DefVarZ() {
	z = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
	z.setName("z");
}


// define the objective function
void SetCoverF::AddObj() {
	IloExpr expr(env);	
	model.add(IloMinimize(env, z));
}

// define the constraints
void SetCoverF::AddCons() {
	SetCover();
	SetZ();
	KnapSack();
}

// set cover constraint
void SetCoverF::SetCover() {
	// read instance
	_g* g = instance;

	// number of generated trees
	int num_trees = g->trees.size();

	// set cover constraint
	for (int v = 0; v < g->num_vertices; v++) {
		IloExpr expr(env);
		for (int T = 0; T < num_trees; T++) {
			for (int i = 0; i < g->trees[T]->num_vertices; i++) {
				if (g->trees[T]->vertices[i] == v) {
					expr += x[T];
					break;
				}
			}
		}
		model.add(expr == 1);
	}
}


// set z constraint
void SetCoverF::SetZ() {
	// read instance
	_g* g = instance;

	// number of generated trees
	int num_trees = g->trees.size();

	for (int v = 0; v < g->num_vertices; v++) {
		IloExpr expr(env);
		for (int T = 0; T < num_trees; T++) {
			for (int i = 0; i < g->trees[T]->num_vertices; i++) {
				if (g->trees[T]->vertices[i] == v) {
					expr += x[T]*g->trees[T]->weight;
					break;
				}
			}
		}
		model.add(expr - z <= 0);
	}
}

// KnapSack constraint
void SetCoverF::KnapSack() {
	// read instance
	_g* g = instance;

	// number of generated trees
	int num_trees = g->trees.size();

	// knapsack constraint
	IloExpr expr(env);
	for (int T = 0; T < num_trees; T++) {
		expr += x[T];
	}
	model.add(expr <= g->num_trees);
}

// force the solution
SetCoverF* SetCoverF::ForceSol() {
	// read instance
	//model.add(x[34] == 1);
	model.add(x[110] == 1);	
	return this;
}

// set integer
SetCoverF* SetCoverF::SetInteger() {
	integer = true;
	return this;
}

// set linear
SetCoverF* SetCoverF::SetLinear() {
	integer = false;
	return this;
}

