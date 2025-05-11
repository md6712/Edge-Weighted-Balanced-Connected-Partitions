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

	PrintSol();

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

// print sol
SetCoverF* SetCoverF::PrintSol() {
	// read instance
	_g* g = instance;
	// number of generated trees
	int num_trees = g->trees.size();
	// print the solution
	for (int T = 0; T < num_trees; T++) {
		double value = cplex.getValue(x[T]);
		if (value > 0.001) {
			printf("X[%3d] = %1.2f\t", T, value);
			printf_s("Tree %d: ", T);
			g->trees[T]->PrintVerticesWeight();
		}
	}
	return this;
}


// set the initial solution
SetCoverF* SetCoverF::SetInitSol() {

	// arrays to store the values of the variables
	IloNumVarArray vars(env);
	IloNumArray vals(env);

	// arrays to store solution values 
	int* sol_x = new int[instance->trees.size()];

	// force sol_x to be zero
	memset(sol_x, 0, sizeof(int) * instance->trees.size());

	int sol_z = 0;

	// for each tree in trees_ub, find the associated tree in the vector trees for that variable set the value of the variable to 1
	for (int i = 0; i < instance->trees_ub.size(); i++) {
		_tree* tree = instance->trees_ub[i];
		// loop over all trees and see if the vertices in two tree are the same 
		for (int j = 0; j < instance->trees.size(); j++) {
			_tree* tree2 = instance->trees[j];
			// compare the two vectors tree->bin_vertices and tree2->bin_vertices
			int k;
			bool result = false;
			if (tree->num_vertices != tree2->num_vertices) {
				continue;
			}

			for (k = 0; k < tree->num_vertices; k++) {
				if (tree->vertices[k] != tree2->vertices[k]) {
					break;
				}
			}

			// if the two trees are the same, set the value of the variable to 1
			if (k == tree->num_vertices) {
				result = true;
			}

			if (result) {
				// set the value of the variable to 1
				sol_x[j] = 1;
				sol_z = max(tree2->weight, sol_z);
				break;
			}
		}
	}

	// add the values to the array
	for (int i = 0; i < instance->num_trees; i++) {		
		vals.add(sol_x[i]);
		vars.add(x[i]);		
	}

	// add the value of z
	vals.add(sol_z);
	vars.add(z);

	// add the initial solution
	
	try {
		cplex.addMIPStart(vars, vals);
	}
	catch (IloException& e) {
		std::cerr << "MIP start rejected: " << e.getMessage() << std::endl;
	}
	vals.end();
	vars.end();

	delete sol_x;
	return this;
}