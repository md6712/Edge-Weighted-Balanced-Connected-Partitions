#include "RepresentativeF.h"
#include "binary.h"

RepresentativeF::RepresentativeF(_g* g, bool redirect) : Cplex(g) {
	this->printCuts = false;
	this->printCycles = false;

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
RepresentativeF::~RepresentativeF() {
}

// run the model
RepresentativeF* RepresentativeF::Run() {
//	cplex.use(callback(env, this));
	Cplex::Run();

	SaveOpt();
	//instance->PrintOptEdges();

	return this;
}

RepresentativeF* RepresentativeF::PrintModel() {
	cplex.exportModel("ModelRepresentative.lp");
	return this;
}

void RepresentativeF::SaveOpt() {

}


void RepresentativeF::DefVar() {
	// Define x variables
	x = NumVarMatrix(env, instance->num_edges);
	for (int i = 0; i < instance->num_edges; i++) {
		x[i] = IloNumVarArray(env, instance->num_edges, 0, 1, ILOBOOL);
	}
}


// objective function
void RepresentativeF::AddObj() {
	
}


// addcons
void RepresentativeF::AddCons() {
	AddConsComputeZ();
	AddConsKeyEdges();
	AddConsIFeeTHENef();
	AddConsBelongToOneTree();
	AddConsTotalNumberOfEdgesSelected();
	AddConsAllCycles();
}

void RepresentativeF::AddConsComputeZ() {	
	for (int e = 0; e < instance->num_edges; e++) {
		IloExpr expr(env);
		for (int f = e; f < instance->num_edges; f++) {
			expr += instance->edges[f][2] * x[e][f];
		}
		model.add(z >= expr);
		expr.end();
	}	
}


// there are at most k = num_trees key edges
void RepresentativeF::AddConsKeyEdges() {
	IloExpr expr(env);
	for (int e = 0; e < instance->num_edges; e++) {
		expr += x[e][e];
	}
	model.add(expr <= instance->num_trees);
	expr.end();
}


// if e is a key edge then f can be in the same tree as e
void RepresentativeF::AddConsIFeeTHENef() {
	for (int e = 0; e < instance->num_edges; e++) {
		for (int f = e; f < instance->num_edges; f++) {
			model.add(x[e][f] <= x[e][e]);
		}
	}
}


// each edge e belongs to at most one tree
void RepresentativeF::AddConsBelongToOneTree() {
	for (int e = 0; e < instance->num_edges; e++) {
		IloExpr expr(env);
		for (int f = 0; f <= e; f++) {
			expr += x[f][e];
		}
		model.add(expr <= 1);
		expr.end();
	}
}

// the total number of edges selected is equal to the number of vertices minus number of trees
void RepresentativeF::AddConsTotalNumberOfEdgesSelected() {
	IloExpr expr(env);
	for (int e = 0; e < instance->num_edges; e++) {
		for (int f = e; f < instance->num_edges; f++) {
			expr += x[e][f];
		}		
	}
	model.add(expr == instance->num_vertices - instance->num_trees);
	expr.end();
}


// all cycles
void RepresentativeF::AddConsAllCycles() {
	// for each edge, for each subset of vertices that is not empty 
	for (int e = 0; e < instance->num_edges; e++) {
		for (int s = 1; s < pow(2, instance->num_vertices); s++) {
			if (instance->nSubsets[s] > 2 && instance->nSubsets[s] < pow(2, instance->num_vertices)) {
				IloExpr cons(env);
				int count = 0;
				for (int f = e; f < instance->num_edges; f++) {
					if (checkbin(instance->subsets[s], instance->edges[f][0]) && checkbin(instance->subsets[s], instance->edges[f][1])) {
						cons += this->x[e][f];
						count++;
					}
				}
				if (count > instance->nSubsets[s] - 1) {
					sprintf_s(name, "ForEachSubset(%d,%d)", e, s);
					model.add(cons <= x[e][e] * (instance->nSubsets[s] - 1)).setName(name);
				}
			}
		}
	}
}