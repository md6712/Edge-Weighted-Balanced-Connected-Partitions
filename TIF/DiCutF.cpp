#include "DiCutF.h"

// constructor 
DiCutF::DiCutF(_g* instance, bool redirect = false) :Cplex(instance) {
	printCuts = false;
	printCycles = false;

	if (!redirect) {
		// Define Variables 
		DefVar();

		// Add Objective function 
		if (second_obj)
			AddObj2();
		else
			AddObj();

		// Add Constraints
		AddCons();
	}
}

// destructor
DiCutF::~DiCutF() {

}

// run the model
DiCutF* DiCutF::Run() {
	//cplex.use(callback(env, this)); // TODO: revert this
	//	cplex.use(callbackuser(env, this));
	Cplex::Run();

	return this;
}

DiCutF* DiCutF::SetPrintCuts(bool printCuts) {
	Cplex::SetPrintCuts(printCuts);
	return this;
}

DiCutF* DiCutF::SetPrintCycles(bool printCycles) {
	Cplex::SetPrintCycles(printCycles);
	return this;
}

DiCutF* DiCutF::PrintModel() {
	cplex.exportModel("ModelDiCut.lp");
	return this;
}

// add obj 
void DiCutF::AddObj() {
	IloExpr exprObj(env);
	// minimize the weight of the first tree
	for (int a = 0; a < instance->num_arcs; a++) {
		exprObj += instance->arcs[a][2] * x[a][0];
	}
	model.add(IloMinimize(env, exprObj));
	exprObj.end();
}

void DiCutF::SaveOpt() {
}

void DiCutF::DefVar() {
	DefVarX();
}

void DiCutF::DefVarX() {
	x = NumVarMatrix(env, instance->num_arcs);
	for (int i = 0; i < instance->num_arcs; i++) {
		x[i] = IloNumVarArray(env, instance->num_trees, 0, 1, ILOBOOL);
	}
}

void DiCutF::AddObj2() {
	IloExpr exprObj(env);
	// minimize the weight of the first tree
	for (int e = 0; e < instance->num_arcs; e++) {
		for (int i = 0; i < instance->num_trees; i++) {
			exprObj += instance->arcs[e][2] * x[e][i];
		}
	}
	model.add(IloMinimize(env, exprObj));
	exprObj.end();
}


void DiCutF::AddCons() {
	AddConsOrderTrees();
	AddConsPicknkArcs();
	AddConsEachToTree();
	AddConsArcsOfVertex();
	
	//AddConsBoundByUB();
}


// AddConsOrderTrees
void DiCutF::AddConsOrderTrees() {
	// add constraints to ensure that the trees are ordered	
	for (int i = 0; i < instance->num_trees - 1; i++) {
		IloExpr consi(env);
		IloExpr consip(env);
		for (int a = 0; a < instance->num_arcs; a++) {
			consi += this->x[a][i] * instance->arcs[a][2];
			consip += this->x[a][i + 1] * instance->arcs[a][2];
		}
		sprintf_s(name, "OrderTrees(%d)", i);
		model.add(consi >= consip).setName(name);
	}
}

// AddConsPicknkArcs
void DiCutF::AddConsPicknkArcs() {
	IloExpr cons(env);
	for (int a = 0; a < instance->num_arcs; a++) {
		for (int i = 0; i < instance->num_trees; i++) {
			cons += x[a][i];
		}		
	}
	sprintf_s(name, "PicknkArcs");
	model.add(cons == instance->num_vertices -  instance->num_trees).setName(name);
	cons.end();
}

// AddConsEachToTree
void DiCutF::AddConsEachToTree() {
	for (int v = 0; v < instance->num_vertices; v++) {
		IloExpr cons(env);
		for (int a = 0; a < instance->num_arcs; a++) {
			// if vertex v is the sink of arc a
			if (instance->arcs[a][1] == v) {
				for (int i = 0; i < instance->num_trees; i++) {
					cons += this->x[a][i];
				}
			}			
		}		
		sprintf_s(name, "EachArcToTree(%d)", v);
		model.add(cons <= 1).setName(name);
	}
}

// AddConsArcsOfVertex
void DiCutF::AddConsArcsOfVertex() {
	// for each vertex
	for (int v = 0; v < instance->num_vertices; v++) {
		// for each source of v 
		for (int b = 0; b < instance->num_arcs; b++) {
			if (instance->arcs[b][0] == v) {
				// for each tree 
				for (int i = 0; i < instance->num_trees; i++) {
					IloExpr cons(env);

					// for each arc a entering v 
					for (int a = 0; a < instance->num_arcs; a++) {
						if (instance->arcs[a][1] == v) {							
							cons += x[a][i];							
						}
					}

					for (int j = 0; j < instance->num_trees; j++) {
						if (i != j) {
							cons += x[b][j];							
						}
					}
					sprintf_s(name, "ArcsOfVertex(%d,%d,%d,%d)", v, b, i);
					model.add(cons <= 1).setName(name);					
				}
			}
		}
	}
}



