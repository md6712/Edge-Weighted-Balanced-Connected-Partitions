#include "Abor.h"

#include <lemon/edmonds_karp.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>

ILOUSERCUTCALLBACK1(callbackuser3, Abor*, abor) {

	// read x values
	for (int a = 0; a < abor->instance->num_arcs; a++) {
		if (abor->cplex.isExtracted(abor->x[a])) {
			double value = getValue(abor->x[a]);
			abor->x_value[a] = value;
		}
		else {
			abor->x_value[a] = 0;
		}
	}

	// read y values
	for (int v = 0; v < abor->instance->num_vertices; v++) {
		if (abor->cplex.isExtracted(abor->y[v])) {
			double value = getValue(abor->y[v]);
			abor->y_value[v] = value;
		}
		else {
			abor->y_value[v] = 0;
		}
	}

	// create a capacity map
	ListDigraph::ArcMap<double> capacity(abor->instance->dg);

	// create a cut map to store the min cut
	ListDigraph::NodeMap<bool> S(abor->instance->dg);

	// set capacity for all arcs (u,v)
	for (ListDigraph::ArcIt ar(abor->instance->dg); ar != INVALID; ++ar) {
		// find associated x value
		for (int a = 0; a < abor->instance->num_arcs; a++) {
			if (abor->instance->dg.id(abor->instance->dg.source(ar)) == abor->instance->arcs[a][0]
				&& abor->instance->dg.id(abor->instance->dg.target(ar)) == abor->instance->arcs[a][1]) {
				capacity[ar] = abor->x_value[a];
			}
		}
	}

	// for each vertex in the graph with y value > 0
	for (int v = 0; v < abor->instance->num_vertices; v++) {
		if (abor->y_value[v] > 0) {
			
			EdmondsKarp<ListDigraph, ListDigraph::ArcMap<double>> ho(
				abor->instance->dg,
				capacity,
				abor->instance->dg.nodeFromId(abor->instance->num_vertices),
				abor->instance->dg.nodeFromId(v)
			);

			ho.run();

			// max flow value 
			double max_flow = ho.flowValue();

			if (max_flow < abor->y_value[v]) {
				
				// map the min cut
				ho.minCutMap(S);

				// create an expression for the cut
				IloExpr cons(abor->env);

				// for each arc in the min cut
				for (ListDigraph::ArcIt a(abor->instance->dg); a != INVALID; ++a) {
					auto u = abor->instance->dg.source(a);
					auto v = abor->instance->dg.target(a);

					if (S[u] && !S[v]) {

						// add the arc to the cut
						for (int e = 0; e < abor->instance->num_arcs; e++) {
							if (abor->instance->arcs[e][0] == abor->instance->dg.id(u) && abor->instance->arcs[e][1] == abor->instance->dg.id(v)) {
								if (abor->cplex.isExtracted(abor->x[e])) {
									cons += abor->x[e];
								}
							}
						}
					}
				}

				// add the cut to the model
				// name
				sprintf(abor->name, "cut(%d)", v);
				// add the cut to the model
				add(cons <= max_flow).setName(abor->name);
			}
		}
	}
}

// constructor
Abor::Abor(_g* g, bool redirect, bool linear) :Cplex(g) {
	printCuts = false;
	printCycles = false;

	// define x_values	
	x_value = new double[this->instance->num_arcs];
	y_value = new double[this->instance->num_vertices];

	// define coef
	upper_bound_weight = IloNum(0);
	vertex_prizes = IloNumArray(env, this->instance->num_vertices);
	costs = IloNumArray(env, this->instance->num_arcs);

	if (linear) {
		SetLinear();
	}
	else {
		SetInteger();
	}
	if (!redirect) {
		// Define Variables 
		DefVars();
		// Add Objective function
		AddObj();
		// Add Constraints
		AddCons();	
	}
}


// destructor
Abor::~Abor() {

	// free coef
	vertex_prizes.end();
	costs.end();

	// free variables
	delete[] x_value;
	delete[] y_value;	
}

Abor* Abor::Run() {

	// set the user callback
	cplex.use(callbackuser3(env, this));

	Cplex::Run();
	return this;
}

void Abor::DefVars() {
	DefVarX();
	DefVarY();	
}

void Abor::DefVarX() {
	// define x variables
	x = IloNumVarArray(env, this->instance->num_arcs);
	for (int a = 0; a < this->instance->num_arcs; a++) {
		sprintf_s(name, "x(%d,%d)", instance->arcs[a][0], instance->arcs[a][1]);
		x[a] = IloNumVar(env, 0, 1, ILOINT, name);
		model.add(x[a]);
	}
}

void Abor::DefVarY() {
	// define y variables
	y = IloNumVarArray(env, this->instance->num_vertices);
	for (int v = 0; v < this->instance->num_vertices; v++) {
		sprintf_s(name, "y(%d)", v);
		y[v] = IloNumVar(env, 0, 1, ILOINT, name);
		model.add(y[v]);
	}
}

void Abor::AddObj() {
	IloExpr obj(env);
	for (int a = 0; a < this->instance->num_arcs; a++) {
		obj += costs[a] * x[a];
	}

	// set the objective function
	objective = IloObjective(env, obj, IloObjective::Maximize, "obj");

	model.add(objective);

	obj.end();
}

void Abor::AddCons() {
	
	AddConsOneIncomingArc();
	AddConsUpperBound();
	AddConsExactlyOneArcFromRoot();
	//AddConsAtLeastOneOutgoingArcTerminal();
	AddConsLoopsOfSizeTwo();	
}

void Abor::AddConsOneIncomingArc() {
	for (int v = 0; v < this->instance->num_vertices; v++) {
		IloExpr cons(env);
		for (int a = 0; a < this->instance->num_arcs; a++) {
			if (this->instance->arcs[a][1] == v) {
				cons += x[a];
			}
		}
		// name
		sprintf(name, "Cons_OneOutgoingArc(%d)", v);
		model.add(cons == y[v]).setName(name);
		cons.end();
	}
}

void Abor::AddConsUpperBound() {
	IloExpr cons(env);
	for (int a = 0; a < this->instance->num_arcs; a++) {
		cons += x[a] * this->instance->arcs[a][2];
	}
	// name
	sprintf(name, "Cons_UpperBound");
	// set the upper bound constraint
	upper_bound = IloRange(env, cons);
	upper_bound.setName(name);
	
	// add to model
	model.add(upper_bound);

	// release the expression
	cons.end();
}

void Abor::AddConsExactlyOneArcFromRoot() {
	IloExpr cons(env);
	for (int a = 0; a < this->instance->num_arcs; a++) {
		if (this->instance->arcs[a][0] == this->instance->num_vertices) {
			cons += x[a];
		}
	}
	// name
	sprintf(name, "Cons_ExactlyOneArcFromRoot");
	model.add(cons == 1).setName(name);
	cons.end();
}

void Abor::AddConsAtLeastOneOutgoingArcTerminal() {
	for (int v = 0; v < this->instance->num_vertices; v++) {
		if (vertex_prizes[v] <= 0) {
			IloExpr cons_outgoing(env);
			for (int a = 0; a < this->instance->num_arcs; a++) {
				if (this->instance->arcs[a][0] == v) {
					cons_outgoing += x[a];
				} 
			}

			IloExpr cons_incoming(env);
			for (int a = 0; a < this->instance->num_arcs; a++) {
				if (this->instance->arcs[a][1] == v) {
					cons_incoming += x[a];
				}
			}
			// name
			sprintf(name, "Cons_AtLeastOneOutgoingArcTerminal(%d)", v);
			model.add(cons_outgoing >= cons_incoming).setName(name);
			cons_outgoing.end();
			cons_incoming.end();			
		}
	}
}

void Abor::AddConsLoopsOfSizeTwo() {

	// [1,num_edges] contains arcs in one direction
	// [num_edges+1,2*num_edges] contains arcs in the other direction


	for (int a = 0; a < this->instance->num_edges; a++) {		
		// find the complement 
		int a1 = a + this->instance->num_edges;

		// check the two vertices
		int u = this->instance->arcs[a][0];
		int v = this->instance->arcs[a][1];

		// one of a and a1 is selected if y[u] = 1 or y[v] = 1, if either is zero the a and a1 are not selected				
		// name
		sprintf(name, "Cons_LoopsOfSizeTwo(%d)1", a);
		model.add(y[u] >= x[a] + x[a1]).setName(name);
		sprintf(name, "Cons_LoopsOfSizeTwo(%d)2", a);
		model.add(y[v] >= x[a] + x[a1]).setName(name);
	}
}

Abor* Abor::Init(_pcst* pcst){
// set the vertex prizes
	for (int v = 0; v < this->instance->num_vertices; v++) {
		vertex_prizes[v] = pcst->vertex_prize[v];
	}
	// set the costs
	costs = IloNumArray(env, this->instance->num_arcs);
	for (int a = 0; a < this->instance->num_arcs; a++) {
		if (pcst->arc_active[a])
			costs[a] = -pcst->arc_cost[a];
		else
			costs[a] = -PCST_LARGE;
	}

	// set the upper bound
	upper_bound_weight = pcst->UB;


	// update objective function with new costs
	IloExpr obj(env);
	for (int a = 0; a < this->instance->num_arcs; a++) {
		obj += costs[a] * x[a];
	}
	// set the objective function
	objective.setExpr(obj);	
	obj.end();


	// update the right hand side of the upper bound constraint
	upper_bound.setUB(upper_bound_weight);

	// update variable bounds for x and y 
	for (int a = 0; a < this->instance->num_arcs; a++) {
		if (pcst->arc_active[a]) {
			x[a].setUB(1);
			x[a].setLB(0);
		}
		else {
			x[a].setUB(0);
			x[a].setLB(0);
		}
	}

	for (int v = 0; v < this->instance->num_vertices; v++) {
		if (pcst->vertex_aborescence_active[v]) {			
			y[v].setUB(1);
			if (pcst->roots[v])
				y[v].setLB(1);
			else
				y[v].setLB(0);
		}
		else {
			y[v].setUB(0);
			y[v].setLB(0);
		}
	}	 

	return this;
}

// print the solution
Abor* Abor::PrintSol() {
	// print the solution
	for (int a = 0; a < this->instance->num_arcs; a++) {
		if (x[a].getUB() > 0) {
			printf("x[%d] = %4.2f\n", a, x[a].getUB());
		}
	}
	return this;
}


// print the model
Abor* Abor::PrintModel() {
	// print the model
	cplex.exportModel("modelABOR.lp");
	return this;
}

