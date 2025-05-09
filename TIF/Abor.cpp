#include "Abor.h"

#include <lemon/edmonds_karp.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>

ILOINCUMBENTCALLBACK2(EarlyAbortCallback, double, threshold, Abor*, abor) {
	if (getObjValue() > threshold) {

		abor->ResetOpt();

		abor->solution_found = true; // solution found
		// save the optimal solution
		for (int a = 0; a < abor->instance->num_arcs; a++) {
			abor->x_value[a] = getValue(abor->x[a]);
		}

		for (int v = 0; v < abor->instance->num_vertices; v++) {			
			abor->y_value[v] = getValue(abor->y[v]);			
		}

		// optimal value 
		abor->opt->value = getObjValue();

		// loop over all vertices
		for (int v = 0; v < abor->instance->num_vertices; v++) {
			if (abor->y_value[v] > 0) {
				addbin(abor->opt->tree->bin_vertices, v);
			}
		}

		// loop over all arcs
		for (int a = 0; a < abor->instance->num_arcs; a++) {
			if (abor->x_value[a] > 0) {
				abor->opt->tree->weight += abor->instance->arcs[a][2];
			}
		}

		abor->opt->tree->print_vertices(abor->instance);

		abort();  // Stop solving early
	}
}

ILOUSERCUTCALLBACK1(callbackuser3, Abor*, abor) {


	if (abor->printCuts) {
		printf_s("user cut\n");
	}

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

	// print x values
	if (abor->printCuts) {
		printf_s("x: ");
		for (int a = 0; a < abor->instance->num_arcs; a++) {
			double value = getValue(abor->x[a]);
			if (value > 0) {
				int u = abor->instance->arcs[a][0];
				int v = abor->instance->arcs[a][1];
				printf("(%d -- %2.2lf --> %d) ", u, value, v);
			}
		}
		printf_s("\n");
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

	//print y values
	if (abor->printCuts) {
		printf_s("y: ");
		for (int v = 0; v < abor->instance->num_vertices; v++) {
			printf_s("%3.2f ", abor->y_value[v]);
		}
		printf_s("\n");
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

			// print the flow value
			if (abor->printCuts) {
				printf_s("max flow: %3.2f\n", max_flow);
			}

			// print y value
			if (abor->printCuts) {
				printf_s("y[%d]: %3.2f\n", v, abor->y_value[v]);
			}

			if (max_flow + 0.001 < abor->y_value[v]) {
				
				// map the min cut
				ho.minCutMap(S);

				// print the min cut
				if (abor->printCuts) {
					printf_s("S: ");
					for (ListDigraph::NodeIt n(abor->instance->dg); n != INVALID; ++n) {
						if (S[n]) {
							printf_s("%d ", abor->instance->dg.id(n));
						}
					}
					printf_s("\n");
				}

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
									if (abor->printCuts) {
										int u = abor->instance->arcs[e][0];
										int v = abor->instance->arcs[e][1];
										printf_s("x[%d,%d] +", u,v);
									}								
								}
							}
						}
					}
				}

				// add the cut to the model
				// name
				sprintf(abor->name, "cut(%d)", v);
				// add the cut to the model
				add(cons >= abor->y[v]).setName(abor->name);


				if (abor->printCuts) {
					printf_s(">= y[%d]\n", v);
				}
			}
		}
	}
}

ILOLAZYCONSTRAINTCALLBACK1(callbacklazy2, Abor*, abor) {	

	if (abor->printCuts) {
		printf_s("lazy cut\n");
	}

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

	// print x values
	if (abor->printCuts) {
		printf_s("x: ");
		for (int a = 0; a < abor->instance->num_arcs; a++) {
			double value = getValue(abor->x[a]);
			if (value > 0) {
				int u = abor->instance->arcs[a][0];
				int v = abor->instance->arcs[a][1];
				printf("(%d -> %d) ", u, v);
			}
		}
		printf_s("\n");
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

	//print y values
	if (abor->printCuts) {
		printf_s("y: ");
		for (int v = 0; v < abor->instance->num_vertices; v++) {
			printf_s("%3.2f ", abor->y_value[v]);
		}
		printf_s("\n");
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

			// print the flow value
			if (abor->printCuts) {
				printf_s("max flow: %3.2f\n", max_flow);
			}

			// print y value
			if (abor->printCuts) {
				printf_s("y[%d]: %3.2f\n", v, abor->y_value[v]);
			}

			if (max_flow + 0.001 < abor->y_value[v]) {

				// map the min cut
				ho.minCutMap(S);

				// print the min cut
				if (abor->printCuts) {
					printf_s("S: ");
					for (ListDigraph::NodeIt n(abor->instance->dg); n != INVALID; ++n) {
						if (S[n]) {
							printf_s("%d ", abor->instance->dg.id(n));
						}
					}
					printf_s("\n");
				}


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
									if (abor->printCuts) {
										printf_s("x[%d] +", e);
									}
								}
							}
						}
					}
				}

				// add the cut to the model
				// name
				sprintf(abor->name, "cut(%d)", v);
				// add the cut to the model
				add(cons >= abor->y[v]).setName(abor->name);

				if (abor->printCuts) {
					printf_s(">= y[%d]\n", v);
				}
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

	// define opt solution
	opt = new AborSol();
	opt->tree = new _small_tree();	

	// define coef
	upper_bound_weight = IloNum(0);
	vertex_prizes = IloNumArray(env, this->instance->num_vertices);
	costs = IloNumArray(env, this->instance->num_arcs);

	// initialize range array
	branch_constraints = IloRangeArray(env);

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

	// free range array
	branch_constraints.end();	

	// free variables
	delete opt->tree;
	delete opt;

	// free variables
	delete[] x_value;
	delete[] y_value;	
}

AborSol* Abor::Run(bool printCuts = false) {
	this->printCuts = printCuts;
	printCycles = false;

	solution_found = false; // solution not found
	ResetOpt(); // reset the optimal solution

	// set the user callback
	cplex.use(callbackuser3(env, this));

	// set the lazy constraint callback	
	cplex.use(callbacklazy2(env, this));

	// set the early abort callback
	//cplex.use(EarlyAbortCallback(env, 0.001, this));

	// run the model
	Cplex::Run();

	// check if solution has not found 
	if (!solution_found) {
		// check if cplex is feasible or optimal	
		if (cplex.getStatus() == IloAlgorithm::Feasible || cplex.getStatus() == IloAlgorithm::Optimal) {
			// save the optimal solution			
			SaveOpt();
		}		
	}


	return opt;
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

	// objective constraint -> obj >= 0
	objective_constraint = IloRange(obj >= 0);

	// add objective
	model.add(objective);

	// add objective constraint
	model.add(objective_constraint);

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

Abor* Abor::AddConstraintsBranching(BP_node* node) {

	for (int i = 0; i < branch_constraints.getSize(); ++i) {
		model.remove(branch_constraints[i]);
		branch_constraints[i].end(); // optional cleanup
	}
	branch_constraints.clear();

	for (int r = 0; r < node->lvl; ++r) {
		int u = node->branch[r].u;
		int v = node->branch[r].v;
		CG_branch_rule rule = node->branch[r].rule;

		IloRange cons;
		if (rule == CG_branch_rule::apart) {
			cons = IloRange(env, -IloInfinity, y[u] + y[v], 1);
		}
		else if (rule == CG_branch_rule::together) {
			cons = IloRange(env, 0, y[u] - y[v], 0);
		}
		else {
			std::cerr << "Error: unknown rule\n";
			continue;
		}

		char name[64];
		sprintf(name, "Cons_Branching(%d,%d)", u, v);
		cons.setName(name);
		branch_constraints.add(cons);  // just collect here
	}

	// Add all at once
	model.add(branch_constraints);

	return this;

}

Abor* Abor::Init(_pcst* pcst, double fixed_cost){
	// set the vertex prizes; these prices are not used directly in the model, but used to dominate some cases
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
		obj += costs[a] * x[a] ;
	}
	// set the objective function
	objective.setExpr(obj - fixed_cost);
	objective_constraint.setExpr(obj - fixed_cost >= 0);
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

void Abor::ResetOpt() {
	// reset vertices in optimal solution
	memset(opt->tree->bin_vertices, 0, sizeof(uint32_t) * SIZE_OF_VERTICES_BINARY);
	opt->tree->weight = 0;
	opt->value = -INFINITY;
}

void Abor::SaveOpt() {
	// save the optimal solution
	for (int a = 0; a < this->instance->num_arcs; a++) {
		if (cplex.isExtracted(x[a])) {
			x_value[a] = cplex.getValue(x[a]);
		}
		else
		{
			x_value[a] = 0;
		}
	}

	for (int v = 0; v < this->instance->num_vertices; v++) {
		if (cplex.isExtracted(y[v])) {
			y_value[v] = cplex.getValue(y[v]);
		}
		else
		{
			y_value[v] = 0;
		}
	}

	// optimal value 
	opt->value = cplex.getObjValue();

	// loop over all vertices
	for (int v = 0; v < this->instance->num_vertices; v++) {
		if (y_value[v] > 0) {
			addbin(opt->tree->bin_vertices, v);
		}
	}

	// loop over all arcs
	for (int a = 0; a < this->instance->num_arcs; a++) {
		if (x_value[a] > 0) {
			opt->tree->weight += this->instance->arcs[a][2];
		}
	}
}

// print the solution
Abor* Abor::PrintSol() {
	// print the solution

	if (cplex.getStatus() != IloAlgorithm::Optimal && cplex.getStatus() != IloAlgorithm::Feasible) {
		std::cerr << "No solution available to print!" << std::endl;
		return this;
	}

	printf("Abor Solution:\n");

	for (int a = 0; a < this->instance->num_arcs; a++) {
		
		if (cplex.isExtracted(x[a])) {
			double value = cplex.getValue(x[a]);
			if (value > 0.001) {
				int v = this->instance->arcs[a][0];
				int u = this->instance->arcs[a][1];
				printf("x(%d,%d) = %4.2f\n", v, u, value);
			}
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


