#include "FlowF.h"
#include <lemon/edmonds_karp.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>

ILOUSERCUTCALLBACK1(callbackuserFlow, FlowF*, flowF) {
	// counter
	flowF->usercallback_count++;

	// create one dimensional array to store the opt x
	double* opt_x = new double[flowF->instance->num_arcs];
	int* S = new int[flowF->instance->num_vertices];


	// store the opt x values in the array opt_x
	for (int a = 0; a < flowF->instance->num_arcs; a++) {
		if (flowF->cplex.isExtracted(flowF->x[a])) {
			double value = getValue(flowF->x[a]);
			opt_x[a] = value;
		}
	}

	// create a capacity map
	lemon::ListDigraph::ArcMap<double> capacity(flowF->instance->dg);
	
	// set capacity for arcs
	for (lemon::ListDigraph::ArcIt a(flowF->instance->dg); a != INVALID; ++a) {	
		double sum_x = 0;
		if (flowF->instance->dg.id(flowF->instance->dg.source(a)) == flowF->instance->num_vertices) {
			int v = flowF->instance->dg.id(flowF->instance->dg.target(a));
			
			for (int aa = 0; aa < flowF->instance->num_arcs; aa++) {
				if (flowF->instance->arcs[aa][0] == v) {
					sum_x += opt_x[aa];
				}
			}
			capacity[a] = sum_x;
		}

		// if v is the target
		else if (flowF->instance->dg.id(flowF->instance->dg.target(a)) == flowF->instance->num_vertices + 1) {
			capacity[a] = 1;
		}

		// otherwise
		else {
			for (int e = 0; e < flowF->instance->num_arcs; e++) {
				if (flowF->instance->arcs[e][0] == flowF->instance->dg.id(flowF->instance->dg.source(a)) && flowF->instance->arcs[e][1] == flowF->instance->dg.id(flowF->instance->dg.target(a))) {
					sum_x += opt_x[e];					
					
				}
				/*else if (flowF->instance->arcs[e][0] == flowF->instance->dg.id(flowF->instance->dg.target(a)) && flowF->instance->arcs[e][1] == flowF->instance->dg.id(flowF->instance->dg.source(a))) {
					sum_x += opt_x[e];										
				}*/
			}
			capacity[a] = sum_x;
		}


	}


	for (ListDigraph::ArcIt a(flowF->instance->dg); a != INVALID; ++a) {
		if (flowF->instance->dg.id(flowF->instance->dg.source(a)) == flowF->instance->num_vertices) {
			double org_capacity = capacity[a];
			capacity[a] = 100;
			// print a
			//if (cutF->printCuts)
				//printf_s("\nArc %d -> %d: %lf  -->", cutF->instance->dg.id(cutF->instance->dg.source(a)), cutF->instance->dg.id(cutF->instance->dg.target(a)), capacity[a]);

			EdmondsKarp<ListDigraph, ListDigraph::ArcMap<double>> ho(
				flowF->instance->dg,
				capacity,
				flowF->instance->dg.nodeFromId(flowF->instance->num_vertices),
				flowF->instance->dg.nodeFromId(flowF->instance->num_vertices + 1)
			);

			ho.run();

			//create a cut map to store the min cut
			ListDigraph::NodeMap<bool> minCut(flowF->instance->dg);

			// map the min cut
			ho.minCutMap(minCut);

			// compute S 
			int nS = 0;
			for (ListDigraph::NodeIt v(flowF->instance->dg); v != INVALID; ++v) {
				if (minCut[v]) {
					if (flowF->instance->dg.id(v) != flowF->instance->num_vertices) {
						S[nS++] = flowF->instance->dg.id(v);
					}
				}
			}

			// compute cost: for each arc that doesnt have a vertex in S,  add the cost to the cut
			double cost = nS;
			for (int e = 0; e < flowF->instance->num_arcs; e++) {
				bool SinS = false; // if source in S
				for (int j = 0; j < nS; j++) {
					int v = S[j];
					if (flowF->instance->arcs[e][0] == v) {
						SinS = true;
						break;
					}
				}
				bool DinS = false;	// if dist in S
				for (int j = 0; j < nS; j++) {
					int v = S[j];
					if (flowF->instance->arcs[e][1] == v) {
						DinS = true;
						break;
					}
				}
				if (!SinS || !DinS) { // if the edge is outside the cut
					for (int i = 0; i < flowF->instance->num_trees; i++) {
						cost += opt_x[e];
					}
				}
			}

			// print cost and cut
			//printf("\nCost: %lf\n", cost);

			// if cost < n - k - 1, add the cut

			if (nS >= 2 && cost < flowF->instance->num_vertices - flowF->instance->num_trees + 1) {
				IloExpr cons(flowF->env);

				if (flowF->printCuts)
					printf_s("\n cut added: ");

				for (int e = 0; e < flowF->instance->num_arcs; e++) {
					bool SinS = false; // if source in S
					for (int j = 0; j < nS; j++) {
						int v = S[j];
						if (flowF->instance->arcs[e][0] == v) {
							SinS = true;
							break;
						}
					}
					bool DinS = false;	// if dist in S
					for (int j = 0; j < nS; j++) {
						int v = S[j];
						if (flowF->instance->arcs[e][1] == v) {
							DinS = true;
							break;
						}
					}

					if (SinS && DinS) {						
						if (flowF->cplex.isExtracted(flowF->x[e])) {
							cons += flowF->x[e];
							if (flowF->printCuts)
								printf_s("x[%d] +", e);
						}
					}
				}

				add(cons <= nS - 1);
				if (flowF->printCuts)
					printf_s("<= %d", nS - 1);

				cons.end();

				// break as soon as one cut is added.
				break;
			}


			// restore the capacity
			capacity[a] = org_capacity;
		}
	}

	delete S;	
	delete opt_x;
}


FlowF::FlowF(_g* instance, bool redirect, bool linear) : Cplex(instance) {
	this->printCycles = false;
	this->printCuts = false;

	if (linear) {
		SetLinear();
	}
	else {
		SetInteger();
	}

	if (!redirect) {
		DefVar();
		AddObj();
		AddCons();
		AddPriorityY();
	}
}

FlowF::~FlowF() {
}

// run
FlowF* FlowF::Run() {

	cplex.use(callbackuserFlow(env, this));

	Cplex::Run();
	SaveOpt();
	PrintSol();
	instance->_opt = cplex.getValue(eta);

	// print cuts added
	printf_s("\nCuts added: %d", usercallback_count);

	return this;
}

// save the optimal solution
void FlowF::SaveOpt() {
	instance->num_opt_edges = 0;	
}

// print the solution
FlowF* FlowF::PrintSol() {
	printf("Objective value: %f\n", cplex.getObjValue());
	for (int i = 0; i < instance->num_arcs; i++) {
		if (cplex.isExtracted(x[i]))
		if (cplex.getValue(x[i]) > 0.00001) {
			printf("x[%d] = %f \n", i, cplex.getValue(x[i]));
		}
	}
	printf("\n");
	for (int i = 0; i < instance->num_vertices; i++) {
		if (cplex.isExtracted(y[i]))
		if (cplex.getValue(y[i]) > 0.00001) {
			printf("y[%d] = %f \n", i, cplex.getValue(y[i]));
		}
	}
	printf("\n");
	for (int i = 0; i < instance->num_arcs; i++) {
		if (cplex.isExtracted(f[i]))
		if (cplex.getValue(f[i]) > 0.00001) {
			printf("f[%d] = %f \n", i, cplex.getValue(f[i]));
		}
	}
	printf("\n");
	for (int i = 0; i < instance->num_vertices; i++) {
		if (cplex.isExtracted(f0[i]))
		if (cplex.getValue(f0[i]) > 0.00001) {
			printf("f0[%d] = %f \n", i, cplex.getValue(f0[i]));
		}
	}
	printf("\n");
	printf("eta = %f\n", cplex.getValue(eta));
	printf("\n");

	// print all theta values for each vertex and for each tree
	for (int v = 0; v < instance->num_vertices; v++) {
		printf("theta[%d]", v);
		for (int i = 0; i < instance->num_trees; i++) {
			if (cplex.isExtracted(theta[v][i]))
				//if (cplex.getValue(theta[v][i]) > 0.00001) {
				printf(" %2.2lf \t", cplex.getValue(theta[v][i]));
			//}		}
		}
		printf("\n");
	}

	return this;
}

// print the model
FlowF* FlowF::PrintModel() {
	// export the model
	cplex.exportModel("FlowF.lp");
	return this;
}

// define variables
void FlowF::DefVar() {
	DefVarX();
	DefVarY();
	DefVarF();
	DefVarF0();
	DefVarEta();

	// theta
	DefVarTheta();
}

// define x
void FlowF::DefVarX() {
	x = IloNumVarArray(env, instance->num_arcs, 0, 1, integer?ILOINT:ILOFLOAT);
	for (int i = 0; i < instance->num_arcs; i++) {
		sprintf_s(name, "x(%d)", i);
		x[i].setName(name);
	}
}

// define y
void FlowF::DefVarY() {
	y = IloNumVarArray(env, instance->num_vertices, 0, 1, integer ? ILOINT : ILOFLOAT);
	for (int i = 0; i < instance->num_vertices; i++) {
		sprintf_s(name, "y(%d)", i);
		y[i].setName(name);
	}
}

// define f
void FlowF::DefVarF() {
	f = IloNumVarArray(env, instance->num_arcs, 0, IloInfinity, ILOFLOAT);
	for (int i = 0; i < instance->num_arcs; i++) {
		sprintf_s(name, "f(%d)", i);
		f[i].setName(name);
	}
}

// define f0
void FlowF::DefVarF0() {
	// for each vertex
	f0 = IloNumVarArray(env, instance->num_vertices, 0, IloInfinity, ILOFLOAT);
	for (int i = 0; i < instance->num_vertices; i++) {
		sprintf_s(name, "f0(%d)", i);
		f0[i].setName(name);
	}
}


// define theta
void FlowF::DefVarTheta() {
	theta = NumVarMatrix(env, instance->num_vertices);
	for (int v = 0; v < instance->num_vertices; v++) {
		theta[v] = IloNumVarArray(env, instance->num_trees, 0, IloInfinity, ILOFLOAT);
		for (int i = 0; i < instance->num_trees; i++) {
			sprintf_s(name, "theta(%d,%d)", v, i);
			theta[v][i].setName(name);
		}
		
	}
}

// define eta
void FlowF::DefVarEta() {
	eta = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
	eta.setName("eta");
}



// add objective function minimize eta
void FlowF::AddObj() {
	IloExpr exprObj(env);
	exprObj += eta;
	model.add(IloMinimize(env, exprObj));
	exprObj.end();
}


// add constraints
void FlowF::AddCons() {
	AddConsComputeEta();
	AddConsSelectKTreeRepresentatives();
	AddConsSelectArcForRepresentatives();
	AddConsFlow();
	AddConsBoundFByX();
	AddConsBoundF0ByY();
	AddConsSymmetryRootVertex();


	// hande suggestion added
	AddHandeSuggestion();

	// constraints for theta
//	AddConsThetaF0();
//	AddConsThetaEta();

	//
	//AddConsCycleXY();
}

// compute eta
void FlowF::AddConsComputeEta() {	
	for (int v = 0; v < instance->num_vertices; v++) {
		model.add(eta >= f0[v]);
	}
}


// select k tree representatives
void FlowF::AddConsSelectKTreeRepresentatives() {
	IloExpr cons(env);
	for (int v = 0; v < instance->num_vertices; v++) {
		cons += y[v];
	}	
	sprintf_s(name, "AddConsSelectKTreeRepresentatives");
	model.add(cons == instance->num_trees).setName(name);
	cons.end();
}


// select arc for representatives
void FlowF::AddConsSelectArcForRepresentatives() {

	// sum of arcs entering to a vertex v + y[v] = 1
	for (int v = 0; v < instance->num_vertices; v++) {
		IloExpr cons(env);
		for (int a = 0; a < instance->num_arcs; a++) {
			if (instance->arcs[a][0] == v) {
				cons += x[a];
			}
		}
		cons += y[v];
		sprintf_s(name, "AddConsSelectArcForRepresentatives(%d)", v);
		model.add(cons == 1).setName(name);
		cons.end();
	}
}

// flow constraints
void FlowF::AddConsFlow() {
	// f0[v] + sum of flow entering v - sum of flow leaving v  = sum of weight of arcs entering v
	for (int v = 0; v < instance->num_vertices; v++) {
		IloExpr cons(env);
		cons += f0[v];
		for (int a = 0; a < instance->num_arcs; a++) {
			if (instance->arcs[a][0] == v) {
				cons += f[a];
			}
			if (instance->arcs[a][1] == v) {
				cons -= f[a];
			}
		}

		IloExpr cons2(env);
		for (int a = 0; a < instance->num_arcs; a++) {
			if (instance->arcs[a][0] == v) {
				cons2 += instance->arcs[a][2] * x[a];
			}
		}

		sprintf_s(name, "AddConsFlow(%d)", v);
		model.add(cons == cons2).setName(name);
		cons.end();
	}
}

// bound f by x times UB
void FlowF::AddConsBoundFByX() {
	for (int a = 0; a < instance->num_arcs; a++) {
		if (instance->arcs[a][0] == instance->num_vertices) {
			break;
		}
		model.add(f[a] >= x[a] * instance->arcs[a][2]);
		sprintf_s(name, "AddConsBoundFByX(%d)", a);
		model.add(f[a] <= x[a] * instance->UB).setName(name);
	}
}

// bound f0 by y times UB
void FlowF::AddConsBoundF0ByY() {
	for (int v = 0; v < instance->num_vertices; v++) {
		sprintf_s(name, "AddConsBoundF0ByY(%d)", v);
		model.add(f0[v] <= y[v] * instance->UB).setName(name);
	}
}	

// constraint for symmetry of the root vertex
void FlowF::AddConsSymmetryRootVertex() {
	// for each vertex and its adjacent outgoing arcs, if y[v], the u is larger than v in the arc (u,v)
	for (int v = 0; v < instance->num_vertices; v++) {
		for (int a = 0; a < instance->num_arcs; a++) {
			if (instance->arcs[a][0] == v) {
				if (instance->arcs[a][0] > instance->arcs[a][1])
					model.add(x[a] + y[v] <= 1);
			}
		}
	}

}

// add constraints for theta
// add constraints for theta and f0
void FlowF::AddConsThetaF0() {
	for (int v = 0; v < instance->num_vertices; v++) {
		IloExpr cons(env);
		for (int i = 0; i < instance->num_trees; i++) {
			cons += theta[v][i];
		}
		model.add(cons == f0[v]);
		cons.end();
	}
}

// add constraints for theta and eta
void FlowF::AddConsThetaEta() {
	for (int i = 0; i < instance->num_trees; i++) {
		IloExpr cons(env);
		for (int v = 0; v < instance->num_vertices; v++) {
			cons += theta[v][i];
		}
		model.add(cons <= eta);
		cons.end();
	}
}

// add priority Y
void FlowF::AddPriorityY() {
	// give branching priority to y
	for (int v = 0; v < instance->num_vertices; v++) {
		cplex.setPriority(y[v], 1);
	}
}

// AddConsCycleXY
void FlowF::AddConsCycleXY() {
	for (int a = 0; a < instance->num_arcs; a++){
		for (int a2 = 0; a2 < instance->num_arcs; a2++) {
			if (instance->arcs[a2][0] == instance->arcs[a][1] && instance->arcs[a2][1] == instance->arcs[a][0]) {
				model.add(x[a] + x[a2] <= 1);
			}
		}
	}
}

// add handle suggestion
void FlowF::AddHandeSuggestion() {
	// eta is greater than the sum of f0 over v divided by k	
	IloExpr cons(env);
	for (int v = 0; v < instance->num_vertices; v++) {
		cons += f0[v];
	}
	model.add(eta >= cons / instance->num_trees);
	cons.end();
}