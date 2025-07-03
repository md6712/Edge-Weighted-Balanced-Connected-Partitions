#include "FlowF.h"
#include <lemon/edmonds_karp.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>

ILOUSERCUTCALLBACK1(callbackuserFlow, FlowF*, flowF) {
	// counter

	if (getNnodes() == 0) {
		flowF->instance->_lb_root = getBestObjValue(); // LP relaxation bound at root		
	}

	if (!flowF->user_cuts_active) {
		return;
	}

	flowF->usercallback_count++;

	// create one dimensional array to store the opt x
	double* opt_x = new double[flowF->instance->num_edges*2];
	int* S = new int[flowF->instance->num_vertices];

	// initialize the array opt_x and S
	memset(opt_x, 0, sizeof(double) * flowF->instance->num_edges * 2);
	memset(S, 0, sizeof(int) * flowF->instance->num_vertices);

	// store the opt x values in the array opt_x
	for (int a = 0; a < flowF->instance->num_edges*2; a++) {
		if (flowF->cplex.isExtracted(flowF->x[a])) {
			double value = getValue(flowF->x[a]);
			opt_x[a] = value;
		}
	}

	//// print the opt x values
	//for (int a = 0; a < flowF->instance->num_edges*2; a++) {
	//	if (opt_x[a] > 0.00001) {
	//		int u = flowF->instance->arcs[a][0];
	//		int v = flowF->instance->arcs[a][1];
	//		printf_s("x[%d,%d] = %lf\n", u, v, opt_x[a]);
	//	}
	//}

	// create a capacity map
	lemon::ListDigraph::ArcMap<double> capacity(flowF->instance->dg);
	
	// set capacity for arcs
	for (lemon::ListDigraph::ArcIt a(flowF->instance->dg); a != INVALID; ++a) {	
		double sum_x = 0;
		int u = flowF->instance->dg.id(flowF->instance->dg.source(a));
		int v = flowF->instance->dg.id(flowF->instance->dg.target(a));
		if (flowF->instance->dg.id(flowF->instance->dg.source(a)) == flowF->instance->num_vertices) {	
			for (int aa = 0; aa < flowF->instance->num_edges*2; aa++) {
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
			for (int e = 0; e < flowF->instance->num_edges*2; e++) {
				if (flowF->instance->arcs[e][0] == flowF->instance->dg.id(flowF->instance->dg.source(a)) && flowF->instance->arcs[e][1] == flowF->instance->dg.id(flowF->instance->dg.target(a))) {
					sum_x += opt_x[e];										
				}
				/*else if (flowF->instance->arcs[e][0] == flowF->instance->dg.id(flowF->instance->dg.target(a)) && flowF->instance->arcs[e][1] == flowF->instance->dg.id(flowF->instance->dg.source(a))) {
					sum_x += opt_x[e];										
				}*/
			}
			capacity[a] = sum_x;
		}
		
		/*if (capacity[a] > 0.00001)			
			printf_s("Arc %d -> %d: %lf\n", u, v, capacity[a]);*/
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

			//// print S 
			//printf("S: ");
			//for (int j = 0; j < nS; j++) {
			//	int v = S[j];
			//	printf("%d ", v);
			//}
			//printf("\n");

			// compute cost: for each arc that doesnt have a vertex in S,  add the cost to the cut
			double cost = nS;
			for (int a = 0; a < flowF->instance->num_edges*2; a++) {
				bool SinS = false; // if source in S
				for (int j = 0; j < nS; j++) {
					int v = S[j];
					if (flowF->instance->arcs[a][0] == v) {
						SinS = true;
						break;
					}
				}
				bool DinS = false;	// if dist in S
				for (int j = 0; j < nS; j++) {
					int v = S[j];
					if (flowF->instance->arcs[a][1] == v) {
						DinS = true;
						break;
					}
				}
				if (!SinS || !DinS) { // if the edge is outside the cut
					
					/*int u = flowF->instance->arcs[a][0];
					int v = flowF->instance->arcs[a][1];
					if (opt_x[a] > 0.00001)
					printf("arc[%d, %d] : cost = %lf + %lf\n", u, v, cost, opt_x[a]);*/

					cost += opt_x[a];
				}
			}

			// print cost and cut
			//printf("Cost: %lf\n", cost);

			// if cost < n - k - 1, add the cut

			if (nS >= 1 && cost + 0.001 < flowF->instance->num_vertices - flowF->instance->num_trees + 1) {
				IloExpr cons(flowF->env);

				if (flowF->printCuts)
					printf_s("\n cut added: ");

				for (int e = 0; e < flowF->instance->num_edges*2; e++) {
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

				flowF->instance->n_user_cuts++;
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


FlowF::FlowF(_g* instance, bool redirect, bool linear, bool mutual_exclusion_cycles) : Cplex(instance) {
	this->printCycles = false;
	this->printCuts = false;
	this->user_cuts_active = false;
	this->mutual_exclusion_cycles = false;
	this->mutual_exclusion_cycles = mutual_exclusion_cycles;

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

	return this;
}

// save the optimal solution
void FlowF::SaveOpt() {
	instance->_opt = opt;
	instance->_gap = gap;	
	instance->n_lazy_cuts = 0;
}

int FlowF::TraveseInitSol(int* sol_x, int* sol_f, _tree* tree, int v, int vp) {

	// total flow 
	int total_flow = 0;

	// for each edge that is in the tree and includes the vertex v
	for (int i = 0; i < tree->num_edges; i++) {
		int v1 = instance->edges[tree->edges[i]][0];
		int v2 = instance->edges[tree->edges[i]][1];
		int w = instance->edges[tree->edges[i]][2];

		int u = -1;
		if (v1 == v) {
			u = v2;
		}
		else if (v2 == v) {
			u = v1;
		}

		// if u is not the previous vertex

		if (u != -1 && u != vp) {
			// find the arc associated (v, u)
			int a = -1;
			for (int j = 0; j < instance->num_edges*2; j++) {
				if (instance->arcs[j][0] == v && instance->arcs[j][1] == u) {
					a = j;
					break;
				}
			}

			// add the arc to the solution
			sol_x[a] = 1;

			// print 
			//printf_s("x[%d,%d] = 1\n", v,u);
			
			// traverse the tree
			int flow_in_childs  = TraveseInitSol(sol_x, sol_f, tree, u, v);

			sol_f[a] = flow_in_childs + w;

			// print 
			//printf_s("f[%d,%d] = %d + %d\n", v, u, flow_in_childs, w);

			// compute the flow
			total_flow += flow_in_childs + w;
		}		
	}

	return total_flow;
}

FlowF* FlowF::SetInitSol() {

	int* sol_x = new int[instance->num_edges*2];
	int* sol_y = new int[instance->num_vertices];
	int* sol_f = new int[instance->num_edges*2];
	int* sol_f0 = new int[instance->num_vertices];
	int sol_eta = 0;

	// initial values all zero using memset
	memset(sol_x, 0, sizeof(int) * instance->num_edges*2);
	memset(sol_y, 0, sizeof(int) * instance->num_vertices);
	memset(sol_f, 0, sizeof(int) * instance->num_edges*2);
	memset(sol_f0, 0, sizeof(int) * instance->num_vertices);


	// print trees
	for (int i = 0; i < instance->num_trees; i++) {
		// get the tree
		_tree* tree = instance->trees_ub[i];
		// print the tree
		tree->PrintVerticesWeight();
	}




	IloNumVarArray vars(env);
	IloNumArray vals(env);	

	// max flow variable
	int max_flow = 0;

	// for each tree, the vertex with lower index is the root
	for (int i = 0; i < instance->num_trees; i++) {
		// get the tree
		_tree* tree = instance->trees_ub[i];
		// get the root
		int root = tree->vertices[0];

		sol_y[root] = 1;

		int flow = TraveseInitSol(sol_x, sol_f, tree, root, -1);

		sol_f0[root] = flow;

		//printf("f0[%d] = %d\n", root, flow);

		max_flow = max(max_flow, flow);
	}

	sol_eta = max_flow;

	//printf("eta = %d\n", max_flow);

	// add variables to the array
	for (int i = 0; i < instance->num_edges*2; i++) {
		if (cplex.isExtracted(x[i])) {
			vals.add(sol_x[i]);
			vars.add(x[i]);
		}
		else {
			printf("x[%d] not extracted\n", i);
		}

		
		if (cplex.isExtracted(f[i])) {
			vals.add(sol_f[i]);
			vars.add(f[i]);
		}
		else {
			printf("f[%d] not extracted\n", i);
		}
	}

	for (int i = 0; i < instance->num_vertices; i++) {
		if (cplex.isExtracted(y[i])) {
			vals.add(sol_y[i]);
			vars.add(y[i]);
		}
		else {
			printf("y[%d] not extracted\n", i);
		}

		if (cplex.isExtracted(f0[i])){
			vals.add(sol_f0[i]);
			vars.add(f0[i]);
		}
		else {
			printf("f0[%d] not extracted\n", i);
		}
	}

	if (cplex.isExtracted(eta)) {
		vals.add(sol_eta);
		vars.add(eta);
	}	
	else {
		printf("eta not extracted\n");
	}

	

	// add the solution to the model
	try {
		cplex.addMIPStart(vars, vals);
	}
	catch (IloException& e) {
		std::cerr << "MIP start rejected: " << e.getMessage() << std::endl;
	}
	vals.end();
	vars.end();


	// let us force solution as constraints

	/*for (int i = 0; i < instance->num_edges*2; i++) {
		if (cplex.isExtracted(x[i])) {
			model.add(x[i] == sol_x[i]);
		}
	}

	for (int i = 0; i < instance->num_vertices; i++) {
		if (cplex.isExtracted(y[i])) {
			model.add(y[i] == sol_y[i]);
		}
	}


	for (int i = 0; i < instance->num_edges*2; i++) {
		if (cplex.isExtracted(f[i])) {
			model.add(f[i] == sol_f[i]);
		}
	}

	for (int i = 0; i < instance->num_vertices; i++) {
		if (cplex.isExtracted(f0[i])) {
			model.add(f0[i] == sol_f0[i]);
		}
	}

	model.add(eta == sol_eta);*/

	delete sol_x;
	delete sol_y;
	delete sol_f;
	delete sol_f0;	

	return this;
}

// print the solution
FlowF* FlowF::PrintSol() {
	printf("Objective value: %f\n", cplex.getObjValue());
	for (int i = 0; i < instance->num_edges*2; i++) {
		if (cplex.isExtracted(x[i]))
		if (cplex.getValue(x[i]) > 0.00001) {
			int u = instance->arcs[i][0];
			int v = instance->arcs[i][1];
			printf("x[%d,%d] = %f \n", u,v, cplex.getValue(x[i]));
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
	for (int i = 0; i < instance->num_edges*2; i++) {
		if (cplex.isExtracted(f[i]))
		if (cplex.getValue(f[i]) > 0.00001) {
			int u = instance->arcs[i][0];
			int v = instance->arcs[i][1];
			printf("f[%d,%d] = %f \n", u, v, cplex.getValue(f[i]));
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
	x = IloNumVarArray(env, instance->num_edges * 2, 0, 1, integer?ILOINT:ILOFLOAT);
	for (int i = 0; i < instance->num_edges * 2; i++) {
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
	f = IloNumVarArray(env, instance->num_edges*2, 0, IloInfinity, ILOFLOAT);
	for (int i = 0; i < instance->num_edges*2; i++) {
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
	if (mutual_exclusion_cycles) {
		AddConsCycleXY();
	}	
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
		for (int a = 0; a < instance->num_edges*2; a++) {
			if (instance->arcs[a][1] == v) {
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
		for (int a = 0; a < instance->num_edges*2; a++) {
			if (instance->arcs[a][1] == v) {
				cons += f[a];
			}
			if (instance->arcs[a][0] == v) {
				cons -= f[a];
			}
		}

		IloExpr cons2(env);
		for (int a = 0; a < instance->num_edges*2; a++) {
			if (instance->arcs[a][1] == v) {
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
	for (int a = 0; a < instance->num_edges*2; a++) {
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
		for (int a = 0; a < instance->num_edges*2; a++) {
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
	for (int a = 0; a < instance->num_edges*2; a++){
		for (int a2 = 0; a2 < instance->num_edges*2; a2++) {
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