#include "CutF.h"
#include "binary.h"

#include <lemon/edmonds_karp.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>


ILOLAZYCONSTRAINTCALLBACK1(callback, CutF*, cutF) {
	cutF->lazycallback_count++;
	cutF->instance->n_lazy_cuts++;

	bool allConflictsAreAdded = true;

	// This part has no effect
	//This code is designed to add a user cut to prevent an infeasible or undesired solution in a graph-based problem, likely involving spanning trees or multi-tree optimization. The cut is based on the following logic:
	//If two edges e and f are incident to the same vertex v and belong to different trees(i and j) :
	//	It checks if both edges are selected(x1 + x2 > 1).
	//	If this condition is true, it adds a constraint to prevent this by ensuring x[e][i] + x[f][j] <= 1.
	//	This effectively prevents two conflicting edges from being included in the solution for different trees at the same time.
	//// Iterate over all vertices
	//for (int v = 0; v < cutF->instance->num_vertices; v++) {
	//	// Iterate over all edges
	//	for (int e = 0; e < cutF->instance->num_edges; e++) {
	//		// Check if the edge 'e' is incident to the vertex 'v'
	//		if (cutF->instance->edges[e][0] == v || cutF->instance->edges[e][1] == v) {
	//			// Iterate over all edges again (to check for pairs of edges)
	//			for (int f = 0; f < cutF->instance->num_edges; f++) {
	//				if (e != f) {
	//					// Check if edge 'f' is also incident to the vertex 'v'
	//					if (cutF->instance->edges[f][0] == v || cutF->instance->edges[f][1] == v) {
	//						// Iterate over all pairs of trees (i, j)
	//						for (int i = 0; i < cutF->instance->num_trees - 1; i++) {
	//							for (int j = i + 1; j < cutF->instance->num_trees; j++) {
	//								// Get the values of the decision variables from the current solution
	//								int x1 = getValue(cutF->x[e][i]);
	//								int x2 = getValue(cutF->x[f][j]);
	//								// Check if the sum of x1 and x2 exceeds 1 (invalid configuration)
	//								if (x1 + x2 > 1) {
	//									// A conflict has been detected, so add a cut
	//									allConflictsAreAdded = false;
	//									IloExpr cons(cutF->env);
	//									// Create the cut: x[e][i] + x[f][j] <= 1
	//									cons += cutF->x[e][i] + cutF->x[f][j];
	//									// Optionally print the cut for debugging
	//									if (cutF->printCuts) {
	//										printf_s("cut added: ");
	//										printf_s("x[%d][%d] + x[%d][%d] <= 1\n", e, i, f, j);
	//									}
	//									// Add the cut to the CPLEX model
	//									add(cons <= 1);
	//									cutF->cplex.addUserCut(cons <= 1);
	//									cons.end();  // Clean up the expression
	//								}
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}


	if (allConflictsAreAdded) {
		//cutF->SaveOpt();  /// there is an bug here, the solution is not saved correctly, the code below is working for now

		// get the values of the variables
		cutF->instance->num_opt_edges = 0;
		for (int i = 0; i < cutF->instance->num_trees; i++) {
			cutF->instance->start_edge_tree[i] = cutF->instance->num_opt_edges;
			for (int e = 0; e < cutF->instance->num_edges; e++) {
				if (cutF->cplex.isExtracted(cutF->x[e][i])) {
					double value = getValue(cutF->x[e][i]);
					if (value > 0.5) {
						cutF->instance->opt_edges[cutF->instance->num_opt_edges++] = e;
					}
				}
			}
		}

		if (cutF->printCuts) {
			// get the objective value
			int opt_value = (int)getObjValue();
			// print objective value
			printf_s("\n**** INTEGER SOL (Callback) ****", opt_value);
			printf_s("\nObjective value: %d", opt_value);
			cutF->instance->PrintOptEdges();
		}

		// check if the solution has cycles
		if (cutF->instance->CheckCyclesInOptEdges())
		{

			// print the cycles
			if (cutF->printCycles) {
				cutF->instance->PrintCycles();
			}


			// add the cycle elimination constraint
			for (int k = 0; k < cutF->instance->num_cycles; k++) {
				IloExpr cons(cutF->env);
				int count = 0;
				if (cutF->printCuts) {
					printf_s("\ncycle cut added: ");
				}

				for (int e = 0; e < cutF->instance->num_edges; e++) {
					if (cutF->instance->cycles[k][e]) {
						count++;
						for (int i = 0; i < cutF->instance->num_trees; i++) {
							cons += cutF->x[e][i];
							if (cutF->printCuts) {
								printf_s("x[%d][%d] +", e, i);
							}
						}
					}
				}
				add(cons <= count - 1);
				//cutF->cplex.addUserCut(cons <= count - 1);
				if (cutF->printCuts) {
					printf_s("<= %d", count - 1);
				}
				cons.end();
			}

		}
	}
};

ILOUSERCUTCALLBACK1(callbackuser, CutF*, cutF) {	
	if (getNnodes() == 0) {
		cutF->instance->_lb_root = getBestObjValue(); // LP relaxation bound at root		
	}

	cutF->usercallback_count++;
	
	// create a two dimensional array to store the opt x

	for (int e = 0; e < cutF->instance->num_edges; e++) {
		for (int i = 0; i < cutF->instance->num_trees; i++) {
			cutF->opt_x[e][i] = 0;
		}
	}

	// store the opt x values in the array opt_x
	for (int i = 0; i < cutF->instance->num_trees; i++) {
		for (int e = 0; e < cutF->instance->num_edges; e++) {
			if (cutF->cplex.isExtracted(cutF->x[e][i])) {
				double value = getValue(cutF->x[e][i]);				
				cutF->opt_x[e][i] = value;
			}
		}
	}

	// print all opt x
	/*for (int i = 0; i < cutF->instance->num_edges; i++) {
		for (int j = 0; j < cutF->instance->num_trees; j++) {
			printf_s("%d,%d,%4.2lf \t", cutF->instance->edges[i][0], cutF->instance->edges[i][1], opt_x[i][j]);
		}
		printf_s("\n");
	}*/

	// create a capacity map
	ListDigraph::ArcMap<double> capacity(cutF->instance->dg);


	// set capacity for all arcs (u,v)
	for (ListDigraph::ArcIt a(cutF->instance->dg); a != INVALID; ++a) {	
		double sum_x = 0;

		// if u is the source
		if (cutF->instance->dg.id(cutF->instance->dg.source(a)) == cutF->instance->num_vertices) {		
				for (int e = 0; e < cutF->instance->num_edges; e++) {
					int v = cutF->instance->dg.id(cutF->instance->dg.target(a));
					for (int i = 0; i < cutF->instance->num_trees; i++) {
						if (cutF->instance->edges[e][0] == v || cutF->instance->edges[e][1] == v) {
							sum_x += cutF->opt_x[e][i];
						}
					}
				}
				capacity[a] = sum_x/2;	
		}
		// if v is the target
		else if (cutF->instance->dg.id(cutF->instance->dg.target(a)) == cutF->instance->num_vertices+1) {			
			capacity[a] = 1;
		}

		


		
		// otherwise
		else {
			for (int e = 0; e < cutF->instance->num_edges; e++) {
				if (cutF->instance->edges[e][0] == cutF->instance->dg.id(cutF->instance->dg.source(a)) && cutF->instance->edges[e][1] == cutF->instance->dg.id(cutF->instance->dg.target(a))) {
					for (int i = 0; i < cutF->instance->num_trees; i++) {
						sum_x += cutF->opt_x[e][i];
					}
				}
				else if (cutF->instance->edges[e][0] == cutF->instance->dg.id(cutF->instance->dg.target(a)) && cutF->instance->edges[e][1] == cutF->instance->dg.id(cutF->instance->dg.source(a))) {
					for (int i = 0; i < cutF->instance->num_trees; i++) {
						sum_x += cutF->opt_x[e][i];
					}					
				}
			}
			capacity[a] = sum_x / 2;
		}	

		// print the capacity of the arc
		//printf("Arc %d -> %d: %lf\n", cutF->instance->dg.id(cutF->instance->dg.source(a)), cutF->instance->dg.id(cutF->instance->dg.target(a)), capacity[a]);
	}

	//// print the capacity 
	//for (ListDigraph::ArcIt a(cutF->instance->dg); a != INVALID; ++a) {
	//	printf("Arc %d -> %d: %lf\n", cutF->instance->dg.id(cutF->instance->dg.source(a)), cutF->instance->dg.id(cutF->instance->dg.target(a)), capacity[a]);
	//}	

	for (ListDigraph::ArcIt a(cutF->instance->dg); a != INVALID; ++a) {
		if (cutF->instance->dg.id(cutF->instance->dg.source(a)) == cutF->instance->num_vertices) {
			double org_capacity = capacity[a];
			capacity[a] = 100;
			// print a
			//if (cutF->printCuts)
				//printf_s("\nArc %d -> %d: %lf  -->", cutF->instance->dg.id(cutF->instance->dg.source(a)), cutF->instance->dg.id(cutF->instance->dg.target(a)), capacity[a]);

			EdmondsKarp<ListDigraph, ListDigraph::ArcMap<double>> ho(
				cutF->instance->dg,
				capacity,
				cutF->instance->dg.nodeFromId(cutF->instance->num_vertices),
				cutF->instance->dg.nodeFromId(cutF->instance->num_vertices + 1)
			);

			ho.run();

			//create a cut map to store the min cut
			ListDigraph::NodeMap<bool> minCut(cutF->instance->dg);

			// map the min cut
			ho.minCutMap(minCut);
			//printf_s("\n S:");

			//// print the min cut
			//for (ListDigraph::NodeIt v(cutF->instance->dg); v != INVALID; ++v) {
			//	if (minCut[v]) {
			//		printf_s("%d ", cutF->instance->dg.id(v));
			//	}
			//}

			// check all arcs that have a source in the cut and a target outside the cut
			//double sum = 0; 
			//for (ListDigraph::ArcIt a(cutF->instance->dg); a != INVALID; ++a) {
			//	if (minCut[cutF->instance->dg.source(a)] && !minCut[cutF->instance->dg.target(a)]) {
			//		// print the arc
			//		if (cutF->printCuts)
			//			if (capacity[a] > 0)
			//			printf_s("\nArc %d -> %d: %lf", cutF->instance->dg.id(cutF->instance->dg.source(a)), cutF->instance->dg.id(cutF->instance->dg.target(a)), capacity[a]);

			//		sum += capacity[a];
			//	}
			//}

			//// print the sum

			//if (cutF->printCuts)
			//	printf_s("\nSum: %lf \n", sum);	


			// compute S 
			int nS = 0;
			for (ListDigraph::NodeIt v(cutF->instance->dg); v != INVALID; ++v) {
				if (minCut[v]) {					
					if (cutF->instance->dg.id(v) != cutF->instance->num_vertices) {
						cutF->S[nS++] = cutF->instance->dg.id(v);
						

					}
				}
			}
			
			// compute cost: for each edge that doesnt have a vertex in S,  add the cost to the cut
			double cost = nS;
			for (int e = 0; e < cutF->instance->num_edges; e++) {
				bool SinS = false; // if source in S
				for (int j = 0; j < nS; j++) {
					int v = cutF->S[j];
					if (cutF->instance->edges[e][0] == v) {
						SinS = true;
						break;
					}
				}
				bool DinS = false;	// if dist in S
				for (int j = 0; j < nS; j++) {
					int v = cutF->S[j];
					if (cutF->instance->edges[e][1] == v) {
						DinS = true;
						break;
					}
				}
				if (!SinS || !DinS) { // if the edge is outside the cut
					for (int i = 0; i < cutF->instance->num_trees; i++) {
						cost += cutF->opt_x[e][i];
					}
				}
			}		

			// print cost and cut
			//printf("\nCost: %lf\n", cost);

			// if cost < n - k - 1, add the cut
			
			if (nS >= 2 && cost < cutF->instance->num_vertices - cutF->instance->num_trees + 1) {
				IloExpr cons(cutF->env);

				if (cutF->printCuts)
					printf_s("\n cut added: ");

				for (int e = 0; e < cutF->instance->num_edges; e++) {
					bool SinS = false; // if source in S
					for (int j = 0; j < nS; j++) {
						int v = cutF->S[j];
						if (cutF->instance->edges[e][0] == v) {
							SinS = true;
							break;
						}
					}
					bool DinS = false;	// if dist in S
					for (int j = 0; j < nS; j++) {
						int v = cutF->S[j];
						if (cutF->instance->edges[e][1] == v) {
							DinS = true;
							break;
						}
					}

					if (SinS && DinS) {
						for (int i = 0; i < cutF->instance->num_trees; i++) {
							if (cutF->cplex.isExtracted(cutF->x[e][i])) {
								cons += cutF->x[e][i];
								if (cutF->printCuts)									
									printf_s("x[%d][%d] +", e, i);
									
							}							
						}
					}
				}

				cutF->instance->n_user_cuts++;
				add(cons <= nS - 1);
				if (cutF->printCuts)
					printf_s("<= %d", nS - 1);												

				cons.end();

				// break as soon as one cut is added.
				break;
			}			
			
			
			// restore the capacity
			capacity[a] = org_capacity;
		}					
	}		
}

void CutF::SaveOpt() {
	this->instance->num_opt_edges = 0;
	for (int i = 0; i < this->instance->num_trees; i++) {		
		this->instance->start_edge_tree[i] = this->instance->num_opt_edges;
		for (int e = 0; e < this->instance->num_edges; e++) {
			if (this->cplex.isExtracted(this->x[e][i])) {
				try
				{		
					double value = this->cplex.getValue(this->x[e][i]);
					if (value > 0.5) {
						this->instance->opt_edges[this->instance->num_opt_edges++] = e;
					}
				}
				catch (const std::exception&)
				{
					printf("Error\n");
				}								
			}
		}
	}

	// get the objective value
	this->instance->_opt = opt;
	this->instance->_gap = gap;
}

CutF::CutF(_g* instance, bool redirect = false, bool linear = false) :Cplex(instance) {
	printCuts = false; 
	printCycles = false;

	// initialize the arrays for callbacks
	opt_x = new double* [instance->num_edges];
	for (int i = 0; i < instance->num_edges; i++) {
		opt_x[i] = new double[instance->num_trees];
	}
	S = new int[instance->num_vertices];

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

CutF::~CutF() {
	// delete x variables
	for (int i = 0; i < instance->num_edges; i++) {
		x[i].end();
	}

	x.end();

	
	// delete callback arrays
	delete S;
	for (int i = 0; i < instance->num_edges; i++) {
		delete opt_x[i];
	}
	delete opt_x;
}

CutF* CutF::Run() {	

	cplex.use(callback(env, this)); // TODO: revert this
	cplex.use(callbackuser(env, this));
	Cplex::Run();

	SaveOpt();
	//instance->PrintOptEdges();

	return this;
}

CutF* CutF::SetPrintCuts(bool printCuts) {	
	Cplex::SetPrintCuts(printCuts);
	return this;
}

CutF* CutF::SetPrintCycles(bool printCycles) {	
	Cplex::SetPrintCycles(printCycles);
	return this;
}

CutF* CutF::PrintModel() {
	cplex.exportModel("ModelCut.lp");
	return this;
}

CutF* CutF::SetInitSol() {
	
	IloNumVarArray vars(env);
	IloNumArray vals(env);

	// create a binary num_edge * num_trees array
	bool** sol = new bool* [instance->num_edges];
	for (int i = 0; i < instance->num_edges; i++) {
		sol[i] = new bool[instance->num_trees];
		memset(sol[i], 0, sizeof(bool) * instance->num_trees);
	}

	for (int i = 0; i < instance->trees_ub.size(); i++) {
		_tree* tree = instance->trees_ub[i];
		for (int j = 0; j < tree->num_edges; j++) {
			int e = tree->edges[j];
			sol[e][i] = 1;
		}
	}

	for (int e = 0; e < instance->num_edges; e++) {		
		for (int i = 0; i < instance->num_trees; i++) {
			vals.add(sol[e][i]);
			vars.add(x[e][i]);
		}
		
	}
	try {
		cplex.addMIPStart(vars, vals);
	}
	catch (IloException& e) {
		std::cerr << "MIP start rejected: " << e.getMessage() << std::endl;
	}	
	vals.end();
	vars.end();

	return this;
}

// Force Solution	
CutF* CutF::ForceSol() {
	/*model.add(x[8][0] == 1);
	model.add(x[4][0] == 1);*/


	for (int i = 0; i < instance->trees_ub.size(); i++) {
		_tree* tree = instance->trees_ub[i];
		for (int j = 0; j < tree->num_edges; j++) {
			int e = tree->edges[j];
			model.add(x[e][i] == 1);
		}
	}

	return this;
}

void CutF::DefVar() {
	DefVarX();
}

void CutF::DefVarX() {
	x = NumVarMatrix(env, instance->num_edges);
	for (int i = 0; i < instance->num_edges; i++) {
		x[i] = IloNumVarArray(env, instance->num_trees, 0, 1, integer?ILOINT:ILOFLOAT);
		for (int j = 0; j < instance->num_trees; j++) {
			sprintf_s(name, "x(%d,%d)", i, j);
			x[i][j].setName(name);			
		}
	}
}

void CutF::AddObj() {
	IloExpr exprObj(env);	
	// minimize the weight of the first tree
	for (int e = 0; e < instance->num_edges; e++) {		
		exprObj += instance->edges[e][2] * x[e][0];		
	}

	model.add(IloMinimize(env, exprObj));
	exprObj.end();
}

void CutF::AddCons() {
	// ordering the trees
	AddConsOrderTrees();
	
	// pick n-k edges
	AddConstraintPicknkV();
	
	// assign each vertex to a tree
	AddConsEachToTree();

	// the next three are alternatives 
	//AddConsTwoEdgesOfVertex();
	AddConsEdgesOfVertex();
	//AddConsEdgesOfVertex_Improved(); // too many combinations
	
	// each vertex has at least one edge
	//AddConsAtLeastOneEdge();  // we cannot add this constraint, because it eliminates singletons 

	// bounded tree by UB
	AddConsBoundByUB();

	// add the forbidden pair edges
	//AddConsForbiddenPairEdges();
}

void CutF::AddConsAtLeastOneEdge() {
	for (int i = 0; i < instance->num_trees; i++) {
		IloExpr cons(env);
		for (int e = 0; e < instance->num_edges; e++) {
			cons += this->x[e][i];
		}
		sprintf_s(name, "AtLeastOneEdge(%d)", i);
		model.add(cons >= 1).setName(name);
		cons.end();
	}
}

void CutF::AddConsOrderTrees()
{
	// add constraints to ensure that the trees are ordered	
	for (int i = 0; i < instance->num_trees-1; i++) {
		IloExpr consi(env);
		IloExpr consip(env);
		for (int e = 0; e < instance->num_edges; e++) {			
			 consi += this->x[e][i] * instance->edges[e][2];			
			 consip += this->x[e][i + 1] * instance->edges[e][2];
		}
		sprintf_s(name, "OrderTrees(%d)", i);
		model.add(consi >= consip).setName(name);
		consi.end();
		consip.end();
	}	
}

void CutF::AddConstraintPicknkV() {
	// add constraints for each subset
	IloExpr cons(env);
	for (int e = 0; e < instance->num_edges; e++) {
		for (int i = 0; i < instance->num_trees; i++) {
			cons += this->x[e][i];
		}
	}
	sprintf_s(name, "PicknkV");
	model.add(cons == instance->num_vertices - instance->num_trees).setName(name);
	cons.end();
}

void CutF::AddConsForEachSubset()
{
	
	this->instance->computeSubsets();
	for (int i = 0; i < instance->num_trees; i++) {
		for (int s = 0; s < pow(2, instance->num_vertices); s++) {
			if (instance->nSubsets[s] > 2 && instance->nSubsets[s] < pow(2, instance->num_vertices)){
				IloExpr cons(env);
				int count = 0;
				for (int e = 0; e < instance->num_edges; e++) {
					if (checkbin(instance->subsets[s], instance->edges[e][0]) && checkbin(instance->subsets[s], instance->edges[e][1])) {
						cons += this->x[e][i];
						count++;
					}
				}
				if (count > instance->nSubsets[s] - 1) {
					sprintf_s(name, "ForEachSubset(%d,%d)", i, s);
					model.add(cons <= instance->nSubsets[s] - 1).setName(name);
				}
				cons.end();
			}
		}
	}

	
}

void CutF::AddConsTwoEdgesOfVertex() {
	for (int v = 0; v< instance->num_vertices; v++) {
		for (int e = 0; e < instance->num_edges; e++) {
			if (instance->edges[e][0] == v || instance->edges[e][1] == v) {
				for (int f=0; f < instance->num_edges; f++) {
					if (e != f) {
						if (instance->edges[f][0] == v || instance->edges[f][1] == v) {						
							for (int i = 0; i < instance->num_trees-1; i++) {
								for (int j = i + 1; j < instance->num_trees; j++) {
									sprintf_s(name, "TwoEdgesOfVertex(%d,%d,%d,%d,%d)", v, e, f, i, j);
									model.add(x[e][i] + x[f][j] <= 1).setName(name);
								}																				
							}
						}
					}
				}
			}
		}
	}	
}

void CutF::AddConsEdgesOfVertex() {
	for (int v = 0; v < instance->num_vertices; v++) {
		for (int e = 0; e < instance->num_edges; e++) {
			if (instance->edges[e][0] == v || instance->edges[e][1] == v) {
				for (int f = e+1; f < instance->num_edges; f++) {
					if (e != f) {
						if (instance->edges[f][0] == v || instance->edges[f][1] == v) {
							for (int i = 0; i < instance->num_trees - 1; i++) {
								IloExpr cons(env); 
								IloExpr cons2(env);
								for (int j = 0; j < instance->num_trees; j++) {
									if (i != j)
									{
										cons += x[f][j];
										cons2 += x[e][j];
									}									
								}
								sprintf_s(name, "EdgesOfVertex(%d,%d,%d,%d)", v, e, f, i);
								model.add(x[e][i] + cons <= 1).setName(name);
								model.add(x[f][i] + cons2 <= 1).setName(name);
								cons.end();
								cons2.end();
							}
						}
					}
				}
			}
		}
	}
}

void CutF::AddConsEdgesOfVertex_Improved() {
	/*for (int v = 0; v < instance->num_vertices; v++) {
		for (int e = 0; e < instance->num_edges; e++) {
			if (instance->edges[e][0] == v || instance->edges[e][1] == v) {
				for (int f = e + 1; f < instance->num_edges; f++) {
					if (e != f) {
						if (instance->edges[f][0] == v || instance->edges[f][1] == v) {
							for (int i = 0; i < instance->num_trees; i++) {
								IloExpr cons(env);
								IloExpr cons2(env);
								for (int j = 0; j < i; j++) {									
									cons += x[e][j];
									cons2 += x[f][j];									
								}
								for (int j = i; j < instance->num_trees; j++) {
									cons += x[e][j];
									cons2 += x[f][j];
								}
								sprintf_s(name, "EdgesOfVertex(%d,%d,%d,%d)", v, e, f, i);
								model.add(x[e][i] + cons <= 1).setName(name);
								model.add(x[f][i] + cons2 <= 1).setName(name);
							}
						}
					}
				}
			}
		}
	}*/
}

void CutF::AddConsEachToTree() {
	for (int e = 0; e < instance->num_edges; e++) {
		IloExpr cons(env);
		for (int i = 0; i < instance->num_trees; i++) {	
			cons += this->x[e][i];			
		}
		sprintf_s(name, "EachEdgeToTree(%d)", e);
		model.add(cons <= 1).setName(name);
		cons.end();
	}
}

void CutF::AddConsBoundByUB() {
	for (int i = 0; i < instance->num_trees; i++) {
		IloExpr exprObj(env);
		for (int e = 0; e < instance->num_edges; e++) {
			exprObj += instance->edges[e][2] * x[e][i];
		}
		model.add(exprObj <= instance->UB + 0.001);
		exprObj.end();
	}
}

void CutF::AddConsForbiddenPairEdges() {
	// for each pair of edges, if their shortest path weight pairs larger than UB, they cannot be in the same tree
	for (int i = 0; i < instance->num_trees; i++) {
		for (int e = 0; e < instance->num_edges; e++) {
			for (int f = e + 1; f < instance->num_edges; f++) {
				if (instance->shortest_path_weight_edges[e][f] > instance->UB) {
					sprintf_s(name, "ForbiddenPairEdges(%d,%d,%d)", e, f, i);
					model.add(x[e][i] + x[f][i] <= 1).setName(name);
				}
			}
		}
	}
}

// set integer
CutF* CutF::SetInteger() {
	Cplex::SetInteger();
	return this;
}

// set linear
CutF* CutF::SetLinear() {
	Cplex::SetLinear();
	return this;
}
