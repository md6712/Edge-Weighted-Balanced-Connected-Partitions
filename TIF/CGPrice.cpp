#include "CGPrice.h"
#include "Cplex.h"
#include "binary.h"	

ILOLAZYCONSTRAINTCALLBACK1(callback_CG_Price, CGPrice*, cg_price) {

	// get current objective value for the last integer solution
	double current_obj = (double)getObjValue();

	// if this is positive, stop the algorithm
	//if (current_obj > 0.001) {
	//	//cg_price->halted = true;
	//	cg_price->ub_positive = current_obj;

	//	// copy x values for each edge
	//	for (int e = 0; e < cg_price->instance->num_edges; e++) {
	//		if (cg_price->cplex.isExtracted(cg_price->x[e])) {
	//			cg_price->x_value[e] = getValue(cg_price->x[e]);
	//		}
	//		else {
	//			cg_price->x_value[e] = 0;
	//		}			
	//	}

	//	// copy y values for each vertex
	//	for (int v = 0; v < cg_price->instance->num_vertices; v++) {
	//		if (cg_price->cplex.isExtracted(cg_price->y[v])) {
	//			cg_price->y_value[v] = getValue(cg_price->y[v]);
	//		}
	//		else {
	//			cg_price->y_value[v] = 0;
	//		}
	//	}		

	//	//abort();
	//	return;
	//}


	bool allConflictsAreAdded = true;

	if (allConflictsAreAdded) {
		cg_price->instance->num_opt_edges = 0;
		
		cg_price->instance->start_edge_tree[0] = 0;
		cg_price->instance->start_edge_tree[1] = 100;
		for (int e = 0; e < cg_price->instance->num_edges; e++) {
			if (cg_price->cplex.isExtracted(cg_price->x[e])) {
				double value = getValue(cg_price->x[e]);
				if (value > 0.5) {
					cg_price->instance->opt_edges[cg_price->instance->num_opt_edges++] = e;					
				}
			}
		}

		

		if (cg_price->printCuts) {
			// get the objective value
			int opt_value = (int)getObjValue();
			// print objective value			
			cg_price->instance->PrintOptEdges();			
			printf_s("Objective value: %d", opt_value);
		}

		// check if the solution has cycles
		if (cg_price->instance->CheckCyclesInOptEdges())
		{

			// print the cycles
			if (cg_price->printCycles) {
				cg_price->instance->PrintCycles();
			}


			// add the cycle elimination constraint
			for (int k = 0; k < cg_price->instance->num_cycles; k++) {
				IloExpr cons(cg_price->env);
				int count = 0;
				if (cg_price->printCuts) {
					printf_s("cycle cut added: ");
				}

				for (int e = 0; e < cg_price->instance->num_edges; e++) {
					if (cg_price->instance->cycles[k][e]) {
						count++;						
						cons += cg_price->x[e];
						if (cg_price->printCuts) {
							printf_s("x[%d] +", e);
						}						
					}
				}
				add(cons <= count - 1);
				cg_price->cplex.addUserCut(cons <= count - 1);
				if (cg_price->printCuts) {
					printf_s("<= %d \n", count - 1);
				}
				cons.end();
			}

		}
	}
};

ILOUSERCUTCALLBACK1(callbackuser_CG_Price, CGPrice*, cg_price) {	
	// create a two dimensional array to store the opt x
	double* opt_x = new double[cg_price->instance->num_edges];
	int* S = new int[cg_price->instance->num_vertices];

	// set opt_x to 0
	memset(opt_x, 0, sizeof(double) * cg_price->instance->num_edges);

	// store the opt x values in the array opt_x	
	for (int e = 0; e < cg_price->instance->num_edges; e++) {
		if (cg_price->cplex.isExtracted(cg_price->x[e])) {
			double value = getValue(cg_price->x[e]);
			opt_x[e] = value;
		}
	}	
	
	// create a capacity map
	ListDigraph::ArcMap<double> capacity(cg_price->instance->dg);


	// set capacity for all arcs (u,v)
	for (ListDigraph::ArcIt a(cg_price->instance->dg); a != INVALID; ++a) {
		double sum_x = 0;

		// if u is the source
		if (cg_price->instance->dg.id(cg_price->instance->dg.source(a)) == cg_price->instance->num_vertices) {
			for (int e = 0; e < cg_price->instance->num_edges; e++) {
				int v = cg_price->instance->dg.id(cg_price->instance->dg.target(a));
				if (cg_price->instance->edges[e][0] == v || cg_price->instance->edges[e][1] == v) {
					sum_x += opt_x[e];
				}
			}
			capacity[a] = sum_x / 2;
		}
		// if v is the target
		else if (cg_price->instance->dg.id(cg_price->instance->dg.target(a)) == cg_price->instance->num_vertices + 1) {
			capacity[a] = 1;
		}

		// otherwise
		else {
			for (int e = 0; e < cg_price->instance->num_edges; e++) {
				if (cg_price->instance->edges[e][0] == cg_price->instance->dg.id(cg_price->instance->dg.source(a)) && cg_price->instance->edges[e][1] == cg_price->instance->dg.id(cg_price->instance->dg.target(a))) {					
					sum_x += opt_x[e];					
				}
				else if (cg_price->instance->edges[e][0] == cg_price->instance->dg.id(cg_price->instance->dg.target(a)) && cg_price->instance->edges[e][1] == cg_price->instance->dg.id(cg_price->instance->dg.source(a))) {					
					sum_x += opt_x[e];					
				}
			}
			capacity[a] = sum_x / 2;
		}		
	}	

	for (ListDigraph::ArcIt a(cg_price->instance->dg); a != INVALID; ++a) {
		if (cg_price->instance->dg.id(cg_price->instance->dg.source(a)) == cg_price->instance->num_vertices) {
			double org_capacity = capacity[a];
			capacity[a] = 100;
			
			EdmondsKarp<ListDigraph, ListDigraph::ArcMap<double>> ho(
				cg_price->instance->dg,
				capacity,
				cg_price->instance->dg.nodeFromId(cg_price->instance->num_vertices),
				cg_price->instance->dg.nodeFromId(cg_price->instance->num_vertices + 1)
			);

			ho.run();

			//create a cut map to store the min cut
			ListDigraph::NodeMap<bool> minCut(cg_price->instance->dg);

			// map the min cut
			ho.minCutMap(minCut);					

			// compute S 
			int nS = 0;
			for (ListDigraph::NodeIt v(cg_price->instance->dg); v != INVALID; ++v) {
				if (minCut[v]) {
					if (cg_price->instance->dg.id(v) != cg_price->instance->num_vertices) {
						S[nS++] = cg_price->instance->dg.id(v);


					}
				}
			}

			// compute cost: for each edge that doesnt have a vertex in S,  add the cost to the cut
			double cost = nS;
			for (int e = 0; e < cg_price->instance->num_edges; e++) {
				bool SinS = false; // if source in S
				for (int j = 0; j < nS; j++) {
					int v = S[j];
					if (cg_price->instance->edges[e][0] == v) {
						SinS = true;
						break;
					}
				}
				bool DinS = false;	// if dist in S
				for (int j = 0; j < nS; j++) {
					int v = S[j];
					if (cg_price->instance->edges[e][1] == v) {
						DinS = true;
						break;
					}
				}
				if (!SinS || !DinS) { // if the edge is outside the cut					
					cost += opt_x[e];					
				}
			}


			// if cost < n - k - 1, add the cut

			if (nS >= 2 && cost < cg_price->instance->num_vertices - cg_price->instance->num_trees + 1) {
				IloExpr cons(cg_price->env);

				if (cg_price->printCuts)
					printf_s("\n cut added: ");

				for (int e = 0; e < cg_price->instance->num_edges; e++) {
					bool SinS = false; // if source in S
					for (int j = 0; j < nS; j++) {
						int v = S[j];
						if (cg_price->instance->edges[e][0] == v) {
							SinS = true;
							break;
						}
					}
					bool DinS = false;	// if dist in S
					for (int j = 0; j < nS; j++) {
						int v = S[j];
						if (cg_price->instance->edges[e][1] == v) {
							DinS = true;
							break;
						}
					}

					if (SinS && DinS) {						
						if (cg_price->cplex.isExtracted(cg_price->x[e])) {
							cons += cg_price->x[e];
							if (cg_price->printCuts)
								printf_s("x[%d] +", e);

						}						
					}
				}

				add(cons <= nS - 1);
				if (cg_price->printCuts)
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

// constructor
CGPrice::CGPrice(_g* g, bool redirect):Cplex(g) {
	printCuts = false;
	printCycles = false;

	// define coef 
	theta = IloNum(0);
	eta = IloNumArray(env, this->instance->num_vertices);
	zeta = IloNumArray(env, this->instance->num_vertices);

	// define aborescence
	aborescence = new Abor(this->instance, false, false);

	if (!redirect) {
		// Define Variables 
		DefVar();

		// Add Objective function 
		AddObj();

		// Add Constraints
		AddCons();
	}

	// init x_value and y_value
	x_value = new double[this->instance->num_edges];
	y_value = new double[this->instance->num_vertices];

	// vertices postive zeta
	vertcies_pos_zeta = new int[this->instance->num_vertices];

	// vertices weight
	vertex_weights = new double[this->instance->mwcs->num_vertices];

	// vertex prizes
	vertex_prizes = new double[this->instance->mwcs->num_vertices];

	// edge costs
	edge_costs = new double[this->instance->mwcs->num_edges];

	// roots pcst
	roots_pcst = new int[this->instance->mwcs->num_vertices];


	
}

// destructor
CGPrice::~CGPrice() {
	eta.end();
	zeta.end();

	delete[] x_value;
	delete[] y_value;
	delete[] vertcies_pos_zeta;
	delete[] vertex_weights;

	delete[] edge_costs;
	delete[] vertex_prizes;
	delete[] roots_pcst;

	delete aborescence;
}


CGPrice* CGPrice::SetPrintCuts(bool printCuts) {
	Cplex::SetPrintCuts(printCuts);
	return this;
}

CGPrice* CGPrice::SetPrintCycles(bool printCycles) {
	Cplex::SetPrintCycles(printCycles);
	return this;
}

// print the model
CGPrice* CGPrice::PrintModel() {
	cplex.exportModel("ModelCGPrice.lp");
	return this;
}

// print the solution
CGPrice* CGPrice::PrintSol() {
	cout << "Objective Value: " << cplex.getObjValue() << endl;
	cout << "x: ";
	for (int e = 0; e < this->instance->num_edges; e++) {
		cout << cplex.getValue(x[e]) << " ";
	}
	cout << endl;

	cout << "y: ";
	for (int v = 0; v < this->instance->num_vertices; v++) {
		cout << cplex.getValue(y[v]) << " ";
	}
	cout << endl;

	cout << "phi: " << endl;
	for (int e = 0; e < this->instance->num_edges; e++) {
		cout << "edge " << e << ": ";
		for (int v = 0; v < this->instance->num_vertices; v++) {
			cout << cplex.getValue(phi[e][v]) << " ";
		}
		cout << endl;
	}
	return this;
}

// run
CGPrice* CGPrice::Run() {
	halted = false;
	// set the callback
	cplex.use(callback_CG_Price(env, this));

	// set the user callback
	cplex.use(callbackuser_CG_Price(env, this));
	
	// run the model
	Cplex::Run();

	if (halted) {
		this->opt = ub_positive;
	}
	return this;
}


// define variables
void CGPrice::DefVar() {
	DefVarX();
	DefVarY();
	DefVarPhi();
}

// define x variables
void CGPrice::DefVarX() {
	x = IloNumVarArray(env, this->instance->num_edges, 0, 1, ILOINT);
	// set names
	for (int e = 0; e < this->instance->num_edges; e++) {
		stringstream ss;
		ss << "x_" << e;
		x[e].setName(ss.str().c_str());
	}
}

// define y variables
void CGPrice::DefVarY() {
	y = IloNumVarArray(env, this->instance->num_vertices, 0, 1, ILOINT);
	// set names
	for (int v = 0; v < this->instance->num_vertices; v++) {
		stringstream ss;
		ss << "y_" << v;
		y[v].setName(ss.str().c_str());
	}
}

// define phi variables
void CGPrice::DefVarPhi() {
	phi = NumVarMatrix(env, this->instance->num_edges);
	for (int v = 0; v < this->instance->num_edges; v++) {
		phi[v] = IloNumVarArray(env, this->instance->num_vertices, 0, 1, ILOINT);
		// set names
		for (int e = 0; e < this->instance->num_vertices; e++) {
			stringstream ss;
			ss << "phi_" << v << "_" << e;
			phi[v][e].setName(ss.str().c_str());
		}
	}
}

// define the objective function
void CGPrice::AddObj() {
	IloExpr obj(env);
	for (int v = 0; v < this->instance->num_vertices; v++) {
		obj += eta[v] * y[v];
		for (int e = 0; e < this->instance->num_edges; e++) {
			obj -= zeta[v] * this->instance->edges[e][2] * phi[e][v];
		}
	}
	
	objective = IloMaximize(env,1000+ obj - theta); 

	model.add(objective);
	obj.end();
}

// fix the solution
CGPrice* CGPrice::FixSol() {
	/*model.add(y[4] == 1);
	model.add(y[2] == 1);
	model.add(y[3] == 1);
	model.add(y[5] == 1);
	model.add(y[9] == 1);
	model.add(x[2] == 1);
	model.add(x[14] == 1);*/

	return this;
}

// define the constraints
void CGPrice::AddCons() {
	AddConsNumEdgesVerticesMinusOne();
	XYRelation();
	PhiRelation();
	AtleastOneEdge();

	// cost of tree extra constraint
	AddConsCostOfTree();
}

// add constraints the number of edges selected = the number of vertices selected minus one
void CGPrice::AddConsNumEdgesVerticesMinusOne() {
	IloExpr sum_x(env);
	IloExpr sum_y(env);
	for (int e = 0; e < this->instance->num_edges; e++) {
		sum_x += x[e];
	}

	for (int v = 0; v < this->instance->num_vertices; v++) {
		sum_y += y[v];
	}
	// name
	stringstream ss;
	ss << "Cons_num_edge_num_vertex";
	model.add(sum_x == sum_y - 1).setName(ss.str().c_str());
	sum_x.end();
}

// add constraints that relate x and y variables
void CGPrice::XYRelation() {
	for (int e = 0; e < this->instance->num_edges; e++) {	
		// name
		stringstream ss;
		ss << "Cons_x_y_" << e;
		model.add(x[e] <= y[this->instance->edges[e][0]]).setName(ss.str().c_str());
		ss.str("");
		ss << "Cons_x_y_" << e;
		model.add(x[e] <= y[this->instance->edges[e][1]]).setName(ss.str().c_str());
	}
}

// add constraints that relate phi variables
void CGPrice::PhiRelation() {
	for (int e = 0; e < this->instance->num_edges; e++) {
		for (int v = 0; v < this->instance->num_vertices; v++) {
			model.add(phi[e][v] <= x[e]);
			model.add(phi[e][v] <= y[v]);
			model.add(phi[e][v] >= x[e] + y[v] - 1);
		}
	}
}

// add constraints that ensure for each v that is selected, at least one edge is selected
void CGPrice::AtleastOneEdge() {
	for (int v = 0; v < this->instance->num_vertices; v++) {
		IloExpr sum(env);
		for (int e = 0; e < this->instance->num_edges; e++) {
			if (this->instance->edges[e][0] == v || this->instance->edges[e][1] == v)
				sum += x[e];			
		}				
		model.add(sum >= y[v]);
		sum.end();
	}
}


// add constraints that relate the cost of the tree
void CGPrice::AddConsCostOfTree() {
	IloExpr sum(env);
	for (int e = 0; e < this->instance->num_edges; e++) {
		sum += this->instance->edges[e][2] * x[e];
	}
	model.add(sum <= this->instance->UB);
	sum.end();
}

//set the dual value of theta
CGPrice* CGPrice::setTheta(IloNum value) {
	theta = value;
	return this;
}

// set the dual values of eta
CGPrice* CGPrice::setEta(IloNumArray values) {
	eta.end();
	eta = values;
	return this;
}

// set the dual values of zeta
CGPrice* CGPrice::setZeta(IloNumArray values) {
	zeta.end();
	zeta = values;
	return this;
}

// update the objective coefficients
CGPrice* CGPrice::UpdateObjectiveCoefficients() {

	//for (int v = 0; v < this->instance->num_vertices; v++) {
	//	objective.setLinearCoef(y[v], eta[v]);
	//}

	//for (int e = 0; e < this->instance->num_edges; e++) {
	//	for (int v = 0; v < this->instance->num_vertices; v++) {
	//		objective.setLinearCoef(phi[e][v], zeta[v] * this->instance->edges[e][2]);
	//	}
	//}

	IloExpr updatedExpr(env);

	// Add linear terms as needed
	for (int v = 0; v < this->instance->num_vertices; v++) {		
		updatedExpr += eta[v] * y[v];
	}

	for (int e = 0; e < this->instance->num_edges; e++) {
		for (int v = 0; v < this->instance->num_vertices; v++) {
			updatedExpr -= zeta[v] * this->instance->edges[e][2] * phi[e][v];
		}
	}

	updatedExpr += -theta ;

	objective.setExpr (updatedExpr);
	return this;
}

// get the tree associated to the optimal solution
_small_tree* CGPrice::GetTree() {

	_small_tree* tree = new _small_tree();	

	memset(tree->bin_vertices, 0, sizeof(uint32_t) * SIZE_OF_VERTICES_BINARY);

	for (int v = 0; v < this->instance->num_vertices; v++) {
		if (!halted) {
			if (cplex.getValue(y[v]) > 0) {
				addbin(tree->bin_vertices, v);
			}
		}
		else {
			if (y_value[v] > 0) {
				addbin(tree->bin_vertices, v);
			}
		}
	}

	tree->weight = 0;

	// loop over all edges
	for (int e = 0; e < this->instance->num_edges; e++) {
		if (!halted) {
			if (cplex.getValue(x[e]) > 0) {
				tree->weight += this->instance->edges[e][2];
			}
		}
		else {
			if (x_value[e] > 0) {
				tree->weight += this->instance->edges[e][2];
			}
		}
	}
	
	return tree;
}

// heuristic 
_tree* CGPrice::heuristic() {
	// set the callback
	
	// create an empty tree 
	_tree* tree = new _tree(this->instance);

	// init sumWeight
	int test_sum_weight = 0;
	int sum_weight = 0;

	// init sumZeta 
	double test_sum_zeta = 0;	
	double sum_zeta = 0;

	// reward 
	double reward = -INFINITY;
	double new_reward = -INFINITY;
	double max_reward = -INFINITY;

	// best
	int best_edge = -1;
	int best_sum_weight = 0;
	double best_sum_zeta = 0;

	// compute the reward of the tree
	double current_reward = ComputeReward(tree);

	while (true) {
		
		// for each edge compute cost 
		for (int e = 0; e < this->instance->num_edges; e++) {
			if (tree->num_vertices == 0) {

				// check both vertices
				int u = this->instance->edges[e][0];
				int v = this->instance->edges[e][1];

				// compute the reward
				test_sum_zeta = zeta[u] + zeta[v];
				test_sum_weight = this->instance->edges[e][2];

				new_reward = eta[u] + eta[v] - (test_sum_zeta)*test_sum_weight;
				
				// check if the reward is better
				if (new_reward > max_reward) {
					max_reward = new_reward;				
					best_edge = e;
					best_sum_weight = test_sum_weight;
					best_sum_zeta = test_sum_zeta;
				}
			}
			else {
				
				// check if the edge is in the tree
				if (tree->IsEdgeInTree(e)) {
					continue;
				}

				// get the vertex that is connecting the tree to the edge
				int v = tree->IncidentVertex(e); // this returns -1 if both vertices are in the tree -> a cycle is created
				if (v == -1) {
					continue;
				}

				// compute the other vertex 
				int u = this->instance->edges[e][0] == v ? this->instance->edges[e][1] : this->instance->edges[e][0];

				// check if the weight of the tree is not greater than the UB
				if (tree->weight + this->instance->edges[e][2] > this->instance->UB) {
					continue;
				}

				// compute the reward
				test_sum_zeta = sum_zeta + zeta[u];
				test_sum_weight = sum_weight + this->instance->edges[e][2];

				
				new_reward = reward + eta[u] + (sum_weight * sum_zeta) - (test_sum_zeta)*test_sum_weight;

				// check if the reward is better
				if (new_reward > max_reward) {					

					max_reward = new_reward;
					best_sum_weight = test_sum_weight;
					best_sum_zeta = test_sum_zeta;
					best_edge = e;
				}
			}
		}	

		if (best_edge == -1) {
			break;
		}
		else {
			tree->AddEdge(best_edge);
			sum_weight = best_sum_weight;
			sum_zeta = best_sum_zeta;
			tree->weight = sum_weight;
			tree->PrintTree	();
			reward = max_reward;
			
			best_edge = -1;

			//// draw the graph
			//instance->DrawGraph(tree->edges, tree->num_edges);

			// print reward
			//printf("reward = %f\n", reward - theta);						

			// print current reward
			//printf("current reward = %f\n", current_reward);
		}		
	}

	//instance->DrawGraph(tree->edges, tree->num_edges);

	// current reward
	current_reward = ComputeReward(tree);
	

	// local search
	bool improved = true; 

	int best_e, best_ep;

	while (improved) {

		improved = false;
		 
		max_reward = current_reward; 

		// reset best edge
		best_e = -1;
		best_ep = -1;

		// for all edges that are in the tree
		for (int e = 0; e < this->instance->num_edges; e++) {
			// check if the edge is in the tree
			if (!tree->IsEdgeInTree(e)) {
				continue;
			}

			// get the vertex that is connecting the tree to the edge
			int v = tree->GetAnyIncidentVertex(e);
			if (v == -1) {
				continue;
			}

			// for all edges that are not in the tree
			for (int ep = 0; ep < this->instance->num_edges; ep++) {
				// check if the edge is in the tree
				if (tree->IsEdgeInTree(ep)) {
					continue;
				}

				// get the vertex that is connecting the tree to the edge
				int vp = tree->GetAnyIncidentVertex(e);
				if (vp == -1) {
					continue;
				}				

				tree->RemoveEdge(e);
				tree->AddEdge(ep);						
							
				// print tree	
				//tree->PrintTree();

				//instance->DrawGraph(tree->edges, tree->num_edges);

				if (tree->IsSpanningTree()) {				
					
					new_reward = ComputeReward(tree);

					if (new_reward > max_reward) {
						max_reward = new_reward;
						improved = true;	
						best_e = e;	
						best_ep = ep;
					}					

					// print the edges removed and added
					//printf("TEST ---- Removed: (%d %d) Added: (%d %d) -- Reward: %4.2lf \n", this->instance->edges[e][0], this->instance->edges[e][1], this->instance->edges[ep][0], this->instance->edges[ep][1], new_reward);
				}
				

				// reverse the change 
				tree->RemoveEdge(ep);
				tree->AddEdge(e);								
			}			
		}

		if (improved) {
			tree->RemoveEdge(best_e);
			tree->AddEdge(best_ep);
			current_reward = max_reward;

			// print the edges removed and added
			printf("DEPLOYED ---- Removed: (%d %d) Added: (%d %d)\n -- Reward: %4.2lf", this->instance->edges[best_e][0], this->instance->edges[best_e][1], this->instance->edges[best_ep][0], this->instance->edges[best_ep][1], current_reward);
		}
	}



	double total_reward = ComputeReward(tree);

	this->opt = total_reward;
	
	tree->reduced_cost = total_reward;	

	return tree;
}

// compute reward
double CGPrice::ComputeReward(_tree* tree) {
	double reward = 0;

	for (int i = 0; i < tree->num_vertices; i++) {
		int v = tree->vertices[i];
		reward += eta[v];
		reward -= zeta[v] * tree->weight;
	}

	reward -= theta;

	return reward;
}


// solve using minimum weight connected subgraph
CGPrice* CGPrice::solve_mwcs() {

	best_mwcs = -INFINITY;

	// get time in ticks
	clock_t start = clock();

	// count and store vertices with positive zeta
	num_vertices_pos_zeta = 0;
	for (int v = 0; v < this->instance->num_vertices; v++) {
		if (zeta[v] > 0.001) {
			vertcies_pos_zeta[num_vertices_pos_zeta++] = v;
		}
	}

	// compute vertex weight from 0 to num_vertex
	for (int v = 0; v < instance->num_vertices; v++) {
		vertex_weights[v] = eta[v];
	}

	double total_zeta = 0;


	// from zero to 2^num_vertices_pos_zeta 
	for (int i = 0; i < pow(2, num_vertices_pos_zeta); i++) {

		total_zeta = 0;
		for (int j = 0; j < num_vertices_pos_zeta; j++) {
			if (checkbinSingle(i, j)) {
				total_zeta += zeta[vertcies_pos_zeta[j]];
				vertex_weights[vertcies_pos_zeta[j]] = eta[vertcies_pos_zeta[j]];
			}
			else {
				vertex_weights[vertcies_pos_zeta[j]] = -1000000;
			}
		}

		// compute vertex weight for each edge
		for (int e = 0; e < instance->num_edges; e++) {
			vertex_weights[instance->num_vertices + e] = -instance->edges[e][2] * total_zeta;
		}

		this->instance->mwcs
			->set_vertex_weights(this->vertex_weights, this->theta)
			->print_in_file()
			->solve()
			->read_solution_from_file();	

		if (this->instance->mwcs->objective_value > best_mwcs) {
			best_mwcs = this->instance->mwcs->objective_value;		
		}

		if (this->instance->mwcs->objective_value > 0) {
			//break;
		}


		
	}

	// print the best mwcs
	printf("Best MWCS: %4.2lf\n", best_mwcs);

	// compute elapsed time
	clock_t end = clock();
	double elapsed_time = double(end - start) / CLOCKS_PER_SEC;

	// print time elapsed
	printf("Time elapsed: %4.3f\n", elapsed_time);

	return this;

}



// solve using prize collecting steiner tree
CGPrice* CGPrice::solve_pcst() {
	best_pcst = -INFINITY;

	// get time in ticks
	clock_t start = clock();

	// count and store vertices with positive zeta
	num_vertices_pos_zeta = 0;
	for (int v = 0; v < this->instance->num_vertices; v++) {
		if (zeta[v] > 0.001) {
			vertcies_pos_zeta[num_vertices_pos_zeta++] = v;
		}
	}
	
	//// min prize
	//double min_prize = 1000000;

	// compute vertex weight from 0 to num_vertex	
	for (int v = 0; v < instance->num_vertices; v++) {
		vertex_prizes[v] = eta[v];  // eta[v] is the prize of vertex v
		//if (vertex_prizes[v] < min_prize) {  // find the minimum prize
		//	min_prize = vertex_prizes[v]; 
		//}
	}

	// compute edge costs 
	for (int e = 0; e < instance->num_edges; e++) {
		edge_costs[e] = instance->edges[e][2]; // edge cost is the negative of the edge weight
	}

	//// set minum prize as zerp
	//if (min_prize > 0) min_prize = 0;

	//// adjust based on the minimum prize
	//if (min_prize < 0) {
	//	// adjust the prizes	
	//	for (int v = 0; v < instance->num_vertices; v++) {
	//		vertex_prizes[v] -= min_prize;
	//	}		
	//}

	double total_zeta = 0;
	int num_roots = 0;
	// from zero to 2^num_vertices_pos_zeta
	for (int i = pow(2, num_vertices_pos_zeta) -1; i >= 0 ; i--) {
		// set the number of roots to zero
		num_roots = 0;
		
		// set the total zeta to zero
		total_zeta = 0;

		// computed adjusted cost of edge
		for (int e = 0; e < instance->num_edges; e++) {
			edge_costs[e] = instance->edges[e][2]; //- min_prize;
		}
		
		for (int j = 0; j < num_vertices_pos_zeta; j++) {
			if (checkbinSingle(i, j)) {
				total_zeta += zeta[vertcies_pos_zeta[j]];
				roots_pcst[num_roots++] = vertcies_pos_zeta[j];
			}
			else {
				// cost of all adjacent edges is very high
				for (int e = 0; e < instance->num_edges; e++) {
					if (instance->edges[e][0] == vertcies_pos_zeta[j] || instance->edges[e][1] == vertcies_pos_zeta[j]) {
						edge_costs[e] = 1000000;
					}
				}
			}
		}

		// compute vertex weight for each edge
		for (int e = 0; e < instance->num_edges; e++) {
			if (edge_costs[e] != 1000000) {
				edge_costs[e] *= total_zeta;
			}
		}

		this->instance->pcst
			->set_upper_bound(this->instance->UB)
			->reset_num_vertices_edges()
			->reset_active_status()
			->set_roots(num_roots, roots_pcst)
			->set_vertex_prizes(this->vertex_prizes)
			->set_edge_costs(this->edge_costs)
			//->print()
			->solve();

		this->aborescence
			->Init(this->instance->pcst)			
			->PrintModel()
			->Run();

	}

	return this;
}
