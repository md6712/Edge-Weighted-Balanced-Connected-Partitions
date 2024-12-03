#include "CGPrice.h"
#include "Cplex.h"
#include "binary.h"	

ILOLAZYCONSTRAINTCALLBACK1(callback_CG_Price, CGPrice*, cg_price) {
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


// constructor
CGPrice::CGPrice(_g* g, bool redirect):Cplex(g) {
	printCuts = false;
	printCycles = false;

	// define coef 
	theta = IloNum(0);
	eta = IloNumArray(env, this->instance->num_vertices);
	zeta = IloNumArray(env, this->instance->num_vertices);

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
CGPrice::~CGPrice() {
	eta.end();
	zeta.end();
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

	// set the callback
	cplex.use(callback_CG_Price(env, this));

	// run the model
	Cplex::Run();
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

	/*model.add(x[13] == 0);
	model.add(x[11] == 0);
	model.add(x[12] == 0);*/
}

// define the constraints
void CGPrice::AddCons() {
	AddConsNumEdgesVerticesMinusOne();
	XYRelation();
	PhiRelation();
	AtleastOneEdge();

	// cost of tree extra constraint
	//AddConsCostOfTree();
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

	model.add(sum_x == sum_y - 1);
	sum_x.end();
}

// add constraints that relate x and y variables
void CGPrice::XYRelation() {
	for (int e = 0; e < this->instance->num_edges; e++) {		
		model.add(x[e] <= y[this->instance->edges[e][0]]);
		model.add(x[e] <= y[this->instance->edges[e][1]]);		
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

	objective.setExpr	(updatedExpr);
	return this;
}


// get the tree associated to the optimal solution
_tree* CGPrice::GetTree(int *vertices, uint32_t* bin_vertices) {

	// set vertices to zero
	memset (vertices, 0, sizeof(int) * this->instance->num_vertices);
	memset (bin_vertices, 0, sizeof(uint32_t) * binaryArrlength(this->instance->num_vertices));

	// init num_vertices
	int num_vertices = 0;
	int num_edges = 0;

	for (int v = 0; v < this->instance->num_vertices; v++) {
		if (cplex.getValue(y[v]) > 0) {			
			addbin(bin_vertices, v);
			vertices[num_vertices++] = v;						
		}
	}

	_tree* tree = new _tree(this->instance);
	tree->CopyVertices(num_vertices, vertices, bin_vertices);

	tree->weight = 0;

	// loop over all edges
	for (int e = 0; e < this->instance->num_edges; e++) {
		if (cplex.getValue(x[e]) > 0) {
			tree->edges[num_edges++] = e;
			tree->weight += this->instance->edges[e][2];
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
				int v = tree->IncidentVertex(e);
				if (v == -1) {
					continue;
				}

				// compute the other vertex 
				int u = this->instance->edges[e][0] == v ? this->instance->edges[e][1] : this->instance->edges[e][0];

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
			//tree->printTree	();
			reward = max_reward;
			
			best_edge = -1;
		}
	}

	// local search

	bool improved = true; 

	while (improved) {

		improved = false;
		 
		// for all edges that are in the tree
		for (int e = 0; e < this->instance->num_edges; e++) {
			if (!tree->IsEdgeInTree(e)) {
				continue;
			}

			// get the vertex that is connecting the tree to the edge
			int v = tree->IncidentVertex(e);
			if (v == -1) {
				continue;
			}

			// for all edges that are not in the tree
			for (int ep = 0; ep < this->instance->num_edges; ep++) {
				if (tree->IsEdgeInTree(ep)) {
					continue;
				}

				// get the vertex that is connecting the tree to the edge
				int vp = tree->IncidentVertex(ep);
				if (vp == -1) {
					continue;
				}				

				
				
			}

			
		}

	}



	double total_reward = -theta + reward;

	this->opt = total_reward;
	
	tree->reduced_cost = total_reward;	

	return tree;
}
