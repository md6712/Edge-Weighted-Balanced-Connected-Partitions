#include "CG.h"
#include <ilcplex/ilocplex.h>
#include "binary.h"
#include "CGCallBack.h"
#include "timer.h"



// constructor
CG::CG(_g* g, bool redirect, bool linear) :Cplex(g) {
	printCuts = false;
	printCycles = false;

	pricing_problem = new CGPrice(this->instance, false);

	// matrix total assignment 
	matrix_total_assignment = new double* [this->instance->num_vertices];
	for (int i = 0; i < this->instance->num_vertices; i++) {
		matrix_total_assignment[i] = new double[this->instance->num_vertices];
	}

	if (linear) {
		SetLinear();
	}
	else {
		SetInteger();
	}

	if (!redirect) {

		// Add Objective function
		AddObj();

		// Add Constraints
		AddCons();

		// Define Variables 
		AddVars(false);				
	}
}

// destructor
CG::~CG() {
	delete pricing_problem;

	// free variables
	for (int i = 0; i < this->instance->num_vertices; i++) {
		delete[] matrix_total_assignment[i];
	}
	delete[] matrix_total_assignment;
}


CG* CG::Run_BP() {
	// create a root node 
	root = new BP_node();	

	root->next = nullptr;
	root->lvl = 0;
	root->LB = 0;
	root->n_trees = instance->select_trees_for_CG.size();
	root->number = 1;

	// column generation timer
	Timer timer_CG;

	// set start time
	timer_CG.setStartTime();
	Run_CG(root, timer_CG);

	this->instance->_lb_root = root->LB; // set the lower bound of the root node

	// create a temp node --- only a pointer 
	BP_node* temp; 

	int number_of_nodes = 1; // number of nodes in the list	

	bool time_limit_reached = false;
	// check if timer has reached the time limit
	timer_CG.setEndTime();
	double elapsedTime = timer_CG.calcElaspedTime_sec();
	if (elapsedTime > TIME_LIMIT) {
		time_limit_reached = true;
	}

	while (!time_limit_reached && root != nullptr && root->LB + 0.001 < instance->UB - 1) {

		// copy root into two child nodes
		BP_node* node_L = new BP_node();
		BP_node* node_R = new BP_node();

		memcpy(node_L, root, sizeof(BP_node));
		memcpy(node_R, root, sizeof(BP_node));

		int u = root->branch_pair[0];
		int v = root->branch_pair[1];

		// set the branch rule u and v apart 
		node_L->branch[root->lvl].u = u;
		node_L->branch[root->lvl].v = v;
		node_L->branch[root->lvl].rule = CG_branch_rule::apart;
		node_L->lvl = root->lvl + 1;
		node_L->number = ++ number_of_nodes;

		// run column generation
		Run_CG(node_L, timer_CG);

		// check if timer has reached the time limit
		timer_CG.setEndTime();
		elapsedTime = timer_CG.calcElaspedTime_sec();
		if (elapsedTime > TIME_LIMIT) {
			break;
		}

		// check if the node LB is not less than the parent lower bound
		if (node_L->LB < root->LB - 0.001) {
			printf("ERROR: Node LB %f < Parent LB %f \n", node_L->LB, root->LB);
			break;
		}

		// add the node to the sorted linked list if the lower bound is less than the upper bound
		if (node_L->LB+0.001 > instance->UB) {
			delete node_L;					
		}
		else {
			AddNodeSorted(node_L);			
		}			

		// set the branch rule u and v together
		node_R->branch[root->lvl].u = u;
		node_R->branch[root->lvl].v = v;
		node_R->branch[root->lvl].rule = CG_branch_rule::together;
		node_R->lvl = root->lvl + 1;
		node_R->number = ++number_of_nodes;

		// run column generation
		Run_CG(node_R, timer_CG);
		// check if timer has reached the time limit
		timer_CG.setEndTime();
		elapsedTime = timer_CG.calcElaspedTime_sec();
		if (elapsedTime > TIME_LIMIT) {
			break;
		}

		// check if the node LB is not less than the parent lower bound
		if (node_R->LB < root->LB - 0.001) {		
			printf("ERROR: Node LB %f < Parent LB %f \n", node_R->LB, root->LB);
			break;
		}

		// add the node to the sorted linked list if the lower bound is less than the upper bound		
		if (node_R->LB+0.001 > instance->UB) {
			delete node_R;
		}
		else {
			AddNodeSorted(node_R);
		}

		// move to next node after root and delete root	
		temp = root;
		root = root->next;
		if (root == nullptr) {			
			printf("No more nodes in the list\n");
			root = temp;
			root->LB = instance->UB; // set the lower bound to upper bound
			break;
		}
		delete temp;
		PrintNodeList(root);

		// check if timer has reached the time limit
		timer_CG.setEndTime();
		elapsedTime = timer_CG.calcElaspedTime_sec();
		if (elapsedTime > TIME_LIMIT) {
			break;
		}
	} 
		
	instance->_opt = instance->UB; 
	if (root != nullptr) {
		instance->_gap = (instance->UB - root->LB+0.001) / instance->UB;
		delete root;
	}
	else {
		instance->_gap = 0;
	}

	return this;
}


// run
CG* CG::Run_CG(BP_node* node, Timer timer_CG) {
	this->instance->num_bp_nodes++;
	Impose_Braching(node);		

	Timer timer;

	double elapsedtimeMaster = 0;
	double elapsedtimePricing = 0;
	double elapsedtimePricingPCST = 0;

	int itr = 0;
	int last_upper_bound_update_num_trees = this->instance->select_trees_for_CG.size();
	int min_tree_added_trigger_upperbound = 20;

	while (true)
	{
		itr++;
		timer.setStartTime();
		Cplex::Run();
		timer.setEndTime();
		elapsedtimeMaster += timer.calcElaspedTime_sec();

		if (opt == -1) {
			// if the model is infeasible, we add artificial columns
			AddArtificialColumns(node);
			continue;
		}

		// print the solution
		//PrintSol();

		printf("\n ***** Master problem itr = %4d  ***** Elapsed  %4.2f %4.2f %4.2f \n", itr, elapsedtimeMaster, elapsedtimePricing, elapsedtimePricingPCST);

		// print the objective value + num_trees
		double obj = cplex.getObjValue();
		printf("obj = %3.2f \t #trees = %d\n", obj, this->instance->select_trees_for_CG.size());
		
		dualCover = IloNumArray(env, this->instance->num_vertices);
		dualZ = IloNumArray(env, this->instance->num_vertices);
		dualKnapsack = IloNum(0);

		// read duals values
		dualKnapsack = cplex.getDual(Knapsack);
		cplex.getDuals(dualCover, Cover);
		cplex.getDuals(dualZ, Z);

		//// print duals
		//printf_s("dualCover\t = ");
		//for (int v = 0; v < this->instance->num_vertices; v++) {
		//	printf_s("%3.2f  ", dualCover[v]);
		//}
		//printf_s("\n");
		//printf_s("dualZ \t\t = ");
		//for (int v = 0; v < this->instance->num_vertices; v++) {
		//	printf_s("%3.2f  ", dualZ[v]);
		//}
		//printf_s("\n");
		//printf_s("dualKnapsack\t = %3.2f\n", dualKnapsack);

		// set dual values in the pricing problem
		timer.setStartTime();
		pricing_problem
			->setTheta(dualKnapsack)
			->setEta(dualCover)
			->setZeta(dualZ)
			->UpdateObjectiveCoefficients()
			->SetQuiet();
		/*->PrintModel()
		->SetPrintCuts(true)
		->SetPrintCycles(true)*/

		// solve the sub-problem heuristically
		//_tree* tree = pricing_problem->heuristic();

		// print pricing optimal
		//printf("\n pricing_problem->opt = %f\n", pricing_problem->opt);

		// solve the sub-problem exact

		//Timer timer_pricing;
		//timer_pricing.setStartTime();
		//// compute time
		//pricing_problem->Run();
		//timer_pricing.setEndTime();
		//double elapsedtimePricingA = timer_pricing.calcElaspedTime_sec();
		//printf("\n ***** Pricing problem itr = %4d  ***** Elapsed  %4.2f \n", itr, elapsedtimePricingA);


		//_small_tree* tree = pricing_problem->GetTree(); // get the tree associated to the optimal solution		

		//// print the solution
		//pricing_problem->PrintSol();

		//// print the tree
		//tree->print_vertices(this->instance);			

		timer.setEndTime();
		elapsedtimePricing += timer.calcElaspedTime_sec();


		timer.setStartTime();
		pricing_problem->solve_pcst(node, false);
		timer.setEndTime();
		elapsedtimePricingPCST += timer.calcElaspedTime_sec();


		// if there is any tree to be added, add them 
		if (pricing_problem->num_trees_to_add > 0 && pricing_problem->best_pcst >= 0.001) {
			for (int i = 0; i < pricing_problem->num_trees_to_add; i++) {
				// get the tree
				_small_tree* tree = pricing_problem->trees_to_add[i];

				// add the tree to the instance
				this->instance->select_trees_for_CG.push_back(tree);

				// print the tree
				//tree->print_vertices(this->instance);

				// add the variable to the model
				AddVar(this->instance->select_trees_for_CG.size() - 1, false);
			}

			// if number of trees added is larger than the last upper bound update, update the upper bound
			if (this->instance->select_trees_for_CG.size() - last_upper_bound_update_num_trees > min_tree_added_trigger_upperbound) {
				
				// if lower_bound is less than UB
				if (obj < instance->UB - 0.999) {

					printf("\n ***** Compute UB: ");

					int previous_UB = instance->UB; // save previous upper bound
					UpdateUB(); // update upper bound				
					printf("%d \n",instance->UB);
					// if upper bound is not updated
					if (instance->UB < previous_UB) {
						printf(" ***** Update UB: %d -> %d\n", previous_UB, instance->UB);
						min_tree_added_trigger_upperbound = 20; // update the trigger for upper bound update
					}
					else {
						min_tree_added_trigger_upperbound = min(50, min_tree_added_trigger_upperbound+20); // update the trigger for upper bound update
					}

					last_upper_bound_update_num_trees = this->instance->select_trees_for_CG.size();
				}
			}
		}
		else {
			// if there is no tree to be added, break the loop
			break;
		}		

		//// print all elapsed time
		//printf("elapsedtimeMaster = %f\n", elapsedtimeMaster);
		//printf("elapsedtimePricing = %f\n", elapsedtimePricing);
		//printf("elapsedtimePricingPCST = %f\n", elapsedtimePricingPCST);


		// check if timer has reached the time limit
		timer_CG.setEndTime();
		double elapsedTime = timer_CG.calcElaspedTime_sec();
		if (elapsedTime > TIME_LIMIT) {
			return this;
		}
	}

	//// print elapsed time
	//printf("elapsedtimeMaster = %f\n", elapsedtimeMaster);
	//printf("elapsedtimePricing = %f\n", elapsedtimePricing);
	//printf("elapsedtimePricingPCST = %f\n", elapsedtimePricingPCST);

	// compute the total assignment
	ComputeTotalAssignment(node);

	//// print the total assignment
	PrintTotalAssignment();

	// lower bound of the node is equal the optimal value of the master
	node->LB = cplex.getObjValue();
	
	if (node->LB < instance->UB) {		
		int previous_UB = instance->UB; // save previous upper bound
		UpdateUB();
		// if upper bound is not updated
		if (instance->UB < previous_UB) {
			printf("\n ***** Update UB: %d -> %d\n", previous_UB, instance->UB);
		}
	}

	// if lower bound is less than upper bound - 1 branch

	if (node->LB + 0.001 < instance->UB - 1) {
		// if the total assignment of the best pair is not integer

		double total_assignment = matrix_total_assignment[best_pair[0]][best_pair[1]];

		if (total_assignment > 0.001 && total_assignment < 0.999) {
			node->branch_pair[0] = best_pair[0];
			node->branch_pair[1] = best_pair[1];
		}		
		else {
			// integer solution has found
			node->branch_pair[0] = -1;
			node->branch_pair[1] = -1;

			// update UB
			instance->UB = node->LB;
		}
	}

	// set the number of trees in the node
	node->n_trees = instance->select_trees_for_CG.size();

	// print node 
	PrintNode(node);

	return this;
}

void CG::AddArtificialColumns(BP_node* node) {

	PrintNode(node);
	
	// create a small tree including all vertices as far as rules allow
	_small_tree** small_tree = new _small_tree*[instance->num_trees];

	int vertices_to_include[5];
	bool* vertex_selected = new bool[this->instance->num_vertices];
	memset(vertex_selected, 0, sizeof(bool) * this->instance->num_vertices);
	
	int current_tree = 0; // current tree index

	while (current_tree < instance->num_trees){
		small_tree[current_tree] = new _small_tree();
		for (int v1 = 0; v1 < this->instance->num_vertices; v1++) {

			if (vertex_selected[v1]) continue;

			int n_vertices_to_include = 0;

			// check if including the vertex doesnt violate all // the branching rules
			bool include_vertex = true;
			for (int i = 0; i < node->lvl; i++) {
				// branch u an v
				int u = node->branch[i].u;
				int v = node->branch[i].v;

				if (u == v1 || v == v1) {
					int u1 = (u == v1) ? v : u; // get the other vertex

					if (node->branch[i].rule == CG_branch_rule::apart) {
						if (checkbin(small_tree[current_tree]->bin_vertices, u1)) {
							// if the rule is apart, we can include the vertex only if the other vertex is not included
							include_vertex = false;
							break;
						}
					}
					else if (node->branch[i].rule == CG_branch_rule::together) {
						// if the rule is together, we can include the vertex only if the other vertex is included
						vertices_to_include[n_vertices_to_include++] = u1; // add the other vertex to the list of vertices to include
					}
				}
			}

			if (include_vertex) {
				// if the vertex can be included, add it to the small tree
				addbin(small_tree[current_tree]->bin_vertices, v1);
				vertex_selected[v1] = true; // mark the vertex as selected

				// add each vertex to include
				for (int i = 0; i < n_vertices_to_include; i++) {
					addbin(small_tree[current_tree]->bin_vertices, vertices_to_include[i]);
					vertex_selected[vertices_to_include[i]] = true; // mark the vertex as selected
				}
			}
		}
		
		current_tree++;
	}



	// print trees
	for (int i = 0; i < current_tree; i++) {
		small_tree[i]->weight = 1000; // reset the weight		
		small_tree[i]->print_vertices(this->instance);

		//add to the select_trees_for_CG
		this->instance->select_trees_for_CG.push_back(small_tree[i]);
		AddVar(this->instance->select_trees_for_CG.size() - 1, false); // add the variable to the model
	}
}

// update upper bound
void CG::UpdateUB() {
	// remove all variables
	RemoveXVars();

	// each tree add x as an integer variable
	for (int i = 0; i < this->instance->select_trees_for_CG.size(); i++) {
		AddVar(i, true);
	}
	model.add(x);

	// solve the master with integer variables	
	Cplex::Run();	

	int UB = (int)(cplex.getObjValue() + 0.001);

	if (UB < instance->UB) {
		instance->UB = UB;
	}

	// reverse the variables
	RemoveXVars();

	// remove all trees that have weight larger than UB
	for (int i = 0; i < this->instance->select_trees_for_CG.size();) {
		if (this->instance->select_trees_for_CG[i]->weight > instance->UB) {
			_small_tree* tree = this->instance->select_trees_for_CG[i]; // get the tree
			this->instance->select_trees_for_CG.erase(this->instance->select_trees_for_CG.begin() + i); // remove the tree from the list
			delete tree; // delete the tree
		}
		else {
			i++;
		}
	}

	for (int i = 0; i < this->instance->select_trees_for_CG.size(); i++) {
		AddVar(i, false);
	}
	model.add(x);

	// reset cplex

	cplex.end();
	cplex = IloCplex(model);
}

// print model
CG* CG::PrintModel() {
	cplex.exportModel("ModelCG.lp");
	return this;
}

// print solution
CG* CG::PrintSol() {
	// loop over all variables
	for (int i = 0; i < x.getSize(); i++) {
		if (cplex.getValue(x[i]) > 0.00001) {
			printf("x[%d] = %1.3f \t", i, cplex.getValue(x[i]));
			this->instance->select_trees_for_CG[i]->print_vertices(this->instance);
		}
	}

	// print the objective value
	printf("obj = %3.2f\n", cplex.getObjValue());

	// print number of trees
	printf("Number of trees: %d\n", this->instance->select_trees_for_CG.size());

	return this;
}

// add objective function
void CG::AddObj() {
	// set the objective function	
	obj = IloObjective(env, 0, IloObjective::Minimize, "obj");
	model.add(obj);
}

// add constraints
void CG::AddCons() {
	// Create arrays	
	Cover = IloRangeArray(env, this->instance->num_vertices);
	Z = IloRangeArray(env, this->instance->num_vertices);

	// set bounds 
	Knapsack = IloRange(env, -this->instance->num_trees, IloInfinity);
	Knapsack.setName("Knapsack");
	for (int v = 0; v < this->instance->num_vertices; v++) {

		// name cover
		sprintf(name, "Cover_%d", v);
		Cover[v] = IloRange(env, 1, 1); 
		Cover[v].setName(name);

		// name Z
		sprintf(name, "Z_%d", v);
		Z[v] = IloRange(env, 0, IloInfinity); 
		Z[v].setName(name);
	}
	
	// add to model
	model.add(Knapsack);
	model.add(Cover);
	model.add(Z);
}

// add variables	
void CG::AddVars(bool integer) {
	// instance g
	_g* g = instance;	

	// define the variable Z
	sprintf(name, "z");
	z = IloNumVar(env, 0, IloInfinity, ILOFLOAT, name);
	
	obj.setLinearCoef(z, 1);

	for (int v = 0; v < g->num_vertices; v++) {
					
		Z[v].setLinearCoef(z, 1);						
	}

	// define the variable X
	x = IloNumVarArray(env);

	for (int i = 0; i < g->select_trees_for_CG.size(); i++) {
		// name 
		AddVar( i, integer);
	}

	// add variables to the model
	model.add(z);
	model.add(x);	
}

// add variable
void CG::AddVar(int i, bool integer = false) {
	// instance g
	_g* g = instance;

	// add variables
//	sprintf(name, "x_%d", i);
IloNumVar Xvar(env, 0, IloInfinity, integer ? ILOINT : ILOFLOAT);

	obj.setLinearCoef(Xvar, 0);
	Knapsack.setLinearCoef(Xvar, -1);

	for (int v = 0; v < g->num_vertices; v++) {
		if (checkbin(g->select_trees_for_CG[i]->bin_vertices, v)) {
			Cover[v].setLinearCoef(Xvar, 1);
		}
		else {
			Cover[v].setLinearCoef(Xvar, 0);
		}

		if (checkbin(g->select_trees_for_CG[i]->bin_vertices, v)) {
			Z[v].setLinearCoef(Xvar, -g->select_trees_for_CG[i]->weight);
		}
		else {
			Z[v].setLinearCoef(Xvar, 0);
		}
	}

	x.add(Xvar);
}

// remove vars 
void CG::RemoveXVars() {
	// remove variables from the model
	for (int i = 0; i < x.getSize(); i++) {
		model.remove(x[i]);
		x[i].end();
	}
	x.clear();
}

// impose branching
CG* CG::Impose_Braching(BP_node* node) {
	// reverse all previous branching rules imposed
	// check all trees and put the bounds to 0 and 1
	for (int i = 0; i < instance->select_trees_for_CG.size(); i++) {
		x[i].setLB(0);
		x[i].setUB(1);
	}

	// impose new rules
	// for each rule in node
	for (int r = 0; r < node->lvl; r++) {
		int u = node->branch[r].u;
		int v = node->branch[r].v;
		CG_branch_rule rule = node->branch[r].rule;

		// check each tree
		for (int i = 0; i < instance->select_trees_for_CG.size(); i++) {
			// get the tree
			_small_tree* tree = instance->select_trees_for_CG[i];

			// check if the tree contains the edge (u,v)
			if (rule == CG_branch_rule::apart) {
				// see if it contains the edge (u,v), if yes, then change the bounds of the variable and put it equal to zero			
				if (checkbin(tree->bin_vertices, u) && checkbin(tree->bin_vertices, v)) {
					x[i].setLB(0);
					x[i].setUB(0);
				}
			}
			else {
				// see if the tree contains one of the vertices but not both, then set the bounds of the variable to 0
				if (checkbin(tree->bin_vertices, u) && !checkbin(tree->bin_vertices, v)) {
					x[i].setLB(0);
					x[i].setUB(0);
			
				}
				else if (!checkbin(tree->bin_vertices, u) && checkbin(tree->bin_vertices, v)) {
					x[i].setLB(0);
					x[i].setUB(0);
				}
			}
		}
	}
	
	// print bounds for variable x
	for (int i = 0; i < x.getSize(); ++i) {
		double lb = x[i].getLB(); // ✅ Valid
		double ub = x[i].getUB(); // ✅ Valid
		//printf("x[%d] = (%.2f, %.2f)\n", i, lb, ub);
	}



	return this;
}


// compute total assignment
void CG::ComputeTotalAssignment(BP_node* node) {
	// instance g
	_g* g = instance;
	// compute the total assignment

	// set best pair
	best_pair[0] = 0;
	best_pair[1] = 1;

	// reset total assignments
	for (int i = 0; i < g->num_vertices; i++) {
		for (int j = 0; j < g->num_vertices; j++) {
			matrix_total_assignment[i][j] = 0;
		}
	}

	for (int i = 0; i < g->num_vertices; i++) {
		for (int j = i+1; j < g->num_vertices; j++) {
			matrix_total_assignment[i][j] = 0;

			// go over all trees that are in the solution and include i and j and sum the x_values
			for (int k = 0; k < g->select_trees_for_CG.size(); k++) {
				// check if the tree is in the solution
				if (cplex.getValue(x[k]) > 0.00001) {
					// check if i and j are in the tree
					if (checkbin(g->select_trees_for_CG[k]->bin_vertices, i) && checkbin(g->select_trees_for_CG[k]->bin_vertices, j)) {
						matrix_total_assignment[i][j] += cplex.getValue(x[k]);
					}
				}
			}

			// if distance of the total assignment from 0.5 is less than 0.5, assign the pair to best pair. 
			
			if (fabs(matrix_total_assignment[i][j] - 0.5) < fabs(matrix_total_assignment[best_pair[0]][best_pair[1]] - 0.5)) {
				// check if it is not part of previous rules

				bool is_part_of_previous_rules = false;
				for (int r = 0; r < node->lvl; r++) {
					if ((node->branch[r].u == i && node->branch[r].v == j) || (node->branch[r].u == j && node->branch[r].v == i)) {
						is_part_of_previous_rules = true;
						break;
					}
				}

				if (!is_part_of_previous_rules) {
					best_pair[0] = i;
					best_pair[1] = j;
				}				
			}
		}
	}	
}

// print node
void CG::PrintNode(BP_node* node) {
	// print the node
	printf("Node %d: Lvl %d \t LB = %3.2f \t #trees = %d \n", node->number, node->lvl, node->LB, node->n_trees);
	for (int i = 0; i < node->lvl; i++) {
		printf("branch[%d] = (%d,%d) \t rule = %d\n", i, node->branch[i].u, node->branch[i].v, node->branch[i].rule);
	}
	printf("***************\n\n");
}

// print total assignment
void CG::PrintTotalAssignment() {
	// instance g
	_g* g = instance;
	// print the total assignment
	for (int i = 0; i < g->num_vertices; i++) {
		for (int j = i+1; j < g->num_vertices; j++) {
			printf("%3.2f ", matrix_total_assignment[i][j]);
		}
		printf("\n");
	}

	// print best pair
	printf("best pair = (%d,%d) \t total assignment = %3.2f\n", best_pair[0], best_pair[1], matrix_total_assignment[best_pair[0]][best_pair[1]]);
}

// add node to linked list
void CG::AddNodeSorted(BP_node* node) {
	// add node to the list
	if (root == nullptr) {
		root = node;
	}
	else {
		BP_node* current = root;
		BP_node* previous = nullptr;
		while (current != nullptr && node->LB+0.001>= current->LB) {
			previous = current;
			current = current->next;
		}
		if (previous == nullptr) {
			node->next = root;
			root = node;
		}
		else {
			node->next = current;
			previous->next = node;
		}
	}
}

void CG::PrintNodeList(BP_node* root) {
	// print the node list
	BP_node* current = root;
	printf("NodeList \n");
	while (current != nullptr) {
		printf("Node %d: %d \t LB = %3.2f \n", current->number, current->lvl, current->LB);
		current = current->next;
	}
}