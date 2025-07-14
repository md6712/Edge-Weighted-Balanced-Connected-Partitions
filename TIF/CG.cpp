#include "CG.h"
#include <ilcplex/ilocplex.h>
#include "binary.h"
#include "CGCallBack.h"
#include "timer.h"



// constructor
CG::CG(_g* g, bool redirect, bool linear) :Cplex(g) {
	int num_trees = this->instance->select_trees_for_CG.size();
	sorted_indices = new int[num_trees + 10000];

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

	delete[] sorted_indices;

	premature_nodes.clear(); // clear the premature nodes vector
}


CG* CG::Run_BP() {


	// recompute the upper bound using the integer master problem
	//UpdateUB();
	
	// create a root node 
	BP_node* root = new BP_node();



	root->lvl = 0;
	root->LB = 0;
	root->n_trees = instance->select_trees_for_CG.size();
	root->number = 1;
	root->premature = false;

	// column generation timer
	Timer timer_CG;

	// set start time
	timer_CG.setStartTime();
	Run_CG(root, timer_CG);

	// add the node
	AddNodeSorted(root);
	
	if (!root->premature)
		this->instance->_lb_root = root->LB; // set the lower bound of the root node

	int number_of_nodes = 1; // number of nodes in the list	

	bool time_limit_reached = false;
	// check if timer has reached the time limit
	timer_CG.setEndTime();
	double elapsedTime = timer_CG.calcElaspedTime_sec();
	if (elapsedTime > TIME_LIMIT) {
		time_limit_reached = true;
	}

	// temp node
	BP_node* first_node = root; // save the root node to be removed later

	while (
			!time_limit_reached 
			&& root != nullptr 
			&& (root->premature || root->LB + 0.001 < instance->UB - 1)
			) 
	{

		// check if there is any premature nodes
		if (premature_nodes.size() > 0) {
			// get the first premature node that is not completed
			for (int i = 0; i < premature_nodes.size(); i++) {
				if (!premature_nodes[i]->premature_branched) {
					root = premature_nodes[i];				
					root->premature_branched = true; // mark the node as branched
					break;
				}
			}
		}
		else if (active_nodes.size() > 0) {
			root = active_nodes[0]; // get the first node in the active nodes vector
			active_nodes.erase(active_nodes.begin()); // remove the first node from the active nodes vector
			nodes_to_remove.push_back(root); // add the node to the nodes to be removed vector		
		}
		else {
			printf("No more nodes in the list\n");
			break;
		}

		// copy root into two child nodes
		BP_node* node_L = new BP_node();
		BP_node* node_R = new BP_node();

		// set them as siblings of each other		
		root->child_L = node_L;
		root->child_R = node_R;

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
		node_L->premature = false; 		
		node_L->premature_branched = false; // mark the node as not branched yet
		node_L->child_L = nullptr; // set the child nodes to nullptr
		node_L->child_R = nullptr; // set the child nodes to nullptr
		node_L->parent = root;		
		node_L->sibling = node_R;
		


		// run column generation
		Run_CG(node_L, timer_CG);

		// check if timer has reached the time limit
		timer_CG.setEndTime();
		elapsedTime = timer_CG.calcElaspedTime_sec();
		if (elapsedTime > TIME_LIMIT) {
			break;
		}

		// add the node to the sorted linked list if the lower bound is less than the upper bound
		AddNodeSorted(node_L);

		// set the branch rule u and v together
		node_R->branch[root->lvl].u = u;
		node_R->branch[root->lvl].v = v;
		node_R->branch[root->lvl].rule = CG_branch_rule::together;
		node_R->lvl = root->lvl + 1;
		node_R->number = ++number_of_nodes;
		node_R->premature = false;		
		node_R->premature_branched = false; // mark the node as not branched yet
		node_R->child_L = nullptr; // set the child nodes to nullptr
		node_R->child_R = nullptr; // set the child nodes to nullptr
		node_R->parent = root;
		node_R->sibling = node_L;

		// run column generation
		Run_CG(node_R, timer_CG);
		// check if timer has reached the time limit
		timer_CG.setEndTime();
		elapsedTime = timer_CG.calcElaspedTime_sec();
		if (elapsedTime > TIME_LIMIT) {
			break;
		}		

		// add node
		AddNodeSorted(node_R);		

		PrintNodeList();

		// check if timer has reached the time limit
		timer_CG.setEndTime();
		elapsedTime = timer_CG.calcElaspedTime_sec();
		if (elapsedTime > TIME_LIMIT) {
			break;
		}
	} 
	
	root = first_node; // restore the root node

	UpdateConfirmedLB(root);	

	PrintNodeList();

	instance->_opt = instance->UB; 	
	instance->_gap = ((double)(instance->UB - root->CLB)) / ((double)instance->UB);
	delete root;

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
	
	int last_upper_bound_update_itr = 0;	 // this is the iteration the UpdateUB is called
	int last_upper_bound_improve_itr = 0;  // this is the iteration the upper bound is improved
	timer_CG.setEndTime();	
	double time_last_upper_bound_update = timer_CG.calcElaspedTime_sec();

	// premature termination of CG
	bool premature_termination = false;

	// array to keep track of last 5 lower bounds 
	double last_5_lower_bounds[5] = { 0, 0, 0, 0, 0 };
	double last_5_pcst_values[5] = { 0, 0, 0, 0, 0 };

	node->completed = false;

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

		// update last 5 lower bounds
		for (int i = 4; i > 0; i--) {
			last_5_lower_bounds[i] = last_5_lower_bounds[i - 1];
		}
		last_5_lower_bounds[0] = cplex.getObjValue();		

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

		//int B = 0;
		////printf_s("dualZ \t\t = ");
		//for (int v = 0; v < this->instance->num_vertices; v++) {
		////	printf_s("%3.2f  ", dualZ[v]);
		//	if (dualZ[v] > 0.001) {
		//		B++;
		//	}
		//}
		//printf_s("\n");
		////printf_s("dualKnapsack\t = %3.2f\n", dualKnapsack);

		//printf_s("Number of vertices with positive dualZ: %d\n", B);

		//// get reduced costs
		//IloNumArray reducedCosts(env);
		//cplex.getReducedCosts(reducedCosts, x);

		//if (this->instance->select_trees_for_CG.size() > 1000 && base_purify_tries < 20) {
		//	if (B > 8 || (itr > 20 && B > 7)) {
		//		// remove variables with large positive reduced cost
		//		base_purify_tries++;
		//		RemoveLargePositiveReducedCost(reducedCosts);
		//		itr--;
		//		continue;
		//	}
		//}
		//base_purify_tries = 0;
						
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

		

		//timer.setStartTime();
		//pricing_problem->Run();	

		//printf("pricing_problem->opt = %f\n", pricing_problem->opt);

		//_small_tree* tree = pricing_problem->GetTree(); // get the tree associated to the optimal solution			

		////// print the tree
		//tree->print_vertices(this->instance);			

		//timer.setEndTime();
		//elapsedtimePricing += timer.calcElaspedTime_sec();


		timer.setStartTime();
		pricing_problem->solve_pcst(node, false);
		timer.setEndTime();
		elapsedtimePricingPCST += timer.calcElaspedTime_sec();


		// salve last five pcst values
		for (int i = 4; i > 0; i--) {
			last_5_pcst_values[i] = last_5_pcst_values[i - 1];
		}
		last_5_pcst_values[0] = pricing_problem->best_pcst;


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

			// check if timer has reached the time limit
			timer_CG.setEndTime();
			double elapsedTime = timer_CG.calcElaspedTime_sec();			

			// if number of trees added is larger than the last upper bound update, update the upper bound
			if (
				this->instance->select_trees_for_CG.size() - last_upper_bound_update_num_trees > min_tree_added_trigger_upperbound
				|| (itr - last_upper_bound_update_itr > 10) // if more than 10 iterations have passed since the last upper bound update and the number of trees is larger than 1000
				|| (elapsedTime - time_last_upper_bound_update > 120) // if more than 60 seconds have passed since the last upper bound update
			) {				
				// if lower_bound is less than UB
				if (obj < instance->UB - 0.999 && node->lvl == 0) {	
					int previous_UB = instance->UB; // save previous upper bound
					UpdateUB(); // update upper bound				
				
					// if upper bound is not updated
					if (instance->UB < previous_UB) {
						min_tree_added_trigger_upperbound = 20; // update the trigger for upper bound update
						last_upper_bound_improve_itr = itr; // update the last upper bound improve iteration
					}
					else {
						min_tree_added_trigger_upperbound = min(50, min_tree_added_trigger_upperbound+20); // update the trigger for upper bound update
					}

					last_upper_bound_update_num_trees = this->instance->select_trees_for_CG.size();
					last_upper_bound_update_itr = itr; // update the last upper bound update iteration
					time_last_upper_bound_update = elapsedTime; // update the last upper bound update time

					// we need to reimpose the branching rules. 
					//Impose_Braching(node);
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


		// print last 5 lower bounds and last 5 pcst values
		printf("Last 5 lower bounds: ");
		for (int i = 0; i < 5; i++) {
			printf("%3.2f ", last_5_lower_bounds[i]);
		}
		printf("\n");
		printf("Last 5 pcst values: ");
		for (int i = 0; i < 5; i++) {
			printf("%3.2f ", last_5_pcst_values[i]);
		}
		printf("\n");


		// if the improvement is less than 0.001, we terminate the CG
		//if (node->lvl < 2)
		//if ((itr - last_upper_bound_improve_itr) > 5 && last_5_lower_bounds[4] - last_5_lower_bounds[0] < 0.01) {
		//	// if pcst is not larger than 1 for the last five iterations, we terminate the CG
		//	// compute maximum pcst 
		//	double max_pcst = 0;
		//	for (int i = 0; i < 5; i++) {
		//		max_pcst = max(max_pcst, last_5_pcst_values[i]);
		//	}
		//	if (max_pcst < 1.0) {
		//		premature_termination = true;				
		//		printf("Premature termination of CG due to no improvement in the last 5 iterations.\n");
		//		break;
		//	}			
		//}
	}	

	node->completed = true;

	// Run ones more
	Cplex::Run();

	//// print elapsed time
	//printf("elapsedtimeMaster = %f\n", elapsedtimeMaster);
	//printf("elapsedtimePricing = %f\n", elapsedtimePricing);
	//printf("elapsedtimePricingPCST = %f\n", elapsedtimePricingPCST);

	// compute the total assignment
	ComputeTotalAssignment(node);

	//// print the total assignment
	//PrintTotalAssignment();

	

	// set if node is premature
	node->premature = premature_termination;

	// lower bound of the node is equal the optimal value of the master
	if (!node->premature)
		node->PLB = node->LB = cplex.getObjValue();
	else {
		node->PLB = cplex.getObjValue(); // save the premature lower bound
		node->LB = 0; // if the node is premature, set the lower bound to zero
	}
		

	// update UB if needed
	if (node->LB < instance->UB || node->premature) {		
		UpdateUB();		
	}

	// if lower bound is less than upper bound - 1 branch
	if (node->LB + 0.001 < instance->UB - 1) {
		// if the total assignment of the best pair is not integer

		double total_assignment = matrix_total_assignment[best_pair[0]][best_pair[1]];

		if (node->premature || (total_assignment > 0.001 && total_assignment < 0.999)) {
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

	UpdatePrematureParent(node); // update premature nodes

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
	printf("\n ***** Update UB *****\n");

	// Timer	
	Timer timer;
	timer.setStartTime();


	// reset x to integer variables
	RemoveXVars();
	AddXVars(true); // add x variables as integer variables

	// solve the master with integer variables	
	Cplex::Run();	
	int UB = (int)(cplex.getObjValue() + 0.001);

	if (UB < instance->UB) {
		int previous_UB = instance->UB; // save previous upper bound	
		instance->UB = UB;
		printf("UB UPDATED : %d -> %d\n", previous_UB, instance->UB);
	}

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

	// reset x as float
	AddXVars(false); // add x variables as integer variables

	// reset cplex
	cplex.end();
	cplex = IloCplex(model);

	timer.setEndTime();
	double elapsedTime = timer.calcElaspedTime_sec();
	printf("Time elapsed UPDATE UB: %4.3f\n", elapsedTime);
}

// update recursive prematurity
void CG::UpdatePrematureParent(BP_node* node) {
	// ignore if node itself is premature
	if (node->premature) {
		return;
	}
	
	// check if node has a premature parent
	if (node->parent != nullptr && node->parent->premature) {

		// since the node is mature, increase the premature nodes counter
		node->parent->n_childs_closed++;	

		// if both children of the parent is closed, parent becomes mature
		if (node->parent->n_childs_closed == 2) {
			node->parent->premature = false; // set the parent as premature
			node->parent->LB = min(node->parent->LB, node->LB);
			
			// recursive call for parent
			UpdatePrematureParent(node->parent);

			// remove parent from the premature nodes vector
			for (auto it = premature_nodes.begin(); it != premature_nodes.end(); ++it) {
				if (*it == node->parent) {
					premature_nodes.erase(it);
					break;
				}
			}


			// add the node to nodes to be removed
			nodes_to_remove.push_back(node->parent); // add the parent to the nodes to be removed vector
			// sort the nodes to be removed vector
			std::sort(nodes_to_remove.begin(), nodes_to_remove.end(), [](BP_node* a, BP_node* b) {
				return a->LB < b->LB; // sort by lower bound
				});
		}		
	}

}

// update confirmed lower bound
void CG::UpdateConfirmedLB(BP_node* node) {
	
	node->CLB = 0;

	if (!node->completed) {
		return;
	}

	if (node->premature) {
		return;
	}

	bool left_complete; 

	node->CLB = ceil(node->LB); // set the confirmed lower bound to the lower bound of the node

	if (node->child_L != nullptr && node->child_L->completed) {
		if (node->child_R != nullptr && node->child_R->completed) {
			UpdateConfirmedLB(node->child_L); // recursive call for left child
			UpdateConfirmedLB(node->child_R); // recursive call for left child
			if (node->child_L->CLB != 0 && node->child_R->CLB != 0) {
				node->CLB = min(node->child_L->CLB, node->child_R->CLB);
			}
		}
	}
}

// RemoveLargePositiveReducedCost 
void CG::RemoveLargePositiveReducedCost(IloNumArray reducedCosts) {
	
	// default number of trees to be removed
	int num_trees_to_remove = 1; // default number of trees to be removed

	// make a array of size equal to the number of trees
	int num_trees = this->instance->select_trees_for_CG.size();

	// sort the indices of the reduced costs in descending order
	for (int i = 0; i < num_trees; i++) {
		sorted_indices[i] = i;
	}

	// sort the indices based on the reduced costs in descending order
	std::sort(sorted_indices, sorted_indices + num_trees, [&](int a, int b) {
		return reducedCosts[a] > reducedCosts[b];
		});

	// sort the first 100 indices in descending order of the index
	std::sort(sorted_indices, sorted_indices + min(num_trees, num_trees_to_remove), [&](int a, int b) {
		return a > b; // sort in descending order of indices
		});

	RemoveXVars();

	// remove the first 100 trees with largest reduced costs from cg_trees 
	for (int i = 0; i < num_trees && i < num_trees_to_remove; i++) {
		int index = sorted_indices[i];
		_small_tree* tree = this->instance->select_trees_for_CG[index]; // get the tree
		if (reducedCosts[index] > 0.0001) { // if the reduced cost is positive
			this->instance->select_trees_for_CG.erase(this->instance->select_trees_for_CG.begin() + index); // remove the tree from the list
		}
	}

	AddXVars(false); // add x variables as float variables
}

// print model
CG* CG::PrintModel() {

	// get current time tick
	long long time_tick = GetCurrentTime();

	// create a string with modelcg-<time_tick>.lp
	sprintf(name, "modelcg-%lld.lp", time_tick);
			
	cplex.exportModel(name);
	return this;
}

// print solution
CG* CG::PrintSol() {

	printf("\n ***** Print MP Solution *****\n");

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

	model.add(z);

	AddXVars(integer); // add x variables
}

// add variable
void CG::AddVar(int i, bool integer = false) {
	// instance g
	_g* g = instance;

	// define xVar
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

	// add it to the x variable array
	x.add(Xvar);
}

// add variable
void CG::AddXVars(bool integer /* = false */)
{
	_g* g = instance;
	const int nCols = g->select_trees_for_CG.size();
	const int nVerts = g->num_vertices;	

	/* ---------- 1. Create all master-problem columns in one shot ---------- */
	x = IloNumVarArray(env, nCols, 0.0, 1.0, integer ? ILOINT : ILOFLOAT);
	model.add(x); // add x variables to the model

	/* ---------- 2. Bulk-set Knapsack row (all −1) ---------- */	
	IloNumArray knCoeff(env, nCols);
	for (int i = 0; i < nCols; ++i)
		knCoeff[i] = -1.0;
	Knapsack.setLinearCoefs(x, knCoeff);
	knCoeff.end();
	

	/* ---------- 3. Build Cover[v] and Z[v] rows in bulk ---------- */
	IloNumArray coverCoeff(env, nCols);                 // reusable buffer
	IloNumArray zCoeff(env, nCols);

	for (int v = 0; v < nVerts; ++v)
	{
		/* fill the two buffers for vertex v */
		for (int i = 0; i < nCols; ++i)
		{
			const auto* tree = g->select_trees_for_CG[i];
			bool inTree = checkbin(tree->bin_vertices, v);

			coverCoeff[i] = inTree ? 1.0 : 0.0;
			zCoeff[i] = inTree ? -tree->weight : 0.0;
		}

		Cover[v].setLinearCoefs(x, coverCoeff);         // one API call
		Z[v].setLinearCoefs(x, zCoeff);             // one API call
	}

	coverCoeff.end();
	zCoeff.end();
}

// remove vars 
void CG::RemoveXVars() {	

	// remove all x variables from the model	
	if (!x.getSize()) {
		printf("No x variables to remove\n");
		return;
	}

	for (int i = 0; i < x.getSize(); ++i) {
		model.remove(x[i]);  // Add the variable to the model to remove it
	}

	model.remove(x);  // Bulk removal	
	x.endElements();     // Bulk memory release
	x.clear();
	x.end(); // End the array, release memory
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
	printf("\n***************\n");
	printf("Node %d: Lvl %d \t LB = %3.2f \t #trees = %d \n", node->number, node->lvl, node->LB, node->n_trees);
	for (int i = 0; i < node->lvl; i++) {
		printf("branch[%d] = (%d,%d) \t rule = %d\n", i, node->branch[i].u, node->branch[i].v, node->branch[i].rule);
	}
	
	// if premature
	if (node->premature) {
		printf("Node is premature\n");
	}
	else {
		printf("Node is not premature\n");
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

	// check if the node LB is not less than the parent lower bound
	if (node->parent != nullptr)
	if (node->LB < node->parent->LB - 0.001) {
		printf("ERROR: Node LB %f < Parent LB %f \n", node->LB, node->parent->LB);
		exit(1);
	}

	// add the node to the sorted linked list if the lower bound is less than the upper bound
	if (node->premature) {
		// add node to the premature list
		premature_nodes.push_back(node); // add the node to the premature nodes vector
	}
	else if (node->LB > instance->UB - 1) {
		// add node to nodes to be deleted
		nodes_to_remove.push_back(node); // add the node to the nodes to be removed vector
		// sort the nodes based on the lower bound in ascending order -- don't touch root
		std::sort(nodes_to_remove.begin(), nodes_to_remove.end(), [](BP_node* a, BP_node* b) {
			return a->LB < b->LB;
			});
	}
	else {
		active_nodes.push_back(node); // add the node to the active nodes vector
		// sort the nodes based on the lower bound in ascending order -- don't touch root
		std::sort(active_nodes.begin(), active_nodes.end(), [](BP_node* a, BP_node* b) {
			return a->LB < b->LB;
			});
	}
}

void CG::PrintNodeList() {
	
	// print active nodes
	printf("Active Nodes: \n");
	for (int i = 0; i < active_nodes.size(); i++) {
		printf("Node %d: %d \t LB = %3.2f \t CLB = %2d \t P: %2d \t L: %2d \t R: %2d \t PM: %1d \n",
			active_nodes[i]->number, active_nodes[i]->lvl, active_nodes[i]->LB,
			active_nodes[i]->CLB,
			active_nodes[i]->parent ? active_nodes[i]->parent->number : -1, // print parent number
			active_nodes[i]->child_L ? active_nodes[i]->child_L->number : -1,
			active_nodes[i]->child_R ? active_nodes[i]->child_R->number : -1,
			active_nodes[i]->premature ? 1 : 0); // print completed status		
	}

	// print premature noeds
	printf("Premature Nodes: \n");
	for (int i = 0; i < premature_nodes.size(); i++) {
		printf("Node %d: %d \t LB = %3.2f \t CLB = %2d \t P: %2d \t L: %2d \t R: %2d \t PM: %1d \t PMB: %1d \t NChC: %1d\n",
			premature_nodes[i]->number, premature_nodes[i]->lvl, premature_nodes[i]->PLB,
			premature_nodes[i]->CLB,
			premature_nodes[i]->parent ? premature_nodes[i]->parent->number : -1, // print parent number
			premature_nodes[i]->child_L ? premature_nodes[i]->child_L->number : -1,
			premature_nodes[i]->child_R ? premature_nodes[i]->child_R->number : -1,
			premature_nodes[i]->premature ? 1 : 0,
			premature_nodes[i]->premature_branched ? 1 : 0, 
			premature_nodes[i]->n_childs_closed); 
	}

	// print nodes to be removed
	printf("Nodes to be removed: \n");
	for (int i = 0; i < nodes_to_remove.size(); i++) {
		printf("Node %d: %d \t LB = %3.2f \t CLB = %2d \t P: %2d \t L: %2d \t R: %2d \t PM: %1d \t C: %1d\n", 
			nodes_to_remove[i]->number, nodes_to_remove[i]->lvl, nodes_to_remove[i]->LB,
			nodes_to_remove[i]->CLB,
			nodes_to_remove[i]->parent ? nodes_to_remove[i]->parent->number : -1, // print parent number
			nodes_to_remove[i]->child_L ? nodes_to_remove[i]->child_L->number : -1,
			nodes_to_remove[i]->child_R ? nodes_to_remove[i]->child_R->number : -1,
			nodes_to_remove[i]->premature ? 1 : 0, 
			nodes_to_remove[i]->completed ? 1 : 0 // print completed status
			);
	}
}