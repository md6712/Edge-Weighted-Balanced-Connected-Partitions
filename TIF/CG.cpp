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
	BP_node* root = new BP_node();	
	root->next = nullptr;
	root->lvl = 0;
	root->LB = 0;
	root->UB = instance->UB;

	Run_CG(root);

	// create a temp node --- only a pointer 
	BP_node* temp; 

	while (root != nullptr && root->LB + 0.001 < instance->UB - 1) {
		// copy root into two child nodes
		BP_node* node_L = new BP_node();
		BP_node* node_R = new BP_node();

		memcpy(node_L, root, sizeof(BP_node));
		memcpy(node_R, root, sizeof(BP_node));

		int u = best_pair[0];
		int v = best_pair[1];

		// set the branch rule u and v apart 
		node_L->branch[node_L->lvl].u = u;
		node_L->branch[node_L->lvl].v = v;
		node_L->branch[node_L->lvl].rule = CG_branch_rule::apart;
		node_L->lvl += 1;

		// run column generation
		Run_CG(node_L);

		// add the node to the sorted linked list if the lower bound is less than the upper bound
		if (node_L->LB+0.001 > instance->UB) {
			delete node_L;					
		}
		else {
			AddNodeSorted(node_L);
		}		

		// set the branch rule u and v together
		node_R->branch[node_R->lvl].u = u;
		node_R->branch[node_R->lvl].v = v;
		node_R->branch[node_R->lvl].rule = CG_branch_rule::together;
		node_R->lvl += 1;

		// run column generation
		Run_CG(node_R);			

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
		delete temp;
	} 
	

	return this;
}


// run
CG* CG::Run_CG(BP_node* node) {

	Impose_Braching(node);	

	// print model
	PrintModel();

	Timer timer;

	double elapsedtimeMaster = 0;
	double elapsedtimePricing = 0;
	double elapsedtimePricingPCST = 0;

	int itr = 0;


	while (true)
	{
		itr++;
		timer.setStartTime();
		Cplex::Run();
		timer.setEndTime();
		elapsedtimeMaster += timer.calcElaspedTime_sec();

		// print the solution
		PrintSol();

		dualCover = IloNumArray(env, this->instance->num_vertices);
		dualZ = IloNumArray(env, this->instance->num_vertices);
		dualKnapsack = IloNum(0);

		// read duals values
		dualKnapsack = cplex.getDual(Knapsack);
		cplex.getDuals(dualCover, Cover);
		cplex.getDuals(dualZ, Z);

		// print duals
		printf_s("dualCover\t = ");
		for (int v = 0; v < this->instance->num_vertices; v++) {
			printf_s("%3.2f  ", dualCover[v]);
		}
		printf_s("\n");
		printf_s("dualZ \t\t = ");
		for (int v = 0; v < this->instance->num_vertices; v++) {
			printf_s("%3.2f  ", dualZ[v]);
		}
		printf_s("\n");
		printf_s("dualKnapsack\t = %3.2f\n", dualKnapsack);

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

		// compute time
		//pricing_problem->Run();

		//_small_tree* tree = pricing_problem->GetTree(); // get the tree associated to the optimal solution		

		// print the solution
		//pricing_problem->PrintSol();

		// print the tree
		//tree->print_vertices(this->instance);			

		timer.setEndTime();
		elapsedtimePricing += timer.calcElaspedTime_sec();


		timer.setStartTime();
		pricing_problem->solve_pcst(node, false);
		timer.setEndTime();
		elapsedtimePricingPCST += timer.calcElaspedTime_sec();


		// if there is any tree to be added, add them 
		if (pricing_problem->num_trees_to_add > 0 && pricing_problem->best_pcst > 0.001) {
			for (int i = 0; i < pricing_problem->num_trees_to_add; i++) {
				// get the tree
				_small_tree* tree = pricing_problem->trees_to_add[i];
				// add the tree to the instance
				this->instance->select_trees_for_CG.push_back(tree);

				// print the tree
				tree->print_vertices(this->instance);

				// add the variable to the model
				AddVar(this->instance->select_trees_for_CG.size() - 1, false);
			}
		}
		else {
			// if there is no tree to be added, break the loop
			break;
		}		

		// print all elapsed time
		printf("elapsedtimeMaster = %f\n", elapsedtimeMaster);
		printf("elapsedtimePricing = %f\n", elapsedtimePricing);
		printf("elapsedtimePricingPCST = %f\n", elapsedtimePricingPCST);

	}

	// print elapsed time
	printf("elapsedtimeMaster = %f\n", elapsedtimeMaster);
	printf("elapsedtimePricing = %f\n", elapsedtimePricing);
	printf("elapsedtimePricingPCST = %f\n", elapsedtimePricingPCST);

	// compute the total assignment
	ComputeTotalAssignment();

	// print the total assignment
	PrintTotalAssignment();

	// lower bound of the node is equal the optimal value of the master
	node->LB = cplex.getObjValue();
	
	//if (node->LB < instance->UB) {
	//	
	//	// remove all variables
	//	RemoveXVars();

	//	// each tree add x as an integer variable
	//	for (int i = 0; i < this->instance->select_trees_for_CG.size(); i++) {
	//		AddVar(i, true);
	//	}
	//	model.add(x);

	//	// solve the master with integer variables
	//	timer.setStartTime();
	//	Cplex::Run();
	//	timer.setEndTime();

	//	int UB = (int) (cplex.getObjValue() + 0.001);

	//	if (UB < instance->UB) {
	//		instance->UB = UB;
	//		node->UB = UB;
	//	}
	//	else {
	//		node->UB = instance->UB;
	//	}

	//	// reverse the variables
	//	RemoveXVars();

	//	for (int i = 0; i < this->instance->select_trees_for_CG.size(); i++) {
	//		AddVar(i, false);
	//	}
	//	model.add(x);
	//}

	// if lower bound is less than upper bound - 1 branch

	if (node->LB + 0.001 < instance->UB - 1) {		
		// copy the best pair of vertices for branching; this is used for Ryan Foster branching
		node->branch_pair[0] = best_pair[0];
		node->branch_pair[1] = best_pair[1];		
	}

	return this;
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
	sprintf(name, "x_%d", i);
	IloNumVar Xvar(env, 0, IloInfinity, integer?ILOINT:ILOFLOAT, name);

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
				else {
					x[i].setLB(0);
					x[i].setUB(1);
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
				else {
					// check if the tree contains both vertices, then set the bounds of the variable to 1
						x[i].setLB(0);
						x[i].setUB(1);
				}
			}
		}
	}
	
	return this;
}


// compute total assignment
void CG::ComputeTotalAssignment() {
	// instance g
	_g* g = instance;
	// compute the total assignment

	// set best pair
	best_pair[0] = 0;
	best_pair[1] = 1;

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
				best_pair[0] = i;
				best_pair[1] = j;
			}
		}
	}	
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
		while (current != nullptr && current->LB < node->LB) {
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