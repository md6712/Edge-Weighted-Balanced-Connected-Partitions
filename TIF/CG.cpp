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
		AddVars();				
	}
}

// destructor
CG::~CG() {
	delete pricing_problem;
}


// run
CG* CG::Run() {

	// create a array of vertices 
	int *vertices = new int[this->instance->num_vertices];

	// create an array of bin vertices
	uint32_t *bin_vertices = new uint32_t[binaryArrlength(this->instance->num_vertices)];

	//// use cplex column genration	
	//cplex.use(callbackCG(env,this));

	// Use an inline lambda function for the callback	
	//cplex.use(new (env) CGColumnGenCallback(this), IloCplex::Callback::Context::Id::ThreadUp);

	// print all trees
	for (int i = 0; i < this->instance->trees.size(); i++) {
		this->instance->trees[i]->PrintTree();
	}

	Timer timer;

	double elapsedtimeMaster = 0;
	double elapsedtimePricing = 0;
	

	while (true)
	{
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
		pricing_problem->Run();
		_tree*  tree = pricing_problem->GetTree(vertices, bin_vertices); // get the tree associated to the optimal solution
		
		printf_s("pricing_problem->opt = %f\n", pricing_problem->opt);

		timer.setEndTime();

		elapsedtimePricing += timer.calcElaspedTime_sec();

		// print the solution
		pricing_problem->PrintSol();
	
		// print the tree
		tree->PrintTree();			

		// add the tree to the instance
		this->instance->trees.push_back(tree);

		// print size of the trees
		//printf("trees.size() = %d\n", this->instance->trees.size());

		if (pricing_problem->opt < 0.001) {
			break;
		}

		// add the variable to the model
		AddVar(this->instance->trees.size() - 1);

	}

	// print elapsed time
	printf("elapsedtimeMaster = %f\n", elapsedtimeMaster);
	printf("elapsedtimePricing = %f\n", elapsedtimePricing);

	delete[] vertices;
	delete[] bin_vertices;

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
			this->instance->trees[i]->PrintVerticesWeight();
		}
	}

	// print the objective value
	printf("obj = %3.2f\n", cplex.getObjValue());

	// print number of trees
	printf("Number of trees: %d\n", this->instance->trees.size());

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
	for (int v = 0; v < this->instance->num_vertices; v++) {
		Cover[v] = IloRange(env, 1, IloInfinity);
		Z[v] = IloRange(env, 0, IloInfinity);
	}
	
	// add to model
	model.add(Knapsack);
	model.add(Cover);
	model.add(Z);
}

// add variables	
void CG::AddVars() {
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

	for (int i = 0; i < g->trees.size(); i++) {
		// name 
		AddVar( i);
	}

	// add variables to the model
	model.add(z);
	model.add(x);	
}

// add variable
void CG::AddVar(int i) {
	// instance g
	_g* g = instance;

	// add variables
	sprintf(name, "x_%d", i);
	IloNumVar Xvar(env, 0, IloInfinity, ILOFLOAT, name);

	obj.setLinearCoef(Xvar, 0);
	Knapsack.setLinearCoef(Xvar, -1);

	for (int v = 0; v < g->num_vertices; v++) {
		if (checkbin(g->trees[i]->bin_vertices, v)) {
			Cover[v].setLinearCoef(Xvar, 1);
		}
		else {
			Cover[v].setLinearCoef(Xvar, 0);
		}

		if (checkbin(g->trees[i]->bin_vertices, v)) {
			Z[v].setLinearCoef(Xvar, -g->trees[i]->weight);
		}
		else {
			Z[v].setLinearCoef(Xvar, 0);
		}
	}

	x.add(Xvar);
}


