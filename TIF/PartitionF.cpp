#include "PartitionF.h"
#include "binary.h"

ILOLAZYCONSTRAINTCALLBACK1(callbackPartition, PartitionF*, partition) {
	// print values of the variables y 
	memset(partition->instance->num_opt_vertices, 0, sizeof(int) * partition->instance->num_trees);

	for (int i = 0; i < partition->instance->num_trees; i++) {
		for (int v = 0; v < partition->instance->num_vertices; v++) {
			if (partition->cplex.isExtracted(partition->y[v][i])) {
				double value = getValue(partition->y[v][i]);
				if (value > 0.5) {
					partition->instance->opt_vertices[i][partition->instance->num_opt_vertices[i]++] = v;
				}
			}
		}
	}


	if (partition->printCuts) {
		// get the objective value
		int opt_value = (int)getObjValue();
		// print objective value
		printf_s("\nObjective value: %d", opt_value);
	}
	
	
	// set the number of edges to zero
	partition->instance->num_opt_edges = 0;

	// apply MST 
	for (int i = 0; i < partition->instance->num_trees; i++) {
		int mst_w = partition->instance->MST(i,partition->instance->opt_vertices[i], partition->instance->num_opt_vertices[i]);
		if (mst_w > 0) {
			for (int i2 = 0; i2 < partition->instance->num_trees; i2++) {
				IloExpr cons(partition->env);
				if (partition->printCuts) printf_s("\nz>= %d", mst_w);
				for (int v = 0; v < partition->instance->num_opt_vertices[i]; v++) {
					cons -= (1 - partition->y[partition->instance->opt_vertices[i][v]][i2]) * partition->instance->sum_adjacent_weight[partition->instance->opt_vertices[i][v]];
					if (partition->printCuts) printf_s(" - %d*y[%d][%d] ", partition->instance->sum_adjacent_weight[partition->instance->opt_vertices[i][v]], partition->instance->opt_vertices[i][v], i2);
				}

				add(partition->z >= mst_w + cons);
				partition->cplex.addUserCut(partition->z >= mst_w - cons);
			}
		}
		else if (mst_w == -1) {
			for (int i2 = 0; i2 < partition->instance->num_trees; i2++) {
				
				// initialize bin_vertices
				memset(partition->instance->bin_vertices, 0, sizeof(uint32_t) * binaryArrlength(partition->instance->num_vertices));
			
				// add vertices to bin_vertices
				for (int v = 0; v < partition->instance->num_opt_vertices[i]; v++) {					
					addbin(partition->instance->bin_vertices, partition->instance->opt_vertices[i][v]);
				}

				// generate cut 

				// introduce cons
				IloExpr cons(partition->env);
				if (partition->printCuts) printf_s("\n");

				// add y[v][i2] to cons if v is in the tree 
				for (int v = 0; v < partition->instance->num_vertices; v++) {
					if (checkbin(partition->instance->bin_vertices, v)) {
						cons += partition->y[v][i2];
						if (partition->printCuts) printf_s(" + y[%d][%d] ", v, i2);
					}				
				}

				// add y[v][i3] to cons if v is not in the tree 
				for (int v = 0; v < partition->instance->num_vertices; v++) {
					if (!checkbin(partition->instance->bin_vertices, v)) {											
						cons += 1 - partition->y[v][i2];
						if (partition->printCuts) printf_s(" + (1 - y[%d][%d]) ", v, i2);						
					}
				}

				if (partition->printCuts) printf_s("<= %d", partition->instance->num_vertices -1);
				add(cons <= partition->instance->num_vertices - 1);
				partition->cplex.addUserCut(cons <= partition->instance->num_vertices - 1);
			}
		}
	}

};

void PartitionF::SaveOpt() {
	memset(instance->num_opt_vertices, 0, sizeof(int) * instance->num_trees);

	for (int i = 0; i < instance->num_trees; i++) {
		for (int v = 0; v < instance->num_vertices; v++) {
			if (cplex.isExtracted(y[v][i])) {
				double value = cplex.getValue(y[v][i]);
				if (value > 0.5) {
					instance->opt_vertices[i][instance->num_opt_vertices[i]++] = v;
				}
			}
		}
	}

	if (printCuts)
	// print vertices in the optimal solution
	for (int i = 0; i < instance->num_trees; i++) {
		printf_s("\nTree %d: ", i);
		for (int v = 0; v < instance->num_opt_vertices[i]; v++) {
			printf_s("%d ", instance->opt_vertices[i][v]);
		}
	}
	


	// set the number of edges to zero
	instance->num_opt_edges = 0;

	// set instance opt values to zero
	instance->_opt = 0;

	// apply MST 
	for (int i = 0; i < instance->num_trees; i++) {
		int mst_w = instance->MST(i, instance->opt_vertices[i], instance->num_opt_vertices[i]);
		if (printCuts) printf_s("\nmst %d: %d", i, mst_w);
		
		if (mst_w > instance->_opt) {
			instance->_opt = mst_w;
		}
	}
		
	if (printCuts) printf_s("\nObjective value: %d",(int) instance->_opt);
}


PartitionF::PartitionF(_g* instance, bool redirect = false) :Cplex(instance) {
	printCuts = false;
	if (!redirect) {
		// Define Variables 
		DefVar();

		// Add Objective function 
		AddObj();

		// Add Constraints
		AddCons();
	}
}

PartitionF* PartitionF::Run() {
	cplex.use(callbackPartition(env, this));
	Cplex::Run();

	SaveOpt();
	//instance->PrintOptEdges();

	return this;
}

PartitionF* PartitionF::PrintModel() {
	cplex.exportModel("submodel_Partition.lp");
	return this;
}

void PartitionF::DefVar() {
//	DefVarX();
	DefVarY();
	DefVarZ();
}

void PartitionF::DefVarX() {
	x = NumVarMatrix(env, instance->num_edges);
	for (int i = 0; i < instance->num_edges; i++) {
		x[i] = IloNumVarArray(env, instance->num_trees, 0, 1, ILOINT);
		for (int j = 0; j < instance->num_trees; j++) {
			sprintf_s(name, "x(%d,%d)", i, j);
			x[i][j].setName(name);
		}
	}
}

void PartitionF::DefVarY() {
	y = NumVarMatrix(env, instance->num_vertices);
	for (int v = 0; v < instance->num_vertices; v++) {
		y[v] = IloNumVarArray(env, instance->num_trees, 0, 1, ILOINT);
		for (int i = 0; i < instance->num_trees; i++) {
			sprintf_s(name, "y(%d,%d)", v, i);
			y[v][i].setName(name);
		}
	}
}

void PartitionF::DefVarZ() {
	z = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
}

void PartitionF::AddObj() {
	IloExpr exprObj(env);	
	model.add(IloMinimize(env, z));
}

void PartitionF::AddCons() {
	AddConsEachVertexInATree();
	AddConsMinSeperators();
	//AddConsOrderTrees();
	//AddConsXYRelation();
}

// add constraints to ensure each vertex is in exactly one tree
void PartitionF::AddConsEachVertexInATree()
{
	for (int v = 0; v < instance->num_vertices; v++) {
		IloExpr cons(env);
		for (int i = 0; i < instance->num_trees; i++) {
			cons += y[v][i];
		}
		sprintf_s(name, "ConsEachVertexInATree(%d)", v);
		model.add(cons == 1).setName(name);
		cons.end();
	}
}

// add constraints to ensure that the trees are ordered
void PartitionF::AddConsOrderTrees()
{
	for (int i = 0; i < instance->num_trees - 1; i++) {
		IloExpr conse(env);		
		for (int e = 0; e < instance->num_edges; e++) {
			conse += (this->x[e][i] - this->x[e][i + 1] )* instance->edges[e][2];			
		}
		sprintf_s(name, "ConsOrderTrees(%d)", i);
		model.add(conse >= 0).setName(name);
		conse.end();
	}
}

// add constraints to make sure that the edges and vertices are in the same tree
void PartitionF::AddConsXYRelation() {	
	for (int i = 0; i < instance->num_trees; i++) {
		IloExpr consv(env); 
		IloExpr conse(env);
		for (int v = 0; v < instance->num_vertices; v++) {
			consv += y[v][i];
		}
		for (int e = 0; e < instance->num_edges; e++) {
			conse += x[e][i];
		}
		sprintf_s(name, "ConsXYRelation(%d)", i);
		model.add(conse <= consv - 1).setName(name); // x[e][i] <= y[v][i] - 1   force number of edges in a tree to be one less than the number of vertices in that tree
		consv.end();
		conse.end();
	}

	for (int i = 0; i < instance->num_trees; i++) {
		for (int v = 0; v < instance->num_vertices; v++) {
			IloExpr conse(env);
			for (int e = 0; e < instance->num_edges; e++) {
				if (instance->edges[e][0] == v || instance->edges[e][1] == v) {
					conse += x[e][i];
					sprintf_s(name, "conseEV(%d,%d,%d)", v, e, i);
					model.add(y[v][i] >= x[e][i]).setName(name); // y[v][i] >= x[e][i] // prevent edges to be in a tree if the vertex is not in that tree
				}
			}
			sprintf_s(name, "ConsVE(%d,%d)", v, i);
			model.add(y[v][i] <= conse).setName(name); // y[v][i] <= conse // prevent vertices to be in a tree if the edges are not in that tree
			conse.end();
		}
	}	
}

// add constraints to ensure that the minimum seperators are satisfied
void PartitionF::AddConsMinSeperators() {
	for (int i = 0; i < instance->num_trees; i++) {
		for (int u = 0; u < instance->num_vertices; u++) {
			for (int v = u + 1; v < instance->num_vertices; v++) {				
				if (instance->UVSeperator(u, v)) {
					IloExpr cons(env);
					cons += y[u][i] + y[v][i];
					//printf("\ny[%d][%d] + y[%d][%d]", u, i, v, i);
					for (int w = 0; w < instance->num_vertices; w++) {
						if (w != u && w != v) {
							if (instance->v_colors[w] == yellow) {
								cons -= y[w][i];
						//		printf("- y[%d][%d]", w, i);
							}
						}
					}
					//printf("<= 1");
					sprintf_s(name, "ConsMinSeperators(%d,%d,%d)", i, u, v);
					model.add(cons <= 1).setName(name);
					cons.end();
				}				
			}
		}
	}
}