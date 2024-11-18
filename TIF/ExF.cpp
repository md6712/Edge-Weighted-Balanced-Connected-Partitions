#include "ExF.h"
#include "Cplex.h"


// constructor

ExF::ExF(_g* g, bool redirect, bool linear = false) :Cplex(g) {
	printCuts = false;
	printCycles = false;
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

// destructor
ExF::~ExF() {

}


// run
ExF* ExF::Run() {

	Cplex::Run();

	//SaveOpt();
	//instance->PrintOptEdges();

	return this;
}

ExF* ExF::PrintModel() {
	cplex.exportModel("ModelExF.lp");
	return this;
}

void ExF::SaveOpt() {
	this->instance->num_opt_edges = 0;
	for (int i = 0; i < this->instance->num_trees; i++) {
		this->instance->start_edge_tree[i] = this->instance->num_opt_edges;
		for (int e = 0; e < this->instance->num_edges; e++) {			
			if (this->cplex.isExtracted(this->x[e][i])) {
				printf("x[%d][%d]", e, i);
				double value = this->cplex.getValue(this->x[e][i]);
				if (value > 0.5) {
					this->instance->opt_edges[this->instance->num_opt_edges++] = e;
				}
			}
		}
	}
}

// define the variables
void ExF::DefVar() {
	DefVarX();
	DefVarXbar();
	DefVarY();
	DefVarZ();
	DefVarEta();
	DefVarF();
}

// define the variables x
void ExF::DefVarX() {
	x = NumVarMatrix(env, instance->num_edges);
	for (int e = 0; e < instance->num_edges; e++) {
		x[e] = IloNumVarArray(env, instance->num_trees, 0, 1, ILOBOOL);		
		for (int i = 0; i < instance->num_trees; i++) {
			char name[50];
			sprintf(name, "x_%d_%d", e, i);
			x[e][i].setName(name);
		}
	}
}

// define the variables xbar
void ExF::DefVarXbar() {
	xbar = NumVarMatrix(env, instance->num_arcs); // 
	for (int a = 0; a < instance->num_arcs; a++) {
		xbar[a] = IloNumVarArray(env, instance->num_trees, 0, 1, ILOBOOL);
		for (int i = 0; i < instance->num_trees; i++) {
			char name[50];
			sprintf(name, "xbar_%d_%d", a, i);
			xbar[a][i].setName(name);
		}
	}
}

// define the variables y
void ExF::DefVarY() {
	y = NumVarMatrix(env, instance->num_vertices);
	for (int v = 0; v < instance->num_vertices; v++) {		
		y[v] = IloNumVarArray(env, instance->num_trees, 0, 1, ILOBOOL);
		for (int i = 0; i < instance->num_trees; i++) {
			char name[50];
			sprintf(name, "y_%d_%d", v, i);
			y[v][i].setName(name);
		}
	}
}

// define the variables z
void ExF::DefVarZ() {
	z = NumVarMatrix(env, instance->num_vertices);
	for (int v = 0; v < instance->num_vertices; v++) {
		z[v] = IloNumVarArray(env, instance->num_trees, 0, 1, ILOBOOL);
		for (int i = 0; i < instance->num_trees; i++) {
			char name[50];
			sprintf(name, "z_%d_%d", v, i);
			z[v][i].setName(name);
		}
	}
}


// define the variable eta	
void ExF::DefVarEta() {
	eta = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
	eta.setName("eta");
}

// define the variable f
void ExF::DefVarF() {
	f = NumVarCube(env, instance->num_arcs);
	for (int a = 0; a < instance->num_arcs; a++) {
		f[a] = NumVarMatrix(env, instance->num_vertices);
		for (int v = 0; v < instance->num_vertices; v++) {
			f[a][v] = IloNumVarArray(env, instance->num_trees, 0, IloInfinity, ILOFLOAT);
			for (int i = 0; i < instance->num_trees; i++) {				
				sprintf(name, "f_%d_%d_%d", a, v, i);
				f[a][v][i].setName(name);
			}
		}
	}
}

// define the objective function
void ExF::AddObj() {
	// for each tree, eta is the minimum sum of selected egdges	
	for (int i = 0; i < instance->num_trees; i++) {
		IloExpr sum_x(env);		
		for (int e = 0; e < instance->num_edges; e++) {
			sum_x += x[e][i] * instance->edges[e][2];
		}

		// generate name for the constraint		
		sprintf(name, "compute_eta_%d", i);
		model.add(eta >= sum_x).setName(name);
	}

	// the objective is to minimize eta 
	model.add(IloMinimize(env, eta));
}

// define the constraints
void ExF::AddCons() {
	AssignVerticesToTrees();
	DecideRootOfTrees();
	SetYZRelation();
	SetXZRelation();
	SumOfFlowFromSToVEqualsZui();
	FlowConservation();
	FlowConservationUV();
	FlowXBarRelation();
	Flow0vYRelation();
	XbarZRelation();
	XXbarRelation();
}


// assign vertices to trees
void ExF::AssignVerticesToTrees() {
	for (int v = 0; v < instance->num_vertices; v++) {
		IloExpr sum_z(env);
		for (int i = 0; i < instance->num_trees; i++) {
			sum_z += z[v][i];
		}

		// generate name for the constraint
		sprintf(name, "assign_vertices_to_trees_%d", v);
		model.add(sum_z == 1).setName(name);
	}
}

// decide the root of trees
void ExF::DecideRootOfTrees() {
	for (int i = 0; i < instance->num_trees; i++) {
		IloExpr sum_y(env);
		for (int v = 0; v < instance->num_vertices; v++) {
			sum_y += y[v][i];
		}

		// generate name for the constraint
		sprintf(name, "decide_root_of_trees_%d", i);
		model.add(sum_y == 1).setName(name);
	}
}

// set the relation between y and z
void ExF::SetYZRelation() {
	for (int i = 0; i < instance->num_trees; i++) {
		for (int v = 0; v < instance->num_vertices; v++) {
			model.add(y[v][i] <= z[v][i]);
		}
	}
}

// set the relation between x and z
void ExF::SetXZRelation() {
	for (int i = 0; i < instance->num_trees; i++) {
		for (int e = 0; e < instance->num_edges; e++) {
			model.add(x[e][i] <= z[instance->edges[e][0]][i]);
			model.add(x[e][i] <= z[instance->edges[e][1]][i]);
		}
	}
}

// sum of flow from s to v equals zui
void ExF::SumOfFlowFromSToVEqualsZui() {
	for (int i = 0; i < instance->num_trees; i++) {
		for (int v = 0; v < instance->num_vertices; v++) {
			IloExpr sum_f(env);
			for (int a = instance->num_edges * 2; a < instance->num_arcs; a++) {
				sum_f += f[a][v][i];
			}
					
			sprintf(name, "sum_of_flow_from_s_to_v_equals_zui_%d_%d", i, v);
			model.add(sum_f == z[v][i]);			
		}
	}
}

// flow conservation
void ExF::FlowConservation() {
	for (int i = 0; i < instance->num_trees; i++) {
		for (int v = 0; v < instance->num_vertices; v++) {
			IloExpr sum_f(env);
			for (int a = 0; a < instance->num_arcs; a++) {
				if (instance->arcs[a][0] == v) {
					sum_f -= f[a][v][i];
				}
				else if (instance->arcs[a][1] == v) {
					sum_f += f[a][v][i];
				}
			}

			sprintf(name, "flow_conservation_%d_%d_%d", i, v);
			model.add(sum_f == z[v][i]);			
		}
	}
}


// flow conservation for any u,v \in V and u != v
void ExF::FlowConservationUV() {
	for (int i = 0; i < instance->num_trees; i++) {
		for (int u = 0; u < instance->num_vertices; u++) {
			for (int v = 0; v < instance->num_vertices; v++) {
				if (u == v) continue;
				IloExpr sum_f(env);
				for (int a = 0; a < instance->num_arcs; a++) {
					if (instance->arcs[a][0] == v) {
						sum_f -= f[a][u][i];
					}
					else if (instance->arcs[a][1] == v) {
						sum_f += f[a][u][i];
					}
				}

				sprintf(name, "flow_conservationUV_%d_%d_%d", i, v);
				model.add(sum_f == 0).setName(name);
			}
		}
	}
}


// flow xbar relation
void ExF::FlowXBarRelation() {
	// for each arc, for each vertex and for each tree, f is less than xbar
	for (int i = 0; i < instance->num_trees; i++) {
		for (int a = 0; a < instance->num_arcs; a++) {
			for (int v = 0; v < instance->num_vertices; v++) {

				sprintf(name, "flow_xbar_relation_%d_%d_%d", a, v, i);
				model.add(f[a][v][i] <= xbar[a][i]).setName(name);
			}
		}
	}	
}


// flow 0v y relation
void ExF::Flow0vYRelation() {
	// for vertices u and v, and tree i, f_0v is less than y_v
	for (int i = 0; i < instance->num_trees; i++) {
		for (int v = 0; v < instance->num_vertices; v++) {
			for (int a = 0; a < instance->num_arcs; a++) {
				if (instance->arcs[a][0] == 0 && instance->arcs[a][1] == v) {
					sprintf(name, "flow_0v_y_relation_%d_%d_%d", v, i);
					model.add(f[a][v][i] <= y[v][i]).setName(name);
				}
			}		
		}
	}
}


// xbar z relation
void ExF::XbarZRelation() {
	// for each vertex, for each tree, xbar is less than z
	for (int i = 0; i < instance->num_trees; i++) {
		// sum xbar 
		IloExpr sum_xbar(env);
		for (int a = 0; a < instance->num_arcs; a++) {
			sum_xbar += xbar[a][i];
		}

		// sum z	
		IloExpr sum_z(env);
		for (int v = 0; v < instance->num_vertices; v++) {
			sum_z += z[v][i];
		}
		
		sprintf(name, "xbar_z_relation_%d", i);
		model.add(sum_xbar == sum_z - 1).setName(name);
	}
}


// x xbar relation
void ExF::XXbarRelation() {
	// for each edge, for each tree, the sum of xbars for the two arcs is equal to x
	for (int i = 0; i < instance->num_trees; i++) {
		for (int e = 0; e < instance->num_edges; e++) {
			// sum xbar
			IloExpr sum_xbar(env);

			// find the two arcs
			for (int a = 0; a < instance->num_arcs; a++) {
				if (instance->arcs[a][0] == instance->edges[e][0] && instance->arcs[a][1] == instance->edges[e][1]) {
					sum_xbar += xbar[a][i];
				}
				else if (instance->arcs[a][0] == instance->edges[e][1] && instance->arcs[a][1] == instance->edges[e][0]) {
					sum_xbar += xbar[a][i];
				}
			}

			sprintf(name, "xxbar_relation_%d_%d", e, i);
			model.add(sum_xbar == x[e][i]).setName(name);
		}
	}	
}