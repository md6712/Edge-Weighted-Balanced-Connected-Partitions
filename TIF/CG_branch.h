#pragma once

#define BP_MAX_DEPTH 10

enum CG_branch_rule {
	together, 
	apart	
};

struct CG_branch {
	int u;
	int v; 
	CG_branch_rule rule;
};

struct BP_node {
	int lvl;
	CG_branch branch[BP_MAX_DEPTH];
	int UB;
	double LB;
	int branch_pair[2];
	BP_node* next;
};