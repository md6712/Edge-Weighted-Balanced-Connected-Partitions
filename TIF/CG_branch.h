#pragma once

#define BP_MAX_DEPTH 20

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
	int number; 
	int lvl;
	CG_branch branch[BP_MAX_DEPTH];
	int n_trees;
	double PLB;
	double LB;
	int CLB; // confirmed LB
	int branch_pair[2];	
	
	bool premature = false; // if the node is premature, it means that CG is not finished for this node
	bool premature_branched = false; // if the prematurity is resolved, it means that the node is not premature anymore, and it can be removed
	bool completed = false;

	BP_node* parent = nullptr; // parent node
	BP_node* sibling = nullptr; // sibling node
	BP_node* child_L = nullptr; // left child
	BP_node* child_R = nullptr; // right child
	
	int n_childs_closed = 0; // number of closed childs -> if 2, then the node is not premature anymore, and it can be removed.
};