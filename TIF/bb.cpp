#include "bb.h"
#include <lemon/kruskal.h>
#include <lemon/list_graph.h>
#include <iostream>
#include <algorithm>
#include "binary.h"


#define COMP_STATE(state1,state2,result)	result = memcmp(state1->bin_vertices,state2->bin_vertices,SIZE_OF_VERTICES_BINARY*sizeof(uint32_t));

bb::bb(_g* instance)
{
	this->instance = instance;
	// initialize the benv
	benv = create_benv(sizeof(_node_t));

	// initialize the benv for trees
	benv_trees = create_benv(sizeof(_small_tree));
	
	// set number of trees
	num_trees_generated = 0;

	// initialize the ssbt
	sbbt = init_multi_hashed_sbbt(SBBT_HASH_SIZE);


	// define graph
	graph = new ListGraph();

	// mst components
	// define nodes 
	nodes = new ListGraph::Node[instance->num_vertices];

	// include the nodes
	for (int i = 0; i < instance->num_vertices; i++) {
		nodes[i] = graph->addNode();
	}

	// define edges
	edges = new ListGraph::Edge[instance->num_edges];

	// define edge_weight
	edge_weights = new ListGraph::EdgeMap<int>(*graph);

	// include the edges
	for (int e = 0; e < instance->num_edges; e++) {
		edges[e] = graph->addEdge(nodes[instance->edges[e][0]], nodes[instance->edges[e][1]]);
		(*edge_weights)[edges[e]] = instance->edges[e][2];
	}

	// define mst 
	mst = new ListGraph::EdgeMap<bool>(*graph, false);

	// define matrix
	matrix = new int* [instance->num_trees - 1];
	for (int i = 0; i < instance->num_trees - 1; i++) {
		matrix[i] = new int[instance->num_vertices - instance->num_trees];
	}

	// define edges to be removed
	memset(edges_to_be_removed, 0, sizeof(uint32_t) * SIZE_OF_EDGES_BINARY);

	// define edges selected
	memset(edges_selected, 0, sizeof(uint32_t) * SIZE_OF_EDGES_BINARY);

	// define tree weights
	tree_weights = new int[instance->num_trees];

	// define num of vertices and edges per tree
	num_vertices = new int[instance->num_trees];
	num_edges = new int[instance->num_trees];

	// define vertices and edges per tree
	vertices = new int*[instance->num_trees];
	edges_ub = new int* [instance->num_trees];
	for (int i = 0; i < instance->num_trees; i++) {
		vertices[i] = new int[instance->num_vertices];
		edges_ub[i] = new int[instance->num_edges];
	}	

	// define bin vertices
	bin_vertices = new uint32_t * [instance->num_trees];
	for (int i = 0; i < instance->num_trees; i++) {
		bin_vertices[i] = new uint32_t[binaryArrlength(instance->num_vertices)];
		memset(bin_vertices[i], 0, sizeof(uint32_t) * binaryArrlength(instance->num_vertices));
	}

	// define decomp
	decomp = new int[instance->num_vertices];

	// init counters
	num_nodes = 0;
	num_kruskal_calls = 0;
	num_ub_calls = 0;

	// init clock ticks
	kruskal_total_time = 0;
	ub_total_time = 0;

}

bb::~bb()
{
	// free the benv
	free_benv(benv);

	// free the benv for trees
	free_benv(benv_trees);

	// free the ssbt
	free_multi_hashed_sbbt(sbbt);

	// delete components mst
	delete[] nodes;
	delete[] edges; 
	delete edge_weights;
	delete mst;
	delete graph;

	// delete components UB
	for (int i = 0; i < instance->num_trees - 1; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
	delete[] tree_weights;
	//delete[] edges_selected;			// this array is fixed size
	//delete[] edges_to_be_removed;		// this array is fixed size
		

	for (int i = 0; i < instance->num_trees; i++) {
		delete[] vertices[i];
		delete[] edges_ub[i];
		delete[] bin_vertices[i];
	}
	delete[] num_vertices;
	delete[] num_edges;
	delete[] vertices;
	delete[] edges_ub;
	delete[] bin_vertices;	

	delete[] decomp;
}

bb* bb::Run() {
	// create temp node
	 _node_t* temp = nullptr;

	// create root node 
	this->root = create_root_node();

	int num_tree_before, num_tree_after; // number of trees before and after branching
	int last_itr_num_tree_increased = 0; // last iteration where the number of trees increased

	int itr = 0;

	// set timer 
	clock_t start = clock();	

	while (root != nullptr) {

		// number of iterations
		itr++;
		
		// check if 1 seconds passed
		clock_t end = clock();
		if ((double)(end - start) / CLOCKS_PER_SEC > 1) {
			break;
		}


		// number of trees before branching
		num_tree_before = num_trees_generated;

		//print_node(root);
		// branch root
		branch(root);

		// number of trees after branching
		num_tree_after = num_trees_generated;

		// check if the number of trees increased
		if (num_tree_after > num_tree_before) {
			last_itr_num_tree_increased = itr;
		}

		// break if num of trees has not increased for 100 iterations
		if (itr - last_itr_num_tree_increased > 100) {
			break;
		}

		if (num_trees_generated > 1000 )
		{
			break;
		}

		if (itr > 50) {
			break;
		}

		if (itr == 1) {
			instance->UB_naive = root->UB;
		}

		// remove root
		_node_t* temp = root->next;
		free_memory(benv, root);
		root = temp;	
		
	}

	// print number of nodes, num krushkal calls, num ub calls
	std::cout << "Number of nodes: " << num_nodes << std::endl;
	std::cout << "Number of Kruskal calls: " << num_kruskal_calls << std::endl;
	std::cout << "Number of UB calls: " << num_ub_calls << std::endl;

	// print time for kruskal in seconds
	std::cout << "Total time for Kruskal: " << (double)kruskal_total_time / CLOCKS_PER_SEC << std::endl;		

	// print time for ub in seconds
	std::cout << "Total time for UB: " << (double)ub_total_time / CLOCKS_PER_SEC << std::endl;

	// print UB naive 
	std::cout << "UB naive: " << instance->UB_naive << std::endl;

	// print UB 
	std::cout << "UB: " << instance->UB << std::endl;

	// add singlton trees
	add_singlton_trees();

	// copy the trees to the vector
	populate_all_trees();

	return this;
}

void bb::branch(_node_t* node) {

	// loop over all edges from the highest to the lowest and see if they are in the mst
	// if they are, then investigate them 
	_node_t* current, *prev;

	for (int e = instance->num_edges -1; e >= 0; e--) {
		if (checkbin(node->forbidden_edges, e)) continue;
		if (checkbin(node->mst_edges, e)) {			
			// create a new node
			_node_t* new_node = (_node_t*)alloc_memory(benv);
			// copy the forbidden edges
			memcpy(new_node->forbidden_edges, node->forbidden_edges, sizeof(uint32_t) * SIZE_OF_EDGES_BINARY);
			// set mst 
			memset(new_node->mst_edges, 0, sizeof(uint32_t) * SIZE_OF_EDGES_BINARY);			
			// add the edge to the forbidden edges
			addbin(new_node->forbidden_edges, e);
			// compute the mst
			compute_mst(new_node);							
			
			// remove the node if LB is greater than UB
			if (new_node->LB > instance->UB) {
				free_memory(benv, new_node);
			}
			else {					
				
				// add the node to linkedlist such that the linked list stay sorted based on LB
				prev = nullptr;
				current = root;
				while (current != nullptr && current->LB <= new_node->LB) {					
					prev = current;
					current = current->next;
				}

				// insert the node
				new_node->next = current;				
				if (prev == nullptr) {
					root = new_node;
				}
				else {
					prev->next = new_node;
				}

				// increase number of nodes
				num_nodes++;
			}
		}
	}	
}

// create root node 
_node_t* bb::create_root_node() {

	// allocate memory
	_node_t* root = (_node_t*)alloc_memory(benv);

	// init the forbidden edges 
	memset(root->forbidden_edges, 0, sizeof(uint32_t) * SIZE_OF_EDGES_BINARY); 
	memset(root->mst_edges, 0, sizeof(uint32_t) * SIZE_OF_EDGES_BINARY);
	root->next = nullptr;

	// compute mst
	compute_mst(root);

	// print node
	//print_node(root);

	// increase number of nodes
	num_nodes++;

	// return root
	return root;
}

void bb::compute_mst(_node_t *node) {

	// update weights for forbidden edges
	for (int e = 0; e < instance->num_edges; e++) {
		if (checkbin(node->forbidden_edges, e)) {
			(*edge_weights)[edges[e]] = 1000000;
		}
	}

	// get clock tick 
	clock_t start = clock();

	// compute kruskal
	kruskal(*graph, *edge_weights, *mst);

	// get clock tick
	clock_t end = clock();

	// update total time
	kruskal_total_time += end - start;

	// add counter for kruskal calls
	num_kruskal_calls++;

	// set weight 0
	int weight = 0;
	int k = 0;

	// loop over edges and see if they are in mst
	for (int e = instance->num_edges - 1; e >= 0 ; e--) {					
		if ((*mst)[edges[e]]) {			
			addbin(node->mst_edges, e);
			if (checkbin(node->forbidden_edges, e)) {
				weight += 1000000;
				break;
			}
			if (++k>instance->num_trees -1)				
				weight += instance->edges[e][2];				
		}
	}

	// check bin array mst and forbidden
	

	if (weight < 1000000) {
		node->LB = ceil((double)weight / instance->num_trees);
		if (node->LB < instance->UB) {

			// get clock tick
			start = clock();

			// compute upper bound
			compute_upper_bound(node);

			// get clock tick
			end = clock();

			// update total time
			ub_total_time += end - start;
		}		
		node->UB = instance->UB;
	}
	else {
		node->LB = 1000000;
		node->UB = 1000000;
	}

	// revert weights
	for (int e = 0; e < instance->num_edges; e++) {
		if (checkbin(node->forbidden_edges, e)) {
			(*edge_weights)[edges[e]] = instance->edges[e][2];
		}
	}
}

void bb::compute_upper_bound(_node_t* node) {
	
	// set edges to be selected to 0
	memset(edges_selected, 0, sizeof(uint32_t) * SIZE_OF_EDGES_BINARY);

	// set edges to be removed to 0
	memset(edges_to_be_removed, 0, sizeof(uint32_t) * SIZE_OF_EDGES_BINARY);

	// select the first k-1 edges to be removed that are in mst
	int k = 0;
	for (int e = instance->num_edges - 1; e >= 0; e--) {
		if (checkbin(node->mst_edges, e)) {
			addbin(edges_to_be_removed, e);	 // to be removed		
			addbin(edges_selected, e); // selected 
			if (++k == instance->num_trees - 1) {
				break;
			}
		}
	}

	// add counter for ub calls	
	num_ub_calls++;

	//// print edges to be removed
	//printf_s("\n Edges to be removed: ");
	//for (int e = 0; e < instance->num_edges; e++) {
	//	if (checkbin(edges_to_be_removed, e)) {
	//		printf_s("(%d,%d)", instance->edges[e][0], instance->edges[e][1]);
	//	}
	//}

	//// print edges to be selected 
	//printf_s("\n Edges selected: ");
	//for (int e = 0; e < instance->num_edges; e++) {
	//	if (checkbin(edges_selected, e)) {
	//		printf_s("(%d,%d)", instance->edges[e][0], instance->edges[e][1]);
	//	}
	//}
	
	int max_forest_weight = 0; 

	// compute spanning k forest
	max_forest_weight = compute_spanning_k_forest(node);
	int min_max_forest_weight = max_forest_weight;

	int a_min, d_min, aa, dd, a, d;

	while (1) {
		a_min = -1;
		d_min = -1;

		aa = 0;
		for (a = 0; a < instance->num_edges; a++) {
			if (!checkbin(edges_to_be_removed, a)) continue;
			dd = 0;
			for (d = 0; d < instance->num_edges; d++) {
				if (checkbin(edges_to_be_removed, d)) continue;
				if (!checkbin(node->mst_edges, d))continue;

				addbin(edges_to_be_removed, d);
				removebin(edges_to_be_removed, a);

				// compute spanning k forest
				int max_forest_weight = compute_spanning_k_forest(node);

				if (max_forest_weight < min_max_forest_weight) {
					min_max_forest_weight = max_forest_weight;
					a_min = a;
					d_min = d;
				}
				else if (max_forest_weight == min_max_forest_weight) {
					if (!checkbin(edges_selected, d)) {
						min_max_forest_weight = max_forest_weight;
						a_min = a;
						d_min = d;
					}
				}

				dd++;

				//reverse 

				addbin(edges_to_be_removed, a);
				removebin(edges_to_be_removed, d);
			}
			aa++; 
		}

		if (a_min == -1 || d_min == -1) {
			break;
		}
		else {
			removebin(edges_to_be_removed, a_min);
			addbin(edges_to_be_removed, d_min);
			addbin(edges_selected, d_min);
		}
	}

	// compute spanning k forest
	node->UB = max_forest_weight = compute_spanning_k_forest(node);

	// update graph upper bound: Min (UB, maxForestWeight)
	bool opt = false; 
	if (node->UB < instance->UB) {
		instance->num_upper_bound_updates++;
		instance->UB = node->UB;
		instance->recomputeLB();
		opt = true;
	}

	// for each tree
	for (int i = 0; i < instance->num_trees; i++) {

		if (tree_weights[i] > instance->UB) {
			continue;
		}

		// create a small tree from benv
		_small_tree* tree = (_small_tree*)alloc_memory(benv_trees);
		memset(tree->bin_vertices, 0, sizeof(uint32_t) * SIZE_OF_VERTICES_BINARY);
		// copy the vertices
		memcpy(tree->bin_vertices, bin_vertices[i], sizeof(uint32_t) * binaryArrlength(instance->num_vertices));		

		// copy weight 
		tree->weight = tree_weights[i];		
		
		uint16_t hash = hash_tree(tree);

		// check if tree exist
		
		sbbt->lvl = 0;

		// check if tree exist
		FIND_STATE_IN_SBBT(sbbt, hash, tree, temp_tree, curr_node, 32, l, _small_tree, COMP_STATE);

		//// print hash value: 
		//printf_s("Hash: %d ", hash);
		//tree->print_vertices(instance);

		if (temp_tree == nullptr){
			// create a new tree
			_small_tree* new_tree = (_small_tree*)alloc_memory(benv_trees);

			// copy the tree
			memcpy(new_tree, tree, sizeof(_small_tree));

			new_tree->part_of_optimal = opt ? instance->num_upper_bound_updates : 0;

			// insert the tree
			ADD_STATE_IN_SBBT(sbbt, hash, new_node, curr_node, new_tree);

			// balance the tree
			BALANCE_SBBT_AFTER_ADD(sbbt, hash, curr_node, mother_node, tmp_node, l);

			// count the number of trees
			num_trees_generated++;
		}
		else {
			// free the tree	
			if (opt)
				temp_tree->part_of_optimal = instance->num_upper_bound_updates;
			free_memory(benv_trees, tree);
		}

		
	}

	//// print edges for ub
	//printf_s("\n Edges for UB: ");
	//for (int i = 0; i < instance->num_trees; i++) {
	//	printf_s("\n i: ", i);
	//	int w = 0;
	//	for (int j = 0; j < num_edges[i]; j++) {
	//		printf_s("(%d,%d)", instance->edges[edges_ub[i][j]][0], instance->edges[edges_ub[i][j]][1]);			
	//		w += instance->edges[edges_ub[i][j]][2];
	//	}
	//	printf_s("\t w: %d ", w);
	//}
}
// spanning k forest
int bb::compute_spanning_k_forest(_node_t* node) {
	
	// traverse the tree to decompose the tree
	traverse_to_decompose(node, 0, -1, 0, 0, edges_to_be_removed);

	// reset num vertices and edges
	memset(num_vertices, 0, sizeof(int) * instance->num_trees);
	memset(num_edges, 0, sizeof(int) * instance->num_trees);
	memset(tree_weights, 0, sizeof(int) * instance->num_trees);

	// reset bin vertices
	for (int i = 0; i < instance->num_trees; i++) {
		memset(bin_vertices[i], 0, sizeof(uint32_t) * binaryArrlength(instance->num_vertices));
	}

	// add vertices to the trees
	for (int v = 0; v < instance->num_vertices; v++) {
		int i = this->decomp[v];		
		vertices[i][num_vertices[i]++] = v;
		addbin(bin_vertices[i], v);
	}



	// compute the edges for each tree
	for (int e = 0; e < instance->num_edges; e++) {
		if (!checkbin(node->mst_edges, e)) continue;

		// get the vertices of the edge
		int u = instance->edges[e][0];
		int v = instance->edges[e][1];
		// get the tree of the vertices
		int treeU = this->decomp[u];
		int treeV = this->decomp[v];
		// if the vertices are in the same tree
		if (treeU == treeV) {
			// add the edge to the tree
			edges_ub[treeU][num_edges[treeU]++] = e;
			// weight of tree
			if (checkbin(node->forbidden_edges, e))
				tree_weights[treeU] += 10000;
			else
				tree_weights[treeU] += instance->edges[e][2];
		}
	}

	// compute max forest weight
	int max_forest_weight = 0;
	for (int i = 0; i < instance->num_trees; i++) {
		if (tree_weights[i] > max_forest_weight) {
			max_forest_weight = tree_weights[i];
		}
	}

	return max_forest_weight;
}

int bb::traverse_to_decompose(_node_t* node, int vertex, int prev_vertex, int i, int largest_i, uint32_t* edge_to_be_removed) {
		
	this->decomp[vertex] = i;

	// print decompsition	
	//std::cout << "Vertex: " << vertex << " Tree: " << i << std::endl;

	int connector; // the vertex to connect to

	int new_tree_vertex[9]; // the vertex to create a new tree
	int new_tree_vertex_count = 0;

	// check all edges
	for (int e = 0; e < instance->num_edges; e++) {
		if (checkbin(node->forbidden_edges, e)) continue;
		if (checkbin(node->mst_edges, e)) {
			// get the vertices of the edge
			int u = instance->edges[e][0];
			int v = instance->edges[e][1];

			connector = -1;

			// if edge is adjacent to the vertex
			if (vertex == u) {
				connector = v;
			}
			else if (vertex == v) {
				connector = u;
			}

			// check if the vertex is connected
			if (connector != -1 && connector != prev_vertex) {

				// if edge is to be excluded
				if (checkbin(edge_to_be_removed,e)) {
					// create a new tree
					new_tree_vertex[new_tree_vertex_count++] = connector;
				}
				else {
					largest_i = traverse_to_decompose(node,connector, vertex, i, largest_i, edge_to_be_removed);

					// print connector vertex i largest i 
					//std::cout << "Connector: " << connector << " Vertex: " << vertex << " i: " << i << " Largest i: " << largest_i << std::endl;

				}

			}
		}

	}

	// if a new tree is to be created
	for (int s = 0; s < new_tree_vertex_count; s++) {
		largest_i = traverse_to_decompose(node, new_tree_vertex[s], vertex, largest_i + 1, largest_i + 1, edge_to_be_removed);
	}

	return largest_i;
}

void bb::print_node(_node_t* node) {
	printf_s("node LB = %5d  UB = %5d  UB* = %5d  \tf ", node->LB, node->UB, instance->UB);

	for (int e = 0; e < instance->num_edges; e++) {
		if (checkbin(node->forbidden_edges, e)) {
			printf("(%d,%d)", instance->edges[e][0], instance->edges[e][1]);
		}
	}
	printf_s("\t m");
	for (int e = 0; e < instance->num_edges; e++) {
		if (checkbin(node->mst_edges, e)) {
			printf("(%d,%d)", instance->edges[e][0], instance->edges[e][1]);
		}
	}
	
	printf_s("\n");
}

__inline uint16_t bb::hash_tree(_small_tree* tree)
{
	uint16_t hash = 0;
	hash = ((uint16_t)tree->bin_vertices[0]+1) % 65521;
	hash = (hash * ((uint32_t)tree->bin_vertices[1]+1)) % 65521;
	hash = (hash * ((uint32_t)tree->bin_vertices[2]+1)) % 65521;
	hash = (hash * ((uint32_t)tree->bin_vertices[3]+1)) % 65521;
	hash = (hash * ((uint32_t)tree->weight+1)) % 65521;
	return hash;
}

// populate all trees
void bb::populate_all_trees() {
	// loop over all possible hash values
	for (int i = 0; i < SBBT_HASH_SIZE; i++) {
		// get the root of the tree
		_sbbt_node* root = sbbt->root[i];
		// if the root is not null
		if (root != nullptr) {			
			
			// create an small tree alias
			_small_tree* t =(_small_tree*) root->item;

			// check if tree size is less than UB 
			if (t->weight <= instance->UB) {				
				// create a small tree not from benv
				_small_tree* tree = new _small_tree();

				// copy the tree
				memcpy(tree, t, sizeof(_small_tree));

				// add the tree to the vector
				instance->select_trees_for_CG.push_back(tree);
			}			

			// investigate the childs
			populate_all_trees_investigate(root);
		}
	}

	// print the trees
	for (int i = 0; i < instance->select_trees_for_CG.size(); i++) {
		instance->select_trees_for_CG[i]->print_vertices(instance);
	}
}

// populate all trees investigate // add the trees in the right and left. 
void bb::populate_all_trees_investigate(_sbbt_node* node) {
	// investigate the left child
	if (node->left != nullptr) {
		
		_small_tree* t = (_small_tree*)((_sbbt_node*)node->left)->item;

		// check if tree size is less than UB
		if (t->weight <= instance->UB) {
			// create a small tree not from benv
			_small_tree* tree = new _small_tree();
			// copy the tree
			memcpy(tree, t, sizeof(_small_tree));
			// add the tree to the vector
			instance->select_trees_for_CG.push_back(tree);
		}
				
		populate_all_trees_investigate((_sbbt_node*)node->left);
	}
	// investigate the right child
	if (node->right != nullptr) {
		_small_tree* t = (_small_tree*)((_sbbt_node*)node->right)->item;

		// check if tree size is less than UB
		if (t->weight < instance->UB) {
			// create a small tree not from benv
			_small_tree* tree = new _small_tree();
			// copy the tree
			memcpy(tree, t, sizeof(_small_tree));
			// add the tree to the vector
			instance->select_trees_for_CG.push_back(tree);
		}

		populate_all_trees_investigate((_sbbt_node*)node->right);
	}
}


// add singleton trees
void bb::add_singlton_trees() {
	// loop over all vertices
	for (int i = 0; i < instance->num_vertices; i++) {
		// create a small tree not from benv
		_small_tree* tree = (_small_tree*)alloc_memory(benv_trees);
		memset(tree->bin_vertices, 0, sizeof(uint32_t) * SIZE_OF_VERTICES_BINARY);
		
		addbin(tree->bin_vertices, i);

		// copy weight 
		tree->weight = 0;

		uint16_t hash = hash_tree(tree);

		// check if tree exist
		sbbt->lvl = 0;

		// check if tree exist
		FIND_STATE_IN_SBBT(sbbt, hash, tree, temp_tree, curr_node, 32, l, _small_tree, COMP_STATE);		

		if (temp_tree == nullptr) {
			// create a new tree
			_small_tree* new_tree = (_small_tree*)alloc_memory(benv_trees);

			// copy the tree
			memcpy(new_tree, tree, sizeof(_small_tree));

			// insert the tree
			ADD_STATE_IN_SBBT(sbbt, hash, new_node, curr_node, new_tree);

			// balance the tree
			BALANCE_SBBT_AFTER_ADD(sbbt, hash, curr_node, mother_node, tmp_node, l);

			// count the number of trees
			num_trees_generated++;
		}
		else {
			// free the tree
			free_memory(benv_trees, tree);
		}
	}
}