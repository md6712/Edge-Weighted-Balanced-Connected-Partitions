/*
* Copyright 2013-2025 Morteza Davari.  All rights reserved.
*
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions and the following
 *      disclaimer in the documentation and/or other materials
 *      provided with the distribution.

 *  $Id: sbbt.cpp $
 *  $Date: 2018/01/28 $
 *  $Author: morteza $
*/

// self balanced binary tree
#include <cstdlib> 
#include "sbbt.h"


_sbbt* init_sbbt(){
	_sbbt *sbbt = new _sbbt;
	sbbt->env_leaves = create_benv(sizeof(_sbbt_leaf));
	sbbt->env_nodes = create_benv(sizeof(_sbbt_node));
	sbbt->root = NULL;
	return sbbt;
}
void free_sbbt(_sbbt *sbbt){
	free_benv(sbbt->env_leaves);
	free_benv(sbbt->env_nodes);
	free(sbbt->root);
	free(sbbt);
}


_multi_hashed_sbbt* init_multi_hashed_sbbt(int size){
	_multi_hashed_sbbt *sbbt = new _multi_hashed_sbbt;	
	sbbt->sizeofHASH = size;
	sbbt->env_nodes = create_benv(sizeof(_sbbt_node));
	sbbt->root =(_sbbt_node **) malloc(sizeof(_sbbt_node*)*size);
	sbbt->branch_nodes = (_sbbt_node **) malloc(sizeof(_sbbt_node*)*40); // max lvl 40
	sbbt->branch_lr = (int8_t *) malloc(sizeof(int8_t)*40); // max lvl 40
	for (int i = 0; i<size; i++){sbbt->root[i] = NULL;}
	return sbbt;
}
void free_multi_hashed_sbbt(_multi_hashed_sbbt *sbbt){
	free_benv(sbbt->env_nodes);
	free(sbbt->root);
	free(sbbt->branch_nodes);
	free(sbbt->branch_lr);
	free(sbbt);
}
void add_state_in_sbbt_func(_multi_hashed_sbbt *sbbt, uint16_t hash_V_16,_sbbt_node* new_node,_sbbt_node* curr_node,void * new_state ) {
	new_node = (_sbbt_node*) alloc_memory(sbbt->env_nodes); 
	new_node->status = 0; new_node->item = new_state; 
	new_node->left = new_node->right = NULL; 
	if (sbbt->lvl > 0){ 
		curr_node = sbbt->branch_nodes[sbbt->lvl-1]; 
		if (sbbt->branch_lr[sbbt->lvl-1] < 0){curr_node->status ++;curr_node->left = new_node;}else{curr_node->status --;curr_node->right = new_node;}
	}else{
		sbbt->root[hash_V_16] = sbbt->curr_root = new_node;
	}
}
void balance_sbbt_after_add_func(_multi_hashed_sbbt *sbbt, uint16_t hash_V_16,_sbbt_node* curr_node,_sbbt_node* mother_node,_sbbt_node *tmp_node,uint8_t l) {
	while (sbbt->lvl > 1 && curr_node->status != 0){
		mother_node = sbbt->branch_nodes[sbbt->lvl-2];
		if (sbbt->branch_lr[sbbt->lvl-2] == -1){
			if (mother_node->status == 1){
				if(curr_node->status == -1){
					l = sbbt->branch_lr[sbbt->lvl-2];ROTATE_LEFT(curr_node,mother_node,tmp_node,l,sbbt->root[hash_V_16]);
					if (curr_node->status == -1){ ((_sbbt_node *)curr_node->left)->status = 1;} else{((_sbbt_node *)curr_node->left)->status = 0;}	curr_node->status = 1;
				}
				if (sbbt->lvl > 2) {l = sbbt->branch_lr[sbbt->lvl-3];ROTATE_RIGHT(mother_node,sbbt->branch_nodes[sbbt->lvl-3],tmp_node,l,sbbt->root[hash_V_16]);}
				else {ROTATE_RIGHT(mother_node,NULL,tmp_node,0,sbbt->root[hash_V_16]);}
				mother_node->status = 0;((_sbbt_node *)mother_node->right)->status = 0;	break;
			}else if (mother_node->status == -1){
				mother_node->status = 0;break;
			}else{
				mother_node->status = 1;
			}
		}else{
			if (mother_node->status == -1){
				if(curr_node->status == 1){
					l = sbbt->branch_lr[sbbt->lvl-2];ROTATE_RIGHT(curr_node,mother_node,tmp_node,l,sbbt->root[hash_V_16]);
					if (curr_node->status == 1){ ((_sbbt_node *)curr_node->right)->status = -1;} else{((_sbbt_node *)curr_node->right)->status = 0;}curr_node->status = -1;
				}
				if (sbbt->lvl > 2) {l = sbbt->branch_lr[sbbt->lvl-3];ROTATE_LEFT(mother_node,sbbt->branch_nodes[sbbt->lvl-3],tmp_node,l,sbbt->root[hash_V_16]);}
				else {ROTATE_LEFT(mother_node,NULL,tmp_node,0,sbbt->root[hash_V_16]);}
				mother_node->status = 0;((_sbbt_node *)mother_node->left)->status = 0;break;
			}else if (mother_node->status == 1){
				mother_node->status = 0;break;
			}else{
				mother_node->status = -1;
			}
		}
		sbbt->lvl--;curr_node = sbbt->branch_nodes[sbbt->lvl-1];
	}
}


void** return_all_items_in_sbbt(_multi_hashed_sbbt* sbbt) {
	void** items = (void**)malloc(sizeof(void*) * 1000000);
	int k = 0;
	for (int i = 0; i < sbbt->sizeofHASH; i++) {
		_sbbt_node* node = sbbt->root[i];
		if (node == NULL) continue;
		_sbbt_node* stack[1000];
		int top = 0;
		while (top > 0 || node != NULL) {
			if (node != NULL) {
				stack[top++] = node;
				node = (_sbbt_node*)node->left;
			}
			else {
				node = stack[--top];
				items[k++] = node->item;
				node = (_sbbt_node*)node->right;
			}
		}
	}
	items[k] = NULL;
	return items;
}