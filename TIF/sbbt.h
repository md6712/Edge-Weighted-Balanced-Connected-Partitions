
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

 *  $Id: sbbt.h $
 *  $Date: 2018/01/28 $
 *  $Author: morteza $
*/

// self balanced binary tree

#pragma once

#include "binary.h"
#include "bmemory.h"

// binary values for the balance factor
#define BF_ZERO 0
#define BF_POSONE 1 
#define BF_NEGONE -1
#define BF_POSTWO 2 
#define BF_NEGTWO -2

#define SBBT_HASH_SIZE 65521

struct _sbbt_leaf{
	int8_t status;
	void *item;
};
struct _sbbt_node{
	int8_t status;
	void *item;
	void *left;
	void *right;
};
struct _sbbt{
	_benv_t *env_nodes; // memory environment for nodes
	_benv_t *env_leaves; // memory environment for leaves
	void* root;  // it is either a node or a leaf
};
struct _multi_hashed_sbbt{
	uint32_t sizeofHASH;
	_benv_t *env_nodes; // memory environment for nodes
	_sbbt_node **root;

	// tracking
	_sbbt_node *curr_root;
	uint8_t lvl;
	_sbbt_node **branch_nodes;
	int8_t *branch_lr;
};


_multi_hashed_sbbt* init_multi_hashed_sbbt(int size);
void free_multi_hashed_sbbt(_multi_hashed_sbbt *sbbt);

_sbbt* init_sbbt();
void free_sbbt(_sbbt *sbbt);

#define ROTATE_RIGHT(node, mother,tmp, LRofmother,root) tmp =(_sbbt_node *)((node)->left);(node)->left = (tmp)->right;(tmp)->right = (node);\
	if((mother)!= NULL){if(LRofmother<0){((_sbbt_node *)mother)->left = tmp;}else{((_sbbt_node *)mother)->right = tmp;}}else{root = tmp;} node = tmp;
#define ROTATE_LEFT(node, mother,tmp, LRofmother,root) tmp =(_sbbt_node *) ((node)->right);(node)->right = (tmp)->left;(tmp)->left = (node);\
	if(mother!= NULL){ if (LRofmother<0){((_sbbt_node *)mother)->left = tmp;}else{((_sbbt_node *)mother)->right = tmp;}}else{root = tmp;} node = tmp;

#define FIND_STATE_IN_SBBT(sbbt,hash_V_16,state,tmp_state,curr_node,nactbin,l,type_item,cmp_func) curr_node = sbbt->curr_root = sbbt->root[hash_V_16];	tmp_state = NULL; sbbt->lvl = 0;\
							if (curr_node != NULL){\
								while(1){ tmp_state = (type_item*)curr_node->item;	cmp_func(state,tmp_state,l);	\
									if (l == 0){break;}	else{sbbt->branch_lr[sbbt->lvl] = l; sbbt->branch_nodes[sbbt->lvl++] = curr_node;if (l > 0){curr_node =(_sbbt_node *) curr_node->right;}else {curr_node = (_sbbt_node *) curr_node->left;}	if (curr_node == NULL){tmp_state = NULL;break;}}\
								}\
							}

#define ADD_STATE_IN_SBBT(sbbt,hash_V_16, new_node, curr_node, new_state ) new_node = (_sbbt_node*) alloc_memory(sbbt->env_nodes); new_node->status = 0; new_node->item = new_state; new_node->left = new_node->right = NULL; \
								if (sbbt->lvl > 0){ curr_node = sbbt->branch_nodes[sbbt->lvl-1]; if (sbbt->branch_lr[sbbt->lvl-1] < 0){curr_node->status ++;curr_node->left = new_node;}else{curr_node->status --;curr_node->right = new_node;}}else{sbbt->root[hash_V_16] = sbbt->curr_root = new_node;}


#define BALANCE_SBBT_AFTER_ADD(sbbt,hash_V_16,curr_node,mother_node,tmp_node,l) \
								while (sbbt->lvl > 1 && curr_node->status != 0){\
									mother_node = sbbt->branch_nodes[sbbt->lvl-2];\
									if (sbbt->branch_lr[sbbt->lvl-2] == -1){\
										if (mother_node->status == 1){\
											if(curr_node->status == -1){\
												l = sbbt->branch_lr[sbbt->lvl-2];ROTATE_LEFT(curr_node,mother_node,tmp_node,l,sbbt->root[hash_V_16]);\
												if (curr_node->status == -1){ ((_sbbt_node *)curr_node->left)->status = 1;} else{((_sbbt_node *)curr_node->left)->status = 0;}	curr_node->status = 1;\
											}\
											if (sbbt->lvl > 2) {l = sbbt->branch_lr[sbbt->lvl-3];ROTATE_RIGHT(mother_node,sbbt->branch_nodes[sbbt->lvl-3],tmp_node,l,sbbt->root[hash_V_16]);}\
											else {ROTATE_RIGHT(mother_node,NULL,tmp_node,0,sbbt->root[hash_V_16]);}\
											mother_node->status = 0;((_sbbt_node *)mother_node->right)->status = 0;	break;\
										}else if (mother_node->status == -1){\
											mother_node->status = 0;break;\
										}else{\
											mother_node->status = 1;\
										}\
									}else{\
										if (mother_node->status == -1){\
											if(curr_node->status == 1){\
												l = sbbt->branch_lr[sbbt->lvl-2];ROTATE_RIGHT(curr_node,mother_node,tmp_node,l,sbbt->root[hash_V_16]);\
												if (curr_node->status == 1){ ((_sbbt_node *)curr_node->right)->status = -1;} else{((_sbbt_node *)curr_node->right)->status = 0;}curr_node->status = -1;\
											}\
											if (sbbt->lvl > 2) {l = sbbt->branch_lr[sbbt->lvl-3];ROTATE_LEFT(mother_node,sbbt->branch_nodes[sbbt->lvl-3],tmp_node,l,sbbt->root[hash_V_16]);}\
											else {ROTATE_LEFT(mother_node,NULL,tmp_node,0,sbbt->root[hash_V_16]);}\
											mother_node->status = 0;((_sbbt_node *)mother_node->left)->status = 0;break;\
										}else if (mother_node->status == 1){\
											mother_node->status = 0;break;\
										}else{\
											mother_node->status = -1;\
										}\
									}\
									sbbt->lvl--;curr_node = sbbt->branch_nodes[sbbt->lvl-1];\
								}\

void add_state_in_sbbt_func(_multi_hashed_sbbt *sbbt, uint16_t hash_V_16,_sbbt_node* new_node,_sbbt_node* curr_node,void * new_state );
void balance_sbbt_after_add_func(_multi_hashed_sbbt *sbbt, uint16_t hash_V_16,_sbbt_node* _curr_node,_sbbt_node* mother_node,_sbbt_node *tmp_node, uint8_t l);

//__forceinline  _sbbt_node * rotate_left(_sbbt *sbbt, _sbbt_node *node, _sbbt_node *mother, int8_t LRofmother,_sbbt_node *tmp, uint16_t &hash_V_16){
//	tmp =(_sbbt_node *) node->right;if (tmp) node->right = tmp->left;tmp->left = node;
//	if(mother!= NULL){ if (LRofmother<0){mother->left = tmp;}else{mother->right = tmp;}}
//	else{*((void**)sbbt->root+hash_V_16) = tmp;}
//	return tmp;
//}
//
//__forceinline  _sbbt_node * rotate_right(_sbbt *sbbt, _sbbt_node *node, _sbbt_node *mother, int8_t LRofmother,_sbbt_node *tmp, uint16_t &hash_V_16){
//	tmp =(_sbbt_node *) node->left;node->left = tmp->right;tmp->right = node;
//	if(mother!= NULL){ if (LRofmother<0){mother->left = tmp;}else{mother->right = tmp;}}
//	else{*((void**)sbbt->root+hash_V_16) = tmp;}
//	return tmp;	
//}