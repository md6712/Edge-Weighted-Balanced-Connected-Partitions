/*
* This an adaptation of the code of Shunji Tanaka -- SEE BELOW.
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
 
 *  $Id: bmemory.h $
 *  $Date: 2018/01/28 $
 *  $Author: morteza $
*/


/*
 * Copyright 2006-2012 Shunji Tanaka.  All rights reserved.
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
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *  $Id: bmemory.h,v 1.7 2012/01/28 03:13:03 tanaka Rel $
 *  $Revision: 1.7 $
 *  $Date: 2012/01/28 03:13:03 $
 *  $Author: tanaka $
 *
 */
#ifndef BMEMORY_H
#define BMEMORY_H

typedef struct {
  int n;
  int nb;
  int cb;
  size_t size;
  int n_per_block;
  void *ptr;
  void *bptr;
  void *cbptr;
  void *mem_env;
} _benv_t;

typedef struct {
	int min_size;
	int max_size;
	int step;
	long long allocated;
	long long used;
	_benv_t **env_s;
	uint32_t* activated;
} _mem_env;


_benv_t *create_benv(size_t);
int _n_per_block(size_t);
void free_benv(_benv_t *);
void *alloc_bmemory(_benv_t *);
void free_unused_bmemory(_benv_t *);
void free_bmemory(_benv_t *);
void *alloc_memory(_benv_t *);
void free_memory(_benv_t *, void *);
long long calc_allocated_mem(_benv_t *benv);
char* calc_and_return_allocated_mem(_benv_t* benv);

_mem_env *create_mem_env(int min_size, int max_size, int step);
void init_mem_env(_mem_env *, int min_size, int max_size, int step); // min size cannot be less than 8
void free_mem_env(_mem_env *);
void *alloc_memory_env(_mem_env *, int size);
void free_memory_env(_mem_env *, int size, void *);
void free_unused_bmemory_env(_mem_env *,int size);
void free_unused_bmemory_env_all(_mem_env *);
void free_bmemory_env(_mem_env *, int size);
void free_bmemory_env_all(_mem_env *);

#endif /* BMEMORY_H */
