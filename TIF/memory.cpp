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
 *  $Id: memory.c,v 1.6 2011/12/28 08:48:56 tanaka Rel $
 *  $Revision: 1.6 $
 *  $Date: 2011/12/28 08:48:56 $
 *  $Author: tanaka $
 *
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "memory.h"

/*
 * xmalloc(size)
 *   allocates memory space of size bytes and returns a pointer to
 *   the allocated memory.
 *      size: bytes to be allocated
 *
 */
void *xmalloc(size_t size)
{
	void *ptr  = NULL;
	try{
		ptr = malloc(size);
	}catch(int a){

	}

  if(ptr == NULL) {
    fprintf(stderr, "xmalloc(): memory allocation failed (%u).\n",
	    (unsigned int) size);
    exit(1007);
  }

  return(ptr);
}

/*
 * xcalloc(nmemb, size)
 *   allocates an array of nmemb elements of size bytes each and
 *   returns a pointer to the allocated memory.  The memory is set
 *   to zero.
 *     nmemb: number of elements
 *      size: size of each element
 *
 */
void *xcalloc(size_t nmemb, size_t size)
{
  void *ptr = calloc(nmemb, size);

  if(ptr == NULL) {
    fprintf(stderr, "xcalloc(): memory allocation failed (%u, %u).\n",
	    (unsigned int) nmemb, (unsigned int) size);
    exit(1007);
  }

  return(ptr);
}

/*
 * xrealloc(ptr, size)
 *   changes the size of the previously allocated memory space and returns
 *   a pointer to the allocated memory.
 *       ptr: a pointer to the previously allocated memory
 *      size: size of each element
 *
 */
void *xrealloc(void *p, size_t size)
{
  void *ptr = realloc(p, size);

  if(ptr == NULL && size != 0) {
    fprintf(stderr, "xrealloc(): memory allocation failed (%p, %u).\n",
	    p, (unsigned int) size);
    exit(1007);
  }

  return(ptr);
}

/*
 * xfree(ptr)
 *   frees the memory space pointed to by ptr.
 *       ptr: pointer to the memory space to be freed.
 *
 */
void xfree(void *ptr)
{
  free(ptr);
}
