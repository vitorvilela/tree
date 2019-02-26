#ifndef PRINT_H
#define PRINT_H

#include <mpi.h>

#include "uthash.h"


/*
 * 
 */
typedef struct nodeKey {
  int i, j, k;
} nodeKey_t;

/*
 * 
 */
typedef struct node {
  nodeKey_t key;
  double x, y, z;
  int id;
  UT_hash_handle hh;
} node_t;

/*
 * 
 */
extern node_t *nodes;

/*
 * 
 */
void print_pointMesh(MPI_Comm comm, leaf_t **ptr_tree, unsigned int tc);

void print_structuredAMR(MPI_Comm comm, leaf_t **ptr_tree, node_t **ptr_nodes, const double spacing, unsigned int tc);
 
#endif