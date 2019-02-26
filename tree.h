#ifndef TREE_H
#define TREE_H

#include <mpi.h>

#include "constants.h"
#include "mesh.h"
#include "uthash.h"


/*
 * finite volume
 * 
 * x, y, z - center coordinate
 * spacing - finite volume lenght (DX=DY=DZ)
 */
typedef struct finiteVolume {
  double x, y, z, spacing;   
  double center_vars[NUMBER_CENTER_VARS];  
  double face_vars[NUMBER_FACE_VARS];
  double fluxes[NUMBER_FLUXES];   
  double rates[NUMBER_RATES];
} finiteVolume_t;

/*
 * idLag. Hash table for Lagrangian elements
 * 
 * id - identifier
 */
typedef struct idLag {
  unsigned int id;
  UT_hash_handle hh;
} idLag_t;

/*
 * hashKey. Key of leaf hash table 
 */
typedef struct hashKey { 
  int index[SIZE_KEY]; 
} hashKey_t;

/*
 * leaf
 *
 * current_rfn - refinment level at present time step 
 * future_rfn - refinment level at future time step
 * fv - finite volume associated to this leaf
 * key - hash table key 
 * id_notional - storage for hash table of Lagrangian identifiers  
 */
typedef struct leaf {  
  int current_rfn, future_rfn;
  finiteVolume_t fv;  
  hashKey_t key;
  idLag_t *id_notional;
  UT_hash_handle hh;
} leaf_t;


/*
 * Global leaf hash table 
 */
extern leaf_t *tree;


/*
 * Check collision in Lagrangian identifier hash table of a leaf. Put further in lagrangian.c 
 */
int lag_checkCollision(MPI_Comm comm, idLag_t **ptr_idLag, unsigned int id);
/*
 * Check collision in leaf hash table
 */
int tree_checkCollision(MPI_Comm comm, leaf_t **ptr_leaf, hashKey_t searched_key);
/*
 * Initialize leaf hash table called tree based on mesh information (see mesh.h)
 */
void tree_initTree(MPI_Comm comm, mesh_t **ptr_mesh, leaf_t **ptr_tree);
/*
 * Search for the host leaf of a Lagrangian element
 */
leaf_t* tree_findLeaf(MPI_Comm comm, leaf_t **ptr_tree, const double x, const double y, const double z, const double *spacing);
/*
 * Print leaf hash table called tree on file
 */
void tree_printTree_onFile(MPI_Comm comm, leaf_t **ptr_tree, const unsigned int time_step);
/*
 * Print host leaf on file
 */
void tree_printLeaf_onFile(MPI_Comm comm, leaf_t **host_leaf, const double x, const double y, const double z);
/*
 * Print host leaf on screen
 */
void tree_printLeaf_onScreen(MPI_Comm comm, leaf_t **host_leaf, const double x, const double y, const double z);
/*
 * Clear leaf hash table called tree 
 */
void tree_clearTree(MPI_Comm comm, leaf_t **ptr_tree);
/*
 * Find the maximum value in a finite volume vector
 */
double tree_fvMaxValue(MPI_Comm comm, leaf_t **ptr_tree, const int idxValue);
/*
 * Find the minimum value in a finite volume vector
 */
double tree_fvMinValue(MPI_Comm comm, leaf_t **ptr_tree, const int idxValue);
/*
 * Set a value for a finite volume vector (which represents a variable)
 */
void tree_fvInit(MPI_Comm comm, finiteVolume_t *fv, const double x, const double y, const double z, const double spacing);
/*
 * Initialize the finite volume of a leaf
 */
void tree_fvSetVariable(MPI_Comm comm, leaf_t **ptr_tree, const int idxValue, const int type);


#endif