#ifndef MESH_H
#define MESH_H

#include <mpi.h>

#include "constants.h"
#include "uthash.h"

/*
 * A SET OF TYPEDEF STRUCT TO READ A MESH
 */

/*
 * cell
 * 
 * x, y, z - center coordinate
 * current_rfn - refinment level of the first level: 0
 * index - identifier. It is used as a hash table key in tree.h 
 */
typedef struct cell {
  double x, y, z;
  int current_rfn;
  int index[SIZE_KEY];
} cell_t;

/*
 * layer - a group of cells with the same refinment level
 * 
 * spacing - cell lenght (DX=DY=DZ) of a layer
 * ncells - number of cells in a layer
 * cells - cells storage in a layer 
 */
typedef struct layer {
  double spacing;
  int ncells;
  cell_t *cells;
} layer_t;

/*
 * mesh - a group of layers
 * 
 * layers - layers storage 
 */
typedef struct mesh {
  layer_t layers[MESH_LEVELS];  
} mesh_t;
  


/*
 * A SET OF READ MESH FUNCTIONS
 */

/*
 * Initialize an uniform mesh automatically based on parameters defined in constants.h
 */
void mesh_uniformMesh(MPI_Comm comm, mesh_t **ptr_mesh, const double *spacing);

#endif