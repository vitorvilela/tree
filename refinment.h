#ifndef REFINMENT_H
#define REFINMENT_H

#include <mpi.h>

#include "constants.h"
#include "tree.h"
#include "uthash.h"

/*
 * Copy some of Lagrangian elements from a coaser leaf (mother)
 * to a finer leaf (which is one of its eight children) - being built
 */
void refinment_inheritLag_toFiner(MPI_Comm comm, leaf_t **finer_leaf);
/*
 * Interpolate finite volume variables from a coarser leaf to a finer leaf - being built
 */
void refinment_fvInterpolate_toFiner(MPI_Comm comm, leaf_t **finer_leaf, leaf_t **coarser_leaf, const int i, const int j, const int k);
/*
 * Create eight finer leaves in the place of the coarser leaf
 */
void refinment_reallocFiner(MPI_Comm comm, leaf_t **ptr_tree, leaf_t *coarser_leaf);
/*
 * Copy all Lagrangian elements from a set of finer leaves (children)
 * to the coaser leaf (mother) - being built
 */
void refinment_inheritLag_toCoarser(MPI_Comm comm, leaf_t **coarser_leaf, leaf_t ***children_leaves);
/*
 * Interpolate finite volume variables from a set of finer leaves (children)
 * to a coarser leaf - being built
 */
void refinment_fvInterpolate_toCoarser(MPI_Comm comm, leaf_t **coarser_leaf, leaf_t ***children_leaves);
/*
 * Create a coarser leaf in the place of eight finer leaves (children)
 */
void refinment_reallocCoarser(MPI_Comm comm, leaf_t **ptr_tree, leaf_t *finer_leaf);
/*
 * Refinment criterion based on stretching following x, y or z direction
 */
void refinment_stretching(MPI_Comm comm, leaf_t **ptr_tree, const char type);
/*
 * Refinment criterion based on a box
 */
void refinment_setBox(MPI_Comm comm, leaf_t **ptr_tree, double *startPoint, double *endPoint);
/*
 * Dynamic refinment criterion based on a zPlaneSlope
 * It is a function of time counter
 */
void refinment_zPlaneSlope(MPI_Comm comm, leaf_t **ptr_tree, int time_counter); 
/*
 * Realloc the whole tree
 */
void refinment_treeRealloc(MPI_Comm comm, leaf_t **ptr_tree);
/*
 * Realloc a leaf
 */
void refinment_leafRealloc(MPI_Comm comm, leaf_t **ptr_tree, leaf_t *ptr_leaf);
/*
 * Check if coarser was already created in refinment process
 */
int refinment_coaserExists(MPI_Comm comm, leaf_t **ptr_tree, hashKey_t coarserKey);
/*
 * Verify aspect ratio of a leaf respect to their neighbors
 * and refine them if necessary (i.e., ratio > 2)
 */
void refinment_verifyAspectRatio(MPI_Comm comm, leaf_t **ptr_tree, const double *mesh_spacing);



#endif