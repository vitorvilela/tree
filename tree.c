#include <stdio.h>    /* fprintf() */
#include <stdlib.h>   /* exit() */
#include <math.h>
#include <mpi.h>

#include "tree.h"

/*
 * Global leaf hash table 
 */
leaf_t *tree = NULL;

/*
 * Check collision in Lagrangian identifier hash table of a leaf. Put further in lagrangian.c 
 * 
 * comm - MPI communicator
 * ptr_idLag - pointer to Lagrangian identifier hash table of type idLag_t
 * id - searched key
 */
int lag_checkCollision(MPI_Comm comm, idLag_t **ptr_idLag, unsigned int id) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  int collision = 1;
  
  idLag_t *ptr_id = NULL;
  
  HASH_FIND(hh, *ptr_idLag, &id, sizeof(unsigned int), ptr_id);

  if(ptr_id == NULL) collision = 0;
  else {    
    fprintf(stderr, "\n\nCOLLISION OCCURS! See lagrangian.c\n\n"); 
    exit(1);
  }
  
  return collision;
  
}

/*
 * Check collision in leaf hash table
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * searched_key - searched key 
 */
int tree_checkCollision(MPI_Comm comm, leaf_t **ptr_tree, hashKey_t searched_key) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  int collision = 1;
  
  leaf_t *ptr_leaf = NULL;
  
  HASH_FIND(hh, *ptr_tree, &searched_key, sizeof(hashKey_t), ptr_leaf);

  if(ptr_leaf == NULL) collision = 0;
  else {    
    fprintf(stderr, "\n\nCOLLISION OCCURS! See tree.c\n\n");     
    exit(1);
  }
  
  return collision;
  
}

/*
 * Initialize leaf hash table called tree based on mesh information (see mesh.h)
 * 
 * comm - MPI communicator
 * ptr_mesh - mesh pointer
 * ptr_tree - pointer to leaf hash table called tree 
 */
void tree_initTree(MPI_Comm comm, mesh_t **ptr_mesh, leaf_t **ptr_tree) {
     
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
     
  for(int l=0; l<MESH_LEVELS; l++) {
    
    int ncells = (*ptr_mesh)->layers[l].ncells;
     
    for(int c=0; c<ncells; c++) {
      
      leaf_t *ptr_leaf = malloc(sizeof(leaf_t));
      
      ptr_leaf->current_rfn = (*ptr_mesh)->layers[l].cells[c].current_rfn;
      ptr_leaf->future_rfn = ptr_leaf->current_rfn;
      ptr_leaf->id_notional = NULL;
      for(int i=0; i<SIZE_KEY; i++) 
        ptr_leaf->key.index[i] = (*ptr_mesh)->layers[l].cells[c].index[i];
      
      /*
       * Init the finite volume of a leaf using center coordinate and fv lenght
       */      
      double x = (*ptr_mesh)->layers[l].cells[c].x;
      double y = (*ptr_mesh)->layers[l].cells[c].y;
      double z = (*ptr_mesh)->layers[l].cells[c].z;
      double spacing = (*ptr_mesh)->layers[l].spacing;
      tree_fvInit(comm, &(ptr_leaf->fv), x, y, z, spacing);

      /*
       * Add a leaf into the hash table
       */      
      int collision = tree_checkCollision(comm, ptr_tree, ptr_leaf->key);
      if(collision == 0)
	HASH_ADD(hh, *ptr_tree, key, sizeof(hashKey_t), ptr_leaf);
      
    }
    
  }
      
}

/*
 * Search for the host leaf of a Lagrangian element  
 * PS. The Lagrangian transport functions must guarantee the particle is not exactly at domain boundary
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * x, y, z - Lagrangian element position or leaf center coordinate
 * spacing - a pointer to spacing of the finite volume of a leaf
 */
leaf_t* tree_findLeaf(MPI_Comm comm, leaf_t **ptr_tree, const double x, const double y, const double z, const double *spacing) {
  
  int process_rank;
  MPI_Comm_rank(comm, &process_rank);
    
  leaf_t *ptr_leaf = NULL;
  
  hashKey_t searched_key;
  for(int i=0; i<SIZE_KEY; i++)
    searched_key.index[i] = 0;
  
  double xOrigin = ORIGIN_X; 
  double yOrigin = ORIGIN_Y;
  double zOrigin = ORIGIN_Z;
  
  int l = 0, i = -1;
  
  while(l < RFN_LEVELS) {
    
    int idx = (int)(1 + (x - xOrigin) / spacing[l]);
    int jdx = (int)(1 + (y - yOrigin) / spacing[l]);
    int kdx = (int)(1 + (z - zOrigin) / spacing[l]);
      
    searched_key.index[++i] = idx;
    searched_key.index[++i] = jdx;
    searched_key.index[++i] = kdx;    
            
    HASH_FIND(hh, *ptr_tree, &searched_key, sizeof(hashKey_t), ptr_leaf);
    
    if(ptr_leaf != NULL) break;
    else {      
      xOrigin += (idx-1)*spacing[l];
      yOrigin += (jdx-1)*spacing[l];
      zOrigin += (kdx-1)*spacing[l];
      l++;
    }
    
  }
  
  return ptr_leaf;
  
}

/*
 * Print leaf hash table called tree on file
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * time_step - time step
 */
void tree_printTree_onFile(MPI_Comm comm, leaf_t **ptr_tree, const unsigned int time_step) {
 
  int process_rank;
  MPI_Comm_rank(comm, &process_rank);
  
  int counter_leaves = 1;
  
  int number_leaves = HASH_COUNT(*ptr_tree);
  
  char fname[100];
  sprintf(fname, "../output/tree_%i_print_%u.tree", process_rank, time_step);
  FILE *ptr_file = fopen(fname, "w+");
  fprintf(ptr_file, "TREE HAS %i LEAVES\n", number_leaves);
  
  leaf_t *iter_leaf, *tmp_leaf;
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {    
    int number_lagrangian = (iter_leaf->id_notional==NULL?0:HASH_COUNT(iter_leaf->id_notional));
    fprintf(ptr_file, "\nLEAF %i\n", counter_leaves);
    fprintf(ptr_file, "Spacing %f\n", iter_leaf->fv.spacing);
    fprintf(ptr_file, "Current %i and future %i refinment\n", iter_leaf->current_rfn, iter_leaf->future_rfn);
    fprintf(ptr_file, "Lagrangian %i\n", number_lagrangian);
    fprintf(ptr_file, "Key ");
    for(int i=0; i<SIZE_KEY; i++) 
      fprintf(ptr_file, "%i ", iter_leaf->key.index[i]);
    fprintf(ptr_file, "\n");
    counter_leaves++;
  }  
  
  fclose(ptr_file);
  
}

/*
 * Print host leaf on file
 * 
 * comm - MPI communicator
 * host_leaf - pointer to host leaf
 * x, y, z - host leaf center coordinate
 */
void tree_printLeaf_onFile(MPI_Comm comm, leaf_t **host_leaf, const double x, const double y, const double z) {
 
  int process_rank;
  MPI_Comm_rank(comm, &process_rank);
    
  char fname[100];
  sprintf(fname, "../output/hostLeaf_%i.leaf", process_rank);
  FILE *ptr_file = fopen(fname, "w+");

  if((*host_leaf) == NULL)
    fprintf(ptr_file, "Particle not found");
  else {
    int number_lagrangian = ((*host_leaf)->id_notional==NULL?0:HASH_COUNT((*host_leaf)->id_notional));
    fprintf(ptr_file, "HOST LEAF FOR LAGRANGIAN %1.3f %1.3f %1.3f\n", x, y, z);
    fprintf(ptr_file, "Spacing %f\n", (*host_leaf)->fv.spacing);
    fprintf(ptr_file, "Current %i and future %i refinment\n", (*host_leaf)->current_rfn, (*host_leaf)->future_rfn);
    fprintf(ptr_file, "Lagrangian %i\n", number_lagrangian);
    fprintf(ptr_file, "Key ");
    for(int i=0; i<SIZE_KEY; i++) 
      fprintf(ptr_file, "%i ", (*host_leaf)->key.index[i]);
    fprintf(ptr_file, "\n");
  }
  
  fclose(ptr_file);
  
}

/*
 * Print host leaf on screen
 * 
 * comm - MPI communicator
 * host_leaf - pointer to host leaf
 * x, y, z - host leaf center coordinate
 */
void tree_printLeaf_onScreen(MPI_Comm comm, leaf_t **host_leaf, const double x, const double y, const double z) {
 
  int process_rank;
  MPI_Comm_rank(comm, &process_rank);
  
  if((*host_leaf) == NULL)
    printf("\n\nParticle not found\n\n");
  else {
    int number_lagrangian = ((*host_leaf)->id_notional==NULL?0:HASH_COUNT((*host_leaf)->id_notional));
    printf("HOST LEAF FOR LAGRANGIAN %1.3f %1.3f %1.3f\n", x, y, z);
    printf("Spacing %f\n", (*host_leaf)->fv.spacing);
    printf("Current %i and future %i refinment\n", (*host_leaf)->current_rfn, (*host_leaf)->future_rfn);
    printf("Lagrangian %i\n", number_lagrangian);
    printf("Key ");
    for(int i=0; i<SIZE_KEY; i++) 
      printf("%i ", (*host_leaf)->key.index[i]);
    printf("\n\n");
  }
    
}

/*
 * Clear leaf hash table called tree 
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 */
void tree_clearTree(MPI_Comm comm, leaf_t **ptr_tree) {
  
  int process_rank;
  MPI_Comm_rank(comm, &process_rank);
  
  leaf_t *iter_leaf, *tmp_leaf;
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {
    HASH_DEL(*ptr_tree, iter_leaf);
    free(iter_leaf);
  }
  
}

/*
 * Find the maximum value in a finite volume vector
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * idxValue - index that identifies which vector will be used 
 */
double tree_fvMaxValue(MPI_Comm comm, leaf_t **ptr_tree, const int idxValue) {
  
  int process_rank;  
  MPI_Comm_rank(comm, &process_rank);
    
  double max = -10000.0;
  
  leaf_t *iter_leaf, *tmp_leaf;
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {    
    double value = iter_leaf->fv.center_vars[idxValue];    
    if(value > max) max = value; 
  }
  
  return max;
  
}

/*
 * Find the minimum value in a finite volume vector
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * idxValue - index that identifies which vector will be used  
 */
double tree_fvMinValue(MPI_Comm comm, leaf_t **ptr_tree, const int idxValue) {
  
  int process_rank;  
  MPI_Comm_rank(comm, &process_rank);
    
  double min = 10000.0;
  
  leaf_t *iter_leaf, *tmp_leaf;
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {    
    double value = iter_leaf->fv.center_vars[idxValue];    
    if(value < min) min = value; 
  }
  
  return min;
  
}

/*
 * Set a value for a finite volume vector (which represents a variable)
 *
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * idxValue - index that identifies which vector will be used 
 * type - type of set value
 */
void tree_fvSetVariable(MPI_Comm comm, leaf_t **ptr_tree, const int idxValue, const int type) {
  
  int process_rank;  
  MPI_Comm_rank(comm, &process_rank);
  
  leaf_t *iter_leaf, *tmp_leaf;
  
  switch(type) {  
    case 0:
      HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) { 
	iter_leaf->fv.center_vars[idxValue] = 0.0;
      }
      break;    
    case 1:      
      HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) { 
	iter_leaf->fv.center_vars[idxValue] = iter_leaf->fv.x;
      }
      break;      
    case 2:      
      HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {
	iter_leaf->fv.center_vars[idxValue] = iter_leaf->fv.y;
      }
      break;	
  }
  
}

/*
 * Initialize the finite volume of a leaf
 * 
 * comm - MPI communicator
 * fv - pointer to finite volume
 * x, y, z - finite volume center coordinate
 * spacing - a pointer to spacing of the finite volume 
 */
void tree_fvInit(MPI_Comm comm, finiteVolume_t *fv, const double x, const double y, const double z, const double spacing) {
  
  int process_rank;  
  MPI_Comm_rank(comm, &process_rank);
  
  
  
  fv->x = x;
  fv->y = y;
  fv->z = z;
  fv->spacing = spacing;
  
  
  for(int i=0; i<NUMBER_CENTER_VARS; i++) fv->center_vars[i] = 0.0;
  for(int i=0; i<NUMBER_FACE_VARS; i++) fv->face_vars[i] = 0.0;
  for(int i=0; i<NUMBER_FLUXES; i++) fv->fluxes[i] = 0.0;
  for(int i=0; i<NUMBER_RATES; i++) fv->rates[i] = 0.0;
  
    
}
