#include <stdio.h>    /* fprintf() */
#include <stdlib.h>   /* exit() */
#include <math.h>
#include <mpi.h>

#include "refinment.h"
#include "tree.h"

/*
 * Copy some of Lagrangian elements from a coaser leaf (mother)
 * to a finer leaf (which is one of its eight children) - being built
 * 
 * comm - MPI communicator
 * finer_leaf - pointer to a finer leaf 
 */
void refinment_inheritLag_toFiner(MPI_Comm comm, leaf_t **finer_leaf) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
 
  /*
   * TODO Temporary inheritance - being built
   */
  (*finer_leaf)->id_notional = NULL;   
  
}

/*
 * Interpolate finite volume variables from a coarser leaf to a finer leaf - being built
 * 
 * comm - MPI communicator
 * finer_leaf - pointer to a finer leaf 
 * coarser_leaf - pointer to a coarser leaf
 * i, j, k - finer leaf index
 */
void refinment_fvInterpolate_toFiner(MPI_Comm comm, leaf_t **finer_leaf, leaf_t **coarser_leaf, const int i, const int j, const int k) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  double xOrigin = (*coarser_leaf)->fv.x - 0.5*((*coarser_leaf)->fv.spacing);
  double yOrigin = (*coarser_leaf)->fv.y - 0.5*((*coarser_leaf)->fv.spacing);
  double zOrigin = (*coarser_leaf)->fv.z - 0.5*((*coarser_leaf)->fv.spacing);

  double spacing = 0.5*((*coarser_leaf)->fv.spacing);
  
  double x = xOrigin + 0.5*spacing + (i-1)*spacing;
  double y = yOrigin + 0.5*spacing + (j-1)*spacing;
  double z = zOrigin + 0.5*spacing + (k-1)*spacing;
  
  /*
   * TODO Temporary interpolation - being built
   */
  tree_fvInit(comm, &((*finer_leaf)->fv), x, y, z, spacing);
    
}

/*
 * Create eight finer leaves in the place of the coarser leaf
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * coarser_leaf - pointer to a coarser leaf
 */
void refinment_reallocFiner(MPI_Comm comm, leaf_t **ptr_tree, leaf_t *coarser_leaf) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
   
  const int start_index = 3*(coarser_leaf->future_rfn);
  
  for(int k=1; k<=2; k++) {
    for(int j=1; j<=2; j++) {
      for(int i=1; i<=2; i++) {
	
	leaf_t *finer_leaf = malloc(sizeof(leaf_t));
	hashKey_t finerKey = coarser_leaf->key;
	finerKey.index[start_index] = i;
	finerKey.index[start_index+1] = j;	
	finerKey.index[start_index+2] = k; 	
	finer_leaf->key = finerKey;  
	HASH_ADD(hh, *ptr_tree, key, sizeof(hashKey_t), finer_leaf);
	
	finer_leaf->current_rfn = coarser_leaf->future_rfn;
	finer_leaf->future_rfn = finer_leaf->current_rfn;
	 	
	refinment_fvInterpolate_toFiner(comm, &finer_leaf, &coarser_leaf, i, j, k);
     
	refinment_inheritLag_toFiner(comm, &finer_leaf);
	
      }
    }
  }
  
  /*
   * TODO Here is the correct position of function call
   * when it is finished - pay attention to &coarser_leaf
   * refinment_inheritLag_toFiner(comm, &coarser_leaf);
   */
  
  HASH_DEL(*ptr_tree, coarser_leaf);
  free(coarser_leaf);   
  
}

/*
 * Copy all Lagrangian elements from a set of finer leaves (children)
 * to the coaser leaf (mother) - being built
 *  
 * comm - MPI communicator
 * coarser_leaf - pointer to a coarser leaf
 * children_leaves - pointer to a set of finer leaves (children) 
 */
void refinment_inheritLag_toCoarser(MPI_Comm comm, leaf_t **coarser_leaf, leaf_t ***children_leaves) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
 
  /* 
   * TODO Temporary inheritance - being built
   */
  (*coarser_leaf)->id_notional = NULL;   
  
}

/*
 * Interpolate finite volume variables from a set of finer leaves (children)
 * to a coarser leaf - being built
 * 
 * comm - MPI communicator
 * coarser_leaf - pointer to a coarser leaf
 * children_leaves - pointer to a set of finer leaves (children) 
 */
void refinment_fvInterpolate_toCoarser(MPI_Comm comm, leaf_t **coarser_leaf, leaf_t ***children_leaves) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
    
  double minX = 10000.0;
  double maxX = -10000.0;
  double minY = 10000.0;
  double maxY = -10000.0;
  double minZ = 10000.0;
  double maxZ = -10000.0;
  
  for(int i=0; i<8; i++) {
    minX = (((*children_leaves)[i]->fv.x)<minX?((*children_leaves)[i]->fv.x):minX);
    maxX = (((*children_leaves)[i]->fv.x)>maxX?((*children_leaves)[i]->fv.x):maxX); 
    minY = (((*children_leaves)[i]->fv.y)<minY?((*children_leaves)[i]->fv.y):minY);
    maxY = (((*children_leaves)[i]->fv.y)>maxY?((*children_leaves)[i]->fv.y):maxY); 
    minZ = (((*children_leaves)[i]->fv.z)<minZ?((*children_leaves)[i]->fv.z):minZ);
    maxZ = (((*children_leaves)[i]->fv.z)>maxZ?((*children_leaves)[i]->fv.z):maxZ);   
  }

  double spacing = 2.0*(*children_leaves)[0]->fv.spacing;
  double x = 0.5*(minX + maxX);
  double y = 0.5*(minY + maxY);
  double z = 0.5*(minZ + maxZ);
   
  /* 
   * TODO Temporary interpolation - being built
   */
  tree_fvInit(comm, &((*coarser_leaf)->fv), x, y, z, spacing);
    
}


/*
 * Create a coarser leaf in the place of eight finer leaves (children)
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * finer_leaf - pointer to a finer leaf
 */
void refinment_reallocCoarser(MPI_Comm comm, leaf_t **ptr_tree, leaf_t *finer_leaf) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
    
        
  leaf_t **children_leaves = malloc(8*sizeof(leaf_t *));
    
  hashKey_t coarserKey = finer_leaf->key;
 
  
  const int start_index = 3*(finer_leaf->current_rfn);
  for(int i=start_index; i<SIZE_KEY; i++) coarserKey.index[i] = 0;  
   
  int exists = refinment_coaserExists(comm, ptr_tree, coarserKey);

  if(exists == 1) {
    HASH_DEL(*ptr_tree, finer_leaf);
    free(finer_leaf);
  }
  else {  
    
    leaf_t *coarser_leaf = malloc(sizeof(leaf_t));
    coarser_leaf->key = coarserKey;  
    HASH_ADD(hh, *ptr_tree, key, sizeof(hashKey_t), coarser_leaf);    
    coarser_leaf->current_rfn = finer_leaf->future_rfn;
    coarser_leaf->future_rfn = coarser_leaf->current_rfn; 
    
    hashKey_t neighborKey = coarserKey;
    int cl = 0;  
    for(int k=1; k<=2; k++) {
      for(int j=1; j<=2; j++) {
	for(int i=1; i<=2; i++) {
	  neighborKey.index[start_index] = i;
	  neighborKey.index[start_index+1] = j;
	  neighborKey.index[start_index+2] = k;	
	  HASH_FIND(hh, *ptr_tree, &neighborKey, sizeof(hashKey_t), children_leaves[cl]);
	  
	  /*
	   * Check if children exists or is still refined
	   */
	  if(children_leaves[cl] != NULL) {
	    cl++;
	  }
	  else {	    
	    hashKey_t finestKey = neighborKey;
	    const int finest_start_index = 3*(coarser_leaf->future_rfn + 1);
	    finestKey.index[finest_start_index] = 1;
	    finestKey.index[finest_start_index+1] = 1;
	    finestKey.index[finest_start_index+2] = 1;	    
	    
	    leaf_t *finest_leaf;
	    HASH_FIND(hh, *ptr_tree, &finestKey, sizeof(hashKey_t), finest_leaf);
	    refinment_reallocCoarser(comm, ptr_tree, finest_leaf);

            HASH_FIND(hh, *ptr_tree, &neighborKey, sizeof(hashKey_t), children_leaves[cl]);
	    cl++;		    
	  }
	}
      }
    }  
    
    refinment_inheritLag_toCoarser(comm, &coarser_leaf, &children_leaves);
    
    refinment_fvInterpolate_toCoarser(comm, &coarser_leaf, &children_leaves);
    
    HASH_DEL(*ptr_tree, finer_leaf);
    free(finer_leaf);
    
  }
  
}

/*
 * Refinment criterion based on stretching following x, y or z direction
 * TODO Maintenance
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * type - stretching type: x, y or z
 */
void refinment_stretching(MPI_Comm comm, leaf_t **ptr_tree, const char type) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  /*
   * TODO Temporary
   */
    
  double valueMax, valueMin;
  
  switch(type) {    
    case 'x':
      valueMax = LX;
      valueMin = ORIGIN_X;
      break;      
    case 'y':
      valueMax = LY;
      valueMin = ORIGIN_Y;
      break;      
    case 'z':
      valueMax = LZ;
      valueMin = ORIGIN_Z;      
      break;      
  }
  
  double dValue = fabs(valueMax - valueMin);
  double dRfn = dValue/RFN_LEVELS;
   
  double *rfnValues = malloc((RFN_LEVELS-1)*sizeof(double));
  rfnValues[0] = valueMin + dRfn; 
  for(int i=1; i<(RFN_LEVELS-1); i++) {
    rfnValues[i] = rfnValues[i-1] + dRfn; 
  }
     
  for(int t=1; t<RFN_LEVELS; t++) {
     
    leaf_t *iter_leaf, *tmp_leaf;
    HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {
    
      double value;
      
      switch(type) {    
	case 'x':
	  value = iter_leaf->fv.x;     
	  break;      
	case 'y':
	  value = iter_leaf->fv.y;     
	  break;      
	case 'z':
	  value = iter_leaf->fv.z;     
	  break;       
      }
     
      int l = 0;     
          
      for(int i=0; i<RFN_LEVELS-1; i++) {    
	if(value > rfnValues[i]) l++;  
	else break;
      }
          
      if(l > (iter_leaf->current_rfn)) {      
	iter_leaf->future_rfn = iter_leaf->current_rfn + 1;
	refinment_leafRealloc(comm, ptr_tree, iter_leaf);
      }
      else if(l < (iter_leaf->current_rfn)) {	
	iter_leaf->future_rfn = iter_leaf->current_rfn - 1;
	refinment_leafRealloc(comm, ptr_tree, iter_leaf);
      }
      
    }      
    
  }  
  
} 

/*
 * Refinment criterion based on a box
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * starPoint - initial x, y and z
 * endPoint - final x, z and z
 */
void refinment_setBox(MPI_Comm comm, leaf_t **ptr_tree, double *startPoint, double *endPoint) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
 
   
  leaf_t *iter_leaf, *tmp_leaf;
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {
    
   
    /*
     * LB - lower bound
     * UB - upper bound 
     */
    double leaf_xLB = (iter_leaf->fv.x) - 0.5*(iter_leaf->fv.spacing);
    double leaf_xUB = (iter_leaf->fv.x) + 0.5*(iter_leaf->fv.spacing);
    double leaf_yLB = (iter_leaf->fv.y) - 0.5*(iter_leaf->fv.spacing);
    double leaf_yUB = (iter_leaf->fv.y) + 0.5*(iter_leaf->fv.spacing);
    double leaf_zLB = (iter_leaf->fv.z) - 0.5*(iter_leaf->fv.spacing);
    double leaf_zUB = (iter_leaf->fv.z) + 0.5*(iter_leaf->fv.spacing);
    
    
    if(leaf_xLB>=startPoint[0] && leaf_yLB>=startPoint[1] && leaf_zLB>=startPoint[2] && leaf_xUB<=endPoint[0] && leaf_yUB<=endPoint[1] && leaf_zUB<=endPoint[2]) {
      if((iter_leaf->current_rfn) < (RFN_LEVELS-1)) {
        iter_leaf->future_rfn = iter_leaf->current_rfn + 1;
      }
      
    }
    
    
  }
  
  refinment_treeRealloc(comm, ptr_tree);
      
}

/*
 * Dynamic refinment criterion based on a zPlaneSlope
 * It is a function of time counter
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * time_counter - integer time counter
 */
void refinment_zPlaneSlope(MPI_Comm comm, leaf_t **ptr_tree, int time_counter) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  const double xStart = ORIGIN_X;
  const double yStart = ORIGIN_Y;
  
  double X = (double)(time_counter)/NX;
  const double slope = 1.0;
  double Y = yStart + slope*(X - xStart);
  const double range = LX/NX;
 
  leaf_t *iter_leaf, *tmp_leaf;
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {
   
    double xLowerBound = (iter_leaf->fv.x) - 0.5*(iter_leaf->fv.spacing);
    double xUpperBound = (iter_leaf->fv.x) + 0.5*(iter_leaf->fv.spacing);
    double yLowerBound = (iter_leaf->fv.y) - 0.5*(iter_leaf->fv.spacing);
    double yUpperBound = (iter_leaf->fv.y) + 0.5*(iter_leaf->fv.spacing);
     
    if(xUpperBound<=X && yUpperBound<=Y && xLowerBound>=(X-range) && yLowerBound>=(Y-range)) {
      if((iter_leaf->current_rfn) < (RFN_LEVELS-1)) 
	iter_leaf->future_rfn = iter_leaf->current_rfn + 1;     
    }
    else {
      if((iter_leaf->current_rfn) >= 1)
	iter_leaf->future_rfn = iter_leaf->current_rfn - 1;       
    }  
       
  }
  
  refinment_treeRealloc(comm, ptr_tree);
  
}
 
/*
 * Realloc the whole tree
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 */
void refinment_treeRealloc(MPI_Comm comm, leaf_t **ptr_tree) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
      
  leaf_t *iter_leaf, *tmp_leaf;
  
  /*
   * TODO This hash_iter also loops over the new created leaves.
   * It can be used to create an algorithm which set a future_rfn=current_rfn+2, 
   * and fragments the refinment process, for example. 
   */
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {     
    refinment_leafRealloc(comm, ptr_tree, iter_leaf);   
  }
  
}

/*
 * Realloc a leaf
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * ptr_leaf - pointer to a leaf
 */
void refinment_leafRealloc(MPI_Comm comm, leaf_t **ptr_tree, leaf_t *ptr_leaf) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
    
  if(ptr_leaf->future_rfn > ptr_leaf->current_rfn) {
    refinment_reallocFiner(comm, ptr_tree, ptr_leaf);
  }  
  else if(ptr_leaf->future_rfn < ptr_leaf->current_rfn) {    
    refinment_reallocCoarser(comm, ptr_tree, ptr_leaf);
  }
  
  
}  

/*
 * Check if coarser was already created in refinment process
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * coarserKey - searched coarser key 
 */
int refinment_coaserExists(MPI_Comm comm, leaf_t **ptr_tree, hashKey_t coarserKey) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  int exists = 1;
  
  leaf_t *ptr_leaf = NULL;
  
  HASH_FIND(hh, *ptr_tree, &coarserKey, sizeof(hashKey_t), ptr_leaf);

  if(ptr_leaf == NULL) exists = 0;
   
  return exists;
  
}

/*
 * Verify aspect ratio of a leaf respect to their neighbors
 * and refine them if necessary (i.e., ratio > 2)
 * 
 * comm - MPI communicator
 * ptr_tree - pointer to leaf hash table called tree 
 * mesh_spacing - layers spacing
 */
void refinment_verifyAspectRatio(MPI_Comm comm, leaf_t **ptr_tree, const double *mesh_spacing) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  leaf_t *iter_leaf, *tmp_leaf, *neighbor_leaf;
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {
    
    if(iter_leaf->current_rfn >= 2) {
      
      double leaf_xi = iter_leaf->fv.x - 0.5*(iter_leaf->fv.spacing);
      double leaf_xf = iter_leaf->fv.x + 0.5*(iter_leaf->fv.spacing);
      double leaf_yi = iter_leaf->fv.y - 0.5*(iter_leaf->fv.spacing);
      double leaf_yf = iter_leaf->fv.y + 0.5*(iter_leaf->fv.spacing);
      double leaf_zi = iter_leaf->fv.z - 0.5*(iter_leaf->fv.spacing);
      double leaf_zf = iter_leaf->fv.z + 0.5*(iter_leaf->fv.spacing);
      
      if(leaf_xi != ORIGIN_X) {	
	double neighbor_x = leaf_xi - 0.1*(iter_leaf->fv.spacing);
	for(int j=-1; j<=1; j++) {
	  for(int k=-1; k<=1; k++) {
	    double neighbor_y = iter_leaf->fv.y + 0.6*j*(iter_leaf->fv.spacing);
	    double neighbor_z = iter_leaf->fv.z + 0.6*k*(iter_leaf->fv.spacing);	  
	    neighbor_leaf = tree_findLeaf(comm, ptr_tree, neighbor_x, neighbor_y, neighbor_z, mesh_spacing);	    
	    if(neighbor_leaf->current_rfn < iter_leaf->current_rfn-1) {
	      neighbor_leaf->future_rfn = neighbor_leaf->future_rfn + 1;
	      refinment_leafRealloc(comm, ptr_tree, neighbor_leaf);
	    }
	  }	   
	}	
      }
      
      if(leaf_xf != (LX-ORIGIN_X)) {
	double neighbor_x = leaf_xf + 0.1*(iter_leaf->fv.spacing);
	for(int j=-1; j<=1; j++) {
	  for(int k=-1; k<=1; k++) {
	    double neighbor_y = iter_leaf->fv.y + 0.6*j*(iter_leaf->fv.spacing);
	    double neighbor_z = iter_leaf->fv.z + 0.6*k*(iter_leaf->fv.spacing); 
	    neighbor_leaf = tree_findLeaf(comm, ptr_tree, neighbor_x, neighbor_y, neighbor_z, mesh_spacing);	    
	    if(neighbor_leaf->current_rfn < iter_leaf->current_rfn-1) {
	      neighbor_leaf->future_rfn = neighbor_leaf->future_rfn + 1;
	      refinment_leafRealloc(comm, ptr_tree, neighbor_leaf);
	    }
	  }	   
	}		
      }
      
      if(leaf_yi != ORIGIN_Y) {
	double neighbor_y = leaf_yi - 0.1*(iter_leaf->fv.spacing);
	double neighbor_x = iter_leaf->fv.x;
	for(int k=-1; k<=1; k++) {
	  double neighbor_z = iter_leaf->fv.z + 0.6*k*(iter_leaf->fv.spacing); 
	  neighbor_leaf = tree_findLeaf(comm, ptr_tree, neighbor_x, neighbor_y, neighbor_z, mesh_spacing);	    
	  if(neighbor_leaf->current_rfn < iter_leaf->current_rfn-1) {
	    neighbor_leaf->future_rfn = neighbor_leaf->future_rfn + 1;
	    refinment_leafRealloc(comm, ptr_tree, neighbor_leaf);
	  }
	}	
      }
      
      if(leaf_yf != (LY-ORIGIN_Y)) {
	double neighbor_y = leaf_yf + 0.1*(iter_leaf->fv.spacing);
	double neighbor_x = iter_leaf->fv.x;
	for(int k=-1; k<=1; k++) {
	  double neighbor_z = iter_leaf->fv.z + 0.6*k*(iter_leaf->fv.spacing); 
	  neighbor_leaf = tree_findLeaf(comm, ptr_tree, neighbor_x, neighbor_y, neighbor_z, mesh_spacing);	    
	  if(neighbor_leaf->current_rfn < iter_leaf->current_rfn-1) {
	    neighbor_leaf->future_rfn = neighbor_leaf->future_rfn + 1;	  
	    refinment_leafRealloc(comm, ptr_tree, neighbor_leaf);
	  }
	}	
      }
      
      if(leaf_zi != ORIGIN_Z) {
	double neighbor_z = leaf_zi - 0.1*(iter_leaf->fv.spacing);
	double neighbor_x = iter_leaf->fv.x;
	for(int j=-1; j<=1; j++) {
	  double neighbor_y = iter_leaf->fv.y + 0.6*j*(iter_leaf->fv.spacing);
	  neighbor_leaf = tree_findLeaf(comm, ptr_tree, neighbor_x, neighbor_y, neighbor_z, mesh_spacing);	    
	  if(neighbor_leaf->current_rfn < iter_leaf->current_rfn-1) {
	    neighbor_leaf->future_rfn = neighbor_leaf->future_rfn + 1;	 
	    refinment_leafRealloc(comm, ptr_tree, neighbor_leaf);
	  }
	}
      }
      
      if(leaf_zf != (LZ-ORIGIN_Z)) {
	double neighbor_z = leaf_zf + 0.1*(iter_leaf->fv.spacing);
	double neighbor_x = iter_leaf->fv.x;
	for(int j=-1; j<=1; j++) {
	  double neighbor_y = iter_leaf->fv.y + 0.6*j*(iter_leaf->fv.spacing);
	  neighbor_leaf = tree_findLeaf(comm, ptr_tree, neighbor_x, neighbor_y, neighbor_z, mesh_spacing);	    
	  if(neighbor_leaf->current_rfn < iter_leaf->current_rfn-1) {
	    neighbor_leaf->future_rfn = neighbor_leaf->future_rfn + 1;
	    refinment_leafRealloc(comm, ptr_tree, neighbor_leaf); 
	  }
	}
      }
                  
    }
        
  }
  
  
  
}



























