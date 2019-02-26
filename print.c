#include <mpi.h>
#include <visit_writer.h>
#include <stdio.h>

#include "constants.h"
#include "tree.h"
#include "uthash.h"
#include "print.h"

node_t *nodes = NULL;

void print_pointMesh(MPI_Comm comm, leaf_t **ptr_tree, unsigned int tc) {

  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  char fname[100];
  sprintf(fname, "../output/pointMesh%i_tc%u", process_rank, tc);
    
  const int NPTS = HASH_COUNT(*ptr_tree);
  
  /* 
   * Create some points and data to save. 
   */

//   float pts[NPTS][3], data[NPTS];
//   int nvars = 2;
//   int vardims[] = {1, 3};
//   const char *varnames[] = {"rfn", "position"};
//   float *vars[] = {data, (float*)pts};

  float pts[NPTS][3];
  float data1[NPTS], data2[NPTS];
  int nvars = 2;
  int vardims[] = {1, 1};
  const char *varnames[] = {"rfn", "rank"};
  float *vars[] = {data1, data2};
  
  int i = 0;
  leaf_t *iter_leaf, *tmp_leaf;
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {
    
    pts[i][0] = iter_leaf->fv.x;
    pts[i][1] = iter_leaf->fv.y;
    pts[i][2] = iter_leaf->fv.z;
    
    
    data1[i] = (float)iter_leaf->current_rfn;
    data2[i] = (float)process_rank;
    
    i++;
    
  }
  
  /* 
   * Pass the mesh and data to visit_writer. 
   */
  write_point_mesh(fname, 1, NPTS, (float*)pts, nvars, vardims, varnames, vars);

}


void print_structuredAMR(MPI_Comm comm, leaf_t **ptr_tree, node_t **ptr_nodes, const double spacing, unsigned int tc) {

  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  char fname[100];
  sprintf(fname, "../output/structuredAMR%i_tc%u", process_rank, tc);
    
  double xOrigin = ORIGIN_X; 
  double yOrigin = ORIGIN_Y;
  double zOrigin = ORIGIN_Z;
  
  int nzones = HASH_COUNT(*ptr_tree);
  
  int zonetypes[nzones];
  for(int i=0; i<nzones; i++) 
    zonetypes[i] = VISIT_HEXAHEDRON;
  
  int nnodes = 0;
  
  const int rfn = 1;
  if(rfn == 1) {
    nnodes = HASH_COUNT(*ptr_nodes);  
    if(nnodes != 0) {
      printf("\n\nnnodes != 0\n\n");
      HASH_CLEAR(hh, *ptr_nodes);
    }  
  }
  nnodes = 0;
  
  float zonal[nzones];
  int connectivity[8*nzones];
    
  int l = 0, c = 0;
  node_t *ptr_node;
    
  leaf_t *iter_leaf, *tmp_leaf;
  HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) { 
    
    zonal[l++] = (float)l;
    
    for(int i=0; i<=1; i++) {
      for(int j=0; j<=1; j++) {
	for(int k=0; k<=1; k++) {  
	  
	  double x = iter_leaf->fv.x + (i-0.5)*(iter_leaf->fv.spacing);
	  double y = iter_leaf->fv.y + (j-0.5)*(iter_leaf->fv.spacing);
	  double z = iter_leaf->fv.z + (k-0.5)*(iter_leaf->fv.spacing);
	  
	  nodeKey_t key;
	  key.i = (int)((x - xOrigin) / spacing);
	  key.j = (int)((y - yOrigin) / spacing);
	  key.k = (int)((z - zOrigin) / spacing);
	  
	  HASH_FIND(hh, *ptr_nodes, &key, sizeof(nodeKey_t), ptr_node);
	    
	  if(ptr_node == NULL) {     
	    ptr_node = malloc(sizeof(node_t));
	    nnodes++;
	    ptr_node->id = nnodes;
	    ptr_node->x = x;
	    ptr_node->y = y;
	    ptr_node->z = z;	    
	    HASH_ADD(hh, *ptr_nodes, key, sizeof(nodeKey_t), ptr_node);
	  }
	      
          connectivity[c++] = ptr_node->id;  
	  printf("\n\nconnectivity[%i] = %i\n\n", c-1, connectivity[c-1]);

	      
	}
      }
    }
        
  }

  
  float pts[3*nnodes], nodal[nnodes];
  
  l = 0;
  c = 0;
  
  node_t *iter_node, *tmp_node;  
  HASH_ITER(hh, *ptr_nodes, iter_node, tmp_node) {
    
     pts[c++] = (float)iter_node->x;
     pts[c++] = (float)iter_node->y;
     pts[c++] = (float)iter_node->z;
    
     nodal[l++] = (float)iter_node->id;
     
  }
  
  int nvars = 2;
  int vardims[] = {1, 1};
  int centering[] = {0, 1};
  const char *varnames[] = {"zonal", "nodal"};
  float *vars[] = {zonal, nodal};
   
  /* 
   * Pass the mesh and data to visit_writer. 
   */
  int useBinary = 1;
  write_unstructured_mesh(fname, useBinary, nnodes, pts, nzones, zonetypes, connectivity, nvars, vardims, centering, varnames, vars);

}
