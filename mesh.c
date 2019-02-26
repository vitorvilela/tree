#include <stdio.h>    /* fprintf() */
#include <stdlib.h>   /* exit() */
#include <math.h>
#include <mpi.h>

#include "mesh.h"

/*
 * Initialize an uniform mesh automatically based on parameters defined at constants.h
 *
 * comm - MPI communicator
 * ptr_mesh - a pointer to mesh storage
 * spacing - a pointer to layers spacing  
 */
void mesh_uniformMesh(MPI_Comm comm, mesh_t **ptr_mesh, const double *spacing) {
  
  int process_rank;
  MPI_Comm_rank(comm, &process_rank);
  
  *ptr_mesh = malloc(sizeof(mesh_t)); 
  int ncells = 0;
    
  (*ptr_mesh)->layers[0].spacing = spacing[0];
  (*ptr_mesh)->layers[0].ncells = NX*NY*NZ;
    
  ncells = (*ptr_mesh)->layers[0].ncells;
  (*ptr_mesh)->layers[0].cells = malloc(ncells*sizeof(cell_t));
  
  int c = 0;
  
  for(int i=1; i<=NX; i++) {
    for(int j=1; j<=NY; j++) {
      for(int k=1; k<=NZ; k++) {
	
	double x = 0.5*(spacing[0]) + (i-1)*(spacing[0]);
	double y = 0.5*(spacing[0]) + (j-1)*(spacing[0]);
	double z = 0.5*(spacing[0]) + (k-1)*(spacing[0]);
  
	(*ptr_mesh)->layers[0].cells[c].current_rfn = 0;	
	(*ptr_mesh)->layers[0].cells[c].x = x;
	(*ptr_mesh)->layers[0].cells[c].y = y;
	(*ptr_mesh)->layers[0].cells[c].z = z;
	
	(*ptr_mesh)->layers[0].cells[c].index[0] = i;
	(*ptr_mesh)->layers[0].cells[c].index[1] = j;
	(*ptr_mesh)->layers[0].cells[c].index[2] = k;
	
	for(int idx=3; idx<SIZE_KEY; idx++)	  
	  (*ptr_mesh)->layers[0].cells[c].index[idx] = 0;
	   
        c++;
	
      }
    }
  } 
  
}