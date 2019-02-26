/*
 * main.c
 *
 *  Created on: Jun 15, 2015
 *      Author: Vitor Vilela
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "distributions.h"
#include "random.h"
#include "mesh.h"
#include "tree.h"
#include "refinment.h"
#include "print.h"

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);  
  
  int process_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  
  /*
   * Create a file per process to monitor the code algorithm
   */   
  monitor_constants(MPI_COMM_WORLD);
  
  /*
   * Creating the mesh and the tree
   */  
  double *mesh_spacing = malloc(RFN_LEVELS*sizeof(double));
  mesh_spacing[0] = LX/NX;
  for(int i=1; i<RFN_LEVELS; i++) {
    mesh_spacing[i] = mesh_spacing[i-1]/2;
  }
       
  mesh_t *mesh = NULL;  
  mesh_uniformMesh(MPI_COMM_WORLD, &mesh, mesh_spacing);  
  
  unsigned int tc = 0; 

  tree_initTree(MPI_COMM_WORLD,  &mesh, &tree);     
  tree_printTree_onFile(MPI_COMM_WORLD, &tree, tc);  
  print_structuredAMR(MPI_COMM_WORLD, &tree, &nodes, 0.5*mesh_spacing[RFN_LEVELS-1], tc);  
  
  
  
  /*
   * Refinment test
   */     
  const int rfn_type = RFN_TYPE;
      
  if(rfn_type == 1) {

    tc = 1; 
    
    double *startPoint = malloc(3*sizeof(double));
    startPoint[0] = 0.375; 
    startPoint[1] = 0.375; 
    startPoint[2] = 0.375;

    double *endPoint = malloc(3*sizeof(double));
    endPoint[0] = 0.625; 
    endPoint[1] = 0.625; 
    endPoint[2] = 0.625;
    
    refinment_setBox(MPI_COMM_WORLD, &tree, startPoint, endPoint);

    tree_printTree_onFile(MPI_COMM_WORLD, &tree, tc);
	
    print_structuredAMR(MPI_COMM_WORLD, &tree, &nodes, 0.5*mesh_spacing[RFN_LEVELS-1], tc);  


    tc = 2;
    
    startPoint[0] = 0.375; 
    startPoint[1] = 0.375; 
    startPoint[2] = 0.375;

    endPoint[0] = 0.5; 
    endPoint[1] = 0.625; 
    endPoint[2] = 0.625;
    
    refinment_setBox(MPI_COMM_WORLD,  &tree, startPoint, endPoint);
    
    tree_printTree_onFile(MPI_COMM_WORLD, &tree, tc);
    
    print_structuredAMR(MPI_COMM_WORLD, &tree, &nodes, 0.5*mesh_spacing[RFN_LEVELS-1], tc);  
    
    
    tc = 3;
    
    refinment_verifyAspectRatio(MPI_COMM_WORLD, &tree, mesh_spacing);
    
    tree_printTree_onFile(MPI_COMM_WORLD, &tree, tc);
    
    print_structuredAMR(MPI_COMM_WORLD, &tree, &nodes, 0.5*mesh_spacing[RFN_LEVELS-1], tc);  

    
  }
  else if(rfn_type == 2) {
    
    const int finalTimeCounter = NX;
    for(int tc=1; tc<=finalTimeCounter; tc++) {      
      refinment_zPlaneSlope(MPI_COMM_WORLD, &tree, tc);
      tree_printTree_onFile(MPI_COMM_WORLD, &tree, tc); 
      print_structuredAMR(MPI_COMM_WORLD, &tree, &nodes, 0.5*mesh_spacing[RFN_LEVELS-1], tc);  
    }   
    
  }
  else if(rfn_type == 3) {    
    
    char type = 'x';
    tc = 1;    
    refinment_stretching(MPI_COMM_WORLD, &tree, type);
    tree_printTree_onFile(MPI_COMM_WORLD, &tree, tc); 
    print_structuredAMR(MPI_COMM_WORLD, &tree, &nodes, 0.5*mesh_spacing[RFN_LEVELS-1], tc);  
    
//      tc = 2;
//      type = 'y';
//      refinment_stretching(MPI_COMM_WORLD, &tree, type);
//      tree_printTree_onFile(MPI_COMM_WORLD, &tree, tc); 
    
  }

    
  /*
   * Search host leaf of Lagrangian
   */
  
  if(FIND_LAG == 1) {
    const double xl = 0.26, yl = 0.26, zl = 0.126;
    leaf_t *host_leaf = tree_findLeaf(MPI_COMM_WORLD, &tree, xl, yl, zl, mesh_spacing);  
    tree_printLeaf_onFile(MPI_COMM_WORLD, &host_leaf, xl, yl, zl);
    tree_printLeaf_onScreen(MPI_COMM_WORLD, &host_leaf, xl, yl, zl);  
  }
    
  tree_clearTree(MPI_COMM_WORLD, &tree);
    

  
  MPI_Finalize();  
  
  return 0;

}
