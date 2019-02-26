#include <stdio.h>    /* fprintf() */
#include <stdlib.h>   /* exit() */
#include <math.h>
#include <mpi.h>

#include "constants.h"

/*
 * Monitoring CONSTANTS. See constants.h
 *
 * comm - MPI communicator
 */
void monitor_constants(MPI_Comm comm) {
  
  int process_rank;
  MPI_Comm_rank(comm, &process_rank);
  
  /*
   * Print on file, per process, to monitor the code algorithm
   */   
  printf("\n\nconstants.h\n");
  printf("\n\n\tGeneral parameters\n");
  printf("\n\tMESH_LEVELS %i\n", MESH_LEVELS);
  printf("\n\tRFN_LEVELS %i\n", RFN_LEVELS);
  printf("\n\tRFN_TYPE %i\n", RFN_TYPE);
  printf("\n\tSIZE_KEY %i\n", SIZE_KEY);
  printf("\n\tFIND_LAG %i\n", FIND_LAG);
  printf("\n\tID_UPPER_LIMIT %u\n", ID_UPPER_LIMIT);
  
  printf("\n\tDomain configurations\n");
  printf("\n\tLX %f\n", LX);
  printf("\n\tLY %f\n", LY);
  printf("\n\tLZ %f\n", LZ);
  printf("\n\tNX %i\n", NX);
  printf("\n\tNY %i\n", NY);
  printf("\n\tNZ %i\n", NZ);
  printf("\n\tORIGIN_X %f\n", ORIGIN_X);
  printf("\n\tORIGIN_Y %f\n", ORIGIN_Y);
  printf("\n\tORIGIN_Z %f\n", ORIGIN_Z);
  
  printf("\n\tMultiprocessing parameters\n");
  printf("\n\tPROCX %i\n", PROCX);
  printf("\n\tPROCY %i\n", PROCY);
  printf("\n\tPROCZ %i\n", PROCZ);
  
  printf("\n\tFinite Volume parameters\n");
  printf("\n\tNUMBER_CENTER_VARS %i\n", NUMBER_CENTER_VARS);
  printf("\n\tNUMBER_FACE_VARS %i\n", NUMBER_FACE_VARS);
  printf("\n\tNUMBER_FLUXES %i\n", NUMBER_FLUXES);
  printf("\n\tNUMBER_RATES %i\n", NUMBER_RATES);

  printf("\n\tFlow data\n");
  printf("\n\tRHO %f\n", RHO);
  
  printf("\n\tMath constants\n");
  printf("\n\tPI %f\n", PI);
  
}