#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <inttypes.h>
#include <limits.h>

/*
 * MESH_LEVELS - Create the first uniform, one level mesh. See mesh.h
 */
#define MESH_LEVELS 1 

/*
 * RNF_LEVELS - Maximum number of levels of the mesh. See refinment.h
 */
#define RFN_LEVELS 3

/*
 * RFN_TYPE - Type of refinment. See refinment.h
 * 
 * 1: setBox. Local refinment - Currently using RFN_LEVELS 3
 * 2: zPlaneSlope. A function of the time counter - Currently using RFN_LEVELS 2
 * 3: stretching. Stretching according to direction x, y or z - Currently using RFN_LEVELS 4  
 */
#define RFN_TYPE 1

/*
 * Size key of leaf hash table. See tree.h
 */
#define SIZE_KEY 3*RFN_LEVELS

/*
 * FIND_LAG - Want to print host leaf data of a Lagrangian element? See main.c
 */
#define FIND_LAG 0

/*
 * Maximum admitted number of Lagrangian elements per process. See lagrangian.c
 * TODO Guarantee unsigned value
 */
#define ID_UPPER_LIMIT UINT_MAX
//UINTMAX_MAX-1

/*
 * Domain configurations 
 * Must be DX=DY=DZ
 * LX, LY, LZ - domain lenght in each direction
 * NX, NY, NZ - first level grid in each direction
 * ORIGIN_X, _Y, _Z - origin coordinate
 */
#define LX 1.0
#define LY 1.0
#define LZ 1.0
#define NX 8
#define NY 8
#define NZ 8
#define ORIGIN_X 0.0
#define ORIGIN_Y 0.0
#define ORIGIN_Z 0.0

/*
 * Multiprocessing parameters
 * 
 * PROCX, PROCY, PROCZ - number of processes in each direction
 */
#define PROCX 1
#define PROCY 1
#define PROCZ 1

/*
 * Finite Volume definitions. Decide latter for velocity at face or center
 * 
 * NUMBER_CENTER_VARS: Pressure, Temperature, (U, V, W or...)
 * NUMBER_FACE_VARS: U, V, W
 * NUMBER_FLUXES: Per face - One velocity component convective and diffusive, Temperature convective and diffusive
 * NUMBER_RATES: Per volume - Velocity convective and diffusive, Temperature convective and diffusive 
 */
#define NUMBER_CENTER_VARS 2 //or 5 to collocated arrangment
#define NUMBER_FACE_VARS 3 //or 0 to collocated arrangment
#define NUMBER_FLUXES 4*NUMBER_FACE_VARS
#define NUMBER_RATES 4

/*
 * Flow data
 * 
 * RHO - fluid density
 */
#define RHO 1.0

/*
 * Math constants
 */
#define PI 4*atan(1.0)


void monitor_constants(MPI_Comm comm);

#endif
