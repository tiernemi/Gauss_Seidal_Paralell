#ifndef MACHINE_H_WH6E1RZG
#define MACHINE_H_WH6E1RZG

/*
 * =====================================================================================
 *
 *       Filename:  machine.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/03/16 10:31:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *   Organization:  
 *
 * =====================================================================================
 */

#include "mpi/mpi.h"

#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

/* 
 * ===  STRUCT  ========================================================================
 *         Name:  Machine
 *       Fields:  int np - Number of processes.
 *                int rank - Rank of process.
 *                int p_nx - Number of rows
 *                int p_ny -
 *                int p_x -
 *                int p_y -
 *  Description:  Structure that stores useful information about the current MPI session.
 * =====================================================================================
 */

typedef struct Machine  {
  int np;
  int rank;
  int neighbour[4];
} Machine ; 

// Host machine visible to all processes. //
extern Machine host; 

void init_machine(int argc, char *argv[]); 
void pprintf(char *format, ...);

#endif /* end of include guard: MACHINE_H_WH6E1RZG */
