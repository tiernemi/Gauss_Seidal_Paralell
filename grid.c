/*
 * =====================================================================================
 *
 *       Filename:  grid.c
 *
 *    Description:  Source code for setting up Grid object.
 *
 *        Version:  1.0
 *        Created:  05/03/16 11:03:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include <stdlib.h> 
#include "machine.h" 
#include "grid.h" 

// -------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  init_grid
 *    Arguments:  char * argv[] - Command Line Arguments passed to program.
 *      Returns:  A pointer to a new Grid object.
 *  Description:  Creates a grid object and checks if it's possible to divide up work
 *                between processes.
 * =====================================================================================
 */

Grid* init_grid(char *argv[]) {
	Grid* g = malloc(sizeof(Grid)); 
	g->nx_local = atoi(argv[1]);
	g->ny_local = atoi(argv[2]);

	// Implement specific geometry. //
    //  5 6 7 
	//  4 . .
	//  2 3 .
	//  0 1 .
	
	if (host.rank < 5) {
		g->x_offset = host.rank%2 * g->nx_local;
		g->y_offset = host.rank/2 * g->ny_local;
	} else if (host.rank == 5) {
		g->x_offset =  0*g->nx_local ;
		g->y_offset = 3*g->ny_local;;
	} else if (host.rank == 6) {
		g->x_offset =  1*g->nx_local ;
		g->y_offset = 3*g->ny_local;;
	} else if (host.rank == 7) {
		g->x_offset =  2*g->nx_local  ;
		g->y_offset = 3*g->ny_local;;
	}

	return g;
}

void free_grid(Grid *g) {
	free(g) ;
}
// =======================================

