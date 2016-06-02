#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include "machine.h" 

/* 
 * ===	FUNCTION  ======================================================================
 *	       Name:  pprintf
 *	  Arguments:  char * format - Format of the print function.
 *	              ... - Unknown further amount of arguments.
 *	Description:  Print function for parallel code. Works like printf such that it takes
 *	              a format string and prints the corresponding arguments on the right
 *	              for each "%".
 * =====================================================================================
 */

void pprintf(char *format, ...) {
	va_list l;
	char* new_format; 
	va_start(l,format);

	// RANK-ZERO ONLY
	if (format[0]=='%' && format[1]=='Z')  {
	  if (host.rank==0) {
	    format += 2; 
	    vprintf(format,l);
	  }
	}
	else  {
	// GLOBAL
	// problem with buffering -- make format and send to vprintf
	  new_format = malloc(sizeof(char)*strlen(format)+6);
	  sprintf(new_format,"P%d:%s",host.rank,format);
	  vprintf(new_format,l);
	  free(new_format); 
	}

	va_end(l);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  init_machine
 *    Arguments:  int argc - Number of arguments in main.
 *                char * argv[] - Argument list of main.
 *  Description:  Initialises MPI Process and host struct.
 * =====================================================================================
 */

void init_machine(int argc, char *argv[]) {
	if (argc != 3) {
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	MPI_Comm_size(MPI_COMM_WORLD,&host.np);
	MPI_Comm_rank(MPI_COMM_WORLD,&host.rank);

	// Implement specific geometry. //
    //  5 6 7 
	//  4 . .
	//  2 3 .
	//  0 1 .

	// Work out neighbouring process ranks
	if(host.rank == 0) {
		host.neighbour[WEST] = -2 ;
		host.neighbour[EAST] = 1 ;
		host.neighbour[NORTH] = 2 ;
		host.neighbour[SOUTH] = -1 ;
	} else if (host.rank == 1) {
		host.neighbour[WEST] = 0 ;
		host.neighbour[EAST] = -2 ;
		host.neighbour[NORTH] = 3 ;
		host.neighbour[SOUTH] = -1 ;
	} else if (host.rank == 2) {
		host.neighbour[WEST] = -2 ;
		host.neighbour[EAST] = 3 ;
		host.neighbour[NORTH] = 4 ;
		host.neighbour[SOUTH] = 0 ;	
	} else if (host.rank == 3) {
		host.neighbour[WEST] = 2 ;
		host.neighbour[EAST] = -2 ;
		host.neighbour[NORTH] = -2 ;
		host.neighbour[SOUTH] = 1 ;
	} else if (host.rank == 4) {
		host.neighbour[WEST] = -2 ;
		host.neighbour[EAST] = -2 ;
		host.neighbour[NORTH] = 5 ;
		host.neighbour[SOUTH] = 2 ;	
	} else if (host.rank == 5) {
		host.neighbour[WEST] = -2 ;
		host.neighbour[EAST] = 6 ;
		host.neighbour[NORTH] = -2 ;
		host.neighbour[SOUTH] = 4 ;
	} else if (host.rank == 6) {
		host.neighbour[WEST] = 5 ;
		host.neighbour[EAST] = 7 ;
		host.neighbour[NORTH] = -2 ;
		host.neighbour[SOUTH] = -2 ;
	} else if (host.rank == 7) {
		host.neighbour[WEST] = 6 ;
		host.neighbour[EAST] = -1 ;
		host.neighbour[NORTH] = -2 ;
		host.neighbour[SOUTH] = -2 ;
	}

	pprintf( "Grid rank is %d , neighbours: N=%d,S=%d,E=%d,W=%d\n",
	   host.rank,
	   host.neighbour[NORTH], host.neighbour[SOUTH],
	   host.neighbour[EAST], host.neighbour[WEST]
	) ;
}
