#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi/mpi.h"

#include "machine.h" 
// Host Machine contains information about the current session. //
Machine host ; 

#include "grid.h" 
#include "field.h" 
#include "comms.h" 

// =======================================

/* 
 * ===	FUNCTION  ======================================================================
 *	       Name:  boundary
 *	  Arguments:  int boundaryID - The number associated with the location of the boundary
 *	    Returns:  Dirichlet data for particular problem.
 *	Description:  Function handling diriclet boundaries.
 * =====================================================================================
 */

// Return Dirichlet data for boundary b (N,S,E,W) at location k
double boundary(int boundaryID, int k, Grid *g) {
	if (boundaryID==NORTH) {
		return 0.0 ;
	}
	if (boundaryID==SOUTH) {
		return 0 ;
	}
	if (boundaryID==EAST) {
		return 1.0 ;
	}
	if (boundaryID==WEST) {
		return 0.0;
	}
	return 0 ;
}

// ------------------------------
double update(Field* phi, int eo, Grid* g) {
	int x,i,j,y,p=eo;
	double diff,diff2_local=0.0,diff2_global;
	double omega=1.9; // Over-relax

	// Update using GS iteration

	for (x=1;x<g->nx_local-1;x++) {
		for (y=p+1;y<g->ny_local-1;y+=2) {

			diff = phi->value[x][y]; 
			phi->value[x][y] = (1.0-omega)*phi->value[x][y] + omega * 0.25 * (
				phi->value[x+1][y] + phi->value[x-1][y]
				+ phi->value[x][y+1] + phi->value[x][y-1]);

			diff -= phi->value[x][y];
			diff2_local += diff*diff;
		}
		p=1-p; 
	}

	int even_odd_state = eo ;
	int ny_local = g->ny_local ;
	int nx_local = g->nx_local ;

	p = 1 - even_odd_state ;
	if (phi->boundary_type[WEST] == BOUNDARY_VN) {
		for (j = p+1 ; j < ny_local-1 ; j+=2) {
			diff = phi->value[0][j] ;
			phi->value[0][j] = 0.125 * (3*phi->value[0][j+1] + 3*phi->value[0][j-1] + 
	    	2*phi->value[1][j]) ;	
			diff -= phi->value[0][j] ;
			diff2_local += diff*diff;
		}
	} else {
		for (j = p+1 ; j < ny_local-1 ; j+=2) {
			diff = phi->value[0][j] ;
			phi->value[0][j] = 0.25 * (phi->value[1][j] + phi->value[-1][j] + 
	    	phi->value[0][j+1] + phi->value[0][j-1]) ;
			diff -= phi->value[0][j];
			diff2_local += diff*diff;
		}
	}
	// Possible VN boundaries. //
	if (phi->boundary_type[EAST] == BOUNDARY_VN) {
		for (j = p+1 ; j < ny_local-1 ; j+=2) {
			diff = phi->value[nx_local-1][j];
			phi->value[nx_local-1][j] = 0.125 * (3*phi->value[nx_local-1][j+1] + 3*phi->value[nx_local-1][j-1] + 
	    	2*phi->value[nx_local-2][j] ) ;			
			diff -= phi->value[nx_local-1][j];
			diff2_local += diff*diff;
		}
	} else {
		for (j = p+1 ; j < ny_local-1 ; j+=2) {
			diff = phi->value[nx_local-1][j];
			phi->value[nx_local-1][j] = 0.25 * (phi->value[nx_local][j] + phi->value[nx_local-2][j] + 
	    	phi->value[nx_local-1][j+1] + phi->value[nx_local-1][j-1]) ;
			diff -= phi->value[nx_local-1][j];
			diff2_local += diff*diff;
		}
	}
	if (phi->boundary_type[NORTH] == BOUNDARY_VN) {
		for (i = p+1 ; i < nx_local-1 ; i+=2) {
			diff = phi->value[i][ny_local-1];
			phi->value[i][ny_local-1] = 0.125 * (3*phi->value[i+1][ny_local-1] + 3*phi->value[i-1][ny_local-1] + 
	    	2*phi->value[i][ny_local-2] ) ;			
			diff -= phi->value[i][ny_local-1];
			diff2_local += diff*diff;
		}
	} else {
		for (i = p+1 ; i < nx_local-1 ; i+=2) {
			diff = phi->value[i][ny_local-1];
			phi->value[i][ny_local-1] = 0.25 * (phi->value[i-1][ny_local-1] + phi->value[i+1][ny_local-1] + 
	    	phi->value[i][ny_local] + phi->value[i][ny_local-2]) ;
			diff -= phi->value[i][ny_local-1];
			diff2_local += diff*diff;
		}
	}
	if (phi->boundary_type[SOUTH] == BOUNDARY_VN) {
		for (i = p+1 ; i < nx_local-1 ; i+=2) {
			diff = phi->value[i][0];
			phi->value[i][0] = 0.125 * (3*phi->value[i+1][0] + 3*phi->value[i-1][0] + 
	    	2*phi->value[i][1]) ;			
			diff -= phi->value[i][0];
			diff2_local += diff*diff;
		}
	} else {
		for (i = p+1 ; i < nx_local-1 ; i+=2) {
			diff = phi->value[i][0];
			phi->value[i][0] = 0.25 * (phi->value[i-1][0] + phi->value[i+1][0] + 
	    	phi->value[i][-1] + phi->value[i][1]) ;
			diff -= phi->value[i][0];
			diff2_local += diff*diff;
		}
	}

	// SOUTH WEST
	if (phi->boundary_type[WEST] == BOUNDARY_VN && phi->boundary_type[SOUTH] == BOUNDARY_VN) {
		phi->value[0][0] = 0.5 * (phi->value[0][1] + phi->value[1][0]) ;
	} else if (phi->boundary_type[WEST] == BOUNDARY_VN && phi->boundary_type[SOUTH] != BOUNDARY_VN) {
		phi->value[0][0] = 0.125 * (3*phi->value[0][1] + 3*phi->value[0][-1] + 
		2*phi->value[1][0]) ;	
	} else if (phi->boundary_type[SOUTH] == BOUNDARY_VN && phi->boundary_type[WEST] != BOUNDARY_VN) {
		phi->value[0][0] = 0.125 * (3*phi->value[1][0] + 3*phi->value[-1][0] + 
		2*phi->value[0][1]) ;
	} else {
		phi->value[0][0] = 0.25 * (phi->value[-1][0] + phi->value[1][0] + 
		phi->value[0][-1] + phi->value[0][1]) ;
	}
	// SOUTH EAST
	if (phi->boundary_type[EAST] == BOUNDARY_VN && phi->boundary_type[SOUTH] == BOUNDARY_VN) {
		phi->value[nx_local-1][0] = 0.5 * (phi->value[nx_local-1][1] + phi->value[nx_local-2][0]) ;
	} else if (phi->boundary_type[EAST] == BOUNDARY_VN && phi->boundary_type[SOUTH] != BOUNDARY_VN) {
		phi->value[nx_local-1][0] = 0.125 * (3*phi->value[nx_local-1][1] + 3*phi->value[nx_local-1][-1] + 
		2*phi->value[nx_local-2][0]) ;	
	} else if (phi->boundary_type[SOUTH] == BOUNDARY_VN && phi->boundary_type[EAST] != BOUNDARY_VN) {
		phi->value[nx_local-1][0] = 0.125 * (3*phi->value[nx_local][0] + 3*phi->value[nx_local-2][0] + 
		2*phi->value[nx_local-1][1]) ;
	} else {
		phi->value[nx_local-1][0] = 0.25 * (phi->value[nx_local-2][0] + phi->value[nx_local][0] + 
		phi->value[nx_local-1][-1] + phi->value[nx_local-1][1]) ;
	}
	// NORTH WEST
	if (phi->boundary_type[WEST] == BOUNDARY_VN && phi->boundary_type[NORTH] == BOUNDARY_VN) {
		phi->value[0][ny_local-1] = 0.5 * (phi->value[0][ny_local-2] + phi->value[1][ny_local-1]) ;
	} else if (phi->boundary_type[WEST] == BOUNDARY_VN && phi->boundary_type[NORTH] != BOUNDARY_VN) {
		phi->value[0][ny_local-1] = 0.125 * (3*phi->value[0][ny_local] + 3*phi->value[0][ny_local-2] + 
		2*phi->value[1][ny_local-1]) ;	
	} else if (phi->boundary_type[NORTH] == BOUNDARY_VN && phi->boundary_type[WEST] != BOUNDARY_VN) {
		phi->value[0][ny_local-1] = 0.125 * (3*phi->value[1][ny_local-1] + 3*phi->value[-1][ny_local-1] + 
		2*phi->value[0][ny_local-2]) ;
	} else {
		phi->value[0][ny_local-1] = 0.25 * (phi->value[-1][ny_local-1] + phi->value[1][ny_local-1] + 
		phi->value[0][ny_local-2] + phi->value[0][ny_local]) ;
	}
	// NORTH EAST
	if (phi->boundary_type[EAST] == BOUNDARY_VN && phi->boundary_type[NORTH] == BOUNDARY_VN) {
		phi->value[nx_local-1][ny_local-1] = 0.5 * (phi->value[nx_local-1][ny_local-2] + phi->value[nx_local-2][ny_local-1]) ;
	} else if (phi->boundary_type[EAST] == BOUNDARY_VN && phi->boundary_type[NORTH] != BOUNDARY_VN) {
		phi->value[nx_local-1][ny_local-1] = 0.125 * (3*phi->value[nx_local-1][ny_local] + 3*phi->value[nx_local-1][ny_local-2] + 
		2*phi->value[nx_local-2][ny_local-1]) ;	
	} else if (phi->boundary_type[NORTH] == BOUNDARY_VN && phi->boundary_type[EAST] != BOUNDARY_VN) {
		phi->value[nx_local-1][ny_local-1] = 0.125 * (3*phi->value[nx_local][ny_local-1] + 3*phi->value[nx_local-2][ny_local-1] + 
		2*phi->value[nx_local-1][ny_local-2]) ;
	} else {
		phi->value[nx_local-1][ny_local-1] = 0.25 * (phi->value[nx_local-2][ny_local-1] + phi->value[nx_local][ny_local-1] + 
		phi->value[nx_local-1][ny_local-2] + phi->value[nx_local-1][ny_local]) ;
	}


	MPI_Allreduce(&diff2_local,&diff2_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	// Communicate boundaries in new
	send_boundary_data(phi,eo,g); 

	return diff2_global;
}

int main(int argc, char *argv[]) {
	Grid* g;
	Field* phi; 
	double diff2;

	MPI_Init(&argc,&argv);
	init_machine(argc,argv); 

	g = init_grid(argv); 

	// Only one field for GS algorithm
	phi = init_field(g,boundary); 

	int counter = 0 ;
	// iterate...
	do {
		diff2 = update(phi,0,g);
		diff2 = update(phi,1,g);
		++counter ;
	} while (diff2 > 1.0e-7) ;


	write_field("grid.out",phi,g);

	if (host.rank == 4) {
		printf("%lf\n", phi->value[g->nx_local/2][g->ny_local/2]);
	}

	free_grid(g);
	free_field(phi);
	MPI_Finalize();

	return EXIT_SUCCESS ;
}
