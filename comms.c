/*
 * =====================================================================================
 *
 *	     Filename:  comms.c
 *
 *	  Description:  Source file containing communication protocol for carrying out
 *	                neighbour communications.
 *
 *	      Version:  1.0
 *	      Created:  05/03/16 10:47:06
 *	     Revision:  none
 *	     Compiler:  gcc
 *
 *	       Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include "mpi/mpi.h"
#include "machine.h" 
#include "grid.h" 
#include "field.h" 
#include "comms.h" 
#include "stdio.h"

/* 
 * ===	STRUCT  ========================================================================
 *	       Name: 
 *	     Fields:  
 *	Description:  
 * =====================================================================================
 */

typedef struct BoundaryComm {
	MPI_Request send[4];
	MPI_Status  send_status[4];
	int nsend;
	MPI_Request recv[4]; 
	MPI_Status  recv_status[4];
	int nrecv;

	double *buffer_send_east;
	double *buffer_recv_east;

	double *buffer_send_west;
	double *buffer_recv_west;

	double *buffer_send_north;
	double *buffer_recv_north;

	double *buffer_send_south;
	double *buffer_recv_south;
} BoundaryComm ; 


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  comms
 *    Arguments:  
 *      Returns:  
 *  Description:  
 * =====================================================================================
 */

static void send_east(Field* f, Grid* g, int eo, BoundaryComm* c) {
	int y ;

	if (f->boundary_type[EAST] == BOUNDARY_ACTIVE) {
		for (y=1-eo;y<g->ny_local;y+=2) {
			c->buffer_send_east[y/2] = f->value[g->nx_local-1][y];
		}
		MPI_Isend(c->buffer_send_east,g->ny_local/2,MPI_DOUBLE,
			host.neighbour[EAST],1000*host.rank,MPI_COMM_WORLD,c->send+c->nsend);
		c->nsend++;
	}
	if (f->boundary_type[WEST] == BOUNDARY_ACTIVE)  {
		MPI_Irecv(c->buffer_recv_east,g->ny_local/2,MPI_DOUBLE,
			host.neighbour[WEST],1000*host.neighbour[WEST],MPI_COMM_WORLD,c->recv+c->nrecv);
	  	c->nrecv++;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  comms
 *    Arguments:  
 *      Returns:  
 *  Description:  
 * =====================================================================================
 */

static void send_west(Field* f, Grid* g, int eo, BoundaryComm* c) {
	int y;

	if (f->boundary_type[WEST] == BOUNDARY_ACTIVE) {
		for (y=eo;y<g->ny_local;y+=2) {
	    	c->buffer_send_west[y/2] = f->value[0][y];
		}
		MPI_Isend(c->buffer_send_west,g->ny_local/2,MPI_DOUBLE,
			host.neighbour[WEST],1000*host.rank+1,MPI_COMM_WORLD,c->send+c->nsend);
		c->nsend++;
	}
	if (f->boundary_type[EAST] == BOUNDARY_ACTIVE) {
		MPI_Irecv(c->buffer_recv_west,g->ny_local/2,MPI_DOUBLE,
			host.neighbour[EAST],1000*host.neighbour[EAST]+1,MPI_COMM_WORLD,c->recv+c->nrecv);
		c->nrecv++;
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  comms
 *    Arguments:  
 *      Returns:  
 *  Description:  
 * =====================================================================================
 */

static void send_north(Field* f, Grid* g, int eo, BoundaryComm* c) {
	int x;
	if (f->boundary_type[NORTH] == BOUNDARY_ACTIVE) {
		for (x=1-eo;x<g->nx_local;x+=2) {
	    	c->buffer_send_north[x/2] = f->value[x][g->ny_local-1];
		}
		MPI_Isend(c->buffer_send_north,g->nx_local/2,MPI_DOUBLE,
			host.neighbour[NORTH],1000*host.rank+2,MPI_COMM_WORLD,c->send+c->nsend);
		c->nsend++;
	}
	if (f->boundary_type[SOUTH] == BOUNDARY_ACTIVE) {
		MPI_Irecv(c->buffer_recv_north,g->nx_local/2,MPI_DOUBLE,
			host.neighbour[SOUTH],1000*host.neighbour[SOUTH]+2,MPI_COMM_WORLD,c->recv+c->nrecv);
		c->nrecv++;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  comms
 *    Arguments:  
 *      Returns:  
 *  Description:  
 * =====================================================================================
 */

static void send_south(Field* f, Grid* g, int eo, BoundaryComm* c) {
	int x; 

	if (f->boundary_type[SOUTH] == BOUNDARY_ACTIVE) {
		for (x=  eo;x<g->nx_local;x+=2) {
			c->buffer_send_south[x/2] = f->value[x][0];
		}
		MPI_Isend(c->buffer_send_south,g->nx_local/2,MPI_DOUBLE,
			host.neighbour[SOUTH],1000*host.rank+3,MPI_COMM_WORLD,c->send+c->nsend);
		c->nsend++;
	}
	if (f->boundary_type[NORTH] == BOUNDARY_ACTIVE) {
		MPI_Irecv(c->buffer_recv_south,g->nx_local/2,MPI_DOUBLE,
			host.neighbour[NORTH],1000*host.neighbour[NORTH]+3,MPI_COMM_WORLD,c->recv+c->nrecv);
		c->nrecv++;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  comms
 *    Arguments:  
 *      Returns:  
 *  Description:  
 * =====================================================================================
 */

static void unpack(Field* f, Grid* g, int eo, BoundaryComm* c) {
	int x,y;

	if (f->boundary_type[SOUTH] == BOUNDARY_ACTIVE) {
		for (x=1-eo;x<g->nx_local;x+=2) {
	    	f->value[x][-1] = c->buffer_recv_north[x/2] ;
		}
	}

	if (f->boundary_type[NORTH] == BOUNDARY_ACTIVE) {
  		for (x=  eo;x<g->nx_local;x+=2) {
	    	f->value[x][g->ny_local] = c->buffer_recv_south[x/2]; 
		}
	}

	if (f->boundary_type[WEST] == BOUNDARY_ACTIVE) {
		for (y=1-eo;y<g->ny_local;y+=2) {
			f->value[-1][y] = c->buffer_recv_east[y/2]; 
		}
	}

	if (f->boundary_type[EAST] == BOUNDARY_ACTIVE) { 
		for (y=  eo;y<g->ny_local;y+=2) {
			f->value[g->nx_local][y] = c->buffer_recv_west[y/2]; 
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  comms
 *    Arguments:  
 *      Returns:  
 *  Description:  
 * =====================================================================================
 */

static BoundaryComm* init_comm(Grid* g, Field * f) {
	BoundaryComm *c = malloc(sizeof(BoundaryComm));

	if (f->boundary_type[NORTH] == BOUNDARY_ACTIVE) 	{
	  c->buffer_send_north = malloc(sizeof(double) * g->nx_local/2); 
	  c->buffer_recv_south = malloc(sizeof(double) * g->nx_local/2); 
	}
	if (f->boundary_type[SOUTH] == BOUNDARY_ACTIVE) {
	  c->buffer_send_south = malloc(sizeof(double) * g->nx_local/2); 
	  c->buffer_recv_north = malloc(sizeof(double) * g->nx_local/2); 
	}
	if (f->boundary_type[EAST] == BOUNDARY_ACTIVE) {
	  c->buffer_send_east  = malloc(sizeof(double) * g->ny_local/2); 
	  c->buffer_recv_west  = malloc(sizeof(double) * g->ny_local/2); 
	}
	if (f->boundary_type[WEST] == BOUNDARY_ACTIVE) {
	  c->buffer_send_west  = malloc(sizeof(double) * g->ny_local/2); 
	  c->buffer_recv_east  = malloc(sizeof(double) * g->ny_local/2); 
	}

	c->nsend = 0;
	c->nsend = 0;
	c->nrecv = 0;
	return c; 
}

static void free_comm(BoundaryComm* c, Field * f) {

	if (f->boundary_type[NORTH] == BOUNDARY_ACTIVE) {
	  free(c->buffer_send_north);
	  free(c->buffer_recv_south);
	}
	if (f->boundary_type[SOUTH] == BOUNDARY_ACTIVE)	{
	  free(c->buffer_send_south);
	  free(c->buffer_recv_north);
	}
	if (f->boundary_type[EAST] == BOUNDARY_ACTIVE) {
	  free(c->buffer_send_east);
	  free(c->buffer_recv_west);
	}
	if (f->boundary_type[WEST] == BOUNDARY_ACTIVE)	{
	  free(c->buffer_send_west);
	  free(c->buffer_recv_east);
	}
	free(c);
}

// ------------------------------------------------------------
//	      Externally visible functions:
// ------------------------------------------------------------

void send_boundary_data(Field* f, int eo, Grid* g) {
	BoundaryComm* comm = init_comm(g,f);

	send_east (f,g,eo,comm);
	send_west (f,g,eo,comm);
	send_north(f,g,eo,comm);
	send_south(f,g,eo,comm);

	MPI_Waitall(comm->nsend,comm->send,comm->send_status);
	MPI_Waitall(comm->nrecv,comm->recv,comm->recv_status);

	unpack(f,g,eo,comm);

	free_comm(comm,f);
}
