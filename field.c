/*
 * =====================================================================================
 *
 *       Filename:  field.c
 *
 *    Description:  Source code for creating field object.
 *
 *        Version:  1.0
 *        Created:  05/03/16 12:19:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 

#include "machine.h" 
#include "grid.h" 
#include "field.h" 

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  init_field
 *    Arguments:  Grid * g - The grid used to divide up the process.
 *                double (*init)(int, int, Grid*) - Function pointer to boundary info.
 *      Returns:  A pointer to a new Field object.
 *  Description:  Initialises the field with boundary information specified by boundary
 *                function pointer with work split among the cores with a geometry
 *                specified by Grid * g.
 * =====================================================================================
 */

Field* init_field(Grid* g, double (*init)(int, int, Grid*) ) {
	int size;
	int i,x,y;
	Field *f = malloc(sizeof(Field));

	// Malloc space for active data //
	size = (g->nx_local+2) * (g->ny_local+2);
	f->active_data = malloc(sizeof(double) * size); 
	for (i=0;i<size;i++) {
		f->active_data[i] = 0.0;
	}

	// Assign pointers to access data //
	f->value = malloc(sizeof(double *) * (g->nx_local+2) );
	for (x=0;x<g->nx_local+2;x++) {
		f->value[x] = f->active_data + x * (g->ny_local+2) + 1;
	}
	f->value++; // shift to start array at -1

    // Create geometry of particular problem. //


	// Now decide if the four boundarys are active (need data updates from other
	// processes) or passive (are actual edges of the grid, so stay unchanged)
	for (i=0;i<4;i++) {
		if (host.neighbour[i]==-1) {
			f->boundary_type[i] = BOUNDARY_DIR ;
		} else if (host.neighbour[i]==-2) {
			f->boundary_type[i] = BOUNDARY_VN ;
		} else {
			f->boundary_type[i] = BOUNDARY_ACTIVE ;
		}
	}

	pprintf("boundary types: N=%d,S=%d,E=%d,W=%d\n",
		 f->boundary_type[NORTH], 
		 f->boundary_type[SOUTH],
		 f->boundary_type[EAST],
		 f->boundary_type[WEST]);

	// Load boundary data from initialiser function
	if (f->boundary_type[NORTH] == BOUNDARY_DIR)
		for (x=0;x<g->nx_local;x++)
			f->value[x][g->ny_local] = init(NORTH,x+g->x_offset,g);
	else 
		for (x=0;x<g->nx_local;x++)
			f->value[x][g->ny_local] = 0.0;

	if (f->boundary_type[SOUTH] == BOUNDARY_DIR)
		for (x=0;x<g->nx_local;x++)
			f->value[x][-1] = init(SOUTH,x+g->x_offset,g);
	else 
		for (x=0;x<g->nx_local;x++)
			f->value[x][-1] = 0.0;

	if (f->boundary_type[EAST ] == BOUNDARY_DIR)
		for (y=0;y<g->ny_local;y++)
			f->value[g->nx_local][y] = init(EAST ,y+g->y_offset,g);
	else 
		for (y=0;y<g->ny_local;y++)
			f->value[g->nx_local][y] = 0.0;

	if (f->boundary_type[WEST ] == BOUNDARY_DIR)
		for (y=0;y<g->ny_local;y++)
			f->value[-1][y] = init(WEST ,y+g->y_offset,g);
	else 
		for (y=0;y<g->ny_local;y++)
			f->value[-1][y] = 0.0;

	return f;
}

/* 
 * ===	FUNCTION ======================================================================
 *         Name: field
 *    Arguments: Field * f - Field object to free.
 *  Description: Frees heap memory associated with field.
 * =====================================================================================
 */

void free_field(Field* f) {
	free(f->active_data);
	// Undo shift of f->value pointer before freeing 
	f->value--;
	free(f->value);
	free(f);
}

/* 
 * ===	FUNCTION ======================================================================
 *	       Name: write_field
 *    Arguments: char * filestub - Name of file.
 *               Field * f - Field to be written.
 *               Grid * g - Grid describing how the field is divided between processes.
 *  Description: Outputs field to file.
 * =====================================================================================
 */

void write_field(char* filestub, Field *f, Grid *g) {
	int x,y;
	FILE* o;
	char *filename = malloc(sizeof(char) * (strlen(filestub)+10));

	sprintf(filename,"%s.%.3d",filestub,host.rank);
	o = fopen(filename,"w");

  // Write the active part, with offsets. Don't write process identities
	for (x=0;x<g->nx_local;x++)
		for (y=0;y<g->ny_local;y++)
			fprintf(o,"%d %d %lf\n",x+g->x_offset,y+g->y_offset,f->value[x][y]);

	fclose(o);
	free(filename);
}
