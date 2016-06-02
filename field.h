#ifndef FIELD_H_NQX3YHOF
#define FIELD_H_NQX3YHOF

/*
 * =====================================================================================
 *
 *       Filename:  field.h
 *
 *    Description:  Header file for field structure.
 *
 *        Version:  1.0
 *        Created:  05/03/16 10:59:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#define BOUNDARY_DIR 0
#define BOUNDARY_VN 1
#define BOUNDARY_ACTIVE 2

/* 
 * ===  STRUCT  ========================================================================
 *         Name:  field_t
 *       Fields:  double * active_data - 1-D Array storing data of points actively being
 *                changed during the simulation.
 *                double ** value - Pointers used to access the active data.
 *                int boundary_type - Array storing whether or not the boundary is 
 *                active (internal) or inactive(external).
 *  Description:  The global field contains all the points being evaluated during the
 *                simulation. Each process contains a different part of the field as
 *                specified by the Grid object.
 * =====================================================================================
 */

typedef struct Field {
  double *active_data;  // Storage for the actively managed data-points
  double **value;	// Access pointers
  // double *boundary[4];
  int boundary_type[4];
} Field ;

Field* init_field(Grid* g, double (*init)(int, int, Grid*));
void free_field(Field* f);
void write_field(char *, Field*, Grid*);

#endif /* end of include guard: FIELD_H_NQX3YHOF */
