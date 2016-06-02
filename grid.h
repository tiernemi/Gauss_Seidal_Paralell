#ifndef GRID_H_HAOVMXGR
#define GRID_H_HAOVMXGR

/*
 * =====================================================================================
 *
 *       Filename:  grid.h
 *
 *    Description:  Header file for Grid object. The grid is the structure of the PDE.
 *
 *        Version:  1.0
 *        Created:  05/03/16 11:00:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

/* 
 * ===  STRUCT  ========================================================================
 *         Name:  grid_t
 *       Fields:  nx - Number of x divisions in grid of dimension nx X ny.
 *                ny - Number of y divisions in grid of dimension nx X ny.
 *                nx_local - Number of rows of grid stored on rank.
 *                ny_local - Number of columns of grid stored on rank.
 *                x_offset - Global x offset of rank's grid section from the
 *                           full grid.
 *                y_offset - Global y offset of rank's grid section from the
 *                           full grid.
 *  Description:  Grid object stores data about the global grid which the program is 
 *                evaluating. The "Field" of the PDE is split into a "Grid" in which 
 *                each grid section is evaluated on a different process. 
 * =====================================================================================
 */

typedef struct Grid {
  int nx_local ;
  int ny_local ; 
  int x_offset ;
  int y_offset ;
} Grid ; 

// -------------------------------------
Grid* init_grid(char *argv[]); 
void free_grid(Grid *g); 

#endif /* end of include guard: GRID_H_HAOVMXGR */
