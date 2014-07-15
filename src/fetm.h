/******************************************************************************
  Module Name : fetm.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description : Finite Element triangle Mesh (FETM) method.

  Others :
      Refers to HYWang's FETM geometry algorithm.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef FETM_H
#define FETM_H

/*------------------------------Defines------------------------------*/

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bulk.h"
#include "target.h"
#include "utils.h"
#include "random.h"
#include "init.h"
#include "material.h"
#include "transport.h"
#include "triangulate.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/* or ------------------------------------*/
//#include "random.h"

/*--------------------------Global variables-------------------------*/
float  box_start_x0, box_start_y0, box_start_z0;  /* fetm box start point */
int    box_count_x, box_count_y, box_count_z;  /* number of box element for a geomentry */
double box_size_x, box_size_y, box_size_z;     /* size of box element for a geometry */
float  box_end_x0, box_end_y0, box_end_z0;
float  box_length_x, box_length_y, box_length_z;
float  z0_max_fetm;
double s0_fetm;  /* substrate section */

/* for makeline */
int    line_flag;
int    di, dj, dk;  /*  */
double ax, ay, az;  /*  */
double rx, ry, rz;  /* relative x, y, z */

struct GRIDBOX  /* geometry */
{
    int   num_triangles;    /* number of triangles in this box */
    int   flag;             /* used for triangulating */
    int   mater_type[800];  /* material type inside this triangle, 100 -> 800 */
    float GBOX[800][9];     /* triangle vertex coordinates */
};
struct GRIDBOX ***GRB;

/*-----------------------------Functions-----------------------------*/
/* read fetm shape parameters form Config file */
int read_fetm_shape (char *file_name);

/* free flight path of ion for fetm geometry method */
void freepath_fetm (double *x, double *y, double *z, double *vx, double *vy, double *vz,
                    int projZ, int projM, float *energy, double flight_length,
                    int *current_material_index, int *ion_left_target, int *in_or_not,
                    float *stopping, float *straggling);

int check_target_type_fetm (double x, double y, double z, double vx, double vy, double vz,
                            double flight_length);

/* triangle interaction for fetm method */
int triangle_intersect_fetm (double trix1, double triy1, double triz1, double trix2,
	                           double triy2, double triz2, double trix3, double triy3,
	                           double triz3, double *temp, double N[], int *flag);

/* step of boxs */
void make_line (double x, double y, double z, double vx, double vy, double vz,
                int *i, int *j, int *k);

/* check up if it still has the point of intersection with geomrtry when ions outside of the target */
int juintersect (double *x, double *y, double *z, double vx, double vy, double vz);

/* determine material types */
void check_material_type (int *in_or_not, double x, double y, double z, int projZ,
                          int *current_material_index, float *current_surf_bind,
                          float *s1, int type);

/* check surface binding energy of ion/atoms */
float check_surf_bind (int current_material_index, int projZ);

/* surface/Interface refraction */
void refraction (double NT[], double *vx, double *vy, double *vz, int *S_flag,
                 float energy, float old_surf_bind, float current_surf_bind);

void transmission (double uz, int idied, float energy, float U0);

#endif
