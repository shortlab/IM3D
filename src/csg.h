/******************************************************************************
  Module Name : csg.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Constructive solid geometry (CSG) method.

  Others :
      Refers to HMLi's CSG geometry algorithm.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef CSG_H
#define CSG_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "target.h"
#include "const.h"
#include "utils.h"
#include "material.h"
#include "random.h"
#include "init.h"
#include "transport.h"
#include "bulk.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/* or ------------------------------------*/
//#include <stdio.h>

//#include "random.h"

/*--------------------------Global variables-------------------------*/
float z0_max_csg;  /* the maximum depth under which the projectile can
                      be traced as in the bulk */

/* shapes information */
int    s_count;       /* number of shape elements */
int    my_shape[20];  /* indexes of shape elements */
int    is_full[20];   /* materials in shape elements */
double s0_csg;        /* substrate section */

/* direction and coordinates */
double LV[3],  LO[3];
double LV1[3], LO1[3];

/* sphere */
double sphere_C[20][3];
double sphere_R2[20];

/* tetrahedron */
double point[20][4][3];

/* cuboid */
double orientation[20][3][3];
double dc[20][3][2];

/* ellipsoid */
double abc2[20][3];

/* taper */
double taper_height[20];
double taper_TgT2[20];

/* paraboloid */
double paraboloid_height[20];
double paraboloid_radius[20];
double paraboloid_RR[20];
double paraboloid_RRH[20];

/* hyperboloid */
double hyperboloid_distance[20];
double hyperboloid_height[20];
double hyperboloid_radius[20];
double hyperboloid_AA[20];
double hyperboloid_BB[20];
double hyperboloid_RR[20];

/* column */
double column_height[20];
double column_R2[20];

/* poly */
int    polygon_Cnt[20];
int    vertex_Cnt[20][20];
double polygon_vertex[20][20][20][3];

/* change_axis */
int    Tr[20];
int    Ro[20];
int    Sh[20];
int    Sc[20];
double rxyz[20][3][3];
double shift[20][3];
double shearing[20][3][3];
double scaling[20][3];
double SH_AA[20];
double R_shearing[20][3][3];

/*-----------------------------Functions-----------------------------*/
/* read shape parameters form Config file */
int read_csg_shape (char* file_name);

/* read parameters of translation */
int read_translation (int i);

/* read parameters of rotation */
int read_rotation (int i, double *ax, double *ay, double *az);

/* read parameters of shearing */
int read_shearing (int i);

/* read parameters of scaling */
int read_scaling (int i);

int write_translation (int i);

int write_rotation (int i, double ax, double ay, double az);

int write_shear_scale (int i);

void change_axis (int Seq);

void trans_TT (int seq, double TT[]);

void freepath_csg (int projZ, int projM, float *energy, double flight_length,
                   double *step, int *current_material_index, int *ion_left_target,
                   float *stopping, float *straggling);

void interface_path (int s_count, double sect[][3], double *s, int projZ,
                     int projM, float *energy, double flight_length,
                     int *current_material_index, int *ion_left_target,
                     float *stopping, float *straggling);

int check_target_type_csg (double flight_length);

int substrate_interface_path (double flight_length, double *s);

void sect_sort (int *n, double sect[]);

int sphere_intersect (int seq, double TT[]);

int tetrahedron_intersect (int seq, double TT[]);

int triangle_intersect_csg (double P[][3], double *temp, int *flag);

int ellipsoid_intersect (int seq, double TT[]);

int plane2_intersect (double N[], double LVO, double D[], double T2[]);

int cuboid_intersect (int seq, double TT[]);

int taper_intersect (int seq, double TT[]);

int paraboloid_intersect (int seq, double TT[]);

int hyperboloid_intersect (int seq, double TT[]);

int column_intersect (int seq, double TT[]);

int polyhedron_intersect (int seq, double TT[]);

int polygon_intersect (double P_temp[][3], int cnt, double P0[]);

#endif
