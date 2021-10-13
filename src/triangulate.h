/******************************************************************************
  Module Name : triangulate.h
  Module Date : 04/08/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : FETM geometry file generation code.

  Others :
      Refers to HYWang's trianglen.f90 code.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef TRIANGULATE_H
#define TRIANGULATE_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fetm.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/*--------------------------Global variables-------------------------*/
/* input */
int file_format;  /* all list in Config.in */
int num_of_material;
int material_type[20];
int num_of_data[20];

/* line */
double xyz_flag;

/* str */
double **triangle;

/*-----------------------------Functions-----------------------------*/
/* generate_triangle */
int generate_fetm_tri_shape (void);

/* input, opengl format */
int read_file_opengl (int n, int num_of_data);

/* input, ply2 format */
int read_file_ply2 (int n);

/* calculate the direction of line */
void cal_direct (double xyz[], double *sin_theta, double *cos_theta,
                 double *sin_phi, double *cos_phi);

/* generates line */
void make_line_tri (int flag_or_not, double sin_theta, double cos_theta,
	                double sin_phi, double cos_phi, double x, double y,
	                double z, int *i, int *j, int *k);

#endif
