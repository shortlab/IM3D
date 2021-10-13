/******************************************************************************
  Module Name : mshwriter.h
  Module Date : 05/02/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Outputfiles with MSH ASCII/binary file format.

  Others :
      Gmsh software's native "MSH" file format.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef MSHWRITER_H
#define MSHWRITER_H

/*------------------------------Includes-----------------------------*/
#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include <stdarg.h>

#include "im3d.h"
#include "target.h"
#include "transport.h"
#include "material.h"
#include "fileio.h"

/*------------------------------Defines------------------------------*/

/*--------------------------Global variables-------------------------*/
int msh_file_type;

int node_x, node_y, node_z;
int node_yz, node_count;

/*-----------------------------Functions-----------------------------*/
int store_results_msh (char *base_name);

int write_int_array_to_msh_file (char *file_name, int mat_i, int *source_array_total,
                                 int *source_array_part[], int count, int n_element);

int write_double_array_to_msh_file (char *file_name, double *source_array1,
                                    double *source_array2, int count);

int get_node_index (int x, int y, int z);

void get_node_XYZ (int index, int *x, int *y, int *z);

void get_node_relative_XYZ (int x, int y, int z, float *rx, float *ry, float *rz);

#endif
