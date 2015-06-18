/******************************************************************************
  Module Name : cfgwriter.h
  Module Date : 09/03/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Outputfiles with .cfg file format.

  Others :
      Revision History:
      Date    Rel Ver.    Notes
      04/01/2014  1.0.0  first coded and integrated in module fileio
      09/03/2014  1.0.4  rewritten in module cfgwriter
******************************************************************************/
#ifndef CFGWRITER_H
#define CFGWRITER_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <string.h>

#include "im3d.h"
#include "target.h"
#include "transport.h"
#include "fileio.h"

/*--------------------------Global variables-------------------------*/

/*-----------------------------Functions-----------------------------*/
/* store the results of the simulation (arrays with distribution of
   implanted ions, defects etc.) in .cfg format */
int store_results_cfg (char *base_name);

/* writes the designated array of Count elements into a .cfg file */
int write_int_array_to_cfg_file (char *file_name, int mat_i, int *source_array_total,
                                 int *source_array_part[], int count, int n_element);
/* writes the double array of energy deposited into a .cfg file */
int write_double_array_to_cfg_file (char *file_name, double *source_array1,
                                    double *source_array2, int count);

#endif
