/******************************************************************************
  Module Name : init.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Contains some initialization functions.

  Others :
      Refer to iradina and Corteo.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef init_H
#define init_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "transport.h"
#include "fileio.h"
#include "im3d.h"
#include "matrix.h"
#include "random.h"
#include "utils.h"
#include "index64.h"
#include "material.h"
#include "target.h"
#include "const.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/* or ------------------------------------*/
//#include "matrix.h"

/*------------------------------Defines------------------------------*/
#define MAXERFLIST 79999

/*--------------------------Global variables-------------------------*/
/* sqrt function related tables, beforehand calculates in order to increase efficiency */
float inv_sqrt_table_exp[256];
float sqrt_table_exp[256];
float inv_sqrt_table[1<<16];
float sqrt_table[1<<16];

float inverse_erf_list[MAXERFLIST];  /* list of evenly distributed but randomly ordered
                          values of inverse erf of -1+1/MAXERFLIST to 1-1/MAXERFLIST */

unsigned int erf_list_pointer;  /* points to next erf element to use */
unsigned int azim_angle;        /* points to next azimuthal angle to choose */
unsigned int ran_list;          /* points to next entry in the random list */
unsigned int ran_log_list;      /* points to next entry in the random sqrt logarithmic list */

/*-----------------------------Functions-----------------------------*/
/* read configuration from file, initialize variables etc. */
int init_configuration (char *ConfigFileName);

void  fill_fast_sqrt_table (void);

/* load tabulated Chu straggling data */
int load_Chu_straggling_values (void);

/* function adapted from corteo.c, loads a list of inverse error function erfinv(x) values */
int load_inverse_Erf (void);

#endif
