/******************************************************************************
  Module Name : mpimod.h
  Module Date : 04/22/2014
  Module Auth : Yonggang Li

  Description : MPI parallel module.

  Others :
      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef MPIMOD_H
#define MPIMOD_H

/*------------------------------Includes-----------------------------*/
#include "mpi.h"

#include "utils.h"
#include "target.h"

/*------------------------------Defines------------------------------*/
#define MPI_PRALLEL  /* serial code when commented out */

#define ROOT 0

#define MAX_SPUTTERED 5000
#define MAX_EL_PER_MAT 10

/*--------------------------Global variables-------------------------*/
int ierr;
int my_node, num_nodes;
int max_no_ions_node;  /* ion num per node */

/* reduced arrays */
int    *reduce_int_data;
double *reduce_double_data;
int     reduce_int_data2[8];
int     reduce_int_data3[MAX_SPUTTERED+1];
int     reduce_int_data4[MAX_EL_PER_MAT*8];

/*-----------------------------Functions-----------------------------*/
/* MPI random */
int  mpi_seed (int seed);

/* arrays used for massage transfer */
int  mpi_init_array (void);

/* distribute jobs to nodes */
void mpi_distribute (void);

/* reduced final output data */
void mpi_reduce_data (void);

#endif
