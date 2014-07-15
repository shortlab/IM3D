/******************************************************************************
  Module Name : bulk.h
  Module Date : 04/24/2014
  Module Auth : Yonggang Li

  Description : Bulk or multi-layer samples.

  Others :
      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef BULK_H
#define BULK_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>

#include "init.h"
#include "target.h"
#include "fetm.h"
#include "material.h"
#include "transport.h"
#include "random.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/*--------------------------Global variables-------------------------*/
FILE *shape_file;  /* shape parameters input file */

int max_no_layers;
int material_of_layer[21];
double thick_of_layer[21];
float sect[21], sect0;

/*-----------------------------Functions-----------------------------*/
/* read bulk shape parameters form Config file */
int read_bulk_shape (char *file_name);

/* free flight path of ion for bulk or multi-layer samples */
void freepath_bulk (int projZ, int projM, float *energy, double z, double vz,
                    double flight_length, double *s, int *current_material_index,
                    int *ion_left_target, float *stopping, float *straggling);

/* check if projectile is still in the bulk/layer */
int check_target_type_bulk (double z, double vz, double flight_length);

/* calculate the distance from ion to the substrate surface */
double dist_to_surf (double z, double z0, double vz);

#endif
