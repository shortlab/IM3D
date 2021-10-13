/******************************************************************************
  Module Name : transport.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Contains the functions to simulate the ion transport, etc.

  Others :
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef TRANSPORT_H
#define TRANSPORT_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <math.h>

#include "im3d.h"
#include "fileio.h"
#include "utils.h"
#include "random.h"
#include "target.h"
#include "bulk.h"
#include "csg.h"
#include "fetm.h"
#include "init.h"
#include "material.h"
#include "matrix.h"
#include "index64.h"
#include "const.h"
#include "magic.h"
#include "cfgwriter.h"
#include "mshwriter.h"
#include "vtkwriter.h"
#include "aivxyz.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/* or ------------------------------------*/
/* 3d geometry modules, ygli */
//#include "csg.h"   /* csg geometry */
//#include "fetm.h"  /* fetm geometry */

/*--------------------------Global variables-------------------------*/
/*======================================================================
  Declare some global variables (not nice, but fast) parameters of the ion beam.
======================================================================*/
int simulation_type;  /* How to do the simulation.
                      Plan: simulation type can be:
                      0 = Full Damage Cascade, follow recoils, free flight path statistically distr;
                      3 = Ions only. Recoils are not followed, No damage profiles stored.
                      NOTE: At the moment only 0 is correctly implemented. */

int detailed_sputtering;  /* Needs to be set to 1, if you want to get better results regarding
                          sputtering. If set to 0, the program doesn't care about surface binding
                          energy at all. However, this takes some calculation time, because in
                          the 3d geometry, it takes some time to find out whether a projectile
                          would be ejected to vacuum or not. */
int store_energy_deposit;  /* if 1, array with deposited energy are created and stored;
                              if 0, not. */

int flight_length_type;    /* Flight lengths between collisions can be selected by three
                              different approaches:
                      0: Poisson distributed flight lengths mit average of inter-atomic spacing;
                      1: constant flight length. Is has be be specified in units of nm.
                      Options 0 ignore the flight_length_constant parameter.
                      NOTE: ygli, 2015.1.8
                            In SRIM, different flight lengths are set for different energy regions,
                            much longer flight lengths are approximately used to increase the
                            efficiency of the calculation, which may introduce some errors due to
                            lower energy loss cost without considering the detailed traces.
                            When use flight_length_type=0, more detailed traces are considered,
                            which should introduce more energy loss and less displacements. */
float flight_length_const;  /* If constant flight length is selected, then this is it (in nm). */

int   max_no_ions;            /* maximum number of ions */
int   ion_Z;                  /* proton number */
float ion_M;                  /* mass of the ion */
float ion_initial_energy;     /* impinging energy */

float ion_vx;                /* vector of ion velocity, normalized to 1 */
float ion_vy;                /* Note, that the ion is NOT described by its actual velocity vector, */
float ion_vz;                /* but by the flying direction vector of length 1 and its energy. */

float min_energy;            /* minimum energy below which all projectiles are stopped */

int   ion_distribution;      /* 0 for random ion entry positions,
                                1 for centered, 2 for specified position,
                                3 for random square around position */
float enter_x, enter_y, enter_z;      /* entry point in nm */
float beam_spread;           /* in nm, only relevant for option 3 */

int    display_interval;        /* display status every so many ions */
int    override_max_ions;       /* for command-line argument override of maximum number of ions */
float  override_energy;         /* for command-line argument override of ion energy */
int    storage_interval;        /* dump the target arrays into the files every such number of ions */
int    status_update_interval;  /* After so many ions, the status is written to the status file
                                (if activated by command line arg). */
int store_transmitted_ions;  /* 1 if true */
int store_path_limit;        /* only for so many ions, the exact paths and recoil cascades are stored */
int transmission_pointer;    /* point to next free index in array */
int store_exiting_recoils;   /* 1 if all exiting recoils should be stored */
int store_exiting_limit;     /* maximum number of exiting recoils to be stored */

struct transmitted_ion       /* describes an ion that has been transmitted (left the target) */
{
    double x;                 /* exit position x */
    double y;                 /* y */
    double z;                 /* z */
    double vx;                /* exit velocity unit vector */
    double vy;
    double vz;
    float  energy;            /* exit energy */
};
struct transmitted_ion *transmit_list;  /* list of transmitted ions */

int leaving_ions[8];        /* number of ions, leaving target in each direction/quadrants, ygli */
int store_ion_paths;  /* 1 if the exact paths should be stored (interesting for debugging stuff */
FILE *ion_paths_fp;         /* the file with the ion paths will be open all the time if the paths
                               should be stored. This points to that file */
int store_recoil_cascades;  /* 1 if the exact recoils cascades should be stored
                               (interesting for debugging stuff...) */
FILE *recoil_cascades_fp;   /* the file with the recoil cascades will be open all the time if the
                               paths should be stored. This points to that file */

float chu_values[98][4];    /* Values to calculate straggling according to Chu's model;
                               fit data from Yang et al. NIMB61(1991)149. */
int straggling_model;       /* how to calc straggling */

int max_annular_coll_volumes;  /* According to W.Eckstein "Computer Simulation if Ion-Solid
                                  Interactions", Springer 1991, p.93, multiple collisions in
                                  annular volumes should be allowed to occur. This number +1
                                  determines the maximum number of collisions. 0 means just 1.
                                  Recommended for sputtering is 2. */

int tracing_recoil_or_not;      /* 1: tracing recoil, 0: only tracing primary ions, KP model, Aug. 05, 2014 */
int scattering_calculation;     /* 0: corteo database, 1: MAGIC  */
int transport_type;             /* 0: accurate, 1: Fast (like corteo) */
int single_ion_sputter_yields;  /* if 1, im3d will store sputter yields for single ions
                  (at the moment those are not separated by type of sputtered particles */
int *sputter_yield_histogram;   /* array that stores single ion sputter yield histogram */
#define MAX_SPUTTERED 5000      /* maximum number of sputtered particles per single ion
                                   possible to store */

/*-----------------------------Functions-----------------------------*/
int irradiate_target (void);  /* Let ions impinge on the target */

void def_ion_enetry_pos (double *x, double *y, double *z);  /* Define ion entry positions */

/* Calculates the path of a projectile with given properties (proton number,
   mass, energy, coordinates, direction) through the target material.
   Can be called recursively to follow recoils.
   The proj_state char contains status bits:
   bit 0: if 1, then the projectile is a recoil which has gained less than
          its displacement energy. That means it may fly around a little and
          kick other atoms, but it probably jumps back to its original locations.
   bit 1: if 1, the projectile replaced another atom, but is still flying
          around a bit and creating damage or possibly might try to escape
          from the solid. */
int projectile_transport (int ion_i, int projZ, float projM, float projE, double proj_x,
	                        double proj_y, double proj_z, double proj_vx, double proj_vy,
                          double proj_vz, int is_ion, int in_or_not, int org_material,
                          int org_element, int org_cell);

int collision (int ion_i, int projZ, float projM, int is_ion, float *energy, double target_x,
               double target_y, double target_z, double *vx, double *vy, double *vz,
               float impact_par, unsigned int *iazim_angle, int target_material_index,
               int org_material, int org_element, int *replaced, int in_or_not);

/* calculates the electronic dedx and energy loss straggling for ion with Z, M and
   energy, in given material [eV/A] */
void e_energy_loss (int proj_Z, int proj_M, float projE, int mater_type, double s,
                    float *stopping, float *straggling);  /* ygli */

/* corteo function to rotate the flying direction as function of scattering angle */
void rotate (double *l, double *m, double *n, unsigned int *iazim_angle,
             float cos_theta, float sin_theta);

/* performs surface refraction of a projectile.
   The direction and energy (passed as refs) are adjusted.
   n is the surface normal vector and E_surf the surface binding energy.
   The surface vector must point to vacuum!
   returns 0 if projectile went into solid.
   returns 1 if projectile tried to leave solid and succeed.
   returns 2 if projectile tried to leave solid but did not succeed. */
int refract_projectile2 (double *vx, double *vy, double *vz, double nx, double ny,
                         double nz, float *energy, float E_surf);

/* update x, y, z */
void update_xyz (double *x, double *y, double *z, double vx, double vy, double vz, double s);

#endif
