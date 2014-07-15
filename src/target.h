/******************************************************************************
  Module Name : target.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description : Main program.

  Others :
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef TARGET_H
#define TARGET_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "iran3d.h"
#include "fileio.h"
#include "utils.h"
#include "transport.h"
#include "material.h"
#include "csg.h"
#include "fetm.h"
#include "msh.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/* or ------------------------------------*/
//#include "material.h"
//#include "csg.h"
//#include "fetm.h"

/*------------------------------Defines------------------------------*/
/* Declare some global variables (not nice, but easier) */

/*======================================================================
  Different "materials" can be defined. A material can be single-elemental
  as well as a compound.
  The target consists of (possibly a large number of) little cells. Each of
  these is filled with one of the defined materials.
  Alternatively, for the dynamic version, a number of seperate elements are
  defined. Each cell has a composition vector describing the fractions of
  each element in that cell. This is called element-based in contrast to
  "material-based".
======================================================================*/
#define MAX_STOPPING_ENTRIES 1000  /* maximum number of values for stopping table */

//#define isnan(x) ((x) != (x))      /* ERROR: VF cannot find isnan, ygli */

/*--------------------------Global variables-------------------------*/
/*======================================================================
  The target is a cuboid-shaped box, which consists of (possibly a large number
  of) of small equal-sized cuboid cells. Several arrays describe the target
  composition, the distributions of implanted ions, recoils, vacancies and so
  on. The dimension of the target is defined by the number and size of the
  small cells in each direction.
  Furthermore, boundary conditions need to be defined in each target direction.
  x is the primary incident ion beam direction. (Though the ion beam may hit
  the target under any angle).
  Note, that integer coordinates always refer to a cell x, y and z number,
  while float coordinates always to refer to nm positions in the target.

  ygli:
  With a half-infinite substrate in the target will make the target more
  realistic and the program more efficient.
  Only the complex structures will be constructed by different geometry
  methods, the other parts can be described by a uniform region with the
  same composition.
  The complex structures can be defined in a special region in term of
  (x0, x0+dx), (y0, y0+dy) and (z0, z0+dz). if substrate = 1, the tracing ion
  is in the substrate region, otherwise it is in a special region.
  z=0 is at the surface of substrate region.
  So, it is convenient for the case of only one special region included(Materials+1).
  For the case of many special structures included in the CSG method, still
  only one special region is included, for the regions without CSG structures
  are also definded as substrate.
  Statistic in substrate region can be counted with a increasing non-uniform mesh.
  Different ions and defects can be saved as their exact coordinates.
======================================================================*/
int geometry_type;     /* 0 - CSG, 1 - FETM */
int no_substrate;      /* 0 - with substrate, 1 - substrate is vaccum */
int gen_shape_or_not;  /* 1 - generate fetm shape file for fetm geometry */

int cell_count_x;      /* number of cells in x-direction (>=1) */
int cell_count_y;      /* number of cells in y-direction (>=1) */
int cell_count_z;      /* number of cells in z-direction (>=1) */
int layer_count_yz;    /* layer_count_xy = cell_count_x * cell_count_y,
                          still use layer_count_yz is in order to make
                          output in style of (x, y, z, value) */
int cell_count;        /* total number of cells */

float cell_size_x;     /* size of cells in x-direction in nm */
float cell_size_y;     /* size of cells in y-direction in nm */
float cell_size_z;     /* size of cells in z-direction in nm */
float cell_volume;     /* product of the above three */

float sub_surf_z;      /* z coordinate of the half-infinite substrate surface in nm */
float target_size_x;   /* size of target in x-direction in nm */
float target_size_y;   /* they are calculated */
float target_size_z;
float target_max_x;    /* maximum allowed position in x-direction in nm */
float target_max_y;    /* these positions are the largest possible values, */
float target_max_z;    /* which are smaller target_size_x etc. */

/* We need to declare some target arrays (or actually pointers to them) */
int   *target_composition;          /* material of each cell */
int   use_density_mult;             /* 0 if multiplicators not used, 1 if used */
int   *target_implanted_ions;       /* implanted ions per cell (interstitials+replacements) */
int   *target_replacing_ions;       /* implanted ions that replaced identical target atoms */
int   *target_total_vacancies;      /* vacancies per cell (of all types) */
int   *target_total_replacements;   /* replacements per cell (of all types) */
int   *target_total_displacements;  /* displacements per cell (of all types) */ /* Disp = Vac + Repl */
int   *target_total_interstitials;  /* sum of interstitials per cell (of all recoil
                                     types, NOT including the implanted ions) */

int  *target_total_sputtered;       /* sum of sputtered atoms from each cell */

int total_sputter_counter[8];       /* number of target atoms leaving the sample in each of the
                                       8 directions/quadrants */
char *TargetCompositionFileName;
int target_composition_file_type;   /* The file which hold the info of what material is in
                                       which cell can have two types: either just one column of
                                       values, which are indexed like the TargetIndexFunction does;
                                       or: four columns with x,y,z values and the material index */

double *target_energy_phonons;      /* all energy deposited into the phononic system */
double *target_energy_electrons;    /* all energy deposited by electronic stopping.
    Note: these two arrays must be of type double, because small values might be added to
          gigantic number, and should not be lost (happens for large number of ions and
          high energys. Doubles reduce this risk as compared to floats) */

/* depth distribution statistics */
int   *target_depth_implanted_ions;
int   *target_depth_replacing_ions;
int   *target_depth_total_vacancies;
int   *target_depth_total_replacements;
int   *target_depth_total_displacements;
int   *target_depth_total_interstitials;

double *target_depth_energy_phonons;
double *target_depth_energy_electrons;

/*-----------------------------Functions-----------------------------*/
/* read target structure (size etc.) from config file and reads target concentration array etc.*/
int init_target_structure (char *file_name);

/* needs to be called from the ini file reader while the traget structure input file is read */
int read_target_structure_data_block (char *block_name);
/* Needs to be called from the ini file reader while the target structure input file is read */
int read_target_structure_data (char *par_name, char *par_value);

/* calculates index for accessing one-dimensional target arrays knowing
   the integer cell coordinates */
int get_target_index (int x, int y, int z);

/* calculates the x,y,z coordinates (cell numbers) for a given index ...
   so the inverse function of TargetIndex() */
void get_target_XYZ (int index, int *x, int *y, int *z);

/* Calculates the relative x,y,z coordinates (unit: traget_size) for a given position */
void get_relative_XYZ (int x, int y, int z, float *rx, float *ry, float *rz);

/* calcualtes index for accessing one-dimensional target arrays from a
   float-value point in space */
int get_cell_index (double x, double y, double z, int *cell_i);

/* this calculates the direction in which the target nucleus is found.
   v is the projectile velocity vector, the IP vector is returned in p components */
void cal_relative_target_atom_position (double vx, double vy, double vz, double *px,
                                        double *py, double *pz, unsigned int iazim_angle);

#endif
