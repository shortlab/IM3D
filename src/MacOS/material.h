/******************************************************************************
  Module Name : material.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Contains the material-related functions, etc.

  Others :
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef MATERIAL_H
#define MATERIAL_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fileio.h"
#include "im3d.h"
#include "transport.h"
#include "const.h"
#include "utils.h"
#include "matrix.h"
#include "target.h"
#include "index.h"

#ifdef MPI_PRALLEL
/* MPI=============================================== */
#include "mpimod.h"
/* MPI=============================================== */
#endif

/* or ------------------------------------*/
//#include "fileio.h"
//#include "init.h"

/*------------------------------Defines------------------------------*/
#define MAX_EL_PER_MAT 10          /* maximum number of elements per material */
#define MAX_NO_MATERIALS 20        /* maximum number of different materials */
#define MAX_ELEMENT_NO 93          /* maximum number of elements in simulation +1 */

/*--------------------------Global variables-------------------------*/
struct material                 /* all properties of a material */
{
    char  name[50];             /* a name can be defined */
    float density;              /* total atomic density in at/cm^3 */
    float density_NM;           /* total atomic density in at/nm^3 */
    float atomic_distance;      /* average inter-atomic distance [nm] */
    float rev_atomic_distance;  /* 1.0/AtomicDistance [nm^-1], ygli */
    float layer_distance;       /* layer distance assuming simple cubic structure[nm] */
    float mean_impact_par;         /* not actually the mean impact parameter but rather
                                      1/sqrt(PI*density*MeanFreePath) */
    float sqrt_rec_fl_density;     /* 1/sqrt(pi*flight_length_constant*density) needed for
                                      constant flight lengths */
    int   elements_Z[MAX_EL_PER_MAT];     /* list of the elements contained in this materials */
    float elements_M[MAX_EL_PER_MAT];     /* list of the element masses in this materials */
    float elements_conc[MAX_EL_PER_MAT];  /* relative concentration of the elements */
    float elements_disp_energy[MAX_EL_PER_MAT];  /* displacement energy for each element in eV */
    float elements_latt_energy[MAX_EL_PER_MAT];  /* lattice energy for each element in eV */
    float elements_surf_energy[MAX_EL_PER_MAT];  /* surface binding energy for each element in eV */
    float ion_surf_energy;                /* Surface binding energy of the ion in this material.
                                             If not provided: use mean */
    int   sputter_counter[MAX_EL_PER_MAT*8];  /* count sputtered atoms leaving the sample in each
                                                 of the possible 8 directions/quadrants
          Note: for all elements all 6 directions in one array, one after another: faster. */
    int   element_count;  /* number of different elements contained in this material */
    int   cell_count;  /* number of cells that consist of this material */
    float mean_Z;      /* average Z, weighted by fraction */
    float mean_M;      /* average M, weighted by fraction */
    float mean_F;      /* average F (reduced energy conversion factor), weighted by fraction */
    float mean_A;      /* average A (screening length), weighted by fraction, */
    float mean_min_red_transfer;  /* average minimum energy transfer in reduced units */
    float mean_Ed;     /* average Ed, displacement energy, Modified Kinchin-Pease Model, Aug. 05, 2014 */

    /* For each element of each material we want to store the recoils, interstitials and so on.
       The equations are:
       Displacements = Vacancies + Replacement Collisions
       Vacancies     = Interstitials + Atoms that left the target */
    /* So will need to make an array of pointers (one pointer for each element), pointing to
       the various target arrays. Further down, these four arrays occur once for the target but
       independently of the elements, and those will hold the sum for all elements */
    int **target_implanted_recoils_int;   /* the implanted recoils stopped as interstitials */
    int **target_implanted_recoils_repl;  /* the implanted recoils stopped as replacement atoms */
    int **target_elemental_vacancies;     /* vacancies of this element left behind */
    int **target_elemental_disp;          /* displacements of this element that took place */

    /* depth distribution statistics */
    int **target_depth_implanted_recoils_int;
    int **target_depth_implanted_recoils_repl;
    int **target_depth_elemental_vacancies;
    int **target_depth_elemental_disp;

    /* radial distribution statistics */
    int **target_radial_implanted_recoils_int;
    int **target_radial_implanted_recoils_repl;
    int **target_radial_elemental_vacancies;
    int **target_radial_elemental_disp;

    int **target_sputtered_atoms;       /* number of atoms sputtered from this cell */
    float **stopping_ZE;                /* points to an array of 92 elements (the Zs) which contain
                                           pointers to logarithmically scaled stopping tables for Z
                                           in this material.*/
    float **straggling_ZE;              /* we need the same for energy loss straggling */
    float *Nd;                          /* Modified Kinchin-Pease Model, Aug. 05, 2014 */
    float **Nd_Z;                       /* Modified Kinchin-Pease Model with different Z, Aug. 13, 2014 */
    struct transmitted_ion** elemental_leaving_recoils;  /* arrays storing recoils that are
                                                            leaving of each element */
    int *leaving_recoils_pointer;       /* array of pointers which for each element of this material
                                           point to next free leaving position */
};
struct material list_of_materials[MAX_NO_MATERIALS];
int    number_of_materials;             /* points to first free index in the ListOfMaterials */
int    number_of_all_target_elements;   /* without ion. */

int existing_elements[MAX_ELEMENT_NO];  /* This array holds a 1 for any element that might exist
                                           in the target, so all these elements might occur as
                                           projectiles. */
int hydrogen_in_target;                 /* Hydrogen is always included in the existing_elements_list,
                                           because we need its stopping powers. However, we do only
                                        need it in the element file, if it is really in the target.*/
int ionZ_in_target;                     /* The ion is always included in the existing_elements_list.
                                           However, we only need it as an extra element in the
                                           element file, if it is really in the target.*/

/*-----------------------------Functions-----------------------------*/
/* read materials from input file and do preparatory calculations */
int init_materials (char *file_name);
/* needs to be called from the init file reader while the materials input file is read */
int read_materials_data_block (char *material_name);
/* needs to be called from the init file reader while the materials input file is read */
int read_materials_data (char *par_name, char *par_value);

/* read stopping data from file and fill arrays etc. */
int prepare_stopping_tables (void);

/* create straggling tables etc. */
/* Note: there is a difference between the material and the elemental version of im3d:
   Omega is stored in the tables in the material version.
   In contrast, Omega^2 is stored in the tables in the element version. */
int prepare_straggling_tables (int model);

/* modified Kinchin-Pease model, Material version */
int prepare_KP_tables1 (void);
int prepare_KP_tables2 (void);

/* returns the number of ones in the provided array */
int count_existing_elements (int *element_array);

/* convert normal material-based input files to element-based input files */
int convert_material_to_element (char *output_file);

#endif
