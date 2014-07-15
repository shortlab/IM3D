/******************************************************************************
  Module Name : const.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description : Global physical constants.

  Others :
      Refers to Corteo.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef CONST_H
#define CONST_H

/*------------------------------Defines------------------------------*/
#define PI 3.1415926535897932384626433832795
#define R_TO_D 57.29577951308233;  /* 180.0 / PI */
#define D_TO_R 0.01745329251994;   /* PI / 180.0 */

#define SCREENCONST 0.46848  /* 0.8853 * 0.5291772108, screening length constant [A] */
#define E2 14.3996445f       /* e^2 / (4 *pi * eps0) = e^2 c^2 in [eV][A] */

#define NUMBERELEMENTS 101   /* maximuum number of elements */

/*--------------------------Global variables-------------------------*/
float most_abundant_isotope[NUMBERELEMENTS];
float atomic_mass[NUMBERELEMENTS];  /* no use now, TODO: use instead of input */
char  atomic_names[NUMBERELEMENTS][3];

//TODO: add float density_of_element[NUMBERELEMENTS];

#endif
