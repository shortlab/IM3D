/******************************************************************************
  Module Name : magic.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : MAGIC approximation for ion scattering angle.

  Others :
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef MAGIC_H
#define MAGIC_H

/*------------------------------Includes-----------------------------*/
#include <math.h>

/*-----------------------------Functions-----------------------------*/
/* returns cos(theta/2) for the given scattering parameters*/
double MAGIC_A (double B, double epsilon);

/* returns the ZBL potential, and via the pointer Vprime its derivative */
double ZBL_and_deri (double R, double *Vprime);

#endif
