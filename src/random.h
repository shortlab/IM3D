/******************************************************************************
  Module Name : random.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description :

  Others :
      Refers to Corteo.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef RANDOM_H
#define RANDOM_H

/*------------------------------Includes-----------------------------*/
#include <math.h>

#include "utils.h"
#include "const.h"

/* or ------------------------------------*/
//#include <math.h>

//#include "utils.h"

/*------------------------------Defines------------------------------*/
/* prime numbers for lists length */
#define MAXLOGLIST 104729
#define MAXAZILIST 72211
#define MAXRANLIST 1000003

/*--------------------------Global variables-------------------------*/
int seed1, seed2;

/* list of evenly distributed but randomly ordered values between 0 and 1 */
float random_list[MAXRANLIST];
float sqrt_random_list[MAXRANLIST];   /* sqrt of randomlist */
float gauss_random_list[MAXRANLIST];  /* Gauss distribution of randomlist */

/* list of evenly distributed but randomly ordered values of sqrt of -log of 1/MAXLOGLIST to 1 */
float log_list[MAXLOGLIST];  /* important modification, sqrt_log_list -> log_list */
float inv_sqrt_log_list[MAXLOGLIST];  /* 1/sqrtloglist */
float sin_azim_angle[MAXAZILIST];  /* list cos and sin components of angles... */

/* ...this angle are evenly distributed but randomly ordered between 0 and 2*PI */
float cos_azim_angle[MAXAZILIST];

/*-----------------------------Functions-----------------------------*/
double randomx (void);

void compute_lists (void);  /* adapted from corteo.h */

void randomize_list (float *list, unsigned int max_list);

#endif
