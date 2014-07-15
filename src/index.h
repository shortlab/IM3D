/******************************************************************************
  Module Name : index.h
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description :

  Others :
      Refers to iradina and Corteo.

      Revision History:
      Date    Rel Ver.    Notes
      04/29/14            add corteo.mat in version_130715;
******************************************************************************/
#ifndef INDEX_H
#define INDEX_H

/*------------------------------Defines------------------------------*/
#define INDEX_BOUND_CHECKING  /* Perform bound cheking during index evaluation
                                 (whole program 10-15% slower). No check may result
                                 out-of-bound matrix access: crash or wrong value */

//#define MATTRIX_VERSION  /* ifdef, Corteo_130715 instead of Corteo_090527: corteo.mat */

#ifdef MATTRIX_VERSION  /* corteo.mat, version_130715 */
    /* matrix index calculation parameters */
    #define MINE   0.0000019073486328125f  /* minimum reduced energy = 2^(-19) */
    #define MAXE   2097152.f               /* maximum reduced energy = 2^21 */
    #define DIME   2560   /* matrix dimension: (21-(-19))*2^6 (8 mantissa bits: *2^6) */
    #define BIASE  6912   /* (127-19)*2^6: exponent bias correction while performing division by 2^(-19) */
    #define SHIFTE 17     /* keep 6 of the 23 mantissa bits (23-6=17) */

    #define MINS   0.00000001490116119384765625f  /* mimimum reduced impact parameter = 2^(-26) */
    #define MAXS   64.f   /* maximum reduced impact parameter = 2^6 */
    #define DIMS   2048   /* (6-(-26))*2^6 */
    #define BIASS  6464   /* (127-26)*2^6: exponent bias correction while performing division by 2^(-26) */
    #define SHIFTS 17     /* keep 6 of the 23 mantissa bits (23-6=17) */

    #define MIND   16.f          /* minimum energy (eV) for stopping power tables = 2^4 */
    #define MAXD   1073741824.f  /* maximum energy (eV) for stopping power tables = 2^30 ~ 1 GeV */
    #define DIMD   416           /* (30-4)*2^4 */
    #define BIASD  2096          /* (127+4)*2^4: exponent bias correction while performing division by 2^10 */
    #define SHIFTD 19            /* keep 4 of the 23 mantissa bits (23-4=19) */
#else  /* corteo.mat, version_090527 */
    /* matrix index calculation parameters */
    #define MINE  0.0000019073486328125f  /* minimum reduced energy = 2^(-19) */
    #define MAXE  2097152.f               /* maximum reduced energy = 2^21 */
    #define DIME  640   /* matrix dimension: (21-(-19))*2^4 (4 mantissa bits: *2^4) */
    #define BIASE 1728  /* (127-19)*2^4: exponent bias correction while performing division by 2^(-4) */
    #define SHIFTE 19   /* keep 4 of the 23 mantissa bits (23-4=19) */

    #define MINS  0.00000001490116119384765625f  /* mimimum reduced impact parameter = 2^(-26) */
    #define MAXS  64.f   /* maximum reduced impact parameter = 2^6 */
    #define DIMS  512    /* (6-(-26))*2^4 */
    #define BIASS 1616   /* (127-26)*2^4: exponent bias correction while performing division by 2^(-16) */
    #define SHIFTS 19    /* keep 4 of the 23 mantissa bits (23-4=19) */

    #define MIND  16.f          /* minimum energy (eV) for stopping power tables = 2^4 */
    #define MAXD  1073741824.f  /* maximum energy (eV) for stopping power tables = 2^30 ~ 1 GeV */
    #define DIMD  416           /* (30-4)*2^4 */
    #define BIASD 2096          /* (127+4)*2^4: exponent bias correction while performing division by 2^10 */
    #define SHIFTD 19           /* keep 4 of the 23 mantissa bits (23-4=19) */
#endif

/*--------------------------Global variables-------------------------*/
/* global variables: number of errors during index fuction calls
   (in work if INDEX_BOUND_CHECKING is defined) */
unsigned long E_min_err, E_max_err, S_min_err, S_max_err, D_min_err, D_max_err;

/* functions that compute an index from a value */
unsigned long E_index (float E_val);
unsigned long S_index (float S_val);
unsigned long D_index (float D_val);

/* functions that return the value corresponding to an index */
float E_val (unsigned long index);
float S_val (unsigned long index);
float D_val (unsigned long index);

#endif
