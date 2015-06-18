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
#include "random.h"

/*=============================================================================
  Function Name : randomx
  Description   : Random function from P. L'Ecuyer,
                  Communications of the ACM 31 (1988) 742.

  Inputs  : no.
  Outputs : no.

  Notes :
      Corteo: Adapted from randomx.c.
=============================================================================*/
double randomx (void) {
    int z, k;

    k = seed1 / 53668;
    seed1 = 40014 * (seed1 - k * 53668) - k * 12211;
    if (seed1 < 0) seed1 = seed1 + 2147483563;

    k = seed2 / 52774;
    seed2 = 40692 * (seed2 - k * 52774) - k * 3791;
    if (seed2 < 0) seed2 = seed2 + 2147483399;

    z = seed1 - seed2;
    if (z < 1) z = z + 2147483562;
    return z * 4.65661305956E-10;
}

/*=============================================================================
  Function Name : compute_lists
  Description   : Generate and randomize lists random_list, sqrt_log_list,
                  sin_azim_angle & cos_azim_angle.

  Inputs  : no.
  Outputs : no.

  Notes :
      Corteo: Functions adapted from corteo.c.

      Function call:
        fromcorteo - void randomize_list (float *list, unsigned int max_list);
                     float sqrtdf (double val);
                     float d2f (double val).
=============================================================================*/
void compute_lists (void) {
    unsigned int irandom_list, ilog_list, iazim_angle, iran_mu;
    float yu, mu;

    /* generate a uniformly spaced list of values between .5/MAXRANLIST and 1-.5/MAXRANLIST */
    for (irandom_list=0; irandom_list<MAXRANLIST; irandom_list++)
        random_list[irandom_list] = (irandom_list + 0.5f) / MAXRANLIST;
    randomize_list (random_list, MAXRANLIST);  /* put the list in random order */
    for (irandom_list=0; irandom_list<MAXRANLIST; irandom_list++)
        /* also compute sqrt of these values */
        sqrt_random_list[irandom_list] = sqrtdf (random_list[irandom_list]);

    /* Gauss distribution randomlist, Box-Muller algorithm,
       p(x) = 1/sqrt(2*pi) * exp(-x^2/2), x -> sigma*x + mean_x */
    for (irandom_list=0; irandom_list<MAXRANLIST; irandom_list++) {
        yu = random_list[irandom_list];
        iran_mu = (int) (yu * (float) MAXRANLIST);
        mu = random_list[iran_mu];
        gauss_random_list[irandom_list] = sqrtdf (-2.0 * log (yu)) * cos (2.8 * PI * mu);
    }


    /* produce a list of MAXLOGLIST sqrt(-log(x)) values for x uniformly distributed
       between x=.5/MAXLOGLIST and x=1-.5/MAXLOGLIST */
    for (ilog_list=0; ilog_list<MAXLOGLIST; ilog_list++)
        /* important modification, sqrt_log_list -> log_list */
        log_list[ilog_list] = - log ((ilog_list + .5) / MAXLOGLIST);
    randomize_list (log_list, MAXLOGLIST);  /* put the list in random order */

    /* pre-compute 1/sqrtloglist[] */
    for (ilog_list=0; ilog_list<MAXLOGLIST; ilog_list++)
        /* important modification, sqrt_log_list -> log_list */
        inv_sqrt_log_list[ilog_list] = 1.0f / sqrtdf (log_list[ilog_list]);

    /* produce a list uniformly distributed but randomly ordered azimuthal angles  */
    for (iazim_angle=0; iazim_angle<MAXAZILIST; iazim_angle++)
        /* cos_azim_angle temporary contains angles */
        cos_azim_angle[iazim_angle] = 2.0f * (float) PI * (float) (iazim_angle)
                                    / (float) (MAXAZILIST);
    randomize_list (cos_azim_angle, MAXAZILIST);  /* put the list in random order */
    for(iazim_angle=0; iazim_angle<MAXAZILIST; iazim_angle++) {
        /* compute the cos and sine of these angles */
        sin_azim_angle[iazim_angle] = d2f (sin (cos_azim_angle[iazim_angle]));
        cos_azim_angle[iazim_angle] = d2f (cos (cos_azim_angle[iazim_angle]));
    }

    return;
}

/*=============================================================================
  Function Name : randomize_list
  Description   : Randomly swap 3 times each elements of list so they come out
                  in random order.

  Inputs  : unsigned int maxlist
  Outputs : float *list

  Notes :
      Corteo: Adapted from corteoutil.c.

      Function call:
        fromcorteo - double randomx ().
=============================================================================*/
void randomize_list (float *list, unsigned int max_list) {
    unsigned int i, ilist, j;
    float x;

    for (i=0; i<max_list*3; i++) {
        ilist = i%max_list;
        j = (unsigned int) floor (randomx () * max_list);
        x = list[ilist];
        list[ilist] = list[j];
        list[j] = x;
    }

    return;
}
