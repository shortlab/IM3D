/******************************************************************************
  Module Name : magic.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li

  Description : MAGIC approximation for ion scattering angle.

  Others :
      Refers to iradina.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "magic.h"

/*=============================================================================
  Function Name : MAGIC_A
  Description   : MAGIC approximation method.

  Inputs  :
            double B - reduced impact par.
            double epsilon - reduced center of mass energy.
  Outputs : returns cos (theta/2) of the scattering event.

  Notes :
      Function call:
        utils - double ZBL_and_deri (double R, double* Vprime).
=============================================================================*/
double MAGIC_A (double B, double epsilon) {
    double cost2;  /* cos(theta/2)*/
    double RoC, delta, R, RR, A, G, alpha, beta, gamma, V, V1, FR, FR1, Q;
    double SQE;
    double C[6];

    C[1] = 0.99229;  /* TRIM 85: */
    C[2] = 0.011615;
    C[3] = 0.0071222;
    C[4] = 14.813;
    C[5] = 9.3066;

    /* initial guess for R: */
    R = B;
    RR = -2.7 * log (epsilon * B);
    if (RR >= B) {
        //if (RR < B) calc potential;
        RR = -2.7 * log (epsilon * RR);
        if (RR >= B) R = RR;
    }
    /* TRIM85: 330 */
    do {
        /* calculate potential and its derivative */
        V = ZBL_and_deri (R, &V1);

        /* excerpt from the TRIM95 code: */
        /*
        EX1 = 0.;
        if (R < 7) EX1 = .18175 * exp (-3.1998 * R);
        EX2 = .50986  * exp (-.94229 * R);
        EX3 = .28022  * exp (-.4029  * R);
        EX4 = .028171 * exp (-.20162 * R);
        V   = (EX1 + EX2 + EX3 + EX4) / R;
        V1  = -(V + 3.1998 * EX1 + .94229 * EX2 + .4029 * EX3 + .20162 * EX4) / R;
        */

        FR  = B * B / R + V * R / epsilon - R;
        FR1 = - B * B / (R * R) + (V + V1 * R) / epsilon - 1.0;
        Q   = FR / FR1;
        R   = R - Q;
    } while (fabs (Q / R) > 0.001);

    RoC = -2.0 * (epsilon - V) / V1;
    SQE = sqrt (epsilon);

    alpha = 1 + C[1] / SQE;
    beta  = (C[2] + SQE) / (C[3] + SQE);  /* TRIM85: CC */
    gamma = (C[4] + epsilon) / (C[5] + epsilon);
    A     = 2 * alpha * epsilon * pow (B, beta);
    G     = gamma / (sqrt ((1.0 + A * A)) - A);  /* TRIM85: 1/FF */
    delta = A * (R - B) / (1 + G);

    /* TRIM85: */
    /*
    CC = (0.011615 + SQE) / (0.0071222 + SQE);
    AA = 2.0 * epsilon * (1.0 + (0.99229 / SQE)) * pow (B, CC);
    FF = (sqrt (AA * AA + 1.0) - AA) * ((9.3066 + epsilon) / (14.813 + epsilon));
    Delta = (R - B) * AA * FF / (FF + 1.0);
    */

    cost2 = (B + RoC + delta) / (R + RoC);

    return cost2;
}

/*=============================================================================
  Function Name : ZBL_and_deri
  Description   : Return ZBL potential, and via the pointer Vprime its derivative,
                  Values are taken from ZBL85.

  Inputs  : double R
  Outputs : double *Vprime

  Notes : no.
=============================================================================*/
double ZBL_and_deri (double R, double *Vprime) {
    double EX1, EX2, EX3, EX4, V;
    /* corteo: return 0.1818 * exp (-3. * x) + 0.5099 * exp (-0.9423 * x)
                      + 0.2802 * exp (-0.4028 * x) + 0.02817 * exp (-0.2016 * x); */
    /*
    EX1 = 0.1818   * exp (-3.0    * R);
    EX2 = 0.5099   * exp (-0.9423 * R);
    EX3 = 0.2802   * exp (-0.4028 * R);
    EX4 = 0.02817  * exp (-0.2016 * R);
    V   = (EX1 + EX2 + EX3 + EX4) / R;
    *Vprime = -(V + 3.0 * EX1 + 0.9423 * EX2 + 0.4028 * EX3 + 0.2016 * EX4) / R;
    return V;
    */

    /* TRIM85: */
    EX1 = 0.18175 * exp (-3.1998 * R);
    /* if (R >= 7) EX1 = 0.0; */  /* According to TRIM95 */
    EX2 = 0.50986  * exp (-0.94229 * R);
    EX3 = 0.28022  * exp (-0.4029  * R);
    EX4 = 0.028171 * exp (-0.20162 * R);
    V = (EX1 + EX2 + EX3 + EX4) / R;
    *Vprime = -(V + 3.1998 * EX1 + 0.94229 * EX2 + 0.4029 * EX3 + 0.20162 * EX4) / R;

    return V;
}
