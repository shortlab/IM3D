/******************************************************************************
  Module Name : index64.c
  Module Date : 02/26/2014
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description :

  Others :
      Refers to iradina and Corteo.

      Revision History:
      Date    Rel Ver.    Notes
      15/06/30            long to int for for 32-bit to 64-bit;
******************************************************************************/
#include "index64.h"

/*=============================================================================
  X_index: compute an integer index using the binary representation of a
  floating point value val/MIN (X is E for reduced energy, S for reduced
  impact parameter or D for energy of stopping).
  These indexes are used to access sin/cos tables of the scattering angle
  and stopping power tables.
  The division by MIN is actually carried out by subtracting the appropriate
  BIAS from the resulting index.

  Input: val, with val >= MIN, (MIN = MINE, MINS, or MIND)

  Returns an integer equal to 16 times the exponent base 2 of val
  + 1 for each 1/16 interval between each power of 2
  (Index of input value = MIN is 0)

  examples
  val/MIN    index
  1           0
  1.0625      1
  1.0626      1
  1.1249      1
  1.1250      2
  1.2         3
  1.25        4
  2           16
  5           40

  BEWARE: architecture dependent!
  Assuming 32-bit 'unsigned long' integer and float
  Assuming IEEE Standard 754 representation for float:
  bit 31: sign
  bit 30-23: exponent (base 2), biased by 127 (i.e. value of exponent for 2^0 is 127)
  bit 22-0: mantissa
=============================================================================*/

/*=============================================================================
  Function Name : E_index
  Description   : Version that compute the index of a reduced energy.

  Inputs  :
            float E_val
  Outputs : no.

  Notes :
      Corteo: Adapted from corteoindex.c.
=============================================================================*/
unsigned int E_index (float E_val) {
    unsigned int ll;

    ll = *(unsigned int *)&E_val;
    ll = (ll >> SHIFTE) - BIASE;

#ifdef INDEX_BOUND_CHECKING
    if (E_val < MINE) {
        E_min_err++;
        return 0;
    }
    if (ll >= DIME) {
        E_max_err++;
        return DIME - 1;
    }
#endif

    return ll;
}

/*=============================================================================
  Function Name : S_index
  Description   : Version that compute the index of a reduced impact parameter.

  Inputs  :
            float S_val
  Outputs : no.

  Notes :
      Corteo: Adapted from corteoindex.c.
=============================================================================*/
unsigned int S_index (float S_val) {
    unsigned int ll;

    ll = *(unsigned int *)&S_val;
    ll = (ll >> SHIFTS) - BIASS;

#ifdef INDEX_BOUND_CHECKING
    if (S_val < MINS) {
        S_min_err++;
        return 0;
    }
    if (ll >= DIMS) {
        S_max_err++;
        return DIMS - 1;
    }
#endif

    return ll;
}

/*=============================================================================
  Function Name : D_index
  Description   : Version that compute the index of the energy of a stopping
                  power and energy straggling.

  Inputs  : float val
  Outputs : no.

  Notes :
      Corteo: Adapted from corteoindex.c.
=============================================================================*/
unsigned int D_index (float D_val) {
    unsigned int ll;

    ll = *(unsigned int *)&D_val;
    ll = (ll >> SHIFTD) - BIASD;

#ifdef INDEX_BOUND_CHECKING
    if (D_val < MIND) {
        D_min_err++;
        return 0;
    }
    if (ll >= DIMD) {
        D_max_err++;
        return DIMD - 1;
    }
#endif

    return ll;
}

/*=============================================================================
  Function Name : E_val
  Description   : Return the float value corresponding to index : reduced energy
                  (Actually returns the average value between the value of this
                  index and the next because when an index is computed, values
                  are truncated to the lower bound of the interval.)

  Inputs  : unsigned long index
  Outputs : no.

  Notes :
      Corteo: Adapted from corteoindex.c.
=============================================================================*/
float E_val (unsigned int index) {
    float temp1, temp2;
    unsigned int ll;

    ll = (index + BIASE) << SHIFTE;
    temp1 = (*(float *)&ll);
    ll = (index + 1 + BIASE) << SHIFTE;
    temp2 = (*(float *)&ll);

    return (temp1 + temp2) * 0.5f;
}

/*=============================================================================
  Function Name : S_val
  Description   : Reduced impact parameter.

  Inputs  : unsigned long index
  Outputs : no.

  Notes :
      Corteo: Adapted from corteoindex.c.
=============================================================================*/
float S_val (unsigned int index) {
    float temp1, temp2;
    unsigned int ll;

    ll = (index + BIASS) << SHIFTS;
    temp1 = (*(float *)&ll);
    ll = (index + 1 + BIASS) << SHIFTS;
    temp2 = (*(float *)&ll);

    return (temp1 + temp2) * 0.5f;
}

/*=============================================================================
  Function Name : D_val
  Description   : The energy of a stopping power and energy straggling.

  Inputs  : unsigned long index
  Outputs : no.

  Notes :
      Corteo: Adapted from corteoindex.c.
=============================================================================*/
float D_val (unsigned int index) {
    float temp1, temp2;
    unsigned int ll;

    ll = (index + BIASD) << SHIFTD;
    temp1 = (*(float *)&ll);
    ll = (index + 1 + BIASD) << SHIFTD;
    temp2 = (*(float *)&ll);

    return (temp1 + temp2) * 0.5f;
}
