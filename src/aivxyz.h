/******************************************************************************
  Module Name : aivxyz.h
  Module Date : 05/12/2015
  Module Auth : Yonggang Li, ygli@theory.issp.ac.cn

  Description : Output the space-positions (x, y, z) of ions (a), interstitials (i)
                and vacancies (v) in the form of .cfg.

  Others :
      Error numbers in this module 13000-13999.

      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#ifndef AIVXYZ_H
#define AIVXYZ_H

/*------------------------------Includes-----------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "im3d.h"

/*------------------------------Defines------------------------------*/
#define aiv_num 10000000

/*--------------------------Global variables-------------------------*/
struct AIV  /* properties of aiv */
{
    int type;       /* type of the defect: 0-ions, 1-interstitial, 2-vacancy. */
    int tab;        /* tab of the ion that a defect (I or V) belong to. */
    int mater;      /* material */
    int element;    /* element in the material */
    float x, y, z;  /* position of the defect */
};

int aiv_i;
//int aiv_num;  /* read from Config.h */
struct AIV *aiv_pos;

/*-----------------------------Functions-----------------------------*/
int init_aivxyz (void);   /* init aiv_pos. */

int store_aivxyz_cfg (char *base_name);  /* Store A-I-V positions in the form of .cfg */

#endif
