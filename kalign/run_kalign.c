/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

/*! \file
 *  \brief 
 *
 *  Algorithms in this file based on [kalign3](https://github.com/TimoLassmann/kalign.git)(GPL-3.0-or-later  &copy;2006, 2019 Timo Lassmann
 */

#include <biomcmc.h>
#include "global.h"
#include "msa.h"
#include "parameters.h"
#include "alignment_parameters.h"
#include "bisectingKmeans.h"
#include "alignment.h"
#include "weave_alignment.h"
#include "misc.h"
#include "alphabet.h"


int run_kalign(struct parameters* param)
{
  struct msa* msa = NULL;
  struct aln_param* ap = NULL;
  int** map = NULL; /* holds all alignment paths  */
  int i;

  /* Step 1: read all input sequences & figure out output  */
  for(i = 0; i < param->num_infiles;i++){ RUNP(msa = read_input(param->infile[i],msa)); }
  LOG_MSG("Detected: %d sequences.", msa->numseq);

  /* allocate aln parameters  */
  RUNP(ap = init_ap(msa->numseq,msa->L ));
  /* Start bi-secting K-means sequence clustering */
  RUN(build_tree_kmeans(msa,ap));
  convert_msa_to_internal(msa, defDNA); // from "detect_alphabet"
  /* Start alignment stuff */
  RUNP(map = hirschberg_alignment(msa, ap));
  /* set to aligned */
  msa->aligned = 1;
  RUN(weave(msa , map, ap->tree));
  /* clean up map */
  for(i = 0; i < msa->num_profiles ;i++){ if(map[i]) MFREE(map[i]); }
  MFREE(map);
  map = NULL;
  /* We are done. */
  RUN(write_msa(msa, param->outfile, param->out_format));

  free_msa(msa);
  free_ap(ap);
  return OK;
  free_parameters(param);
ERROR:
  return FAIL;
  free_parameters(param);
}
