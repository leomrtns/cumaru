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
#include "alignment_parameters.h"
#include "bisectingKmeans.h"
#include "alignment.h"
#include "misc.h"

char_vector kalign3_from_char_vector (char_vector dna)
{
  struct msa* msa = NULL;
  struct aln_param* ap = NULL;
  int** map = NULL; /* holds all alignment paths  */
  int i;
  char_vector aligned_charvector = NULL;

  /* Step 1: read all input sequences & figure out output  */
  msa = read_char_vector_to_msa (dna);
  ap = init_ap (msa->numseq);/* allocate aln parameters  */
  /* Start bi-secting K-means sequence clustering */
  build_tree_kmeans(msa,ap);
  /* Start alignment stuff */
  map = hirschberg_alignment(msa, ap);
  weave(msa , map, ap->tree); /* it's aligned already */
  /* clean up map */
  for (i = 0; i < msa->num_profiles ;i++) if(map[i]) free (map[i]); 
  if (map) free (map);
  /* We are done. */
  aligned_charvector = aligned_msa_to_charvector (msa);
  free_msa(msa);
  free_ap(ap);
  return aligned_charvector;
}
