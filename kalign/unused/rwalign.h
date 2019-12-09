/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

#ifndef RWALIGN_H
#define RWALIGN_H

#include <biomcmc.h>
//#include <unistd.h>
//#include <stdint.h>
#include "parameters.h"
#include "global.h"
#include "alignment_parameters.h"

#define MSA_NAME_LEN 128
#define FORMAT_FA 1
#define FORMAT_MSF 2
#define FORMAT_CLU 3

struct msa_seq {
  char* name;
  char* seq;
  uint8_t* s;
  int* gaps;
  int len;
  int alloc_len;
};

struct msa {
  struct msa_seq** sequences;
  int** sip;
  int* nsip;
  int* plen;
  int numseq;
  int num_profiles;
  int alloc_numseq;
  int aligned;
  int letter_freq[128]; // leo:: used only to detect dna/prot/aligned (can be removed)
  int L;
};

extern struct alignment* detect_and_read_sequences(struct parameters* param);
extern int make_dna(struct alignment* aln);
extern void free_aln(struct alignment* aln);
extern int output (struct alignment* aln,struct parameters* param);
extern int make_aliged_seq (uint8_t* aligned, uint8_t* unaligned, int* gaps,int len);

int convert_msa_to_internal(struct msa* msa);/* convert */
struct msa* read_input(char* infile,struct msa* msa);/* rw functions */
int write_msa(struct msa* msa, char* outfile, int type);
void free_msa(struct msa* msa);

#endif
