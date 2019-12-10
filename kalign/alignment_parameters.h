/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

#ifndef ALIGNMENT_PARAMETERS_H
#define ALIGNMENT_PARAMETERS_H

#include <biomcmc.h>
#include "global.h"
#include "misc.h"

#define MSA_NAME_LEN 128
#define FORMAT_FA 1
#define FORMAT_MSF 2
#define FORMAT_CLU 3
#define defPROTEIN 21
#define redPROTEIN 13
#define defDNA 5

struct msa_seq {
  char* seq;
  uint32_t id; // leo: we don't store seq names here, external char_vector should have it
  uint8_t* s;
  int* gaps;
  int len;
};

struct msa {
  struct msa_seq** sequences;
  int** sip;
  int* nsip;
  int* plen;
  int numseq;
  int num_profiles;
};

struct aln_param {
  struct rng_state* rng;
  float** subm;
  float gpo;
  float gpe;
  float tgpe;
  int* tree;
};

struct alphabet {
  int8_t to_internal[128];
  int8_t to_external[32];
  int type;
};

/**< main input function **/
struct msa* read_char_vector_to_msa (char_vector dna);
/**< main output function **/
char_vector aligned_msa_to_charvector (struct msa* msa);

struct aln_param* init_ap (int numseq);
void free_ap (struct aln_param* ap);
struct alphabet* create_dna_alphabet (void);
int* pick_anchor (struct msa* msa, int* n);
int make_aliged_seq (uint8_t* aligned, uint8_t* unaligned, int* gaps,int len);

void convert_msa_to_internal(struct msa* msa);/* convert */
int write_msa(struct msa* msa, char* outfile, int type);
void free_msa(struct msa* msa);

int weave(struct msa* msa, int** map, int* tree);
int clean_aln(struct msa* msa);
// euclidean_distance.h
extern  int edist_256(const float* a,const float* b, const int len, float* ret);
extern int edist_serial(const float* a,const float* b,const int len, float* ret);
extern int edist_serial_d(const double* a,const double* b,const int len, double* ret);

#endif
