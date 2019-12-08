/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

#ifndef RWALIGN_H
#define RWALIGN_H

#include <biomcmc.h>
#include <unistd.h>
#include "parameters.h"
#include "global.h"
#include "msa.h"
#include "alphabet.h"

extern struct alignment* read_alignment(char* infile);
extern struct alignment* detect_and_read_sequences(struct parameters* param);
extern int make_dna(struct alignment* aln);
extern void free_aln(struct alignment* aln);
extern int convert_alignment_to_internal(struct alignment* aln, int type);
extern int dealign(struct alignment* aln);
extern int output(struct alignment* aln,struct parameters* param);
extern int make_aliged_seq(uint8_t* aligned, uint8_t* unaligned, int* gaps,int len);

#endif
