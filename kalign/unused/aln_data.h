/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

#ifndef ALN_DATA_H
#define ALN_DATA_H

#include "tldevel.h"

struct kalign_sequence {
        int* s;
        char* seq;
        int alloc_seq_len;
        int len;
        uint8_t* name;
        int id;
};

struct kalign_alignmment{
        struct kalign_sequence** s_arr;
        int numseq; /**< number of sequences */
        int alloc_numseq;
};

extern struct kalign_alignmment* kalign_aln_alloc(void);
extern int kalign_alignment_resize(struct kalign_alignmment* aln);
extern void kalign_alignmment_free(struct kalign_alignmment* aln);

extern int kalign_seq_resize(struct kalign_sequence* ks);

#endif
