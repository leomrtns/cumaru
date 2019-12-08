/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

#ifndef ALPHABET_H
#define ALPHABET_H

#include "global.h"

struct alphabet {
        int8_t to_internal[128];
        int8_t to_external[32];
        int type;
        int L;
};

extern struct alphabet* create_dna_alphabet (void);

#endif
