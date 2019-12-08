/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

#include "alphabet.h"

int clean_and_set_to_extern(struct alphabet* a);
int merge_codes(struct alphabet*a,const int X, const int Y);

struct alphabet* create_dna_alphabet(void)
{
  struct alphabet* a = NULL;
  int i;
  char dnacode[16] = "ACGTUNRYSWKMBDHV";
  int code = 0;
  MMALLOC(a, sizeof(struct alphabet));

  for(i = 0; i < 128;i++) a->to_internal[i] = -1;
  for(i = 0; i < 32;i++)  a->to_external[i] = -1;
  for(i = 0; i < 16;i++) {
    a->to_internal[(int) dnacode[i]] = code;
    code++;
  }

  merge_codes(a,'U','T');
  merge_codes(a,'N','R');  /* R.................A or G */
  merge_codes(a,'N','Y');  /* Y.................C or T */
  merge_codes(a,'N','S');  /* S.................G or C */
  merge_codes(a,'N','W');  /* W.................A or T */
  merge_codes(a,'N','K');  /* K.................G or T */
  merge_codes(a,'N','M');  /* M.................A or C */
  merge_codes(a,'N','B');  /* B.................C or G or T */
  merge_codes(a,'N','D');  /* D.................A or G or T */
  merge_codes(a,'N','H');  /* H.................A or C or T */
  merge_codes(a,'N','V');  /* V.................A or C or G */

  clean_and_set_to_extern(a);
  return a;
}

int merge_codes(struct alphabet*a,const int X, const int Y)
{
  int min;
  min = MACRO_MIN(a->to_internal[X],a->to_internal[Y]);
  ASSERT(min != -1, "code not set!");
  a->to_internal[X] = min;
  a->to_internal[Y] = min;
  return OK;
ERROR:
  return FAIL;
}

int clean_and_set_to_extern(struct alphabet* a)
{
  int i;
  int code = 0;
  int8_t trans[32];
  for(i = 0; i < 32;i++) trans[i] = -1;

  for(i = 64; i < 96;i++) if(a->to_internal[i] != -1) trans[a->to_internal[i]] = 1;
  code = 0;
  for(i = 0; i < 32;i++) if(trans[i] == 1) {
    trans[i] = code;
    code++;
  }
  a->L = code;
  for(i = 64; i < 96;i++) if(a->to_internal[i] != -1) {
    a->to_internal[i] = trans[a->to_internal[i]];//a->to_internal[i]];
    a->to_internal[i+32] = a->to_internal[i];
  }
  for(i = 64;i < 96;i++) if(a->to_internal[i] != -1) a->to_external[a->to_internal[i]] = i;
  return OK;
}
