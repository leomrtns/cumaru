/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

#include "alignment_parameters.h"

struct sort_struct{
  int len;
  int id;
};

/* only local; */
struct line_buffer {
  struct out_line** lines;
  int max_line_len;
  int alloc_num_lines;
  int num_line;
};

struct out_line {
  char* line;
  int block;
  int seq_id;
};

struct msa_seq* msa_seq_from_char_vector_string (char *string, size_t nchars, uint32_t id);
struct msa* read_char_vector_to_msa (char_vector dna);

void set_subm_gaps_DNA(struct aln_param* ap);
int clean_and_set_to_extern(struct alphabet* a);
void merge_codes (struct alphabet*a,const int X, const int Y);
int sort_by_len (const void *a, const void *b);
int* select_seqs (struct msa* msa, int num_anchor);

/* from rwalign.c */
void free_msa_seq(struct msa_seq* seq);
struct line_buffer* alloc_line_buffer(int max_line_len);
int resize_line_buffer(struct line_buffer* lb);
void free_line_buffer(struct line_buffer* lb);
static int set_sip_nsip(struct msa* msa);/* local helper functions  rwalign */
static int sort_out_lines(const void *a, const void *b);
static int make_linear_sequence(struct msa_seq* seq, char* linear_seq);
/* end rwalign.c */

struct aln_param* init_ap (int numseq)
{
  struct aln_param* ap =  (struct aln_param*) biomcmc_malloc (sizeof(struct aln_param));
  int i,j;

  ap->tree = NULL;
  ap->tree = (int*) biomcmc_malloc (sizeof(int) * (numseq*3+1));

  for(i = 0;i < (numseq*3+1);i++) ap->tree[i] = 0;

  ap->rng = init_rng(42);
  ap->subm = (float**) biomcmc_malloc (sizeof(float*) * 21);

  for (i = 21;i--;){
    ap->subm[i] = NULL;
    ap->subm[i] = (float*) biomcmc_malloc (sizeof(float) * 21);
    for (j = 21;j--;) ap->subm[i][j] = 0.0f;
  }
  set_subm_gaps_DNA (ap);
  return ap;
}

void free_ap (struct aln_param* ap)
{
  int i;
  if(!ap) return;
  if(ap->subm){
    for (i = 21;i--;) MFREE(ap->subm[i]);
    MFREE(ap->subm);
  }
  if(ap->rng) free_rng(ap->rng);
  if(ap->tree) MFREE(ap->tree);
  MFREE(ap);
}

/* Timo: These are old parameters from kalign 2 */
void set_subm_gaps_DNA (struct aln_param* ap)
{
  int i,j;
  for(i = 0; i < 5; i++) for(j =0; j < 5;j++) ap->subm[i][j] = 283;
  //	A   91 -114  -31 -123    0  -43
  ap->subm[0][0] += 91;
  ap->subm[0][1] += -114;
  ap->subm[0][2] += -31;
  ap->subm[0][3] += -123;
  //	C -114  100 -125  -31    0  -43
  ap->subm[1][0] += -114;
  ap->subm[1][1] += 100;
  ap->subm[1][2] += -125;
  ap->subm[1][3] += -31;
  //	G  -31 -125  100 -114    0  -43
  ap->subm[2][0] += -31;
  ap->subm[2][1] += -125;
  ap->subm[2][2] += 100;
  ap->subm[2][3] += -114;
  //	T -123  -31 -114   91    0  -43
  ap->subm[3][0] += -123;
  ap->subm[3][1] += -31;
  ap->subm[3][2] += -114;
  ap->subm[3][3] += 91;

  ap->gpo = 217;
  ap->gpe = 39.4;
  ap->tgpe =  292.6;
  //param->secret = 28.3;
}

// alphabet.c

struct alphabet* create_dna_alphabet (void)
{
  struct alphabet* a = (struct alphabet*) biomcmc_malloc (sizeof(struct alphabet));
  int i;
  char dnacode[16] = "ACGTUNRYSWKMBDHV";
  int code = 0;

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

void merge_codes (struct alphabet*a, const int X, const int Y)
{
  int min;
  min = MACRO_MIN(a->to_internal[X],a->to_internal[Y]);
  ASSERT(min != -1, "code not set!");
  a->to_internal[X] = min;
  a->to_internal[Y] = min;
}

int clean_and_set_to_extern (struct alphabet* a)
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
  for(i = 64; i < 96;i++) if(a->to_internal[i] != -1) {
    a->to_internal[i] = trans[a->to_internal[i]];//a->to_internal[i]];
    a->to_internal[i+32] = a->to_internal[i];
  }
  for(i = 64;i < 96;i++) if(a->to_internal[i] != -1) a->to_external[a->to_internal[i]] = i;
  return OK;
}

// pick_anchor.c

int* pick_anchor(struct msa* msa, int* n)
{
  int* anchors = NULL;
  int num_anchor = 0, powlog2;
  ASSERT(msa != NULL, "No alignment.");
  powlog2 = (int) pow(log2((double) msa->numseq), 2.0);
  num_anchor = MAX(MIN(32, msa->numseq), powlog2);
  anchors = select_seqs (msa, num_anchor);
  *n = num_anchor;
  return anchors;
}

int* select_seqs(struct msa* msa, int num_anchor)
{
  struct sort_struct** seq_sort = (struct sort_struct**) biomcmc_malloc (sizeof(struct sort_struct*) * msa->numseq);
  int* anchors = NULL;
  int i,stride;
  for(i = 0; i < msa->numseq;i++) {
    seq_sort[i] = (struct sort_struct*) biomcmc_malloc (sizeof(struct sort_struct));
    seq_sort[i]->id = i;
    seq_sort[i]->len = msa->sequences[i]->len;//  aln->sl[i];
  }

  qsort (seq_sort, msa->numseq, sizeof(struct sort_struct*),sort_by_len);
  anchors = (int*) biomcmc_malloc (sizeof(int) * num_anchor);
  stride = msa->numseq / num_anchor;
  for(i = 0; i < num_anchor;i++) anchors[i] = seq_sort[i*stride]->id;
  ASSERT(i == num_anchor,"Cound not select all anchors\tnum_anchor:%d\t numseq:%d",num_anchor,msa->numseq);
  for(i = 0; i < msa->numseq;i++) MFREE(seq_sort[i]);
  MFREE(seq_sort);
  return anchors;
}

int sort_by_len(const void *a, const void *b)
{
  struct sort_struct* const *one = a;
  struct sort_struct* const *two = b;
  if((*one)->len > (*two)->len) return -1;
  else  return 1;
}

// rwalign.c 

/**< only DNA sequences and integer IDs, no names */
struct msa* read_char_vector_to_msa (char_vector dna) 
{
  struct msa* msal = (struct msa*) biomcmc_malloc (sizeof (struct msa));
  int i;
  msal->sequences = NULL;
  msal->numseq = dna->nstrings;
  msal->plen = NULL;
  msal->sip = NULL;
  msal->nsip = NULL;
  msal->sequences = (struct msa_seq**) biomcmc_malloc (sizeof(struct msa_seq*) * msal->numseq);
  for (i = 0; i < msal->numseq; i++) msal->sequences[i] = msa_seq_from_char_vector_string (dna->string[i], dna->nchars[i], (uint32_t) i);
  convert_msa_to_internal (msal);
  set_sip_nsip (msal);
  return msal;
}

int set_sip_nsip(struct msa* msa)
{
  int i;
  ASSERT(msa!= NULL, "No msa");
  if(msa->plen){
    for (i = msa->num_profiles;i--;) if(msa->sip[i]) MFREE(msa->sip[i]);
    if(msa->plen) MFREE(msa->plen);
    if(msa->sip)  MFREE(msa->sip);
    if(msa->nsip) MFREE(msa->nsip);
    msa->plen = NULL;
    msa->sip = NULL;
    msa->nsip = NULL;
  }

  msa->num_profiles = (msa->numseq << 1 ) - 1;
  MMALLOC(msa->sip,sizeof(int*)* msa->num_profiles);
  MMALLOC(msa->nsip,sizeof(int)* msa->num_profiles);
  MMALLOC(msa->plen,sizeof(int)* msa->num_profiles);

  for (i =0;i < msa->num_profiles;i++){
    msa->sip[i] = NULL;
    msa->nsip[i] = 0;
  }

  for(i = 0;i < msa->numseq;i++) {
    MMALLOC(msa->sip[i],sizeof(int));
    msa->nsip[i] = 1;
    msa->sip[i][0] = i;
    msa->plen[i] = 0;
  }
  return OK;
ERROR:
  return FAIL;
}

void convert_msa_to_internal (struct msa* msa)
{
  struct alphabet* a = NULL;
  struct msa_seq* seq = NULL;
  int8_t* t = NULL;
  int i,j;

  a = create_dna_alphabet();
  t = a->to_internal;
  for(i = 0; i <  msa->numseq;i++){
    seq = msa->sequences[i];
    for(j =0 ; j < seq->len;j++){
      if(t[(int) seq->seq[j]] == -1){
        WARNING_MSG("there should be no character not matching the alphabet");
        WARNING_MSG("offending character: >>>%c<<<", seq->seq[j]);
      }else{
        seq->s[j] = t[(int) seq->seq[j]];
      }
    }
  }
  if (a) free (a);
}

int make_linear_sequence(struct msa_seq* seq, char* linear_seq)
{
  int c,j,f;
  f = 0;
  for(j = 0;j < seq->len;j++){
    for(c = 0;c < seq->gaps[j];c++){
      linear_seq[f] = '-';
      f++;

    }
    linear_seq[f] = seq->seq[j];
    f++;
  }
  for(c = 0;c < seq->gaps[ seq->len];c++){
    linear_seq[f] = '-';
    f++;
  }
  linear_seq[f] = 0;
  return OK;
}

char_vector aligned_msa_to_charvector (struct msa* msa)
{
  uint32_t c, f, i, j, id;
  size_t aln_len = msa->sequences[0]->len;
  char_vector aligned;

  for (j = 0; j <= msa->sequences[0]->len; j++) aln_len += msa->sequences[0]->gaps[j];

  aligned = new_char_vector_fixed_length (msa->numseq, aln_len);

  for(i = 0; i < msa->numseq;i++) {
    id = msa->sequences[i]->id; /* I don't know if the order can change so better safe than sorry */
    f = 0;
    for (j = 0; j < msa->sequences[i]->len; j++) {
      for (c = 0; c < msa->sequences[i]->gaps[j]; c++)  aligned->string[id][f++] = '-';
      aligned->string[id][f++] = msa->sequences[i]->seq[j];
    }
    for (c = 0; c < msa->sequences[i]->gaps[j]; c++) aligned->string[id][f++] = '-';
  }
}

void free_msa (struct msa* msa)
{
  int i;
  if(!msa) return;
  for(i = msa->numseq - 1; i>=0; --i) free_msa_seq (msa->sequences[i]);
  for (i = msa->num_profiles;i--;) if(msa->sip[i]) free(msa->sip[i]);
  if (msa->plen) free (msa->plen);
  if (msa->sip)  free (msa->sip);
  if (msa->nsip) free (msa->nsip);
  if (msa->sequences) free (msa->sequences);
}

struct msa_seq* msa_seq_from_char_vector_string (char *string, size_t nchars, uint32_t id)
{
  struct msa_seq* seq = (struct msa_seq*) biomcmc_malloc (sizeof (struct msa_seq));
  int i;
  seq->id = id;
  seq->len = 0;
  nchars++; // we incorporate null_terminate_sequences() whereby we add final 0

  seq->seq  = (char*) biomcmc_malloc (sizeof(char) * nchars);
  seq->gaps = (int*) biomcmc_malloc (sizeof(int) * (nchars + 1)); // in kalign3 it's one more than the dna (to allow trailing indels?)
  for(i = 0;i < nchars + 1; i++) seq->gaps[i] = 0;
  for(i = 0;i < nchars - 1; i++) { // minus one since we added one 
    if( isalpha((int) string[i])) seq->seq[seq->len++] = string[i];
    if( ispunct((int) string[i])) seq->gaps[seq->len]++; // no increment len++ here
  }
  seq->seq[seq->len] = 0; // null_terminate_sequences() in kalign3
  if (seq->len < nchars - 1) {
    seq->seq  = (char*) biomcmc_realloc ((char*)seq->seq, sizeof(char) * seq->len + 1);
    seq->gaps = (int*)  biomcmc_realloc ((int*)seq->gaps, sizeof(int) * (nchars + 2)); // in kalign3 it's one more than the dna (to allow trailing indels?)
  }
  seq->s  = (uint8_t*) biomcmc_malloc (sizeof(uint8_t) * seq->len + 1);
  return seq;
}

void free_msa_seq(struct msa_seq* seq)
{
  if (!seq) return;
  if (seq->seq)  free (seq->seq);
  if (seq->s)    free (seq->s);
  if (seq->gaps) free (seq->gaps);
  if (seq) free (seq);
}

struct line_buffer* alloc_line_buffer(int max_line_len)
{
  struct line_buffer* lb = NULL;
  int i;
  ASSERT(max_line_len > 60, "max_line_len:%d too small", max_line_len);
  MMALLOC(lb, sizeof(struct line_buffer));
  lb->alloc_num_lines = 1024;

  lb->num_line = 0;
  lb->lines = NULL;
  lb->max_line_len = max_line_len;
  MMALLOC(lb->lines, sizeof(struct out_line*) * lb->alloc_num_lines);

  for(i = 0; i < lb->alloc_num_lines;i++){
    lb->lines[i] = NULL;
    MMALLOC(lb->lines[i], sizeof(struct out_line));
    lb->lines[i]->block = 0;
    lb->lines[i]->seq_id =0;
    lb->lines[i]->line = NULL;
    MMALLOC(lb->lines[i]->line,sizeof(char) * lb->max_line_len);
  }

  return lb;
ERROR:
  return NULL;
}


int resize_line_buffer(struct line_buffer* lb)
{
  int old_len = 0;
  int i;
  old_len = lb->alloc_num_lines;
  lb->alloc_num_lines = lb->alloc_num_lines + 1024;

  MREALLOC(lb->lines, sizeof(struct out_line*) * lb->alloc_num_lines);

  for(i = old_len; i < lb->alloc_num_lines;i++){
    lb->lines[i] = NULL;
    MMALLOC(lb->lines[i], sizeof(struct out_line));
    lb->lines[i]->block = 0;
    lb->lines[i]->seq_id =0;
    lb->lines[i]->line = NULL;
    MMALLOC(lb->lines[i]->line,sizeof(char) * lb->max_line_len);
  }
  return OK;
ERROR:
  return FAIL;
}


void free_line_buffer(struct line_buffer* lb)
{
  int i;

  if(lb){
    for(i = 0; i < lb->alloc_num_lines;i++){
      MFREE(lb->lines[i]->line);
      MFREE(lb->lines[i]);
    }
    MFREE(lb->lines);
    MFREE(lb);
  }
}


int sort_out_lines(const void *a, const void *b)
{
  struct out_line* const *one = a;
  struct out_line* const *two = b;
  if((*one)->block > (*two)->block) {
    return 1;
  } else if((*one)->block == (*two)->block) {
    if((*one)->seq_id > (*two)->seq_id) {
      return 1;
    } else if((*one)->seq_id == (*two)->seq_id) {
      return 0;
    } else{ return -1; }
  } else{ return -1;  }
}

// weave_alignment.c

void make_seq (struct msa* msa,int a,int b,int* path);
int update_gaps (int old_len,int*gis,int new_len,int *newgaps);

int weave(struct msa* msa, int** map, int* tree)
{
  int i, a,b;
  for (i = 0; i < (msa->numseq-1)*3;i +=3){
    a = tree[i];
    b = tree[i+1];
    make_seq(msa,a,b,map[tree[i+2]]);
  }
}

void make_seq (struct msa* msa,int a,int b,int* path)
{
  int *gap_a, *gap_b, i, c, posa = 0, posb = 0;

  gap_a = (int*) biomcmc_malloc ((path[0]+1)*sizeof(int));
  gap_b = (int*) biomcmc_malloc ((path[0]+1)*sizeof(int));

  for (i = path[0]+1;i--;) gap_a[i] = gap_b[i] = 0;
  c = 1;
  while(path[c] != 3) {
    if (!path[c]) {
      posa++;
      posb++;
    } else {
      if (path[c] & 1){
        gap_a[posa] += 1;
        posb++;
      } else if (path[c] & 2){
        gap_b[posb] += 1;
        posa++;
      }
    }
    c++;
  }
  for (i = msa->nsip[a];i--;) update_gaps(msa->sequences[msa->sip[a][i]]->len, msa->sequences[msa->sip[a][i]]->gaps, path[0], gap_a);
  for (i = msa->nsip[b];i--;) update_gaps(msa->sequences[msa->sip[b][i]]->len, msa->sequences[msa->sip[b][i]]->gaps, path[0], gap_b);
  if (gap_a) free (gap_a);
  if (gap_b) free (gap_b);
}

int update_gaps (int old_len,int*gis,int new_len,int *newgaps)
{
  unsigned int i,j;
  int add = 0;
  int rel_pos = 0;
  for (i = 0; i <= old_len;i++){
    add = 0;
    for (j = rel_pos;j <= rel_pos + gis[i];j++) if (newgaps[j] != 0) add += newgaps[j];
    rel_pos += gis[i]+1;
    gis[i] += add;
  }
  return OK;
}
