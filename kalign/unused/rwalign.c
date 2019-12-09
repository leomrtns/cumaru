/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

#include "rwalign.h"

/* only local; */
struct line_buffer{
  struct out_line** lines;
  int max_line_len;
  int alloc_num_lines;
  int num_line;
};

struct out_line{
  char* line;
  int block;
  int seq_id;
};

struct msa* read_fasta(char* infile, struct msa* msa);

int write_msa_fasta(struct msa* msa,char* outfile);

/* memory functions  */
struct msa* alloc_msa(void);
int resize_msa(struct msa* msa);

struct msa_seq* alloc_msa_seq(void);
int resize_msa_seq(struct msa_seq* seq);
void free_msa_seq(struct msa_seq* seq);

struct line_buffer* alloc_line_buffer(int max_line_len);
int resize_line_buffer(struct line_buffer* lb);
void free_line_buffer(struct line_buffer* lb);

/* local helper functions  */
static int set_sip_nsip(struct msa* msa);
static int null_terminate_sequences(struct msa* msa);
static int sort_out_lines(const void *a, const void *b);
static int make_linear_sequence(struct msa_seq* seq, char* linear_seq);
int GCGMultchecksum(struct msa* msa);
/* Taken from squid library by Sean Eddy  */
int GCGchecksum(char *seq, int len);

struct msa* read_input(char* infile, struct msa* msa)
{
  int type;
  ASSERT(infile != NULL,"No input file");
  /* sanity checks  */
  if(!my_file_exists(infile)){
    ERROR_MSG("File: %s does not exist.",infile);
  }

  RUNP(msa = read_fasta(infile,msa));
  RUN(convert_msa_to_internal(msa));
  msa->aligned = 0;
  RUN(set_sip_nsip(msa));

  return msa;
ERROR:
  return NULL;
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

int convert_msa_to_internal(struct msa* msa)
{
  struct alphabet* a = NULL;
  struct msa_seq* seq = NULL;
  int8_t* t = NULL;
  int i,j;

  RUNP(a = create_dna_alphabet());

  t = a->to_internal;
  msa->L = a->L;
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
  MFREE(a);
  return OK;
ERROR:
  if(a){
    MFREE(a);
  }
  return FAIL;
}

/* rw functions; I wand fasta, msf and clustal */

int write_msa_fasta(struct msa* msa,char* outfile)
{
  FILE* f_ptr = NULL;
  int i,j,c,f;

  if(!outfile){
    f_ptr = stdout;
  }else{
    RUNP(f_ptr = fopen(outfile, "w"));
  }

  for(i = 0; i < msa->numseq;i++){
    fprintf(f_ptr,">%s\n", msa->sequences[i]->name);
    f = 0;
    for(j = 0;j < msa->sequences[i]->len;j++){

      for(c = 0;c < msa->sequences[i]->gaps[j];c++){
        fprintf(f_ptr,"-");
        f++;
        if(f == 60){
          fprintf(f_ptr, "\n");
          f = 0;
        }
      }
      fprintf(f_ptr,"%c", msa->sequences[i]->seq[j]);
      f++;
      if(f == 60){
        fprintf(f_ptr, "\n");
        f = 0;
      }
    }
    for(c = 0;c < msa->sequences[i]->gaps[ msa->sequences[i]->len];c++){
      fprintf(f_ptr,"-");
      f++;
      if(f == 60){
        fprintf(f_ptr, "\n");
        f = 0;
      }
    }
    if(f){
      fprintf(f_ptr,"\n");
    }
  }
  if(outfile){
    fclose(f_ptr);
  }

  return OK;
ERROR:
  return FAIL;
}

struct msa* read_fasta(char* infile,struct msa* msa)
{
  struct msa_seq* seq_ptr = NULL;
  FILE* f_ptr = NULL;
  char line[BUFFER_LEN];
  int line_len;
  int i;

  /* sanity checks  */
  if(!my_file_exists(infile)){
    ERROR_MSG("File: %s does not exist.",infile);
  }
  if(msa == NULL){
    msa = alloc_msa();
  }
  RUNP(f_ptr = fopen(infile, "r"));

  while(fgets(line, BUFFER_LEN, f_ptr)){
    line_len = strnlen(line, BUFFER_LEN);
    if(line[0] == '>'){
      /* alloc seq if buffer is full */
      if(msa->alloc_numseq == msa->numseq){
        RUN(resize_msa(msa));
      }

      line[line_len-1] = 0;
      for(i =0 ; i < line_len;i++){
        if(isspace(line[i])){
          line[i] = 0;
        }
      }
      seq_ptr = msa->sequences[msa->numseq];
      snprintf(seq_ptr->name ,MSA_NAME_LEN ,"%s",line+1);
      msa->numseq++;
    }else{
      for(i = 0;i < line_len;i++){
        msa->letter_freq[(int)line[i]]++;
        if(isalpha((int)line[i])){

          if(seq_ptr->alloc_len == seq_ptr->len){
            resize_msa_seq(seq_ptr);
          }

          seq_ptr->seq[seq_ptr->len] = line[i];
          seq_ptr->len++;
        }
        if(ispunct((int)line[i])){
          seq_ptr->gaps[seq_ptr->len]++;
        }
      }
    }
  }
  RUN(null_terminate_sequences(msa));
  fclose(f_ptr);
  return msa;
ERROR:
  free_msa(msa);
  return NULL;
}

int null_terminate_sequences(struct msa* msa)
{
  int i;
  struct msa_seq* seq_ptr = NULL;
  /* 0 terminate sequences  */
  for(i = 0; i < msa->numseq;i++){
    seq_ptr = msa->sequences[i];
    if(seq_ptr->alloc_len == seq_ptr->len){
      resize_msa_seq(seq_ptr);
    }
    seq_ptr->seq[seq_ptr->len] = 0;

  }
  return OK;
}

int make_linear_sequence(struct msa_seq* seq, char* linear_seq)
{
  int c,j,f;
  f = 0;
  for(j = 0;j < seq->len;j++){
    //LOG_MSG("%d %d",j,seq->gaps[j]);
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
  //fprintf(stdout,"LINEAR:%s\n",linear_seq);
  return OK;
}

/* memory functions */
struct msa* alloc_msa(void)
{
  struct msa* msa = NULL;
  int i;
  MMALLOC(msa, sizeof(struct msa));
  msa->sequences = NULL;
  msa->alloc_numseq = 512;
  msa->numseq = 0;
  msa->L = 0;
  msa->aligned = 0;
  msa->plen = NULL;
  msa->sip = NULL;
  msa->nsip = NULL;

  MMALLOC(msa->sequences, sizeof(struct msa_seq*) * msa->alloc_numseq);

  for(i = 0; i < msa->alloc_numseq;i++){
    msa->sequences[i] = NULL;
    RUNP(msa->sequences[i] = alloc_msa_seq());
  }
  for(i = 0; i < 128; i++) msa->letter_freq[i] = 0;
  return msa;
ERROR:
  free_msa(msa);
  return NULL;
}

int resize_msa(struct msa* msa)
{
  int i;
  int old_size;

  old_size = msa->alloc_numseq;
  msa->alloc_numseq = msa->alloc_numseq + 512;

  MREALLOC(msa->sequences, sizeof(struct msa_seq*) * msa->alloc_numseq);

  for(i = old_size; i < msa->alloc_numseq;i++){
    msa->sequences[i] = NULL;
    RUNP(msa->sequences[i] = alloc_msa_seq());
  }
  return OK;
ERROR:
  return FAIL;
}

void free_msa(struct msa* msa)
{
  int i;
  if(msa){
    for(i = 0; i < msa->alloc_numseq;i++){
      free_msa_seq(msa->sequences[i]);
    }

    for (i = msa->num_profiles;i--;){
      if(msa->sip[i]){
        MFREE(msa->sip[i]);
      }
    }
    MFREE(msa->plen);
    MFREE(msa->sip);
    MFREE(msa->nsip);

    MFREE(msa->sequences);
    MFREE(msa);
  }
}

struct msa_seq* alloc_msa_seq(void)
{
  struct msa_seq* seq = NULL;
  int i;
  MMALLOC(seq, sizeof(struct msa_seq));
  seq->name = NULL;
  seq->seq = NULL;
  seq->s = NULL;
  seq->gaps = NULL;
  seq->len = 0;
  seq->alloc_len = 512;

  MMALLOC(seq->name, sizeof(char)* MSA_NAME_LEN);

  MMALLOC(seq->seq, sizeof(char) * seq->alloc_len);
  MMALLOC(seq->s, sizeof(uint8_t) * seq->alloc_len);
  MMALLOC(seq->gaps, sizeof(int) * (seq->alloc_len+1));
  for(i =0;i < seq->alloc_len+1;i++){
    seq->gaps[i] = 0;
  }
  return seq;
ERROR:
  free_msa_seq(seq);
  return NULL;
}


int resize_msa_seq(struct msa_seq* seq)
{
  int old_len;
  int i;
  old_len = seq->alloc_len;
  seq->alloc_len = seq->alloc_len + 512;

  MREALLOC(seq->seq, sizeof(char) * seq->alloc_len);
  MREALLOC(seq->s, sizeof(uint8_t) * seq->alloc_len);
  MREALLOC(seq->gaps, sizeof(int) * (seq->alloc_len+1));

  for(i = old_len;i < seq->alloc_len+1;i++){
    seq->gaps[i] = 0;
  }

  return OK;
ERROR:
  return FAIL;
}

void free_msa_seq(struct msa_seq* seq)
{
  if(seq){
    MFREE(seq->name);
    MFREE(seq->seq);
    MFREE(seq->s);
    MFREE(seq->gaps);
    MFREE(seq);
  }
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

  if((*one)->block > (*two)->block){
    return 1;
  }else if((*one)->block == (*two)->block){
    if((*one)->seq_id > (*two)->seq_id){
      return 1;
    }else if((*one)->seq_id == (*two)->seq_id){
      return 0;
    }else{
      return -1;
    }
  }else{
    return -1;
  }
}

/* alignment checksum  */
int GCGMultchecksum(struct msa* msa)
{
  int chk = 0;
  for (int idx = 0; idx < msa->numseq; idx++) chk = (chk + GCGchecksum(msa->sequences[idx]->seq,  msa->sequences[idx]->len)) % 10000;
  return chk;
}

/* Taken from squid library by Sean Eddy  */
int GCGchecksum(char *seq, int len)
{
  int chk = 0;			/* calculated checksum  */
  for (int i = 0; i < len; i++) chk = (chk + (i % 57 + 1) * (toupper((int) seq[i]))) % 10000;
  return chk;
}
