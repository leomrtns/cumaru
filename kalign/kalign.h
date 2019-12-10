/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This file is part of cumaru, a multiple sequence alignment.
 */

/*! \file
 *  \brief library based on [kalign3 multiple sequence alignment program](https://github.com/TimoLassmann/kalign.git)
 */

#ifndef _cumaru_kalign_h
#define _cumaru_kalign_h

#include "alignment.h"
#include "alignment_parameters.h"
#include "bisectingKmeans.h"
#include "bpm.h"
#include "global.h"
#include "kmeans.h"
#include "misc.h"
#include "rng.h"
#include "sequence_distance.h"
#include "tldevel.h"

char_vector kalign3_from_char_vector (char_vector dna);
#endif
