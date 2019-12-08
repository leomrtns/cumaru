AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  -I$(srcdir)/../@BIOMCMCLIB@/lib  @OPENMP_CPPFLAGS@ @ZLIB_LIBS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @ZLIB_LIBS@ ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)

common_headers = 
tldevel.h \
rng.h \
global.h \
msa.h \
parameters.h \
aln_data.h \
alphabet.h \
alignment_parameters.h \
pick_anchor.h \
bisectingKmeans.h \
sequence_distance.h \
euclidean_dist.h \
alignment.h \
weave_alignment.h \
misc.h \
bpm.h

common_src = run_kalign.c \
tldevel.c \
rng.c \
parameters.c \
aln_data.c \
alphabet.c \
alignment_parameters.c \
pick_anchor.c \
bisectingKmeans.c \
sequence_distance.c \
euclidean_dist.c \
alignment.c \
weave_alignment.c \
rwalign.c \
misc.c \
bpm.c

noinst_LTLIBRARIES = libkalign.la   ## change to noinst_LT... if should be linked statically (not installed globally)
libkalign_la_SOURCES = config.h $(common_headers) $(common_src) 

