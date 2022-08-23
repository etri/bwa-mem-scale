##/*************************************************************************************
##                           The MIT License
##
##   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
##   Copyright (C) 2019  Intel Corporation, Heng Li.
##
##   Permission is hereby granted, free of charge, to any person obtaining
##   a copy of this software and associated documentation files (the
##   "Software"), to deal in the Software without restriction, including
##   without limitation the rights to use, copy, modify, merge, publish,
##   distribute, sublicense, and/or sell copies of the Software, and to
##   permit persons to whom the Software is furnished to do so, subject to
##   the following conditions:
##
##   The above copyright notice and this permission notice shall be
##   included in all copies or substantial portions of the Software.
##
##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
##   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
##   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
##   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
##   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
##   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
##   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##   SOFTWARE.
##
##Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
##                                Heng Li <hli@jimmy.harvard.edu> 
##*****************************************************************************************/

ifneq ($(portable),)
	STATIC_GCC=-static-libgcc -static-libstdc++
	STATIC_SAFESTR=1
else
	STATIC_GCC=
	STATIC_SAFESTR=0
endif

#CXX=		icpc
ifeq ($(CXX), icpc)
	CC= icc
else ifeq ($(CXX), g++)
	CC=gcc
endif		
ARCH_FLAGS=	-msse4.1
MEM_FLAGS=	-DSAIS=1
CPPFLAGS+=	-DENABLE_PREFETCH -DV17=1 $(MEM_FLAGS)
INCLUDES=   -Isrc -Iext/safestringlib/include
LIBS=		-lpthread -lm -lz -L. -lbwa -Lext/safestringlib -lsafestring $(STATIC_GCC)
OBJS=		src/fastmap.o src/main.o src/utils.o src/memcpy_bwamem.o src/kthread.o \
			src/kstring.o src/ksw.o src/bwt.o src/ertindex.o src/bntseq.o src/bwamem.o src/ertseeding.o src/profiling.o src/bandedSWA.o \
			src/FMI_search.o src/read_index_ele.o src/bwamem_pair.o src/kswv.o src/bwa.o \
			src/bwamem_extra.o src/bwtbuild.o src/QSufSort.o src/bwt_gen.o src/rope.o src/rle.o src/is.o src/kopen.o src/bwtindex.o \
			src/perfect_index.o src/perfect_map.o src/bwa_shm.o
BWA_LIB=    libbwa.a
SAFE_STR_LIB=    ext/safestringlib/libsafestring.a

EXE_BASE=bwa-mem2
DEPEND_CPPFLAGS=

ifeq ($(memscale),1)
CPPFLAGS+= -DMEMSCALE
shm=1
perfect=1
accel=1
ert=1
rwopt=1
else
DEPEND_CPPFLAGS+= -DMEMSCALE
endif

ifeq ($(shm),1)
CPPFLAGS+= -DUSE_SHM
LIBS+= -lrt
EXE_SHM=$(addsuffix .shm,$(EXE_BASE))
else
DEPEND_CPPFLAGS+= -DUSE_SHM
EXE_SHM=$(EXE_BASE)
endif

ifeq ($(perfect),1)
CPPFLAGS+= -DPERFECT_MATCH
EXE_PERFECT=$(addsuffix .perfect,$(EXE_SHM))
else
DEPEND_CPPFLAGS+= -DPERFECT_MATCH
EXE_PERFECT=$(EXE_SHM)
endif

ifeq ($(accel),1)
CPPFLAGS+= -DSMEM_ACCEL
EXE_ACCEL=$(addsuffix .accel,$(EXE_PERFECT))
else
DEPEND_CPPFLAGS+= -DSMEM_ACCEL
EXE_ACCEL=$(EXE_PERFECT)
endif

ifeq ($(ert),1)
CPPFLAGS+= -DDEFAULT_USE_ERT=1
EXE_ERT=$(addsuffix .ert,$(EXE_ACCEL))
else
DEPEND_CPPFLAGS+= -DDEFAULT_USE_ERT=1
EXE_ERT=$(EXE_ACCEL)
endif

ifeq ($(rwopt),1)
EXE_RWOPT=$(addsuffix .rwopt,$(EXE_ERT))
CPPFLAGS+= -DOPT_RW
else
EXE_RWOPT=$(EXE_ERT)
DEPEND_CPPFLAGS+= -DOPT_RW
endif

ifeq ($(memscale),1)
EXE_MEMSCALE=$(addsuffix .memscale,$(EXE_BASE))
else
EXE_MEMSCALE=$(EXE_RWOPT)
endif

ifeq ($(aff),1)
CPPFLAGS+= -DAFF=1
EXE_AFF=$(addsuffix .aff,$(EXE_MEMSCALE))
else
DEPEND_CPPFLAGS+= -DAFF=1
EXE_AFF=$(EXE_MEMSCALE)
endif

EXE_NOARCH=$(EXE_AFF)

ifdef batch_size
CPPFLAGS+= -DCONFIG_BATCH_SIZE=$(batch_size)
endif

ifeq ($(arch),sse)
	ARCH_FLAGS=-msse4.1
	EXE=$(addsuffix .sse41,$(EXE_NOARCH))
else ifeq ($(arch),avx2)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-march=core-avx2 #-xCORE-AVX2
	else	
		ARCH_FLAGS=-mavx2
	endif
	EXE=$(addsuffix .avx2,$(EXE_NOARCH))
else ifeq ($(arch),avx512)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512
	else	
		ARCH_FLAGS=-mavx512bw
	endif
	EXE=$(addsuffix .avx512bw,$(EXE_NOARCH))
else ifeq ($(arch),native)
	ARCH_FLAGS=-march=native
	EXE=$(EXE_NOARCH)
else ifneq ($(arch),)
# To provide a different architecture flag like -march=core-avx2.
	ARCH_FLAGS=$(arch)
	EXE=$(EXE_NOARCH)
else
EXE=$(EXE_NOARCH)

myall:multi
endif

CXXFLAGS+=	-g -O3 -fpermissive $(ARCH_FLAGS) #-Wall ##-xSSE2
#CXXFLAGS+=	-g -O0 -fpermissive $(ARCH_FLAGS) #-Wall ##-xSSE2

.PHONY:all clean depend multi
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:exec

multi:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse  CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx2   CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx512 CXX=$(CXX) all
	$(CXX) -Wall -O3 src/runsimd.cpp -Iext/safestringlib/include -Lext/safestringlib/ -lsafestring $(STATIC_GCC) -o $(EXE)


exec:$(BWA_LIB) $(SAFE_STR_LIB) src/main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/main.o $(BWA_LIB) $(LIBS) -o $(EXE)

$(BWA_LIB):$(OBJS)
	ar rcs $(BWA_LIB) $(OBJS)

$(SAFE_STR_LIB):
	cd ext/safestringlib/ && $(MAKE) clean && $(MAKE) CC=$(CC) portable=$(STATIC_SAFESTR) directories libsafestring.a

clean_bwa:
	rm -fr src/*.o 
	rm -fr $(BWA_LIB) 

clean:
	rm -fr src/*.o 
	#rm -fr $(BWA_LIB) $(EXE) $(EXE).sse41 $(EXE).avx2 $(EXE).avx512bw
	rm -fr $(BWA_LIB) $(EXE)
	cd ext/safestringlib/ && $(MAKE) clean

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) $(DEPEND_CPPFLAGS) -I. -- src/*.cpp)

# DO NOT DELETE

src/FMI_search.o: src/sais.h src/FMI_search.h src/read_index_ele.h
src/FMI_search.o: src/utils.h src/bntseq.h src/macro.h src/bwa.h src/bwt.h
src/FMI_search.o: src/perfect.h src/memcpy_bwamem.h src/profiling.h
src/FMI_search.o: src/bwa_shm.h
src/bandedSWA.o: src/bandedSWA.h src/macro.h
src/bntseq.o: src/bntseq.h src/utils.h src/macro.h src/kseq.h
src/bntseq.o: src/memcpy_bwamem.h src/khash.h
src/bwa.o: src/bntseq.h src/bwa.h src/bwt.h src/macro.h src/perfect.h
src/bwa.o: src/ksw.h src/utils.h src/kstring.h src/memcpy_bwamem.h src/kvec.h
src/bwa.o: src/kseq.h
src/bwa_shm.o: src/bwa_shm.h src/perfect.h src/FMI_search.h
src/bwa_shm.o: src/read_index_ele.h src/utils.h src/bntseq.h src/macro.h
src/bwa_shm.o: src/bwa.h src/bwt.h src/fastmap.h src/bwamem.h src/kthread.h
src/bwa_shm.o: src/bandedSWA.h src/kstring.h src/memcpy_bwamem.h src/ksw.h
src/bwa_shm.o: src/kvec.h src/ksort.h src/profiling.h src/kseq.h
src/bwamem.o: src/bwamem.h src/bwt.h src/bntseq.h src/bwa.h src/macro.h
src/bwamem.o: src/perfect.h src/kthread.h src/bandedSWA.h src/kstring.h
src/bwamem.o: src/memcpy_bwamem.h src/ksw.h src/kvec.h src/ksort.h
src/bwamem.o: src/utils.h src/profiling.h src/FMI_search.h
src/bwamem.o: src/read_index_ele.h src/kbtree.h
src/bwamem_extra.o: src/bwa.h src/bntseq.h src/bwt.h src/macro.h
src/bwamem_extra.o: src/perfect.h src/bwamem.h src/kthread.h src/bandedSWA.h
src/bwamem_extra.o: src/kstring.h src/memcpy_bwamem.h src/ksw.h src/kvec.h
src/bwamem_extra.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/bwamem_extra.o: src/read_index_ele.h
src/bwamem_pair.o: src/kstring.h src/memcpy_bwamem.h src/bwamem.h src/bwt.h
src/bwamem_pair.o: src/bntseq.h src/bwa.h src/macro.h src/perfect.h
src/bwamem_pair.o: src/kthread.h src/bandedSWA.h src/ksw.h src/kvec.h
src/bwamem_pair.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/bwamem_pair.o: src/read_index_ele.h src/kswv.h
src/bwt.o: src/utils.h src/bwt.h src/kvec.h src/malloc_wrap.h
src/bwt_gen.o: src/QSufSort.h src/malloc_wrap.h
src/bwtbuild.o: src/sais.h src/utils.h src/bntseq.h
src/bwtindex.o: src/bntseq.h src/bwa.h src/bwt.h src/macro.h src/perfect.h
src/bwtindex.o: src/utils.h src/rle.h src/rope.h src/malloc_wrap.h
src/bwtindex.o: src/FMI_search.h src/read_index_ele.h src/bwtbuild.h
src/fastmap.o: src/fastmap.h src/bwa.h src/bntseq.h src/bwt.h src/macro.h
src/fastmap.o: src/perfect.h src/bwamem.h src/kthread.h src/bandedSWA.h
src/fastmap.o: src/kstring.h src/memcpy_bwamem.h src/ksw.h src/kvec.h
src/fastmap.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/fastmap.o: src/read_index_ele.h src/kseq.h src/bwa_shm.h
src/kopen.o: src/memcpy_bwamem.h
src/kstring.o: src/kstring.h src/memcpy_bwamem.h
src/ksw.o: src/ksw.h src/macro.h
src/kswv.o: src/kswv.h src/macro.h src/ksw.h src/bandedSWA.h
src/kthread.o: src/kthread.h src/macro.h src/bwamem.h src/bwt.h src/bntseq.h
src/kthread.o: src/bwa.h src/perfect.h src/bandedSWA.h src/kstring.h
src/kthread.o: src/memcpy_bwamem.h src/ksw.h src/kvec.h src/ksort.h
src/kthread.o: src/utils.h src/profiling.h src/FMI_search.h
src/kthread.o: src/read_index_ele.h
src/main.o: src/main.h src/kstring.h src/memcpy_bwamem.h src/utils.h
src/main.o: src/macro.h src/bandedSWA.h src/profiling.h src/fastmap.h
src/main.o: src/bwa.h src/bntseq.h src/bwt.h src/perfect.h src/bwamem.h
src/main.o: src/kthread.h src/ksw.h src/kvec.h src/ksort.h src/FMI_search.h
src/main.o: src/read_index_ele.h src/kseq.h
src/malloc_wrap.o: src/malloc_wrap.h
src/memcpy_bwamem.o: src/memcpy_bwamem.h
src/perfect_index.o: src/bwa.h src/bntseq.h src/bwt.h src/macro.h
src/perfect_index.o: src/perfect.h src/utils.h src/fastmap.h src/bwamem.h
src/perfect_index.o: src/kthread.h src/bandedSWA.h src/kstring.h
src/perfect_index.o: src/memcpy_bwamem.h src/ksw.h src/kvec.h src/ksort.h
src/perfect_index.o: src/profiling.h src/FMI_search.h src/read_index_ele.h
src/perfect_index.o: src/kseq.h
src/perfect_map.o: src/bntseq.h src/bwa.h src/bwt.h src/macro.h src/perfect.h
src/perfect_map.o: src/ksw.h src/utils.h src/kstring.h src/memcpy_bwamem.h
src/perfect_map.o: src/kvec.h src/bwa_shm.h src/kseq.h
src/profiling.o: src/macro.h src/profiling.h
src/read_index_ele.o: src/read_index_ele.h src/utils.h src/bntseq.h
src/read_index_ele.o: src/macro.h src/bwa_shm.h src/perfect.h
src/utils.o: src/utils.h src/ksort.h src/kseq.h src/memcpy_bwamem.h
src/rle.o: src/rle.h
src/rope.o: src/rle.h src/rope.h
src/is.o: src/malloc_wrap.h
src/QSufSort.o: src/QSufSort.h
src/ertindex.o: src/ertindex.h src/bwt.h src/kvec.h src/macro.h
src/ertseeding.o: src/ertseeding.h src/bwamem.h src/bwt.h src/bntseq.h src/bwa.h src/macro.h 
src/ertseeding.o: src/kthread.h src/bandedSWA.h src/kstring.h src/ksw.h
src/ertseeding.o: src/kvec.h src/ksort.h src/utils.h src/profiling.h
src/ertseeding.o: src/FMI_search.h src/read_index_ele.h src/kbtree.h
