/*************************************************************************************
                           The MIT License

   BWA-MEM-SCALE (Memory-Scalable Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2022 Electronics and Telecommunications Research Institute (ETRI), Changdae Kim.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   Contacts: Changdae Kim <cdkim@etri.re.kr>

** This software builds upon BWA-MEM2, and includes several performance optimization techniques.
   For BWA-MEM2, refer to the follows.

   BWA-MEM2 (Sequence alignment using Burrows-Wheeler Transform)
   Copyright ⓒ 2019 Intel Corporation, Heng Li
   The MIT License
   Website: https://github.com/bwa-mem2/bwa-mem2

*****************************************************************************************/

// ----------------------------------
#include "main.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "2.0"
#endif


// ----------------------------------
uint64_t proc_freq, tprof[LIM_R][LIM_C], prof[LIM_R];
// ----------------------------------

int usage()
{
    fprintf(stderr, "Usage: bwa-mem2.scale <command> <arguments>\n");
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "  index         create index\n");
    fprintf(stderr, "  perfect-index create index for perfect match\n");
    fprintf(stderr, "  smem-table    create index for FM-index accelerator\n");
    fprintf(stderr, "  mem           alignment\n");
    fprintf(stderr, "  load-shm      load index on process shared memory\n");
    fprintf(stderr, "  remove-shm    remove index from process shared memory\n");
    fprintf(stderr, "  version       print version number\n");
    return 1;
}

int main(int argc, char* argv[])
{
        
    // ---------------------------------    
    uint64_t tim = __rdtsc();
    sleep(1);
    proc_freq = __rdtsc() - tim;

    int ret = -1;
    if (argc < 2) return usage();

    if (strcmp(argv[1], "index") == 0)
    {
         uint64_t tim = __rdtsc();
         ret = bwa_index(argc-1, argv+1);
         fprintf(stderr, "Total time taken: %0.4lf\n", (__rdtsc() - tim)*1.0/proc_freq);
         return ret;
    }
    else if (strcmp(argv[1], "mem") == 0)
    {
        tprof[MEM][0] = __rdtsc();
        kstring_t pg = {0,0,0};
        extern char *bwa_pg;

        fprintf(stderr, "-----------------------------\n");
#if __AVX512BW__
        fprintf(stderr, "Executing in AVX512 mode!!\n");
#endif
#if ((!__AVX512BW__) & (__AVX2__))
        fprintf(stderr, "Executing in AVX2 mode!!\n");
#endif
#if ((!__AVX512BW__) && (!__AVX2__) && (__SSE2__))
        fprintf(stderr, "Executing in SSE4.1 mode!!\n");
#endif
#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
        fprintf(stderr, "Executing in Scalar mode!!\n");
#endif
        fprintf(stderr, "-----------------------------\n");

        #if SA_COMPRESSION
        fprintf(stderr, "SA compression enable with xfactor (2^): %d !!!\n", SA_COMPX);
        #endif
        
        ksprintf(&pg, "@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);

        for (int i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
        ksprintf(&pg, "\n");
        bwa_pg = pg.s;
        ret = main_mem(argc-1, argv+1);
        free(bwa_pg);
        
        /** Enable this return to avoid printing of the runtime profiling **/
        //return ret;
    }
#ifdef PERFECT_MATCH
	else if (strcmp(argv[1], "perfect-index") == 0)
	{
		// build perfect table
		uint64_t tim = __rdtsc();
		ret = perfect_index(argc-1, argv+1);
        fprintf(stderr, "Total time taken: %0.4lf\n", (__rdtsc() - tim)*1.0/proc_freq);
        return ret;
	}
#endif
#ifdef SMEM_ACCEL
	else if (strcmp(argv[1], "smem-table") == 0)
	{
		// build two tables for FM-index walking
		if (argc < 3) {
			printf("usage: %s accel <idxbase>\n"
				   "       build two smem tables for FM-index walking acceleration.\n",
				   argv[0]);
			return -1;
		}
		uint64_t tim = __rdtsc();
		if (build_smem_tables(argv[2]) != 0) 
			fprintf(stderr, "Failed to build tables for smem acceleration\n");
        fprintf(stderr, "Total time taken: %0.4lf\n", (__rdtsc() - tim)*1.0/proc_freq);
        return 0;
	}
#endif
#ifdef USE_SHM
	else if (strcmp(argv[1], "load-shm") == 0) {
		uint64_t tim = __rdtsc();
		int ret = bwa_shm_load(argc-1, argv+1);
        fprintf(stderr, "Total time taken: %0.4lf\n", (__rdtsc() - tim)*1.0/proc_freq);
		return ret;
	}
	else if (strcmp(argv[1], "remove-shm") == 0)
	{
		uint64_t tim = __rdtsc();
		int ret = bwa_shm_remove();
        fprintf(stderr, "Total time taken: %0.4lf\n", (__rdtsc() - tim)*1.0/proc_freq);
		return ret;
	}
#endif
    else if (strcmp(argv[1], "version") == 0)
    {
        puts(PACKAGE_VERSION);
        return 0;
    } else {
        fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
        return 1;
    }
        
    fprintf(stderr, "\nImportant parameter settings: \n");
    fprintf(stderr, "\tBATCH_SIZE: %d\n", BATCH_SIZE);
    fprintf(stderr, "\tMAX_SEQ_LEN_REF: %d\n", MAX_SEQ_LEN_REF);
    fprintf(stderr, "\tMAX_SEQ_LEN_QER: %d\n", MAX_SEQ_LEN_QER);
    fprintf(stderr, "\tMAX_SEQ_LEN8: %d\n", MAX_SEQ_LEN8);
    fprintf(stderr, "\tSEEDS_PER_READ: %d\n", SEEDS_PER_READ);
    fprintf(stderr, "\tSIMD_WIDTH8 X: %d\n", SIMD_WIDTH8);
    fprintf(stderr, "\tSIMD_WIDTH16 X: %d\n", SIMD_WIDTH16);
    fprintf(stderr, "\tAVG_SEEDS_PER_READ: %d\n", AVG_SEEDS_PER_READ);
    
    return 0;
}
