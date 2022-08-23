/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

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

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
*****************************************************************************************/

#ifndef _FMI_SEARCH_H
#define _FMI_SEARCH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>
#include <limits.h>
#include <fstream>
#include "macro.h"

#include "read_index_ele.h"
#include "bwa.h"

#define DUMMY_CHAR 6

#define assert_not_null(x, size, cur_alloc) \
        if (x == NULL) { fprintf(stderr, "Allocation of %0.2lf GB for " #x " failed.\nCurrent Allocation = %0.2lf GB\n", size * 1.0 /(1024*1024*1024), cur_alloc * 1.0 /(1024*1024*1024)); exit(EXIT_FAILURE); }

#define CP_BLOCK_SIZE 64
#define CP_FILENAME_SUFFIX ".bwt.2bit.64"
#define CP_MASK 63
#define CP_SHIFT 6

typedef struct checkpoint_occ_scalar
{
    int64_t cp_count[4];
    uint64_t one_hot_bwt_str[4];
}CP_OCC;

#if defined(__clang__) || defined(__GNUC__)
static inline int _mm_countbits_64(unsigned long x) {
    return __builtin_popcountl(x);
}
#endif

#define \
GET_OCC(pp, c, occ_id_pp, y_pp, occ_pp, one_hot_bwt_str_c_pp, match_mask_pp) \
                int64_t occ_id_pp = pp >> CP_SHIFT; \
                int64_t y_pp = pp & CP_MASK; \
                int64_t occ_pp = cp_occ[occ_id_pp].cp_count[c]; \
                uint64_t one_hot_bwt_str_c_pp = cp_occ[occ_id_pp].one_hot_bwt_str[c]; \
                uint64_t match_mask_pp = one_hot_bwt_str_c_pp & one_hot_mask_array[y_pp]; \
                occ_pp += _mm_countbits_64(match_mask_pp);

typedef struct smem_struct
{
#ifdef DEBUG
    uint64_t info; // for debug
#endif
    uint32_t rid;
    uint32_t m, n;
    int64_t k, l, s;
}SMEM;

#define SAL_PFD 16

#ifdef SMEM_ACCEL

#define __num_smem_table_entry(len) (1 << ((len) * 2))

/* ALL SMEM TABLE
 * we skip the first bp, since it is easily computed with count[].
 * An entry with 2 cache lines (128-byte) stores 11-bp. => 2^22 entries required.
 * 128-byte * 2^22 = 2^7 * 2^22 = 2^29 = 512MB 
 */
#define ALL_SMEM_MAX_BP 11
#define ALL_SMEM_TABLE_SIZE (__aligned_size((__num_smem_table_entry(ALL_SMEM_MAX_BP) * sizeof(all_smem_t)), 64))
typedef struct __attribute__ (( __packed__)) {
	uint32_t last_avail; /* the last elem with s > 0 */
	struct {
		// Note that this is the backward extension of 3 - a
		uint32_t k32; // k = prev_k + k32
		uint32_t l32; // l = count[3 - b] + l32
		uint32_t s32; // s = s32
	} list[ALL_SMEM_MAX_BP - 1]; /* 12 * 10 = 120-byte */
	uint8_t __pad[__pad_size((sizeof(uint32_t) + (sizeof(uint32_t) * 3)*(ALL_SMEM_MAX_BP - 1)), 64)];
} all_smem_t;

/* LAST SMEM TABLE
 * each element takes 16-byte. With N-bp, 2^(2*N) elements are required.
 * Totally, 2^(4 + 2*N) bytes are required.
 * N=13 => 2^30 = 1GB
 * N=14 => 2^32 = 4GB
 * N=15 => 2^34 = 16GB
 * N=16 => 2^36 = 64GB
 */
typedef struct __attribute__ ((__packed__)) {
	uint8_t bp; /* the number of bp this entry actually indicates. */
	int8_t kms, lms, sms;
	uint32_t kls, lls, sls;
} last_smem_t;
#define LAST_SMEM_MAX_BP 13
#define LAST_SMEM_TABLE_SIZE (__aligned_size((__num_smem_table_entry(LAST_SMEM_MAX_BP) * sizeof(last_smem_t)), 64))

#define __combine_ms_ls(ms, ls) ((((int64_t) (ms)) << 32) | ((int64_t) ls))

int build_smem_tables(char *prefix);

#endif /* SMEM_ACCEL */

#ifdef USE_SHM
size_t ____size_mlt(const char *prefix, const char *ref_file_name);
#endif

class FMI_search: public indexEle
{
    public:
    FMI_search(const char *fname);
    ~FMI_search();
    //int64_t beCalls;
    
    int build_index();
    void load_index();
    void load_index_other_elements(int which);
#ifdef SMEM_ACCEL
	all_smem_t *build_all_smem_table(int len);
	last_smem_t *build_last_smem_table(int len);
#endif

    void getSMEMs(uint8_t *enc_qdb,
                  int32_t numReads,
                  int32_t batch_size,
                  int32_t readlength,
                  int32_t minSeedLengh,
                  int32_t numthreads,
                  SMEM *matchArray,
                  int64_t *numTotalSmem);
    
    void getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                 int16_t *query_pos_array,
                                 int32_t *min_intv_array,
                                 int32_t *rid_array,
                                 int32_t numReads,
                                 int32_t batch_size,
                                 const bseq1_t *seq_,
                                 int32_t *query_cum_len_ar,
                                 int32_t  max_readlength,
                                 int32_t minSeedLen,
                                 SMEM *matchArray,
                                 int64_t *__numTotalSmem);
    
    void getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                 int32_t *min_intv_array,
                                 int32_t *rid_array,
                                 int32_t numReads,
                                 int32_t batch_size,
                                 const bseq1_t *seq_,
                                 int32_t *query_cum_len_ar,
                                 int32_t max_readlength,
                                 int32_t minSeedLen,
                                 SMEM *matchArray,
                                 int64_t *__numTotalSmem);
        
    
    int64_t bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                           int32_t *max_intv_array,
                                           int32_t numReads,
                                           const bseq1_t *seq_,
                                           int32_t *query_cum_len_ar,
                                           int32_t minSeedLen,
                                           SMEM *matchArray);
        
    void sortSMEMs(SMEM *matchArray,
                   int64_t numTotalSmem[],
                   int32_t numReads,
                   int32_t readlength,
                   int nthreads);
    int64_t get_sa_entry(int64_t pos);
    void get_sa_entries(int64_t *posArray,
                        int64_t *coordArray,
                        uint32_t count,
                        int32_t nthreads);
    void get_sa_entries(SMEM *smemArray,
                        int64_t *coordArray,
                        int32_t *coordCountArray,
                        uint32_t count,
                        int32_t max_occ);
    int64_t get_sa_entry_compressed(int64_t pos, int tid=0);
    void get_sa_entries(SMEM *smemArray,
                        int64_t *coordArray,
                        int32_t *coordCountArray,
                        uint32_t count,
                        int32_t max_occ,
                        int tid);
    int64_t call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset);
    void get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                 int64_t *coordCountArray, int64_t count,
                                 const int32_t max_occ, int tid, int64_t &id_);
    
    int64_t reference_seq_len;
    int64_t sentinel_index;
#ifdef PERFECT_MATCH
	perfect_table_t *perfect_table;
#endif
 	
	int useErt;
    uint64_t         *kmer_offsets;
    uint8_t          *mlt_table;
	void load_ert_index();

private:
        char file_name[PATH_MAX];
        int64_t index_alloc;
        int64_t count[5];
        uint32_t *sa_ls_word;
        int8_t *sa_ms_byte;
        CP_OCC *cp_occ;


#ifdef SMEM_ACCEL
		all_smem_t *all_smem_table;
		last_smem_t *last_smem_table;
#endif
        uint64_t *one_hot_mask_array;
   

		void init(const char *fname);

        int64_t pac_seq_len(const char *fn_pac);
        void pac2nt(const char *fn_pac,
                    std::string &reference_seq);
        int build_fm_index(const char *ref_file_name,
                               char *binary_seq,
                               int64_t ref_seq_len,
                               int64_t *sa_bwt,
                               int64_t *count);
        SMEM backwardExt(SMEM smem, uint8_t a);

#ifdef SMEM_ACCEL
		void __build_all_smem_table(uint8_t *seq, int len, all_smem_t *ent);
		void __build_last_smem_table(uint8_t *seq, int len, last_smem_t *ent);
		void load_smem_table();
#endif /* SMEM_ACCEL */
};

#endif
