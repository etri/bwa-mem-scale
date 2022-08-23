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

#ifndef BWA_SHM_HPP
#define BWA_SHM_HPP
#include <stdio.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <linux/mman.h>
#include "macro.h"
#include "bwamem.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "safe_mem_lib.h"
#ifdef __cplusplus
}
#endif

#ifdef PERFECT_MATCH
#include "perfect.h"
#endif
#include <errno.h>

#ifdef SMEM_ACCEL
#include "FMI_search.h"
#endif

#ifdef USE_SHM

#ifdef PERFECT_MATCH
#define DEFAULT_MMAP_PERFECT 0
#endif

enum bwa_shm_type {
	BWA_SHM_INFO = 0,
	BWA_SHM_BWT,
	BWA_SHM_REF,
	BWA_SHM_PAC,
	BWA_SHM_KMER, /* for ERT */
	BWA_SHM_MLT, /* for ERT */
#ifdef PERFECT_MATCH
	BWA_SHM_PERFECT,
#endif
#ifdef SMEM_ACCEL
	BWA_SHM_SALL, /* SMEM ALL */
	BWA_SHM_SLAST, /* SMEM LAST */
#endif
	NUM_BWA_SHM,
};

enum bwa_shm_mode {
	BWA_SHM_MATCHED = 0,
	BWA_SHM_DISABLE = 1,
	BWA_SHM_RENEWAL = 2
};

enum hugetlb_mode {
	BWA_SHM_NORMAL_PAGE = 0,
	BWA_SHM_HUGE_PAGE = 1,
	BWA_SHM_HUGE_2MB = 2,
	BWA_SHM_HUGE_1GB = 3,
};

extern enum bwa_shm_mode bwa_shm_mode;
extern int shm_fd[NUM_BWA_SHM];
extern void *shm_ptr[NUM_BWA_SHM];

#define BWA_SHM_STATE_NOT_INIT 0
#define BWA_SHM_STATE_MODIFY 1
#define BWA_SHM_STATE_WAIT 2
#define BWA_SHM_STATE_AVAIL 3
typedef struct {
	int lock;
	int state; /* 0: not initialized.
				  1: index manager is modifying the index data.
				  2: index manager waits for finishing current jobs.
	              3: available for mapping jobs
				*/
	int num_map_read; /* the number of mem processes using the process shared memory. */ 
	int num_map_manager; /* the number of manager processes using the process shared memory. */ 

	int hugetlb_flags;
	int useErt;

#ifdef MEMSCALE
	int bwt_on;
	int pac_on;
	int ref_on;
	int kmer_on;
	int mlt_on;
	int perfect_on;
	int smem_all_on;
	int smem_last_on;

	uint32_t pt_num_seed_entry_loaded;
#endif
#ifdef PERFECT_MATCH
	uint32_t pt_num_loc_entry;
	uint32_t pt_num_seed_entry; 
	int pt_seed_len;
	int pt_mmap;
#endif

	/* to distinguish the loaded index */
	int64_t reference_len; /* size(in bytes) + 1 of prefix.0123 file */
	struct timespec mtim_ref; /* last modification time of prefix.0123 file */
	int ref_file_name_len;
	char ref_file_name[0]; /* absolute path of bwt file */
} bwa_shm_info_t;

extern bwa_shm_info_t *bwa_shm_info;
extern bwa_shm_info_t *loading_info; /* for bwa_shm_load */

#define bwa_shm_rlen() (bwa_shm_info ? bwa_shm_info->reference_len : 0)
static inline int bwa_shm_hugetlb_flags() {
	if (loading_info)
		return loading_info->hugetlb_flags;
	else if (bwa_shm_info)
		return bwa_shm_info->hugetlb_flags;
	else
		return 0;
}
#define use_hugetlb(m) ((m) != BWA_SHM_INFO && bwa_shm_hugetlb_flags() != 0)

#define bwa_shm_size_info(abs_path_len) (sizeof(bwa_shm_info_t) + (abs_path_len) + 1)

typedef struct {
	int64_t reference_len;
	int64_t count[5];
	int64_t sentinel_index;
} shm_bwt_header_t;

#define bwa_shm_size_bwt_header() __aligned_size(sizeof(shm_bwt_header_t), 64) 
#define bwa_shm_size_bwt_cp_occ(rlen) __aligned_size((sizeof(CP_OCC) * (((rlen) >> CP_SHIFT) + 1)), 64)
#if SA_COMPRESSION
#define bwa_shm_size_bwt_sa_ms_byte(rlen) __aligned_size(((((rlen) >> SA_COMPX) + 1) * sizeof(int8_t)), 64)
#define bwa_shm_size_bwt_sa_ls_word(rlen) __aligned_size(((((rlen) >> SA_COMPX) + 1) * sizeof(uint32_t)), 64)
#else
#define bwa_shm_size_bwt_sa_ms_byte(rlen) __aligned_size(((rlen) * sizeof(int8_t)), 64)
#define bwa_shm_size_bwt_sa_ls_word(rlen) __aligned_size(((rlen) * sizeof(uint32_t)), 64)
#endif
#define bwa_shm_size_bwt(rlen) \
							bwa_shm_size_bwt_header() \
							+ bwa_shm_size_bwt_cp_occ(rlen) \
							+ bwa_shm_size_bwt_sa_ms_byte(rlen) \
							+ bwa_shm_size_bwt_sa_ls_word(rlen)

#define bwa_shm_size_ref(rlen) __aligned_size((rlen) - 1, 64)

#define bwa_shm_size_pac(rlen) ((((rlen) - 1) >> 3) + 1) // == ((((rlen) - 1)/2)/4 + 1)

#define bwa_shm_size_kmer() (numKmers * sizeof(uint64_t))

static inline size_t bwa_shm_size_mlt(bwa_shm_info_t *info) {
	if (info->ref_file_name == NULL) {
		fprintf(stderr, "ERROR: cannot get the size of mlt since there is no reference file name in bwa_shm_info.\n");
		return 0;
	}
	return ____size_mlt(NULL, info->ref_file_name);
}

#ifdef SMEM_ACCEL
#define bwa_shm_size_accel() \
				(ALL_SMEM_TABLE_SIZE + LAST_SMEM_TABLE_SIZE)
#define bwa_shm_size_sall() (ALL_SMEM_TABLE_SIZE)
#define bwa_shm_size_slast() (LAST_SMEM_TABLE_SIZE)
#endif

#define bwa_shm_create_flags (O_RDWR | O_CREAT | O_TRUNC)
//#define bwa_shm_create_flags (O_RDWR | O_CREAT)
#define bwa_shm_create_mode (0666)
//#define bwa_shm_create_mode (0755)

static inline size_t __get_shm_size(int shmfd) {
	struct stat buf;
	int ret = fstat(shmfd, &buf);

	if (ret < 0) return ret;

	return buf.st_size;
}

int bwa_shm_create(int m, size_t s);
int bwa_shm_open(int m);
void *bwa_shm_map(int m);
int bwa_shm_unmap(int m);
int __bwa_shm_remove(int m);
int bwa_shm_remove(void);

enum bwa_shm_init_mode {
	BWA_SHM_INIT_NEW,
	BWA_SHM_INIT_READ,
#ifdef MEMSCALE
	BWA_SHM_INIT_MODIFY
#endif
};


void bwa_shm_init(const char *ref_file_name, int *useErt, enum bwa_shm_init_mode mode);
void bwa_shm_complete(enum bwa_shm_init_mode mode);
void bwa_shm_final(enum bwa_shm_init_mode mode);

int use_mmap(int m);
int __bwa_shm_load_file(const char *prefix, const char *postfix, int m, void **ret_ptr);

#endif /* USE_SHM */

#endif /* BWA_SHM_HPP */
