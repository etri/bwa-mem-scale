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

#ifdef PERFECT_MATCH
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <assert.h>
#include "bntseq.h"
#include "bwa.h"
#include "ksw.h"
#include "utils.h"
#include "kstring.h"
#include "kvec.h"
#include <string>
#include "bwa_shm.h"
#include "safe_lib.h"
#include <pthread.h>
#include "profiling.h"

#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "perfect.h"
#define BYTE_TO_GIGABYTE(x) ((x + (1 << 30) - 1) / (1 << 30))

//#define PRINT_DETAIL
#ifdef PRINT_DETAIL
#define detail(x...) printf(x)
#else
#define detail(x...) do{}while(0)
#endif

#ifdef PERFECT_PROFILE
int64_t *perfect_profile;
uint32_t perfect_profile_num_seed_entry;
uint32_t perfect_profile_seq_len;

int64_t *perfect_profile_rid;
int64_t *perfect_profile_rid_multi_loc;
int perfect_profile_num_rid;

void set_perfect_profile_rid(int n_seq) {
	perfect_profile_num_rid = n_seq;
	perfect_profile_rid = (int64_t *) calloc(n_seq, sizeof(int64_t));
	perfect_profile_rid_multi_loc = (int64_t *) calloc(n_seq, sizeof(int64_t));
}
#endif

perfect_table_t *perfect_table;
int perfect_table_seed_len; /* PT_SEED_LEN_NO_TABLE means perfect_match is off.
                               PT_SEED_LEN_AUTO_TABLE means auto detection of seedlen.
							   This variable is set before the table is loaded. */
static pthread_mutex_t perfect_auto_load_lock = PTHREAD_MUTEX_INITIALIZER;
static char *perfect_auto_load_prefix;
static uint8_t *perfect_auto_load_reference;
static FMI_search *perfect_auto_load_fmi;

#ifdef USE_SHM
int __shm_remove(int m);

static int ____load_perfect_table_from_shm(char *file_name, int len, perfect_table_t **ret_ptr) {
	perfect_table_t *pt;
	perfect_table_t *shm_ptr = NULL;
	uint8_t *ptr = NULL;

#ifdef MEMSCALE
	if (bwa_shm_info->perfect_on == 0) {
		*ret_ptr = NULL;
		return 0;
	}
#endif

	pt = (perfect_table_t *)_mm_malloc(sizeof(perfect_table_t), 64);
	if (!pt)
		return -1;

	if (bwa_shm_open(BWA_SHM_PERFECT) >= 0) {
		shm_ptr = (perfect_table_t *) bwa_shm_map(BWA_SHM_PERFECT);
	} else
		fprintf(stderr, "[bwa_shm] failed to open BWA_SHM_PERFECT\n");

	if (!shm_ptr) return -1;
	
	if (shm_ptr->seed_len != len) {
		fprintf(stderr, "[bwa_shm] perfect_table for different seed length is on memory. (%d != %d)\n",
							shm_ptr->seed_len, len);
		__bwa_shm_remove(BWA_SHM_PERFECT); /* try to remove */	
		return -1;
	}

	__lpt_link_shm_to_pt(pt, shm_ptr);

	*ret_ptr = pt;
	return 0;
}

int ____load_perfect_table_on_shm(char *file_name, int len, 
									uint32_t num_seed_load __maybe_unused, 
									perfect_table_t **ret_ptr) 
{
	perfect_table_t *pt, *shm_ptr = NULL;
	FILE *fp;
	int pct;
	size_t shm_size;
	uint32_t num_loaded;
	int fd;

	pt = (perfect_table_t *)_mm_malloc(sizeof(perfect_table_t), 64);
	if (!pt)
		return -1;

	fp = xopen(file_name, "rb");
	if (!fp) {
		fprintf(stderr, "%s: failed to open %s\n", __func__, file_name);
		goto err_struct_alloc;
	}

	__lpt_load_head(pt, fp);
#ifdef MEMSCALE
	__lpt_set_num_seed_load(pt, num_seed_load);
#endif
	bwa_shm_info->pt_seed_len = pt->seed_len;
	bwa_shm_info->pt_num_loc_entry = pt->num_loc_entry;
	bwa_shm_info->pt_num_seed_entry = pt->num_seed_entry;
#ifdef MEMSCALE
	bwa_shm_info->pt_num_seed_entry_loaded = pt->num_seed_load; 
#endif

	shm_size = __lpt_shm_size(pt);

	fprintf(stderr, "INFO: shm_create for perfect table. size: %ld hugetlb_flag: %x\n", 
					shm_size, bwa_shm_hugetlb_flags());
	if (bwa_shm_create(BWA_SHM_PERFECT, shm_size) >= 0)
		shm_ptr = (perfect_table_t *) bwa_shm_map(BWA_SHM_PERFECT);

	if (!shm_ptr) {
		goto err_file_open;
	}
	
	__lpt_show_info(pt);

	if (use_mmap(BWA_SHM_PERFECT) == 0) {
		memcpy(shm_ptr, pt, sizeof(perfect_table_t));
		__lpt_set_table_ptr(pt, shm_ptr);

		__lpt_load_loc_table(pt, fp);
		__lpt_load_seed_table(pt, fp);
	} else {
		__lpt_set_table_ptr(pt, shm_ptr);
	}

	err_fclose(fp);	

	fprintf(stderr, "Reading perfect table: Done\n");
	fflush(stderr);

	if (ret_ptr) *ret_ptr = pt;
	return 0;

err_file_open:
	fclose(fp);
err_struct_alloc:
	_mm_free(pt);
	if (ret_ptr) *ret_ptr = NULL;
	return -1;
}

static int ____load_perfect_table_without_shm(char *file_name, int len, 
												perfect_table_t **ret_ptr) 
{
	perfect_table_t *pt;
	FILE *fp;
	int pct;
	uint32_t num_loaded;
	
	pt = (perfect_table_t *)_mm_malloc(sizeof(perfect_table_t), 64);
	if (!pt)
		return -1;
	fp = xopen(file_name, "rb");
	if (!fp) {
		fprintf(stderr, "%s: failed to open %s\n", __func__, file_name);
		goto err_struct_alloc;
	}

	__lpt_load_head(pt, fp);

	__lpt_show_info(pt);
	
	pt->loc_table = (uint32_t *)_mm_malloc(pt->num_loc_entry * sizeof(uint32_t), 64);
	if (!pt->loc_table) {
		fprintf(stderr, "%s: failed to memory allocation (%.2fGB) for seed_len %d\n", 
						__func__, 
						(double) BYTE_TO_GIGABYTE(pt->num_loc_entry * sizeof(uint32_t)), 
						len);
		goto err_file_open;
	}

	__lpt_load_loc_table(pt, fp);

	pt->seed_table = (seed_entry_t *)_mm_malloc(pt->num_seed_entry * sizeof(seed_entry_t), 64);
	if (!pt->seed_table) {
		fprintf(stderr, "%s: failed to memory allocation (%.2fGB) for seed_len %u\n", 
						__func__, 
						(double) BYTE_TO_GIGABYTE(pt->num_seed_entry * sizeof(seed_entry_t)),
						len);
		goto err_table_alloc;
	}

	__lpt_load_seed_table(pt, fp);

	err_fclose(fp);	

	fprintf(stderr, "Reading perfect table: Done\n");
	fflush(stderr);

	*ret_ptr = pt;

	return 0;

err_table_alloc:
	_mm_free(pt->loc_table);
err_file_open:
	fclose(fp);
err_struct_alloc:
	_mm_free(pt);
	*ret_ptr = NULL;
	return -1;
}

int __load_perfect_table(char *file_name, int len, perfect_table_t **ret_ptr) {
	int ret;
	
	if (bwa_shm_mode == BWA_SHM_MATCHED) {
		ret = ____load_perfect_table_from_shm(file_name, len, ret_ptr);
		if (ret == 0) return 0;
	}

	if (bwa_shm_mode != BWA_SHM_DISABLE) {
		ret = ____load_perfect_table_on_shm(file_name, len, 0, ret_ptr);
		if (ret == 0) return 0;
	}

	/* ERROR or BWA_SHM_DISABLE */
	return ____load_perfect_table_without_shm(file_name, len, ret_ptr); 
}
#else /* !USE_SHM */
int __load_perfect_table(char *file_name, int len, perfect_table_t **ret_ptr) {
	perfect_table_t *pt;
	FILE *fp;
	int pct;
	uint32_t num_loaded;
	
	pt = (perfect_table_t *)_mm_malloc(sizeof(perfect_table_t), 64);
	if (!pt)
		return -1;
	fp = xopen(file_name, "rb");
	if (!fp) {
		fprintf(stderr, "%s: failed to open %s\n", __func__, file_name);
		goto err_struct_alloc;
	}

	__lpt_load_head(pt, fp);

	__lpt_show_info(pt);

	pt->loc_table = (uint32_t *)_mm_malloc(pt->num_loc_entry * sizeof(uint32_t), 64);
	if (!pt->loc_table) {
		fprintf(stderr, "%s: failed to memory allocation (%.2fGB) for seed_len %d\n", 
						__func__, 
						(double) (pt->num_loc_entry * sizeof(uint32_t)) / (1L << 30), 
						len);
		goto err_table_alloc;
	}
	__lpt_load_loc_table(pt, fp);

	pt->seed_table = (seed_entry_t *)_mm_malloc(pt->num_seed_entry * sizeof(seed_entry_t), 64);
	if (!pt->seed_table) {
		fprintf(stderr, "%s: failed to memory allocation (%.2fGB) for seed_len %d\n", 
						__func__, 
						(double) (pt->num_seed_entry * sizeof(seed_entry_t)) / (1L << 30),
						len);
		goto err_file_open;
	}

	__lpt_load_seed_table(pt, fp);
	
	err_fclose(fp);	

	fprintf(stderr, "Reading perfect table: Done\n");
	fflush(stderr);

	*ret_ptr = pt;
	return 0;

err_table_alloc:
	_mm_free(pt->seed_table);
err_file_open:
	fclose(fp);
err_struct_alloc:
	_mm_free(pt);
	*ret_ptr = NULL;
	return -1;
}
#endif /* !USE_SHM */

void init_auto_load_perfect_table(const char *prefix, uint8_t **reference, FMI_search *fmi); 

int load_perfect_table(const char *prefix, int len, 
						uint8_t **reference, FMI_search *fmi)
{
	char file_name[PATH_MAX];
	int ret;

	if (reference == NULL || *reference == NULL) {
		fprintf(stderr, "ERROR: reference for perfect table is not given.\n");
		return -1;
	}
	
	perfect_table_seed_len = len;

	if (len <= 0) {
		perfect_table = NULL;
		if (len < 0) 
			init_auto_load_perfect_table(prefix, reference, fmi);
		return -1;	
	}
	
	snprintf_s_si(file_name, PATH_MAX, "%s.perfect.%d", prefix, len);


	if (__load_perfect_table(file_name, len, &perfect_table)) {
		perfect_table = NULL;
		perfect_table_seed_len = PT_SEED_LEN_NO_TABLE;
		return -1;
	}

	if (perfect_table)
		perfect_table->ref_string = *reference;
	if (fmi) fmi->perfect_table = perfect_table;

#ifdef PERFECT_PROFILE
	perfect_profile_num_seed_entry = perfect_table->num_seed_entry;
	perfect_profile_seq_len = perfect_table->seq_len;

	perfect_profile = (int64_t *) calloc(perfect_profile_num_seed_entry, sizeof(int64_t));
	if (perfect_profile == NULL)
		return -1;
#endif

	return 0;
}

void init_auto_load_perfect_table(const char *prefix, uint8_t **reference, FMI_search *fmi) 
{
	int prefix_len = strnlen_s(prefix, PATH_MAX);
	perfect_auto_load_prefix = (char *) malloc(__aligned_size(prefix_len + 1, 64));
	strncpy_s(perfect_auto_load_prefix, PATH_MAX, prefix, prefix_len);
	perfect_auto_load_reference = *reference;
	perfect_auto_load_fmi = fmi;
}

void auto_load_perfect_table(int len) {
	uint64_t t_beg;	
	pthread_mutex_lock(&perfect_auto_load_lock);

	if (perfect_table) goto unlock;
	
	t_beg = __rdtsc();
	load_perfect_table(perfect_auto_load_prefix, len, 
					   &perfect_auto_load_reference, 
					   perfect_auto_load_fmi);
	tprof[PERFECT_TABLE_READ][0] += __rdtsc() - t_beg;
unlock:
	pthread_mutex_unlock(&perfect_auto_load_lock);
}

void free_perfect_table() {
	if (perfect_auto_load_prefix) free(perfect_auto_load_prefix);
#ifdef PERFECT_PROFILE
#define MAX_NUM_LOC_DIST 30
//#define MAX_LOC_DIST 50
#define MAX_IDX_DIST 100
	int64_t idx;
	int64_t num_matched = 0;
	int64_t num_twice_matched = 0;
	uint32_t num_matched_entry = 0;
	uint32_t num_twice_matched_entry = 0;
	seed_entry_t *ent;
	int num_loc;
	int64_t num_loc_dist[MAX_NUM_LOC_DIST];
	int64_t num_loc_twice_dist[MAX_NUM_LOC_DIST];
	//int64_t loc_dist[MAX_LOC_DIST];
	int64_t idx_dist[MAX_IDX_DIST];
	for (idx = 0; idx < MAX_NUM_LOC_DIST; ++idx)
		num_loc_dist[idx] = 0;
	for (idx = 0; idx < MAX_NUM_LOC_DIST; ++idx)
		num_loc_twice_dist[idx] = 0;
	//for (idx = 0; idx < MAX_LOC_DIST; ++idx)
	// 	loc_dist[idx] = 0;
	for (idx = 0; idx < MAX_IDX_DIST; ++idx)
	 	idx_dist[idx] = 0;
	for (idx = 0; idx < perfect_profile_num_seed_entry; ++idx) {
		if (perfect_profile[idx] == 0)
			continue;
		ent = get_seed_entry(perfect_table, idx);
		num_loc = __get_num_location(ent->flags, perfect_table->loc_table);
		num_matched += perfect_profile[idx];
		num_matched_entry++;
		//fprintf(stderr, "[PERFECT_PROFILE_EACH] idx: %10ld num: %6ld loc: %10u num_loc: %6u\n",
		//				idx, perfect_profile[idx], ent->location, num_loc); 
		if (num_loc < MAX_NUM_LOC_DIST)
			num_loc_dist[num_loc - 1] += perfect_profile[idx];
		else
			num_loc_dist[MAX_NUM_LOC_DIST - 1] += perfect_profile[idx];
		if (perfect_profile[idx] > 1) {
			num_twice_matched += perfect_profile[idx];
			num_twice_matched_entry++;
			if (num_loc < MAX_NUM_LOC_DIST)
				num_loc_twice_dist[num_loc - 1] += perfect_profile[idx];
			else
				num_loc_twice_dist[MAX_NUM_LOC_DIST - 1] += perfect_profile[idx];
		}

		//loc_dist[((int64_t) ent->location) * MAX_LOC_DIST / perfect_profile_seq_len] += perfect_profile[idx];
		idx_dist[((int64_t) idx) * MAX_IDX_DIST / perfect_profile_num_seed_entry] += perfect_profile[idx];
	}

	fprintf(stderr, "[PERFECT_PROFILE_TOTAL] seq_len: %12u num_matched: %12ld num_twice_matched: %12ld %8.2f%%\n"
			"[PERFECT_PROFILE_TOTAL] num_matched_entry: %10u / %10u %8.4f%% matched_twice_and_more_entry: %10u / %10u %8.4f%%\n",
						perfect_profile_seq_len, 
						num_matched, num_twice_matched,
						((float) num_twice_matched) * 100 / num_matched,
						num_matched_entry, perfect_profile_num_seed_entry,
						((float) num_matched_entry) * 100 / perfect_profile_num_seed_entry,
						num_twice_matched_entry, num_matched_entry,
						((float) num_twice_matched_entry) * 100 / num_matched_entry
						);
	for (idx = 0; idx < MAX_NUM_LOC_DIST; ++idx)
		fprintf(stderr, "[PERFECT_PROFILE_NUM_LOC_DIST] num_loc: %2ld num: %10ld / %10ld %8.4f%% twice: %10ld / %10ld %8.4f%%\n",
				idx + 1, 
				num_loc_dist[idx], num_matched, 
				((float) num_loc_dist[idx]) * 100 / num_matched,
				num_loc_twice_dist[idx], num_twice_matched, 
				((float) num_loc_twice_dist[idx]) * 100 / num_twice_matched
				);
	/*for (idx = 0; idx < MAX_LOC_DIST; ++idx)
		fprintf(stderr, "[perfect_profile_loc_dist] idx: %3ld loc: %10ld ~ %10ld num: %10ld / %10ld %6.2f%%\n",
				idx, 
				((int64_t) perfect_profile_seq_len) * idx / max_loc_dist,
				((int64_t) perfect_profile_seq_len) * (idx + 1) / max_loc_dist - 1,
				loc_dist[idx], num_matched, 
				((float) loc_dist[idx]) * 100 / num_matched); */
	for (idx = 0; idx < MAX_IDX_DIST; ++idx)
		fprintf(stderr, "[PERFECT_PROFILE_IDX_DIST] dist: %3ld idx: %10ld ~ %10ld num: %10ld / %10ld %6.2f%%\n",
				idx, 
				((int64_t) perfect_profile_num_seed_entry) * idx / MAX_IDX_DIST,
				((int64_t) perfect_profile_num_seed_entry) * (idx + 1) / MAX_IDX_DIST - 1,
				idx_dist[idx], num_matched, 
				((float) idx_dist[idx]) * 100 / num_matched);
	for (idx = 0; idx < perfect_profile_num_rid; ++idx)
		fprintf(stderr, "[PERFECT_PROFILE_RID_DIST] rid: %3ld num: %10ld / %10ld %8.4f%% multi_loc_num: %10ld %8.4f%%\n",
				idx, 
				perfect_profile_rid[idx],
				num_matched, 
				((float) perfect_profile_rid[idx]) * 100 / num_matched,
				perfect_profile_rid_multi_loc[idx],
				((float) perfect_profile_rid_multi_loc[idx]) * 100 / num_matched
				);

#endif
	if (perfect_table == NULL) return;
#ifdef USE_SHM
	if (bwa_shm_unmap(BWA_SHM_PERFECT)) {
		_mm_free(perfect_table->loc_table);
		_mm_free(perfect_table->seed_table);
	}
#else
	_mm_free(perfect_table->loc_table);
	_mm_free(perfect_table->seed_table);
#endif
	_mm_free(perfect_table);
}

#define seedcmp_find(pt, ent, seed, fw_less) \
				__seedcmp(pt->ref_string + ent->location, is_fw_less_entry(ent), \
						  seed,	fw_less, \
						  pt->seed_len)

/*#define seedmatch_find(pt, ent, seed) \
				__seedmatch(pt, pt->ref_string + ent->location, seed) */

static int seedmatch_further(perfect_table_t *pt, seed_entry_t *ent, 
								uint8_t *seed, int fw_less, int len, 
								bseq1_perfect_t *ret) {
	
	int is_rev = is_fw_less_entry(ent) == fw_less ? 0 : 1;
	uint32_t location = UINT32_MAX;
	uint32_t multi_loc;
	if (__seedmatch_further(pt, ent->location, seed,
						is_rev, len)) {
		location = ent->location;
		goto out;
	}

	multi_loc = get_multi_location(ent);
	if (multi_loc == 0)
		goto out;
	else {
		uint32_t i, num_fw, num_rc, *loc_fw, *loc_rc;

		GET_MULTI_FW_AND_RC(pt->loc_table, multi_loc,
							num_fw, loc_fw,
							num_rc, loc_rc);
	
		for (i = 0; i < num_fw; ++i) {
			if (__seedmatch_further(pt, loc_fw[i], seed, is_rev, len)) {
				location = loc_fw[i];
				goto out;
			}
		}
	
		is_rev = !is_rev;
		for (i = 0; i < num_rc; ++i) {
			if (__seedmatch_further(pt, loc_rc[i], seed, is_rev, len)) {
				location = loc_rc[i];
				goto out;
			}
		}
	}

out:
	if (location == UINT32_MAX)
		//return FIND_PERFECT_NOT_MATCHED;
		return FIND_PERFECT_SEED_ONLY_MATCHED;
	else {
		ret->location = location;
		if (is_rev == 0) {
			ret->flags = ent->flags & (~(FLAG_RC)) | (FLAG_VALID);
			return FIND_PERFECT_FW_MATCHED;
		} else {
			ret->flags = ent->flags | (FLAG_RC) | (FLAG_VALID);
			return FIND_PERFECT_RC_MATCHED;
		}
	}
}

static int __find_perfect_match_entry(perfect_table_t *pt, 
									 uint8_t *seed, int len, bseq1_perfect_t *ret) {
	int64_t idx;
	seed_entry_t *ent;
	int pt_len = pt->seed_len;
	int fw_less = __compare_fw_rc(seed, pt->seed_len);
	int cmp;

	idx = __get_hash_idx_seed(pt, seed, fw_less);
	ent = get_seed_entry(pt, idx);

	if (!is_hash_matched_entry(ent))
		return FIND_PERFECT_NOT_MATCHED;

	do {
		cmp = seedcmp_find(pt, ent, seed, fw_less);
		if (cmp == 0) { /* found */
			if (len == pt_len) {
				ret->location = ent->location;
				ret->flags = is_fw_less_entry(ent) == fw_less
								? (ent->flags & (~(FLAG_RC)) | (FLAG_VALID))
								: (ent->flags | (FLAG_RC) | (FLAG_VALID));
#ifdef PERFECT_PROFILE
				perfect_profile[idx]++;
#endif
				return is_fw_less_entry(ent) == fw_less
							? FIND_PERFECT_FW_MATCHED
							: FIND_PERFECT_RC_MATCHED;
			} else {
				int retval = seedmatch_further(pt, ent, seed, fw_less, len, ret);
#ifdef PERFECT_PROFILE
				if (retval != FIND_PERFECT_NOT_MATCHED)
					perfect_profile[idx]++;
#endif
				return retval;
			}
		} else if (cmp > 0) { // ent->seed > seed ==> go to left
			idx = ent->left;
			ent = get_seed_entry(pt, idx);
		} else { // ent->seed > seed ==> go to right
			idx = ent->right;
			ent = get_seed_entry(pt, idx);
		}
	} while (ent);

	return FIND_PERFECT_NOT_MATCHED;
}

static int seed_with_N(uint8_t *seed, int len) {
	int ret = 0, i;
	for (i = 0; i < len; ++i)
		ret = ret | (seed[i] & 0xC);
	return ret;
}

int find_perfect_match_entry(perfect_table_t *pt, bseq1_t *seq, int len) {
	uint8_t *seed = (uint8_t *) seq->seq;

	/* NOTE: initial value of seq->perfect.exist is 0. */
	if (len < perfect_table_seed_len) {
		if (perfect_table_seed_len != PT_SEED_LEN_AUTO_TABLE)
			return FIND_PERFECT_NO_TABLE;
		else  {
			auto_load_perfect_table(len);
			if (len >= perfect_table_seed_len)
				pt = perfect_table;
		} 
	}
	
	if (pt == NULL)
		return FIND_PERFECT_NO_TABLE;

	if (seed_with_N(seed, len))
		return FIND_PERFECT_WITH_N;

	return __find_perfect_match_entry(pt, seed, len, &seq->perfect);
}

void init_mem_aln_perfect(mem_aln_perfect_t *a, int64_t pos, int len, int is_rev, const bntseq_t *bns, int seed_len) {
	bntann1_t *ann;
	/* find_perfect_match_entry() finds the exact locations for both of FW and RC matched cases.
	 * Thus, we do not need to do things like bns_depos() */ 
	a->loc = pos; // for pair-end
	a->rid = bns_pos2rid(bns, pos);
	ann = &bns->anns[a->rid];
	if (len != seed_len && is_rev)
		pos = pos - (len - seed_len);
	a->pos = pos - ann->offset;
	a->flag = 0;
	a->is_rev = is_rev;
	a->is_alt = ann->is_alt;
	//a->XA = 0;
	a->sub = 0;
}

static inline void init_mem_aln_perfect_multi_loc(mem_aln_perfect_v *av, uint32_t num_loc, uint32_t *locs, 
											bseq1_t *s, int is_rev, 
											const bntseq_t *bns, perfect_table_t *pt) {
	uint32_t i, idx;
	int seed_len = pt->seed_len;
	uint32_t loc;
	uint32_t matched_loc = s->perfect.location;
	for (i = 0; i < num_loc; ++i) {
		idx = is_rev ? num_loc - 1 - i : i; /* sort the vector by rb */
		loc = locs[idx];
		if (loc == matched_loc)
			continue;
		if (seed_len == s->l_seq 
				|| __seedmatch_further(pt, loc, (uint8_t *) s->seq, is_rev, s->l_seq))
			init_mem_aln_perfect(&av->a[av->n++], (int64_t) loc,
								s->l_seq, is_rev, bns, pt->seed_len);
	}
}

/* av.a include FW matched entries first */
mem_aln_perfect_v get_perfect_locations(bseq1_t *s, const bntseq_t *bns, perfect_table_t *pt) {
	mem_aln_perfect_v av = {0, 0, 0};
	uint32_t flags = s->perfect.flags;
	int rc_matched = __is_rc_matched(flags);
	uint32_t multi_loc = __get_multi_location(flags);
	uint32_t location = s->perfect.location;
	
	assert(s->perfect.exist != 0);

	av.m += __get_num_location(flags, pt->loc_table);
	assert(av.m > 0);

	av.a = (mem_aln_perfect_t *) calloc(av.m, sizeof(mem_aln_perfect_t));

	if (multi_loc == 0) {
		init_mem_aln_perfect(&av.a[av.n++], 
							(int64_t) s->perfect.location, 
							s->l_seq, 
							rc_matched ? 1 : 0,
							bns, pt->seed_len);
	} else {
		uint32_t i, num_fw, num_rc, *loc_fw, *loc_rc;

		GET_MULTI_FW_AND_RC(pt->loc_table, multi_loc,
							num_fw, loc_fw,
							num_rc, loc_rc);

		if (!rc_matched) {
			init_mem_aln_perfect(&av.a[av.n++], (int64_t) s->perfect.location, 
											s->l_seq, 0, bns, pt->seed_len);
			init_mem_aln_perfect_multi_loc(&av, num_fw, loc_fw,
											s, 0, bns, pt);
			init_mem_aln_perfect_multi_loc(&av, num_rc, loc_rc,
											s, 1, bns, pt);
		} else { /* rc_matched */
			init_mem_aln_perfect_multi_loc(&av, num_rc, loc_rc,
											s, 0, bns, pt);
			init_mem_aln_perfect(&av.a[av.n++], (int64_t) s->perfect.location, 
											s->l_seq, 1, bns, pt->seed_len);
			init_mem_aln_perfect_multi_loc(&av, num_fw, loc_fw,
											s, 1, bns, pt);
		}
	}

#ifdef PERFECT_PROFILE
	if (multi_loc == 0) {
		int rid = av.a[av.n-1].rid;
		__sync_add_and_fetch(&perfect_profile_rid[rid], 1);
		__sync_add_and_fetch(&perfect_profile_rid_multi_loc[rid], 1);
	} else {
		int i;
		int got_1st = 0; /* found the entry that a->pos == s->perfect.location */
		int profile_rid[bns->n_seqs];
		for (i = 0; i < bns->n_seqs; ++i)
			profile_rid[i] = 0;

		for (i = 0; i < av.n; i++) {
			rid = av.a[i].rid;
			profile_rid[rid] = 1;

			if (got_1st == 0 && rc_matched == av.a[i].is_rev) {
				__sync_add_and_fetch(&perfect_profile_rid[rid], 1);
				got_1st = 1;
			}
		}

		for (i = 0; i < bns->n_seqs; ++i) {
			if (profile_rid[i] == 0)
				continue;
			__sync_add_and_fetch(&perfect_profile_rid_multi_loc[i], 1);
		}
	}
#endif
	
	return av;
}

#define check_excluded(a) ((a)->rid < 0)
#define set_excluded(a) do{(a)->rid = -1;}while(0)
int perfect_dedup_patch(const mem_opt_t *opt, int n, int l_seq, 
						mem_aln_perfect_t *a)
{
	int m, i, j;
    if (n <= 1) return n;
	// all lengths are same and already sorted from get_perfect_locations
    
    for (i = 1; i < n; ++i)
    {
        mem_aln_perfect_t *p = &a[i];
		// we can compare a->pos and b->pos instead of a->rb and b->rb,
		// since we compare them only when a->rid == b->rid
        if (p->rid != a[i-1].rid || p->is_rev != a[i-1].is_rev || p->pos >= a[i-1].pos + l_seq + opt->max_chain_gap)
            continue; // then no need to go into the loop below
       
        for (j = i - 1; j >= 0 && p->rid == a[j].rid && p->is_rev == a[j].is_rev && p->pos < a[j].pos + l_seq + opt->max_chain_gap; --j) {
            mem_aln_perfect_t *q = &a[j];
            if (check_excluded(q)) continue; // a[j] has been excluded

			// default value of mask_level_redun = 0.95
            if (q->pos + l_seq - p->pos  > opt->mask_level_redun * l_seq)
				set_excluded(q);
			// mem_patch_reg() always return 0 since q->qb == p->qb == 0
        }
    }
    for (i = 0, m = 0; i < n; ++i) // exclude identical hits
        if (!check_excluded(&a[i])) {
            if (m != i) a[m++] = a[i];
            else ++m;
        }
    n = m;
    return m;

}

int mem_perfect2reg(const mem_opt_t *opt, perfect_table_t *pt, const bntseq_t *bns, 
						bseq1_t *s, mem_alnreg_v *reg) 
{
	mem_aln_perfect_v av = {0, 0, 0};
	int ret;
	int i;
	int l_seq = s->l_seq;

	av = get_perfect_locations(s, bns, pt);
	av.n = perfect_dedup_patch(opt, av.n, l_seq, av.a);
	ret = av.a[0].is_rev;

	reg->n = av.n;
	reg->m = av.n;
	reg->a = (mem_alnreg_t *) calloc(av.n, sizeof(mem_alnreg_t)); 

	for (i = 0; i < av.n; ++i) {
		mem_aln_perfect_t *p = &av.a[i];
		mem_alnreg_t *r = &reg->a[i];

		// FW: [ p->loc , p->loc + s->l_seq )
		// RC: [ (bns->l_pac<<1) + 1 - (p->loc) , (bns->l_pac<<1) + 1 - (p->loc + s->l_seq) )
		//     => ( (bns->l_pac<<1) + 1 - (p->loc - 1) , (bns->l_pac<<1) + 1 - (p->loc + s->l_seq - 1) ]
		//     => [ (bns->l_pac<<1) + 1 - (p->loc + s->l_seq - 1) , (bns->l_pac<<1) + 1 - (p->loc - 1) ) 
		//     => [ (bns->l_pac<<1) - (p->loc + s->l_seq) , (bns->l_pac<<1) - (p->loc) ) 
		if (!p->is_rev) {
			r->rb = p->loc;
			r->re = p->loc + l_seq;
		} else {
			r->rb = (bns->l_pac<<1) - (p->loc + l_seq);
			r->re = (bns->l_pac<<1) - p->loc;
		}
		r->qb = 0;
		r->qe = l_seq;
		r->rid = p->rid;
		r->score = l_seq * opt->a;
		r->truesc = l_seq * opt->a;

		/* due to calloc(), we can skip 0 values */
		//r->sub = 0;
		//r->alt_sc = 0;
		//r->csub = 0;
		//r->sub_n = 0;
		r->w = opt->w;
		//r->seedcov = 0;
		//r->secondary = 0;
		//r->secondary_all = 0;
		r->seedlen0 = l_seq;
		r->n_comp = 1;
		r->is_alt = p->is_alt;
		//r->frac_rep = 0;
		//r->hash = 0;
		//r->flg = 0;
	}

	free(av.a);
	return ret;
}

#endif /* PERFECT_MATCH */
