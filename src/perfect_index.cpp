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
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include <emmintrin.h> // _mm_free
#include <cassert>
#include <errno.h>
#include "bwa.h"
#include "utils.h"
#include "fastmap.h"
#include "read_index_ele.h"
#include "perfect.h"
#include "safe_lib.h"
#include "bwa_shm.h"
#include <semaphore.h>
#include <pthread.h>
#include <sched.h>
#include <sys/time.h>

//#define DEBUG_MESSAGE
//#define ASSERT_IN_DETAIL

#ifdef ASSERT_IN_DETAIL
#define d_assert(x, pt, key) do { \
			if (x) break; \
			show_perfect_table_related(pt, key); \
			fflush(stdout); fflush(stderr); \
			assert(x); \
		} while(0)
#else
#define d_assert(x, y, z) assert(x)
#endif
#ifdef DEBUG_MESSAGE
#define dbg_printf(x...) printf(x)
#define dbg_show_seed_entry(x,y) show_seed_entry(x,y)
#define dbg_show_perfect_table_related(x,y) show_perfect_table_related(x,y)
#else
#define dbg_printf(x...) do{}while(0)
#define dbg_show_seed_entry(x,y) do{}while(0)
#define dbg_show_perfect_table_related(x,y) do{}while(0)
#endif

#define free_ptr(ptr) do { if (ptr) { free(ptr); ptr = NULL; } } while(0)

#define PREFETCH_DISTANCE 10

#define NUM_ENTRY_PER_PRINT 1000000
//#define NUM_ENTRY_PER_PRINT 100
static uint32_t total_added_entry = 0;
static uint32_t total_hole_entry = 0;
static uint32_t total_moved_entry = 0;


int mode_build = 0; /* use this variable only for debugging. now, affect to "show_seed_entry()" */

#define get_num_seed(start, end, len) ((end) - (start) >= (len) ? ((end) - (start) - (len) + 1) : 0)

static inline void *recallocarray(void *ptr, size_t old_n, size_t new_n, size_t size)
{
	ptr = reallocarray(ptr, new_n, size);
	if (!ptr) return ptr;

	if (old_n < new_n) {
		uint8_t *_ptr = (uint8_t *) ptr;
		memset(_ptr + (old_n * size), 0, (new_n - old_n) * size);
	}

	return ptr;
}

static bntann1_t *ann_restore(const char *prefix, int64_t *_seq_len, int32_t *_n_seqs) {
	char str[8193];
	char fn_ann[PATH_MAX];
	bntann1_t *anns = NULL;
	int32_t n_seqs;
	FILE *fp;
	int scanres;
	int i;
	long long xx;
	int64_t seq_len;
	uint32_t seed;
	
	strcpy_s(fn_ann, PATH_MAX, prefix);
	strcat_s(fn_ann, PATH_MAX, ".ann");
		
	fp = xopen(fn_ann, "r");
	scanres = fscanf(fp, "%lld%d%u", &xx, &n_seqs, &seed);
	assert(n_seqs >= 0 && n_seqs <= INT_MAX);
	if (scanres != 3) goto badread;

	seq_len = xx;
	anns = (bntann1_t*)calloc(n_seqs, sizeof(bntann1_t));
	assert(anns != NULL);
	for (i = 0; i < n_seqs; ++i) {
		bntann1_t *p = anns + i;
		char *q = str;
		int c;
		// read gi and sequence name
		scanres = fscanf(fp, "%u%8192s", &p->gi, str);
		if (scanres != 2) goto badread;
		p->name = strdup(str);
		// read fasta comments 
		while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
		while (c != '\n' && c != EOF) c = fgetc(fp);
		if (c == EOF) {
			scanres = EOF;
			goto badread;
		}
		*q = 0;
		assert(strnlen_s(str, 8192) < 8192);
		if (q - str > 1 && strcmp(str, " (null)") != 0) p->anno = strdup(str + 1); // skip leading space
		else p->anno = strdup("");
		// read the rest
		scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
		if (scanres != 3) goto badread;
		p->offset = xx;
	}
	err_fclose(fp);

	*_seq_len = seq_len;
	*_n_seqs = n_seqs;

	return anns;

badread:
	if (anns) free(anns);
	if (EOF == scanres) {
		err_fatal(__func__, "Error reading %s : %s\n", fn_ann, ferror(fp) ? strerror(errno) : "Unexpected end of file");
	}
	err_fatal(__func__, "Parse error reading %s\n", fn_ann);
}

static bntamb1_t *amb_restore(const char *prefix, int64_t *_seq_len, int32_t *_n_holes) {
	char str[8193];
	char fn_amb[PATH_MAX];
	FILE *fp;
	bntamb1_t *ambs;
	long long xx;
	int i;
	int scanres;
	int64_t seq_len;
	int32_t n_seqs, n_holes;

	strcpy_s(fn_amb, PATH_MAX, prefix);
	strcat_s(fn_amb, PATH_MAX, ".amb");
	fp = xopen(fn_amb, "r");
	scanres = fscanf(fp, "%ld%d%d", &seq_len, &n_seqs, &n_holes);
	assert(n_holes >= 0 && n_holes <= INT_MAX);
	if (scanres != 3) goto badread;
	ambs = n_holes ? (bntamb1_t *)calloc(n_holes, sizeof(bntamb1_t)) : 0;
	for (i = 0; i < n_holes; ++i) {
		bntamb1_t *p = ambs + i;
		scanres = fscanf(fp, "%lld%d%8192s", &xx, &p->len, str);
		if (scanres != 3) goto badread;
		p->offset = xx;
		p->amb = str[0];
	}
	err_fclose(fp);

	*_seq_len = seq_len;
	*_n_holes = n_holes;

	return ambs;

badread:
	if (ambs) free(ambs);
	if (EOF == scanres) {
		err_fatal(__func__, "Error reading %s : %s\n", fn_amb, ferror(fp) ? strerror(errno) : "Unexpected end of file");
	}
	err_fatal(__func__, "Parse error reading %s\n", fn_amb);
}

uint32_t get_empty_idx(perfect_table_t *pt, uint32_t key) {
	uint32_t idx;
	if (!is_valid_entry(get_seed_entry(pt, key)))
		return key;

	// simple linear search to find the nearest cache line
	for (idx = (key + 1) % pt->num_seed_entry;
			idx != key;
				idx = (idx + 1) % pt->num_seed_entry) {
		if (!is_valid_entry(get_seed_entry(pt, idx)))
			return idx;
	}

	return NO_ENTRY;
}

/* build_loc and related functions */
typedef struct {
	int fw_n; // number of FW locations
	int fw_m; // allocated size
	uint32_t *fw;
	int rc_n; // number of RC locations
	int rc_m; // allocated size
	uint32_t *rc;
} build_loc_t;

int build_loc_n = 0;
int build_loc_m = 0;
build_loc_t *build_loc = NULL;

static void build_loc_init() {
	/* actually, nothing to do, but... */
	build_loc_n = 0;
	build_loc_m = 0;
	build_loc = NULL;
}

static uint32_t build_loc_new() {
	/* initially, build_loc_n == build_loc_m == 0 */	
	uint32_t ret = build_loc_n;
	
	if (build_loc_n >= build_loc_m) {
		/* the number of entry to add should be larger than 2, since we don't use 0th entry */
		build_loc = (build_loc_t *) recallocarray(build_loc, build_loc_m, build_loc_m + 256, sizeof(build_loc_t));
		build_loc_m += 256;
		assert(build_loc != NULL);
	
		/* we don't use 0th entry. return 1 for the very first allocation. */
		if (build_loc_n == 0) {
			build_loc_n = 2;
			return 1;
		}
	}

	build_loc_n++;
	return ret;
}

static build_loc_t *build_loc_get_entry(uint32_t multi_loc) {
	if (multi_loc == 0 || multi_loc >= build_loc_n)
		return NULL;
	return &build_loc[multi_loc];
}

static int build_loc_add_loc(uint32_t multi_loc, uint32_t loc, int is_rev) {
	build_loc_t *bloc;	
	if (multi_loc == 0)
		return -1;
	
	bloc = &build_loc[multi_loc];

	if (!is_rev) {
		if (bloc->fw_n >= bloc->fw_m) {
			int old_m = bloc->fw_m;
			bloc->fw_m = old_m == 0 ? 4 : old_m * 2;
			bloc->fw = (uint32_t *) recallocarray(bloc->fw, old_m, bloc->fw_m, sizeof(uint32_t));
			assert(bloc->fw != NULL);	
		}

		bloc->fw[bloc->fw_n++] = loc;
	} else {
		if (bloc->rc_n >= bloc->rc_m) {
			int old_m = bloc->rc_m;
			bloc->rc_m = old_m == 0 ? 4 : old_m * 2;
			bloc->rc = (uint32_t *) recallocarray(bloc->rc, old_m, bloc->rc_m, sizeof(uint32_t));
			assert(bloc->rc != NULL);	
		}

		bloc->rc[bloc->rc_n++] = loc;
	}
	return 0;
}

static inline uint32_t *build_loc_to_loc_table(uint32_t *loc_n, 
						uint32_t **__multi_loc_map, uint32_t *map_n) {
	uint32_t b; /* index for build_loc */
	uint32_t *loc_table;
	uint32_t *multi_loc_map;
	uint32_t n = 0; /* number of entries for loc_table */
	uint32_t i, j;
	uint32_t _n, i_many;
	build_loc_t *bloc;
	uint32_t num_fw_loc = 0, num_rc_loc = 0, num_many = 0, num_fw_loc_many = 0, num_rc_loc_many = 0;

	if (build_loc_n >= (FLAG_MULTI_LOC_MAX - 1) / 2)
		return NULL;

	multi_loc_map = (uint32_t *) malloc(build_loc_n * sizeof(uint32_t));
	*map_n = build_loc_n;
	if (multi_loc_map == NULL)
		return NULL;

	/* count loc_n */
	n = 1;
	i_many = 1;
	for (b = 1; b < build_loc_n; ++b) {
		bloc = build_loc_get_entry(b);
		_n = bloc->fw_n + bloc->rc_n;
		if (bloc->fw_n < LOC_MANY && bloc->rc_n < LOC_MANY) {
			n += 1 + _n; /* encoded num_fw and num_rc */
			i_many += 1 + _n;
		} else {
			n += 3 + _n; /* index of second entry, num_fw, and num_rc */
			i_many += 1;
		}
	}

	/* increase n to align the loc_table in the cache line size,
		for mmap()ed perfect table */ 
	n = n + ((64 / sizeof(uint32_t)) - (n % (64 / sizeof(uint32_t))));
	loc_table = (uint32_t *) calloc(n, sizeof(uint32_t));
	*loc_n = n;
	
	if (loc_table == NULL)
		return NULL;

	/* fill up loc_table */
	loc_table[0] = 0; /* NULL entry */
	multi_loc_map[0] = 0;
	n = 1; /* starting index */

	for (b = 1; b < build_loc_n; ++b) {
		bloc = build_loc_get_entry(b);
		_n = bloc->fw_n + bloc->rc_n;
		num_fw_loc += bloc->fw_n;
		num_rc_loc += bloc->rc_n;
		
		assert(n <= FLAG_MULTI_LOC_MAX);
		multi_loc_map[b] = n;
		
		if (bloc->fw_n < LOC_MANY && bloc->rc_n < LOC_MANY) {
			loc_table[n++] = (bloc->fw_n << 16) + bloc->rc_n;
			i = n;
			n += _n;
		} else {
			loc_table[n++] = 0x80000000 | i_many;
			loc_table[i_many++] = bloc->fw_n;
			loc_table[i_many++] = bloc->rc_n;
			i = i_many;
			i_many += _n;

			num_fw_loc_many += bloc->fw_n;
			num_rc_loc_many += bloc->rc_n;
			num_many++;
		}

		for (j = 0; j < bloc->fw_n; ++j) {
			loc_table[i++] = bloc->fw[j];
		}

		for (j = 0; j < bloc->rc_n; ++j) {
			loc_table[i++] = bloc->rc[j];
		}
	}

	printf("%s: num_loc_entry: %u num_seed: %u num_fw_loc: %u num_rc_loc: %u"
			" num_seed_many: %u num_fw_loc_many: %u num_rc_loc_many: %u\n", 
				__func__, *loc_n, build_loc_n - 1, num_fw_loc, num_rc_loc,
				num_many, num_fw_loc_many, num_rc_loc_many);

	*__multi_loc_map = multi_loc_map;
	return loc_table;	
}

static void __build_loc_free(build_loc_t *bloc) {
	bloc->fw_n = 0;
	bloc->fw_m = 0;
	free_ptr(bloc->fw);
	bloc->rc_n = 0;
	bloc->rc_m = 0;
	free_ptr(bloc->rc);
}

static void build_loc_free() {
	uint32_t i;
	for (i = 0; i < build_loc_n; ++i)
		__build_loc_free(&build_loc[i]);
	build_loc_n = 0;
	build_loc_m = 0;
	free_ptr(build_loc);
}

#define seedcmp_entries(pt, a, b) \
				__seedcmp((pt)->ref_string + (a)->location, is_fw_less_entry(a), \
						  (pt)->ref_string + (b)->location, is_fw_less_entry(b), \
						  (pt)->seed_len)

static inline int seedmatch_loc_to_loc(perfect_table_t *pt, uint32_t aloc, uint32_t bloc) {
	return __seedmatch(pt, pt->ref_string + aloc, pt->ref_string + bloc);
}

void show_perfect_table_stat(perfect_table_t *pt, uint32_t loc) {

	printf("HASH_TABLE: [%4.1f%%] seed_len: %u seq_len: %u #seq: %u "
		   "#moved: %u (%.1f%%) #seed_entry: %u #used_seed: %u (%.1f%%) "
		   "#seed_key: %u collision: %5.2f%% #loc_entry: %u (%.2f%%)\n",
				((float) loc) * 100 / (float) pt->seq_len,	
				pt->seed_len, pt->seq_len, 
				total_added_entry,
				total_moved_entry,
				(float) total_moved_entry * 100 / (float) total_added_entry,
				pt->num_seed_entry,
				pt->num_seed_used,
				(float) pt->num_seed_used * 100 / (float) pt->num_seed_entry,
				pt->num_seed_key,
				(float) (pt->num_seed_used - pt->num_seed_key) * 100 / (float) pt->num_seed_used,
				mode_build ? build_loc_n : pt->num_loc_entry, (float) (mode_build ? build_loc_n : pt->num_loc_entry) * 100 / (float) pt->num_seed_used);
	fflush(stdout);
}

/* DO NOT use this function in parallel or twice in a line with str == NULL*/
static char *__get_seed_str(uint32_t loc, int len, perfect_table_t *pt, char *str) {
	static char *__str = NULL;
	static int __len = 0;
	uint32_t i;
	uint32_t beg = loc;
	uint32_t end = loc + len;
	uint8_t s;
	int pos = 0;

	if (end > pt->seq_len) {
		end = pt->seq_len;
		len = end - beg;
	}

	if (str == NULL) {
		if (__len < len) {
			__str = (char *)realloc(__str, len + 1);
			assert(__str != NULL);
			__len = len;
		}
		str = __str;
	}
	
	for (i = beg; i < end; ++i) {
		s = get_seed_loc(pt, i);
		str[pos++] = "ACGTN"[s];
	}
	str[pos] = '\0';

	return str;
}

inline void print_seed_str(const char *head, seed_entry_t *ent, perfect_table_t *pt) {
	char str[pt->seed_len + 1];
	printf("%s%s", head, __get_seed_str(ent->location, pt->seed_len, pt, str));
}

void show_seed_entry(perfect_table_t *pt, uint32_t key) {
	int i, collision, num_multi_fw, num_multi_rc;
	uint32_t multi_loc;
	uint32_t *loc_fw, *loc_rc;
	seed_entry_t *ent = get_seed_entry(pt, key);

	printf("SEED_ENTRY[%08x] ", key);
	if (ent == NULL) {
		printf("NULL\n");
		return;
	}
	
	if (!is_valid_entry(ent)) {
		printf("invalid flags: %8x location: %8x left: %8x right: %8x\n",
					ent->flags, ent->location, ent->left, ent->right);
		return;
	}

	multi_loc = get_multi_location(ent);

	if (multi_loc == 0) {
		num_multi_fw = 0;
		num_multi_rc = 0;
		loc_fw = NULL;
		loc_rc = NULL;
	} else {
		if (mode_build) { /* mode: build */
			build_loc_t *bloc = build_loc_get_entry(multi_loc);
			assert(bloc);
			num_multi_fw = bloc->fw_n;
			loc_fw = bloc->fw;
			num_multi_rc = bloc->rc_n;
			loc_rc = bloc->rc;
		} else { /* mode: mapping */
			GET_MULTI_FW_AND_RC(pt->loc_table, multi_loc,
								num_multi_fw, loc_fw,
								num_multi_rc, loc_rc);
		}
	}

	printf("%7s %9s #loc: %4d %8s left: %8x right: %8x ",
			is_fw_less_entry(ent) ? "fw_less" : "rc_less",
			is_collision_entry(ent) ? "collision" : "matched",
			num_multi_fw + num_multi_rc + 1,
			multi_loc ? "(multi)" : "(single)",
			ent->left,
			ent->right);

	print_seed_str("seed: ", ent, pt);

	printf(" location: %8x", ent->location);
	for (i = 0; i < num_multi_fw; ++i)
		printf(" %8x", loc_fw[i]);
	if (num_multi_rc > 0) {
		printf(" RC:");
		for (i = 0; i < num_multi_rc; ++i)
			printf(" %8x", loc_rc[i]);
	}
	printf("\n");
	return;	
}

void show_perfect_table(perfect_table_t *pt) {
	uint32_t i;
	seed_entry_t *ent;
	printf("=================================================================\n");
	show_perfect_table_stat(pt, pt->seq_len);	
	printf("=================================================================\n");
	for (i = 0; i < pt->num_seed_entry; i++) {
		ent = get_seed_entry(pt, i);
		if (!is_valid_entry(ent))
			continue;
		show_seed_entry(pt, i);
	}
	printf("=================================================================\n");
}

void __show_perfect_table_related(perfect_table_t *pt, uint32_t start) {
	seed_entry_t *ent = get_seed_entry(pt, start);	
	show_seed_entry(pt, start);	
	if (!is_valid_entry(ent))
		return;
	
	if (!mode_build && ent->left != NO_ENTRY)
		__show_perfect_table_related(pt, ent->left);

	if (start == ent->right) {
		printf("a circle is detected!\n");
		return;
	}
	if (ent->right != NO_ENTRY)
		__show_perfect_table_related(pt, ent->right);
}

void show_perfect_table_related(perfect_table_t *pt, uint32_t start) {
	printf("[SHOW_PERFECT_TABLE_RELATED] START:%08x\n"
		   "=================================================================\n", 
		   start);
	__show_perfect_table_related(pt, start);
	printf("=================================================================\n");
}

static inline void set_fw_less_entry(seed_entry_t *entry, int is_fw_less) {
	if (is_fw_less)
		entry->flags = entry->flags | FLAG_FW_LESS;
	else
		entry->flags = entry->flags & (~FLAG_FW_LESS);
}

static inline void set_collision_entry(seed_entry_t *entry, int is_collision) {
	if (is_collision)
		entry->flags = entry->flags | FLAG_COLLISION;
	else
		entry->flags = entry->flags & (~FLAG_COLLISION);
}

static inline int set_multi_location(seed_entry_t *entry, uint32_t multi_loc) {
	if (multi_loc >= FLAG_MULTI_LOC_MAX)
		return -1;
	entry->flags = ((entry->flags & (~FLAG_MULTI_LOC_MASK)) | (multi_loc << FLAG_MULTI_LOC_SHIFT)); 
	return 0;
}

#define INIT_SEED_ENTRY(ent, loc, fw_less, collision) do { \
		(ent)->flags = ((fw_less) ? FLAG_FW_LESS : 0) \
					| ((collision) ? FLAG_COLLISION : 0); \
		(ent)->location = (loc); \
		(ent)->left = NO_ENTRY; \
		(ent)->right = NO_ENTRY; \
} while (0)

void __add_to_hash(perfect_table_t *pt, uint32_t loc, uint32_t key, int fw_less, int len) {
	uint32_t key_idx, new_idx, prev_idx;
	seed_entry_t *key_ent, *new_ent, *prev_ent;

	// find entry
	key_idx = key;
	key_ent = get_seed_entry(pt, key_idx);
	
	dbg_printf("%s: START seed: %s location: %8x key: %8x fw_less: %d\n", __func__,
		__get_seed_str(loc, len, pt, NULL), loc, key_idx, fw_less);
	
	// if tne entry is a collision entry (seed is not matched for the key), move the entry away
	if (is_collision_entry(key_ent)) {
		dbg_printf("%s: START seed: %s location: %8x key: %8x\n", __func__, 
			__get_seed_str(loc, len, pt, NULL), loc, key_idx);
		dbg_show_seed_entry(pt, key_idx);
		dbg_show_seed_entry(pt, new_idx);
		// MOVE ENTRY: do not increment statistics
		new_idx = get_empty_idx(pt, key_idx); // get the nearby empty entry.
		if (new_idx == NO_ENTRY) goto no_empty_entry;
		
		dbg_printf("%s: MOVE %08x -> %08x\n", __func__, key_idx, new_idx);
		dbg_show_perfect_table_related(pt, key_idx);
		
		new_ent = get_seed_entry(pt, new_idx);	
		memcpy(new_ent, key_ent, sizeof(seed_entry_t)); // copy values 
	
		/* check the chain */
		dbg_printf("%s: CHECK THE CHAIN OF COLLISION ENTRY[%08x]\n", __func__, key_idx);
		dbg_show_seed_entry(pt, key_idx);

		prev_idx = get_hash_idx_ent(pt, key_ent);
		prev_ent = get_seed_entry(pt, prev_idx);
		dbg_show_seed_entry(pt, prev_idx);
		fflush(stdout);	
		d_assert(!is_collision_entry(prev_ent), pt, prev_idx);
		while (prev_ent->right != key_idx && prev_ent->right != NO_ENTRY) {
			prev_idx = prev_ent->right;
			prev_ent = get_seed_entry(pt, prev_idx);
			dbg_show_seed_entry(pt, prev_idx);
		}
		d_assert(prev_ent->right == key_idx, pt, prev_idx);
		prev_ent->right = new_idx;
		d_assert(prev_idx != new_idx, pt, prev_idx);
		dbg_printf("%s: CHAIN [%08x].right -> [%08x]\n", __func__, prev_idx, new_idx);
	
		dbg_show_perfect_table_related(pt, new_idx);
		INIT_SEED_ENTRY(key_ent, NO_ENTRY, 0, 0);
		dbg_show_perfect_table_related(pt, key_idx);
		d_assert(!is_valid_entry(key_ent), pt, key_idx);
		total_moved_entry++;
	}

	// add entry
	//     if this is an invalid entry, just add
	//     if this is already valid, 
	//            if the chain already have the seed, add location (in the new or previous location entry)
	//            if the chain don't have the seed, add new seed entry to the chain

	if (!is_valid_entry(key_ent)) { // new entry (seed is matched for the key)
		// NEW ENTRY CASE#1: a seed entry without collision yet
		dbg_printf("%s: NEW MATCHED SEED ENTRY[%08x]\n", __func__, key_idx);
		INIT_SEED_ENTRY(key_ent, loc, fw_less, 0);
		dbg_show_seed_entry(pt, key_idx);

		pt->num_seed_used++;
		pt->num_seed_key++;
	} else {
		// find seed matching entry
		int matched = 0;
		dbg_printf("%s: FIND SEED MACHTING ENTRY from [%08x]\n", __func__, key_idx);
		dbg_show_perfect_table_related(pt, key_idx);	
		prev_idx = NO_ENTRY;
		prev_ent = NULL;
		new_idx = key_idx;
		new_ent = key_ent;
		while (new_idx != NO_ENTRY) {
			dbg_show_seed_entry(pt, new_idx);
			if (new_ent->right != NO_ENTRY)
				__builtin_prefetch(pt->ref_string + get_seed_entry(pt, new_ent->right)->location);
			matched = seedmatch_loc_to_loc(pt, loc, new_ent->location);
			if (matched != 0)
				break;
			prev_idx = new_idx;
			prev_ent = new_ent;
			new_idx = new_ent->right;
			new_ent = get_seed_entry(pt, new_idx);
		}
			
		dbg_printf("%s: FIND SEED MATCHING ENTRY from [%08x]: %s (%08x)\n", __func__, key_idx, new_idx != NO_ENTRY ? "succeed" : "fail", new_idx);
		dbg_show_seed_entry(pt, prev_idx);
		dbg_show_seed_entry(pt, new_idx);

		// if seed matching entry is not found, allocate new one, and exit
		if (matched == 0) {
			// note that entries between the original key_idx and prev_idx are used, since we use linear search for empty entry.
			// NEW ENTRY CASE#2: a seed entry with collision
			assert(new_idx == NO_ENTRY);
			new_idx = get_empty_idx(pt, prev_idx);
			if (new_idx == NO_ENTRY) goto no_empty_entry;
			d_assert(prev_idx != new_idx, pt, prev_idx);
			new_ent = get_seed_entry(pt, new_idx);
			dbg_printf("%s: NEW COLLISION SEED ENTRY[%08x]\n", __func__, new_idx);

			INIT_SEED_ENTRY(new_ent, loc, fw_less, 1);
			
			prev_ent->right = new_idx;
			d_assert(prev_idx != new_idx, pt, prev_idx);
			dbg_printf("%s: CHAIN [%08x].right -> [%08x]\n", __func__, prev_idx, new_idx);

			pt->num_seed_used++;
			dbg_show_perfect_table_related(pt, new_idx);
		} else { // if seed matching entry is found, add location
			uint32_t multi_loc;
			uint32_t n;
			build_loc_t *bloc;
			dbg_printf("%s: ADD LOCATION %8x to SEED ENTRY[%08x]\n", __func__, loc, new_idx);
			dbg_show_seed_entry(pt, new_idx);

			multi_loc = get_multi_location(new_ent);
			if (multi_loc == 0) {
				dbg_printf("%s: SEED ENTRY[%08x] is CHANGED to MULTI-LOCATION ENTRY\n", __func__, new_idx);
				
				multi_loc = build_loc_new();
				assert(multi_loc != 0);
				assert(is_valid_entry(new_ent));
				set_multi_location(new_ent, multi_loc); 
				assert(is_valid_entry(new_ent));
			} 
			assert(multi_loc != 0);
		
			build_loc_add_loc(multi_loc, loc, matched == 1 ? 0 : 1);
			dbg_show_seed_entry(pt, new_idx);
		}
	}

	total_added_entry++;

	return;

no_empty_entry:
	fprintf(stderr, "ERROR: cannot allocate a seed entry of perfect table. Is something wrong? or slack < 1?\n"
					"       seed_len: %u seq_len: %u #seed_entry: %u\n"
					"       #used_seed: %u #seed_key: %u #build_loc_entry: %u\n",
					pt->seed_len, pt->seq_len, pt->num_seed_entry,
					pt->num_seed_used, pt->num_seed_key, build_loc_n);
	exit(EXIT_FAILURE);
}

static inline void __swap(seed_entry_t *node_list, int a, int b) {
	seed_entry_t tmp;
	memcpy(&tmp, &node_list[a], sizeof(seed_entry_t));
	memcpy(&node_list[a], &node_list[b], sizeof(seed_entry_t));
	memcpy(&node_list[b], &tmp, sizeof(seed_entry_t));
}

/* refer to https://www.geeksforgeeks.org/quick-sort/ */
static void __quick_sort(seed_entry_t *node_list, int low, int high, perfect_table_t *pt) {
	int pivot, i, j;	
	seed_entry_t tmp;

	//printf("%s: low: %u high: %u\n", __func__, low, high);
	
	if (low >= high)
		return;

	/* partition */
	pivot = high;
	i = low - 1;
	
	for (j = low; j <= high -1; j++) {
		if (seedcmp_entries(pt, &node_list[j], &node_list[pivot]) < 0) {
			i++;
			__swap(node_list, i, j);
		}
	}
	__swap(node_list, i + 1, pivot);

	pivot = i + 1;
	
	//printf("%s: new pivot: %u\n", __func__, pivot);

	/* recursive call */
	__quick_sort(node_list, low, pivot - 1, pt);
	__quick_sort(node_list, pivot + 1, high, pt);
}

inline void quick_sort(seed_entry_t *node_list, int n, perfect_table_t *pt) {
	__quick_sort(node_list, 0, n - 1, pt);
}

/* Rebuild Perfect Table */
static inline void update_multi_loc(seed_entry_t *ent, uint32_t *multi_loc_map) {
	uint32_t multi_loc = get_multi_location(ent);
	if (multi_loc)
		set_multi_location(ent, multi_loc_map[multi_loc]);
}

static inline void update_multi_loc_list(seed_entry_t *ent, uint32_t *multi_loc_map, int n) {
	int i;
	for (i = 0; i < n; ++i)
		update_multi_loc(&ent[i], multi_loc_map);
}

static inline int get_children_list(perfect_table_t *pt, uint32_t root_idx, 
					uint32_t **__idx_list, seed_entry_t **__node_list, int *__num_list)
{
	int i, n;
	uint32_t idx;
	seed_entry_t *ent;
	uint32_t *idx_list;
	seed_entry_t *node_list;

	/* count the #children */
	n = 0;
	idx = root_idx; // start idx

	while (idx != NO_ENTRY) {
		ent = get_seed_entry(pt, idx);
		n++;
		idx = ent->right;
	}

	if (n > *__num_list) {
		*__idx_list = (uint32_t *) reallocarray(*__idx_list, n, sizeof(uint32_t));
		*__node_list = (seed_entry_t *) reallocarray(*__node_list, n, sizeof(seed_entry_t));
		*__num_list = n;
		assert((*__idx_list) != NULL && (*__node_list) != NULL);
	}

	idx_list = *__idx_list;
	node_list = *__node_list;

	/* copy the entries */
	n = 0;
	idx = root_idx;
	while (idx != NO_ENTRY) {
		ent = get_seed_entry(pt, idx);
		idx_list[n] = idx;
		memcpy(&node_list[n], ent, sizeof(seed_entry_t));
		n++;
		idx = ent->right;
	}

	return n;
}

void __convert_to_bst(perfect_table_t *pt, uint32_t *idx_list, int *idx_next, 
						seed_entry_t *node_list, uint32_t root_idx, int low, int high) {
	int mid;
	seed_entry_t *ent;

	if (low > high)
		return;

	/* find middle entry */
	mid = (low + high) / 2;

	/* copy the middle entry to root idx */
	ent = get_seed_entry(pt, root_idx);
	memcpy(ent, &node_list[mid], sizeof(seed_entry_t));

	/* we pop up the indexes first, since we want the children of an entry exists on the adjacent cache lines. */
	/* if left child exists, pop idx for it */
	if (mid > low) {
		ent->left = idx_list[*idx_next];
		(*idx_next)++;
	} else
		ent->left = NO_ENTRY;

	/* if right child exists, pop idx for it */
	if (mid < high) {
		ent->right = idx_list[*idx_next];
		(*idx_next)++;
	} else
		ent->right = NO_ENTRY;

	/* call for left child */
	if (ent->left != NO_ENTRY)
		__convert_to_bst(pt, idx_list, idx_next, node_list, ent->left, low, mid - 1);

	/* call for right child */
	if (ent->right != NO_ENTRY)
		__convert_to_bst(pt, idx_list, idx_next, node_list, ent->right, mid + 1, high);
}

void _convert_to_bst(perfect_table_t *pt, uint32_t *idx_list, seed_entry_t *node_list, int n) {
	int idx_next = 1; /* idx_list[0] is given as the first root_idx */
	__convert_to_bst(pt, idx_list, &idx_next, node_list, idx_list[0], 0, n - 1);
}

static inline void convert_to_bst(perfect_table_t *pt, uint32_t *idx_list, 
						seed_entry_t *node_list, int n) 
{
	int i;	
	seed_entry_t *ent;

	assert(n > 1);
	
	quick_sort(node_list, n, pt);

	_convert_to_bst(pt, idx_list, node_list, n);

	/* set collision flags properly */
	ent = get_seed_entry(pt, idx_list[0]);
	set_collision_entry(ent, 0);
	for (i = 1; i < n; ++i) {
		ent = get_seed_entry(pt, idx_list[i]);
		set_collision_entry(ent, 1);
	}
}

void rebuild_perfect_table_for_mapping(perfect_table_t *pt) {
	uint32_t *loc_table = NULL;
	uint32_t *multi_loc_map = NULL;
	uint32_t loc_n = 0, multi_loc_map_n = 0;
	uint32_t idx, next_idx, i;
	uint32_t pf_idx;
	build_loc_t *bloc;
	seed_entry_t *ent;
	uint32_t multi_loc;
	/* to boost rebalancing collision entries, 
	   while first path, we make chain of root entries with children using ent->left_idx.
	   Note that ent->left_idx is unused on building */
	uint32_t *idx_list = NULL;
	seed_entry_t *node_list = NULL;
	int num_children, num_list = 0;


	printf("[Rebuilding#1] build loc_table for %u seed entries\n", build_loc_n);
	fflush(stdout);
	loc_table = build_loc_to_loc_table(&loc_n, &multi_loc_map, &multi_loc_map_n);
	build_loc_free();
	pt->num_loc_entry = loc_n;
	pt->loc_table = loc_table;
	printf("[Rebuilding#1] done\n");
	printf("[Rebuilding#2] scan %u entries: set multi_loc and convert collision entries to BST\n", pt->num_seed_entry);
	fflush(stdout);

	for (pf_idx = 0; pf_idx < PREFETCH_DISTANCE && pf_idx < pt->num_seed_entry; ++pf_idx)
		__builtin_prefetch(get_seed_entry(pt, pf_idx));

	for (idx = 0; idx < pt->num_seed_entry; idx++) {
		if ((idx + 1) % 100000000 == 0) {
			printf("[Rebuilding#2] (%.1f%%) %u/%u entries\n",
						(float) (idx + 1) * 100 / pt->num_seed_entry, 
						idx + 1, pt->num_seed_entry);
			fflush(stdout);
		}

		if (pf_idx < pt->num_seed_entry)
			__builtin_prefetch(get_seed_entry(pt, pf_idx++));

		ent = get_seed_entry(pt, idx);
		if (!is_valid_entry(ent))
			continue;
		
		if (is_collision_entry(ent)) 
			/* the root entry does everything for its children */
			continue;

		if (ent->right == NO_ENTRY) {
			update_multi_loc(ent, multi_loc_map);
			continue;
		}

		num_children = get_children_list(pt, idx, &idx_list, &node_list, &num_list);

		update_multi_loc_list(node_list, multi_loc_map, num_children);
		
		/* make collision entries as balanced tree */
		convert_to_bst(pt, idx_list, node_list, num_children);
	}
	
	printf("[Rebuilding#2] done. #seed_entry: %u #loc_entry: %u\n", pt->num_seed_entry, pt->num_loc_entry);
	fflush(stdout);

	mode_build = 0; /* mode_build is related to multi_location */
	free(idx_list);
	free(node_list);
}

#if 0
static inline void add_region_to_hash(perfect_table_t *pt, int seed_len,
									  uint32_t start_loc, uint32_t end_loc) {
	uint32_t loc;
	if (end_loc - start_loc < seed_len)
		return;
	end_loc -= seed_len;
	for (loc = start_loc; loc < end_loc; ++loc)
		add_to_hash(pt, loc, seed_len);	
}
#endif


#define NUM_KEY_THREAD 8
#define NUM_LOC_PER_STEP 3000000
	
typedef struct loc_key_data_s {
		uint32_t key;
		int fw_less;
} loc_key_data_t;

typedef struct loc_key_s {
	int tid, num_thread;
	sem_t read_sem, write_sem;
	uint32_t start, end; // location
	int last;
	perfect_table_t *pt;
	bntann1_t *anns;
	int n_seqs;
#ifndef PERFECT_MATCH_IGNORE_HOLE
	bntamb1_t *ambs;
	int n_holes;
#endif
	loc_key_data_t data[NUM_LOC_PER_STEP];
} loc_key_t;


loc_key_t *loc_key[NUM_KEY_THREAD];

static inline void __calc_loc_key_set(perfect_table_t *pt, uint32_t loc, int len, loc_key_data_t *d) {
	d->fw_less = __compare_fw_rc(pt->ref_string + loc, len);
	d->key = __get_hash_idx_seed(pt, pt->ref_string + loc, d->fw_less);
}

static inline void __calc_loc_key_hole(loc_key_data_t *d) {
	d->fw_less = 0;
	d->key = NO_ENTRY;
}

static void *calc_loc_key(void *arg) {
	loc_key_t *loc_key = (loc_key_t *) arg;
	uint32_t loc, idx;
	uint32_t next = loc_key->tid * NUM_LOC_PER_STEP;
	perfect_table_t *pt = loc_key->pt;
	int seed_len = pt->seed_len;
	uint32_t seq_len = pt->seq_len;
	bntann1_t *anns = loc_key->anns;
	int n_seqs = loc_key->n_seqs;
	int seq_id;
#ifndef PERFECT_MATCH_IGNORE_HOLE
	bntamb1_t *ambs = loc_key->ambs;
	int n_holes = loc_key->n_holes;
	int hole_id;
#endif
	uint32_t end;

	seq_id = 0;
	hole_id = 0;
	while (next < seq_len) {
		sem_wait(&loc_key->write_sem);
		
		loc_key->start = next;
		loc_key->end = next + NUM_LOC_PER_STEP;
		if (loc_key->end >= seq_len) {
			loc_key->end = seq_len;
			loc_key->last = 1;
		}

		idx = 0;
#if DONT_INCLUDE
#ifdef PERFECT_MATCH_IGNORE_HOLE
		if (loc_key->last) {
			uint32_t end = loc_key->end - seed_len;
			for (loc = loc_key->start; loc < end; ++loc)
				__calc_loc_key_set(pt, loc, seed_len, &(loc_key->data[idx++]));
			for ( ; loc < loc_key->end; ++loc)
				__calc_loc_key_hole(&(loc_key->data[idx++]));
		} else {
			for (loc = loc_key->start; loc < loc_key->end; ++loc)
				__calc_loc_key_set(pt, loc, seed_len, &(loc_key->data[idx++]));
		}
#else
		hole_id = 0;
		while (hole_id < n_holes && loc_key->start < ambs[hole_id].offset + ambs[hole_id].len)
			hole_id++;

		for (loc = loc_key->start; loc < loc_key->end; ++loc) {
			if (hole_id < n_holes && loc >= ambs[hole_id].offset + ambs[hole_id].len)
				hole_id++;

			if (hole_id < n_holes && loc > ambs[hole_id].offset - seed_len)
				/* NOTE: if hole_id < n_holes, loc < ambs[hole_id].offset + ambs[hole_id].len by above loop and if-stmt */
				__calc_loc_key_hole(&(loc_key->data[idx++]));
			else if (loc_key->last && loc > loc_key->end - seed_len)
				__calc_loc_key_hole(&(loc_key->data[idx++]));
			else
				__calc_loc_key_set(pt, loc, seed_len, &(loc_key->data[idx++]));
		}
#endif
#endif
		loc = loc_key->start;
		/* find the appropriate ann entry */
		while (seq_id < n_seqs 
					&& loc >= (anns[seq_id].offset + anns[seq_id].len))
			seq_id++;

		/* find the appropriate amb entry */
		while (hole_id < n_holes 
					&& loc >= ambs[hole_id].offset + ambs[hole_id].len)
			hole_id++;

		/* calculate loc keys. */
		while (loc < loc_key->end) {
			if (seq_id < n_seqs
						&& loc >= (anns[seq_id].offset + anns[seq_id].len))
				seq_id++;

			if (hole_id < n_holes 
						&& loc >= ambs[hole_id].offset + ambs[hole_id].len)
				hole_id++;

			if (hole_id < n_holes 
					&& loc > ambs[hole_id].offset - seed_len) {
				end = ambs[hole_id].offset + ambs[hole_id].len;
				if (end > loc_key->end)
					end = loc_key->end;
				while (loc < end) {
					__calc_loc_key_hole(&(loc_key->data[idx++]));
					loc++;
				}
			} else if (seq_id < n_seqs
							&& loc > anns[seq_id].offset + anns[seq_id].len - seed_len) {
				end = anns[seq_id].offset + anns[seq_id].len;
				if (end > loc_key->end)
					end = loc_key->end;
				while (loc < end) {
					__calc_loc_key_hole(&(loc_key->data[idx++]));
					loc++;
				}
			} else {
				end = anns[seq_id].offset + anns[seq_id].len - seed_len + 1;
				if (hole_id < n_holes && end > ambs[hole_id].offset - seed_len)
					end = ambs[hole_id].offset - seed_len + 1;
				if (end > loc_key->end)
					end = loc_key->end;
				while (loc < end) {
					__calc_loc_key_set(pt, loc, seed_len, &(loc_key->data[idx++]));
					loc++;
				}
			}
		}

		sem_post(&loc_key->read_sem);

		next += loc_key->num_thread * NUM_LOC_PER_STEP;
	}

	return (void *) 0;
}

static void add_to_hash(perfect_table_t *pt, loc_key_t *loc_key_list[], int list_len) {
	loc_key_t *loc_key;
	loc_key_data_t *data;
	int i = 0;
	int done = 0;
	int seed_len = pt->seed_len;
	uint32_t loc, idx, key;
	uint32_t pf_idx, idx_len;
	uint32_t pf_key;
	int i_next;
	uint32_t pf_next_idx;
	//uint32_t pf_val = 0;
	//seed_entry_t *pf_ent;
	
	while (!done) {
		loc_key = loc_key_list[i];
		sem_wait(&loc_key->read_sem);

		pf_next_idx = 0;
		i_next = (i + 1) % NUM_KEY_THREAD;

		/* prefetch */
		idx_len = loc_key->end - loc_key->start;
		for (pf_idx = 0; pf_idx < PREFETCH_DISTANCE; ++pf_idx) {
			pf_key = loc_key->data[pf_idx].key;
			if (pf_key != NO_ENTRY) {
				__builtin_prefetch(get_seed_entry(pt, pf_key));
				//pf_ent = get_seed_entry(pt, pf_key);
				//if (is_valid_entry(pf_ent)) /* prefetch the seed_entry */
				//	pf_val ^= *(pt->ref_string + pf_ent->location); /* prefetch the seed string */
			}
		}

		for (idx = 0, loc = loc_key->start; loc < loc_key->end; ++idx, ++loc) {
			key = loc_key->data[idx].key;
			if (key != NO_ENTRY)
				__add_to_hash(pt, loc, key, loc_key->data[idx].fw_less, seed_len);
			
			/* prefetch */
			if (pf_idx < idx_len) { 
				pf_key = loc_key->data[pf_idx++].key;
				if (pf_key != NO_ENTRY) {
					__builtin_prefetch(get_seed_entry(pt, pf_key));
					//pf_ent = get_seed_entry(pt, pf_key);
					//if (is_valid_entry(pf_ent)) /* prefetch the seed_entry */
					//	pf_val ^= *(pt->ref_string + pf_ent->location); /* prefetch the seed string */
				}
			} else {
				__builtin_prefetch(&loc_key_list[i_next]->data[pf_next_idx++]);
				/* after few prefetches, next-line-prefetcher works in the main loop */
			}
		}

		done = loc_key->last;
		sem_post(&loc_key->write_sem);
		
		show_perfect_table_stat(pt, loc);

		i = (i + 1) % list_len;
	}

	return;	
}


#define NUM_NEW_SEED_TABLE_THREAD 8

typedef struct new_seed_table_s {
	uint32_t start, end;
	seed_entry_t *table;
} new_seed_table_t;

static void *__new_seed_table(void *arg) {
	new_seed_table_t *nst = (new_seed_table_t *) arg;
	seed_entry_t *table = nst->table;
	size_t i, step, next;
	size_t start = nst->start, end = nst->end;

	step = 4096 / sizeof(seed_entry_t);
	assert(start % step  == 0);

	printf("new_seed_table: start: %12ld end: %12ld\n", start, end);

	if (start == end)
		return NULL;

	/* init the first entry */
	INIT_SEED_ENTRY(&table[start], NO_ENTRY, 0, 0);

	if (end - start < step) {
		for (i = start + 1; i < end; ++i) 
			memcpy(&table[i], &table[start], sizeof(seed_entry_t));
		return NULL;
	}

	/* init a page by copying the first entry */
	next = start + step;
	for (i = start + 1; i < next; ++i) 
		memcpy(&table[i], &table[start], sizeof(seed_entry_t));

	/* init pages as much as possible */
	next = end - (end % step);
	for ( ; i < next; i += step) 
		memcpy(&table[i], &table[start], sizeof(seed_entry_t) * step);

	/* init remaining entries */
	if (i < end)
		memcpy(&table[i], &table[start], sizeof(seed_entry_t) * (end - i));

	return NULL;
}

seed_entry_t *new_seed_table(size_t nelem) {
	seed_entry_t *table = (seed_entry_t *) malloc(nelem * sizeof(seed_entry_t));
	new_seed_table_t nst[NUM_NEW_SEED_TABLE_THREAD];
	pthread_t nst_thread[NUM_NEW_SEED_TABLE_THREAD];
	size_t start, end, per_thread, per_page, num_page;
	int i, num_thread = 1;

	if (!table) return NULL;

	per_page = (4096 + sizeof(seed_entry_t) - 1) / sizeof(seed_entry_t);
	num_page = (nelem + per_page - 1) / per_page; 
	per_thread = ((num_page + NUM_NEW_SEED_TABLE_THREAD - 1) / NUM_NEW_SEED_TABLE_THREAD) * per_page;

	start = 0;
	end = 0;
	for (i = 0; i < NUM_NEW_SEED_TABLE_THREAD; ++i) {
		start = end;
		end = start + per_thread;
		if (end > nelem)
			end = nelem;
		assert(end >= start);
		nst[i].start = start;
		nst[i].end = end;
		nst[i].table = table;
		if (end == nelem) { 
			num_thread = i + 1;
			break;
		}
	}

	for (i = 0; i < num_thread; ++i)
		pthread_create(&nst_thread[i], NULL, __new_seed_table, &nst[i]);

	for (i = 0; i < num_thread; ++i)
		pthread_join(nst_thread[i], NULL);

	return table;
}

int __perfect_build_index(const char *pt_fn, uint8_t *ref_string,
							int64_t seq_len, double slack, int seed_len,
							bntann1_t *anns, int32_t n_seqs,
							bntamb1_t *ambs, int32_t n_holes) {

	int64_t num_seed_entry;
	perfect_table_t pt;
	int32_t hole;
	int64_t start_loc, end_loc;
	bntamb1_t *amb;
	FILE *fp;
	int i;
	pthread_t key_thread[NUM_KEY_THREAD];
	cpu_set_t cpumask;
	struct timeval t_beg, t_end;
	uint32_t *loc_table;
	seed_entry_t *seed_table;
	
	assert(sizeof(perfect_table_t) % 64 == 0);

	/* initialize global statistics */
	total_added_entry = 0;
	total_hole_entry = 0;
	total_moved_entry = 0;

	CPU_ZERO(&cpumask);
	/* super set */
	for (i = 0; i < NUM_KEY_THREAD; ++i)
		CPU_SET(i, &cpumask);
	for (i = 0; i < NUM_NEW_SEED_TABLE_THREAD; ++i)
		CPU_SET(i, &cpumask);
	sched_setaffinity(0, sizeof(cpumask), &cpumask);	

	pt.seed_len = seed_len;
	if (seq_len >= UINT32_MAX) {
		fprintf(stderr, "ERROR: perfect match does not support genome reference whose sequence length exceeds %u\n", UINT32_MAX);
		exit(EXIT_FAILURE);
	}

	num_seed_entry = (uint64_t) ((double)seq_len * slack);
	if (num_seed_entry > UINT32_MAX) {
		fprintf(stderr, "ERROR: the number of seed entry should be less than %u. The slack should be decreased. (the maximum slack is %f)\n", UINT32_MAX, (double) UINT32_MAX / (double) seq_len);
		exit(EXIT_FAILURE);
	}

	pt.num_loc_entry = 0;
	pt.num_seed_entry = (uint32_t) num_seed_entry;
#ifdef MEMSCALE
	pt.num_seed_load = (uint32_t) num_seed_entry;
#endif
	pt.ref_string = ref_string;
	pt.loc_table = NULL;
	build_loc_init();
	printf("Allocate memory for seed entries of perfect table (%.3fGB)\n",
			(double) pt.num_seed_entry * sizeof(seed_entry_t) / (1024*1024*1024));
	fflush(stdout);
	gettimeofday(&t_beg, NULL);
	pt.seed_table = new_seed_table(pt.num_seed_entry);
	gettimeofday(&t_end, NULL);
	printf("allocation_time: %.3fs\n", t_end.tv_sec - t_beg.tv_sec + (t_end.tv_usec - t_beg.tv_usec) / 1e6);
	pt.seq_len = (uint32_t) seq_len;
	pt.num_seed_used = 0;
	pt.num_seed_key = 0;
	printf("Build perfect table seq_len: %ld\n", seq_len);
	fflush(stdout);
	
	for (i = 0; i < NUM_KEY_THREAD; ++i) {
		loc_key[i] = (loc_key_t *) calloc(1, sizeof(loc_key_t));
		loc_key[i]->tid = i;
		loc_key[i]->num_thread = NUM_KEY_THREAD;
		sem_init(&loc_key[i]->read_sem, 0, 0);
		sem_init(&loc_key[i]->write_sem, 0, 1);
		loc_key[i]->pt = &pt;
		loc_key[i]->anns = anns;
		loc_key[i]->n_seqs = n_seqs;
#ifndef PERFECT_MATCH_IGNORE_HOLE
		loc_key[i]->ambs = ambs;
		loc_key[i]->n_holes = n_holes;
#endif
	}

	for (i = 0; i < NUM_KEY_THREAD; ++i)
		pthread_create(&key_thread[i], NULL, calc_loc_key, loc_key[i]);

	/* main_thread */
	add_to_hash(&pt, loc_key, NUM_KEY_THREAD);

	for (i = 0; i < NUM_KEY_THREAD; ++i)
		pthread_join(key_thread[i], NULL);
	
	printf("Re-build perfect table for mapping\n");
	fflush(stdout);
	
	rebuild_perfect_table_for_mapping(&pt);
	
	printf("Write perfect table to %s\n", pt_fn);
	fflush(stdout);

	/* reset some runtime specific values */
#ifdef MEMSCALE
	pt.num_seed_load = 0;
#endif
	pt.ref_string = NULL;
	loc_table = pt.loc_table;
	pt.loc_table = NULL;
	seed_table = pt.seed_table;
	pt.seed_table = 0;

	fp = xopen(pt_fn, "wb");
	err_fwrite(&pt, sizeof(perfect_table_t), 1, fp);
	err_fwrite(loc_table, sizeof(uint32_t), pt.num_loc_entry, fp);
	err_fwrite(seed_table, sizeof(seed_entry_t), pt.num_seed_entry, fp);
	err_fflush(fp);
	err_fclose(fp);
	free(seed_table);
	free(loc_table);
	printf("Done\n");
	fflush(stdout);
	
	return 0;
}

int perfect_build_index(const char *prefix, int seed_len, double slack)
{
	clock_t t;
	int64_t seq_len;
	int32_t n_holes;
	bntamb1_t *ambs;
	int32_t n_seqs;
	bntann1_t *anns;

  	uint8_t *ref_string;

	char file_name[PATH_MAX];

	anns = ann_restore(prefix, &seq_len, &n_seqs);
	ambs = amb_restore(prefix, &seq_len, &n_holes);
	load_ref_string(prefix, &ref_string);
	
	snprintf_s_si(file_name, PATH_MAX, "%s.perfect.%d", prefix, seed_len);
	__perfect_build_index(file_name, ref_string, seq_len, slack, seed_len, 
							anns, n_seqs, ambs, n_holes);
	_mm_free(ref_string);
	free(ambs);

	return 0;
}

static inline int64_t *__array64_fit(int64_t *ar, ssize_t *size, ssize_t target) {
	/* Assume target > *size */
	int64_t after, i;
	int64_t *ret;

	if (*size == 0) {
		ssize_t s = 16;
		while (target > s)
			s <<= 1;
		*size = s;
		return (int64_t *)calloc(s, sizeof(int64_t));
	}

	after = *size;
	while (target > after)
		after <<= 1;

	ret = (int64_t *)calloc(after, sizeof(int64_t));
	if (!ret)
		return NULL;
	for (i = 0; i < *size; ++i) 
		ret[i] = ar[i];
	*size = after;
	return ret;
}

/* increase the size of 64-bit entry array to include @target number of entries */
static inline void array64_fit(int64_t **ar, ssize_t *size, ssize_t target) {
	if (target < *size)
		return;

	*ar = __array64_fit(*ar, size, target);
	return;
}

void __show_distribution(int64_t *ar, ssize_t size, const char *head, const char *x_name, const char *y_name, int64_t total) {
	int64_t i;

	for (i = 0; i < size; ++i) {
		if (ar[i] == 0) continue;
		printf("[%s] %s: %6ld %s: %16ld (%5.2f)\n",
			head, x_name, i, y_name, ar[i], (float) ar[i] * 100 / total);
	}
}

struct perfect_table_stat {
	int64_t total_key; /* # key entries */
	int64_t total_seed; /* # seed entries */
	int64_t total_loc; /* # location entries */
	int64_t total_valid; /* # valid entries */
	int64_t total_unique; /* no-next-seed and no-multiple-location entries */
								// for no-next-seed only, see @num_seed_dist_key[1] 
								// for no-multiplie-location, see @num_loc_dist_seed[1] 
	
	int64_t *cont_valid_dist;
	int64_t *cont_invalid_dist;
	int64_t *depth_dist_seed;    /* # cache lines to retrieve information related a seed */
	int64_t *max_depth_dist_key; /* max depth per key.... the longest hop entry from the key */
	int64_t *num_loc_dist_seed;  /* # locations per seed */
	int64_t *num_loc_dist_key;   /* # locations per key */
	int64_t *num_seed_dist_key;  /* # seed per key */
	int64_t *range_dist_key;         /* hash table entry range width per key */

	/* memory allocation for statistics above */
	int64_t m_cont_valid_dist;
	int64_t m_cont_invalid_dist;
	int64_t m_depth_dist_seed;
	int64_t m_max_depth_dist_key;
	int64_t m_num_loc_dist_seed; 
	int64_t m_num_loc_dist_key; 
	int64_t m_num_seed_dist_key;    
	int64_t m_range_dist_key;        

	/* temporal variables */
	int64_t cont_valid, cont_invalid, cont_prev;
	int64_t num_seed; 
	int64_t depth;    /* num_seed + num location entry */
	int64_t max_depth;	
	int64_t num_loc;  /* num locations for a seed */
	int64_t num_loc_key; /* num location for seeds attached to a key */
	int64_t min_idx, max_idx; /* for range distribution */
	int64_t range;
};

void __stat_perfect_table(struct perfect_table *pt, int64_t idx, struct perfect_table_stat *s /* stat */) {
	seed_entry_t *ent = get_seed_entry(pt, idx);
	int64_t depth;
	int32_t multi_loc;

	s->depth++;

	if (idx > s->max_idx || (s->max_idx > s->min_idx && idx < s->min_idx))
		s->max_idx = idx;

	s->num_seed++;

	multi_loc = get_multi_location(ent);
	if (!multi_loc) {
		s->num_loc = 1;
		depth = 1;
	} else {
		uint64_t start = (uint64_t) &pt->loc_table[multi_loc];
		uint64_t end;

		s->num_loc = pt->loc_table[multi_loc];
		end = start + sizeof(int64_t) * (s->num_loc + 1);
		
		/* cache-line aligned region */
		start = start & (~0x3fULL);
		end = (end & 0x3fULL) ? ((end & (~0x3fULL)) + 0x40) : end;
		depth = 1 + (int64_t) (end - start) / 64;
	}
	
	depth += s->depth;
	if (depth > s->max_depth)
		s->max_depth = depth;
		
	s->num_loc_key += s->num_loc;

	/* update distribution for this seed */
	array64_fit(&s->depth_dist_seed, &s->m_depth_dist_seed, depth);
	s->depth_dist_seed[depth]++;

	array64_fit(&s->num_loc_dist_seed, &s->m_num_loc_dist_seed, s->num_loc);
	s->num_loc_dist_seed[s->num_loc]++;

	/* go to children */
	if (ent->left != NO_ENTRY)
		__stat_perfect_table(pt, ent->left, s);
	if (ent->right != NO_ENTRY)
		__stat_perfect_table(pt, ent->right, s);

	s->depth--;
}

void stat_perfect_table(struct perfect_table *pt) {
	struct perfect_table_stat s;
	memset(&s, 0, sizeof(perfect_table_stat));

	int64_t idx = 0, seed_idx, loc_idx;
	seed_entry_t *ent, *seed_ent, *loc_ent;

	int64_t i, prev;

	for (idx = 0; idx < pt->num_seed_entry; idx++) {
		if (idx % 10000000 == 0) {
			fprintf(stderr, "[progress] (%.2f%%) idx: %ld total_valid: %ld\n", (float) idx * 100 / pt->num_seed_entry, idx, s.total_valid);
			fflush(stderr);
		}

		ent = get_seed_entry(pt, idx);
		//printf("[IDX] "); show_seed_entry(pt, idx);
		if (!is_valid_entry(ent)) {
			if (s.cont_invalid == 0 && s.cont_valid > 0) {
				/* end of continued valid entries */
				//printf("CONTINUOUS VALID   ENTRIES %16ld ~ %16ld => %16ld\n", s.cont_prev, idx - 1, s.cont_valid);
				//fflush(stdout);
				array64_fit(&s.cont_valid_dist, &s.m_cont_valid_dist, s.cont_valid);
				s.cont_valid_dist[s.cont_valid]++;
				s.cont_valid = 0;
				s.cont_prev = idx;
			}
			s.cont_invalid++;
			continue;
		}
			
		if (s.cont_valid == 0 && s.cont_invalid > 0) {
			/* end of continued invalid entries */
			//printf("CONTINUOUS INVALID ENTRIES %16ld ~ %16ld => %16ld\n", s.cont_prev, idx - 1, s.cont_invalid);
			//fflush(stdout);
			array64_fit(&s.cont_invalid_dist, &s.m_cont_invalid_dist, s.cont_invalid);
			s.cont_invalid_dist[s.cont_invalid]++;
			s.cont_invalid = 0;
			s.cont_prev = idx;
		}
		s.cont_valid++;

		s.total_valid++;
		if (!is_collision_entry(ent)) {
			s.total_key++;
			s.total_seed++;
			if (!get_multi_location(ent) && ent->left == NO_ENTRY && ent->right == NO_ENTRY)
				s.total_unique++;
		} else { /* is_collision_entry(ent) == true */
			s.total_seed++;
			continue; /* further analysis on collision entry had done at the root entry */
		}

		/* init temporal values */
		s.min_idx = idx;
		s.max_idx = idx;
		s.num_seed = 0;
		s.num_loc_key = 0;
		s.max_depth = 0;
		s.depth = 0;

		/* traverse seed entries */
		__stat_perfect_table(pt, idx, &s);

		/* update stat for a key */
		array64_fit(&s.num_seed_dist_key, &s.m_num_seed_dist_key, s.num_seed);
		s.num_seed_dist_key[s.num_seed]++;
		
		array64_fit(&s.max_depth_dist_key, &s.m_max_depth_dist_key, s.max_depth);
		s.max_depth_dist_key[s.max_depth]++;
		
		array64_fit(&s.num_loc_dist_key, &s.m_num_loc_dist_key, s.num_loc_key);
		s.num_loc_dist_key[s.num_loc_key]++;
		
		s.range = s.max_idx >= s.min_idx ? s.max_idx + 1 - s.min_idx
					: (pt->num_seed_entry - s.min_idx + s.max_idx);
		array64_fit(&s.range_dist_key, &s.m_range_dist_key, s.range);
		s.range_dist_key[s.range]++;
	}

	if (s.cont_valid > 0) {
		//printf("CONTINUOUS VALID   ENTRIES %16ld ~ %16ld => %16ld\n", s.cont_prev, idx - 1, s.cont_valid);
		array64_fit(&s.cont_valid_dist, &s.m_cont_valid_dist, s.cont_valid);
		s.cont_valid_dist[s.cont_valid]++;
	} else if (s.cont_invalid > 0) {
		//printf("CONTINUOUS INVALID ENTRIES %16ld ~ %16ld => %16ld\n", s.cont_prev, idx - 1, s.cont_invalid);
		array64_fit(&s.cont_invalid_dist, &s.m_cont_invalid_dist, s.cont_invalid);
		s.cont_invalid_dist[s.cont_invalid]++;
	}

	printf("STATISTICS OF PERFECT TABLE\n");
	fflush(stdout);
	printf("total_valid: %16ld (%.2f%%)\n", s.total_valid, (float) s.total_valid * 100 / pt->num_seed_entry);
	printf("total_key: %16ld (%.2f%%) (%.2f%%)\n", s.total_key, (float) s.total_key * 100 / pt->num_seed_entry, (float) s.total_key * 100 / s.total_valid);
	printf("total_seed: %16ld (%.2f%%) (%.2f%%)\n", s.total_seed, (float) s.total_seed * 100 / pt->num_seed_entry, (float) s.total_seed * 100 / s.total_valid);
	printf("total_unique: %16ld (%.2f%%) (%.2f%%)\n", s.total_unique, (float) s.total_unique * 100 / pt->num_seed_entry, (float) s.total_unique * 100 / s.total_valid);
	printf("SEED STATISTICS=============================================================\n");
	printf("\n");	
	__show_distribution(s.depth_dist_seed, s.m_depth_dist_seed, "DEPTH_SEED", "depth", "seed", s.total_seed);
	printf("\n");	
	__show_distribution(s.num_loc_dist_seed, s.m_num_loc_dist_seed, "NUM_LOCATIONS", "num_loc", "seed", s.total_seed);
	printf("\n");	
	printf("KEY STATISTICS==============================================================\n");
	__show_distribution(s.max_depth_dist_key, s.m_max_depth_dist_key, "MAX_DEPTH", "max_depth", "key", s.total_key);
	printf("\n");	
	__show_distribution(s.num_loc_dist_key, s.m_num_loc_dist_key, "NUM_LOCATIONS", "num_loc", "key", s.total_key);
	printf("\n");	
	__show_distribution(s.num_seed_dist_key, s.m_num_seed_dist_key, "NUM_SEED", "width", "key", s.total_key);
	printf("\n");	
	__show_distribution(s.range_dist_key, s.m_range_dist_key, "RANGE", "width", "key", s.total_key);
	printf("\n");
	printf("TABLE STATISTICS==============================================================\n");
	__show_distribution(s.cont_valid_dist, s.m_cont_valid_dist, "CONT_VALID", "width", "count", pt->num_seed_entry);
	printf("\n");
	__show_distribution(s.cont_invalid_dist, s.m_cont_invalid_dist, "CONT_INVALID", "width", "count", pt->num_seed_entry);
	fflush(stdout);
}

perfect_table_t *load_perfect_table(const char *prefix, int len, uint8_t **ref_string, FMI_search *fmi);

void display_perfect_table_stat(char *prefix, int seed_len) {
  	uint8_t *ref_string;
	perfect_table_t *pt;

	load_ref_string(prefix, &ref_string);
	load_perfect_table(prefix, seed_len, &ref_string, NULL); 
	pt = perfect_table;
	show_perfect_table(pt);

	printf("Statistics of perfect table: start\n");
	fflush(stdout);
	stat_perfect_table(pt);
	printf("Statistics of perfect table: done\n");
	fflush(stdout);
}

void usage_perfect_index() {
	fprintf(stderr, "Usage: bwa-mem2 perfect-index [-l seed_length] [-s slack] <prefix>\n");
	fprintf(stderr, "       -s (float) ==> the hash table will have (slack) * (length of reference sequence) entries\n");
}

int perfect_index(int argc, char *argv[]) // the "perfect-index" command
{
	int c;
	int seed_len = -1;
	double slack = 1.1;
	int opt_display_stat = 0;
	char *prefix = 0, *str;
	while ((c = getopt(argc, argv, "l:s:d")) >= 0) {
		if (c == 'l') {
			seed_len = atoi(optarg);
			if (seed_len <= 0) {
				fprintf(stderr, "ERROR: the seed length should be larger than 0, but %d is given.\n", seed_len);
				return -1;
			}
		} else if (c == 's') slack = atof(optarg);
		else if (c == 'd') opt_display_stat = 1;
		else {
			usage_perfect_index();
			return -1;
		}	
	}
			
	if (seed_len <= 0) {
		fprintf(stderr, "ERROR: the seed length must be given.\n");
		usage_perfect_index();
		return -1;
	}

	if (optind == argc || optind + 1 > argc) {
		usage_perfect_index();
		return -1;
	}
	
	if (opt_display_stat) {
		display_perfect_table_stat(argv[optind], seed_len);
		return 0;
	}

	mode_build = 1;
	perfect_build_index(argv[optind], seed_len, slack);
	mode_build = 0;
	return 0;
}
#endif
