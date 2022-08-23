/* Contact: Changdae Kim <cdkim@etri.re.kr> */

#ifndef BWT_PERFECT_H
#define BWT_PERFECT_H

#ifdef PERFECT_MATCH

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include "utils.h"

/* for development, do normal operations even if the perfect match is found */
//#define DO_NORMAL

/* for development, MUST USE WITH DO_NORMAL */
/* MUST COMPILE USING $ make perfect=1 arch=sse42 rwopt=1 -j */
//#define PRINT_PERFECT_AND_REG

/* for debug, ignore hole is needed */
//#define PERFECT_MATCH_IGNORE_HOLE

/* Note that the exact value requires information of all non-perfect matches, which we don't want to find.
 * FYI, for unknown mapq, mapq is 255 by definition.
 */
#define MAPQ_PERFECT_MATCH 60 

/* The original bwa-mem2 gives $(l_seq)M as cigar string for perfectly matched reads.
 * If you want to distinguish perfect matched reads to others, you can use 'm' instead.
 */
#define PERFECT_MATCH_CIGAR 'M'
//#define PERFECT_MATCH_CIGAR 'm'

/* seed_entry_t */
#define FLAG_FW_LESS		0x1
#define FLAG_COLLISION      0x2
#define NO_ENTRY            UINT32_MAX

/* bseq1_perfect_t */ 
#define FLAG_VALID			0x1
#define FLAG_RC				0x2 /* RC matched entry */

/* both of seed_entry_t and bseq1_perfect_t */
#define FLAG_MULTI_LOC_SHIFT 2
#define FLAG_MULTI_LOC_MASK (0xFFFFFFFF ^ ((0x1 << FLAG_MULTI_LOC_SHIFT) - 1))
#define FLAG_MULTI_LOC_MAX  ((1 << (32 - FLAG_MULTI_LOC_SHIFT)) - 1)

/* SEED TABLE 
 *
 * for perfect matching, we do not concatenate RC-mapping. Thus, 4-byte location and #entries are enough.
 *
 * a seed_entry includes both of FW and RC locations.
 * the key of a seed_entry is the hashed value of alphabetically smaller string of FW and RC.
 */
typedef struct {
	uint32_t flags; /* [0]: is_fw_less: the seed at @location is alphabetically less 
										than its reverse complemented one.
					   [1]: is_collision: hash value of the seed at @location 
					   					  differs from the index of this entry.
					   [2:31]: if multi-location entry, start_index in loc_table. 
					           Otherwise, 0. Don't use 0th element in loc_table.
					   NOTE: for human genome with 150-bp, 99.09% of seed has one location.
					         30-bit start index is likely to enough. */
	uint32_t location; /* If NO_ENTRY (UINT32_MAX), this is invalid entry */
	/* binary search tree for collision entries.
	   UINT32_MAX indicates no child.
	   While index-building, the chain is a singly-linked list, and only @right is used. */
	uint32_t left;
	uint32_t right;
} seed_entry_t;

#define LOC_MANY 256 /* this should be less than (1<<15) */
/* LOCATION TABLE
 * a location table is an array of uint32_t.
 *
 * loc[0]: unused. use as a NULL pointer.
 *
 * in location table, FW means the forward direction of seed_entry_t->location.
 *
 * CASE1) if MSB of the first entry (pointed by flags[2:13] in the seed entry) is 0,
 *        this is the case that #FW entries < LOC_MANY and #RC entries < LOC_MANY.
 *        - first[31:16] and first[15:0] is the number of FW and RC entries respectively.
 *        - the following entires are the locations, first FW entries, and then RC entries.
 *
 * CASE2) if MSB of the first entry (pointed by flags[2:13] in the seed entry) is 1,
 *        this is the case that #FW entries >= LOC_MANY or #RC entries >= LOC_MANY.
 *        - first[30:0] indicates the second entry in the location table.
 *        - the second entry is the number of FW entries.
 *        - the next of second entry is the number of RC entries.
 *        - the following entires are the locations, first FW entries, and then RC entries.
 */

typedef union {
	struct {
		uint32_t flags;
		uint32_t location; /* if len != seed_len, the first matched location */
	};
	uint64_t exist; /* if perfect_matching exists,
						at least, valid bit in flags is 1 */
} bseq1_perfect_t;

#define is_valid_entry(ent) ((ent)->location != NO_ENTRY)
#define is_fw_less_entry(ent)  (((ent)->flags & FLAG_FW_LESS) != 0)
#define is_collision_entry(ent)  (((ent)->flags & FLAG_COLLISION) != 0)
#define __is_rc_matched(flags)  ((flags & FLAG_RC) != 0)
#define is_rc_matched(ent)  (((ent)->flags & FLAG_RC) != 0)
#define is_hash_matched_entry(ent) ((ent) ? (is_valid_entry(ent) && (!is_collision_entry(ent))) \
										  : 0)

static inline int __get_num_location(uint32_t flags, uint32_t *loc_table) {
	uint32_t multi_loc = flags >> FLAG_MULTI_LOC_SHIFT;
	uint32_t first;
	if (multi_loc == 0)
		return 1;
	first = loc_table[multi_loc];

	if (first & 0x80000000) {
		uint32_t second = first & 0x7FFFFFFF;
		return 1 + loc_table[second] + loc_table[second + 1];
	} else {
		return 1 + ((first >> 16) & 0xFFFF) + (first & 0xFFFF);
	}
}

static inline uint32_t __get_multi_location(uint32_t flags) {
	return flags >> FLAG_MULTI_LOC_SHIFT; 
}

static inline uint32_t get_multi_location(seed_entry_t *entry) { 
	return entry->flags >> FLAG_MULTI_LOC_SHIFT; 
}

#define GET_MULTI_FW_AND_RC(loc_table, multi_loc, num_fw, loc_fw, num_rc, loc_rc) \
	do { \
		int ____many = loc_table[multi_loc] & 0x80000000 ? 1 : 0; \
		uint32_t ____start = ____many == 0 ? multi_loc : loc_table[multi_loc] & 0x7FFFFFFF; \
		if (____many == 0) { \
			num_fw = (loc_table[____start] >> 16) & 0xFFFF; \
			num_rc = loc_table[____start] & 0xFFFF; \
			loc_fw = &loc_table[____start + 1]; \
			loc_rc = &loc_table[____start + 1 + num_fw]; \
		} else { \
			num_fw = loc_table[____start]; \
			num_rc = loc_table[____start + 1]; \
			loc_fw = &loc_table[____start + 2]; \
			loc_rc = &loc_table[____start + 2 + num_fw]; \
		} \
	} while (0)

typedef struct __attribute__ ((__packed__)) perfect_table {
	int seed_len; // # alphabets in a seed
	uint32_t num_loc_entry; /* the size of location table */
	uint32_t num_seed_entry; // # entries of seed table
#ifdef MEMSCALE
	uint32_t num_seed_load;
#else
	uint32_t __dummy_memscale; /* for alignment */
#endif

	uint8_t *ref_string;

	/* location table for multi-location entries.
	   the start index entry has the number of locations,
	   and the followings have the locations. 
	   Don't use 0th entry of loc_table. It indicates no more locations. */
	uint32_t *loc_table; 
	// seed table
	seed_entry_t *seed_table;

	uint32_t seq_len;
	uint32_t num_seed_used; // # seed entries in use (including collision entries)
	uint32_t num_seed_key; // # non-collision entries in use (distinguished hash key values)
	
	uint8_t __pad[__pad_size(sizeof(int) + sizeof(uint32_t) * 6 + sizeof(void *) * 3, 64)];
} perfect_table_t;


#define get_seed_loc(pt, loc) (*((pt)->ref_string + (loc)))

static inline uint64_t __get_fw8(uint8_t *fw) {
	return *((uint64_t *) fw);
}

static inline uint64_t __get_rc8(uint8_t *fw) {
#ifdef __GNUC__
	return (~(__builtin_bswap64(*((uint64_t *)fw)))) & 0x0303030303030303LL;
#else
	return (~(    (((uint64_t) fw[0]) << 56)
				| (((uint64_t) fw[1]) << 48)
				| (((uint64_t) fw[2]) << 40)
				| (((uint64_t) fw[3]) << 32)
				| (((uint64_t) fw[4]) << 24)
				| (((uint64_t) fw[5]) << 16)
				| (((uint64_t) fw[6]) << 8)
				| (((uint64_t) fw[7])) ) && 0x0303030303030303LL);
#endif
}

static inline uint64_t __get_fw1(uint8_t *fw) {
	return *fw;
}

static inline uint64_t __get_rc1(uint8_t *fw) {
	return 3 - (*fw);
}

#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
#define __get_ordered_fw8(fw) __get_fw8(fw)
#define __get_ordered_rc8(fw) __get_rc8(fw)
#else
static inline uint64_t __get_ordered_fw8(uint8_t *fw) {
#ifdef __GNUC__
	return __builtin_bswap64(*((uint64_t *)fw));
#else
	return    (((uint64_t) fw[0]) << 56)
			| (((uint64_t) fw[1]) << 48)
			| (((uint64_t) fw[2]) << 40)
			| (((uint64_t) fw[3]) << 32)
			| (((uint64_t) fw[4]) << 24)
			| (((uint64_t) fw[5]) << 16)
			| (((uint64_t) fw[6]) << 8)
			| (((uint64_t) fw[7]));
#endif
}

static inline uint64_t __get_ordered_rc8(uint8_t *fw) {
	return (~(*((uint64_t *)fw))) & 0x0303030303030303LL;
}
#endif

#define __get_ordered_fw1(fw) __get_fw1(fw)
#define __get_ordered_rc1(fw) __get_rc1(fw)

/* it'd be better to use constants for a_fw and b_fw. */
#define ____seedcmp(a, a_fw, b, b_fw, len, ret) do { \
	uint8_t *____a = a_fw ? a : a + len; \
	uint8_t *____b = b_fw ? b : b + len; \
	uint64_t ____a8, ____b8; \
	uint8_t ____a1, ____b1; \
	\
	ret = 0; \
	\
	while (len >= sizeof(uint64_t)) { \
		if (a_fw) { \
			____a8 = __get_ordered_fw8(____a); \
			____a += sizeof(uint64_t); \
		} else { \
			____a -= sizeof(uint64_t); \
			____a8 = __get_ordered_rc8(____a); \
		} \
	\
		if (b_fw) { \
			____b8 = __get_ordered_fw8(____b); \
			____b += sizeof(uint64_t); \
		} else { \
			____b -= sizeof(uint64_t); \
			____b8 = __get_ordered_rc8(____b); \
		} \
	\
		if (____a8 != ____b8) { \
			ret = ____a8 > ____b8 ? 1 : -1; \
			len = 0; /* to skip the next loop */ \
			break; \
		} \
	\
		len -= sizeof(uint64_t); \
	} \
	\
	while (len > 0) { \
		if (a_fw) { \
			____a1 = __get_ordered_fw1(____a); \
			____a++; \
		} else { \
			____a--; \
			____a1 = __get_ordered_rc1(____a); \
		} \
	\
		if (b_fw) { \
			____b1 = __get_ordered_fw1(____b); \
			____b++; \
		} else { \
			____b--; \
			____b1 = __get_ordered_rc1(____b); \
		} \
	\
		if (____a1 != ____b1) { \
			ret = ____a1 > ____b1 ? 1 : -1; \
			break; \
		} \
	\
		len--; \
	} \
} while(0)

static inline int __seedcmp(uint8_t *a, int afl, uint8_t *b, int bfl, int len) {
	int ret;
	
	switch (afl * 2 + bfl) {
	case 0: /* afl == 0 and bfl == 0 */
			____seedcmp(a, 0, b, 0, len, ret);
			return ret;

	case 1: /* afl == 0 and bfl == 1 */
			____seedcmp(a, 0, b, 1, len, ret);
			return ret;

	case 2: /* afl == 1 and bfl == 0 */
			____seedcmp(a, 1, b, 0, len, ret);
			return ret;

	case 3: /* afl == 1 and bfl == 1 */
			____seedcmp(a, 1, b, 1, len, ret);
			return ret;

	default:
			return 0;
	}
}

/* compare FW and RC string in alphabetical order.
 * return 1 if FW <= RC.
 * return 0 if RC < FW.
 */
static inline int __compare_fw_rc(uint8_t *seed, int len) {
	int ret;
	int half = (len + 1) / 2;
	____seedcmp(seed, 1, seed + (len - half), 0, half, ret);
	//____seedcmp(seed, 1, seed, 0, len, ret);
	return ret <= 0 ? 1 : 0;
}

static inline uint64_t __seedmatch_rc(uint8_t *a, uint8_t *b, int len) {
	b = b + len;

	while (len >= sizeof(uint64_t)) {
		b -= sizeof(uint64_t);
		if (__get_ordered_fw8(a) != __get_ordered_rc8(b))
			return 0;
		a += sizeof(uint64_t);
		len -= sizeof(uint64_t);
	}

	while (len > 0) {
		b--;
		if (__get_ordered_fw1(a) != __get_ordered_rc1(b))
			return 0;
		a++;
		len--;
	}

	return 1;
}

static inline uint64_t __seedmatch_fw(uint8_t *a, uint8_t *b, int len) {
	while (len >= sizeof(uint64_t)) {
		if (__get_ordered_fw8(a) != __get_ordered_fw8(b))
			return 0;
		a += sizeof(uint64_t);
		b += sizeof(uint64_t);
		len -= sizeof(uint64_t);
	}
	
	while (len > 0) {
		if (__get_ordered_fw1(a) != __get_ordered_fw1(b))
			return 0;
		a++;
		b++;
		len--;
	}

	return 1;
}

/* return 0 if aloc and bloc are not matched 
 * return 1 if aloc and bloc are same
 * return 2 if bloc is a reverse complement of aloc.
 */
static inline int __seedmatch(perfect_table_t *pt, uint8_t *a, uint8_t *b) {
	int len = pt->seed_len;

	if (__seedmatch_fw(a, b, len))
		return 1;
	else if (__seedmatch_rc(a, b, len))
		return 2;
	else
		return 0;
}

/* r: reference, s: seed */
static inline uint64_t __seedmatch_further_rc(uint8_t *r, uint8_t *s, int offset, int len) {
	s = s + offset + len;
	r = r - len;

	while (len >= sizeof(uint64_t)) {
		s -= sizeof(uint64_t);
		if (__get_ordered_fw8(r) != __get_ordered_rc8(s))
			return 0;
		r += sizeof(uint64_t);
		len -= sizeof(uint64_t);
	}

	while (len > 0) {
		s--;
		if (__get_ordered_fw1(r) != __get_ordered_rc1(s))
			return 0;
		r++;
		len--;
	}

	return 1;
}

/* r: reference, s: seed */
static inline uint64_t __seedmatch_further_fw(uint8_t *r, uint8_t *s, int offset, int len) {
	r = r + offset;
	s = s + offset;

	while (len >= sizeof(uint64_t)) {
		if (__get_ordered_fw8(r) != __get_ordered_fw8(s))
			return 0;
		r += sizeof(uint64_t);
		s += sizeof(uint64_t);
		len -= sizeof(uint64_t);
	}
	
	while (len > 0) {
		if (__get_ordered_fw1(r) != __get_ordered_fw1(s))
			return 0;
		r++;
		s++;
		len--;
	}

	return 1;
}
static inline int __seedmatch_further(perfect_table_t *pt, uint32_t loc,
									uint8_t *seed, int is_rev, int len) {
	len = len - pt->seed_len;	
	assert(len > 0);
	if (is_rev == 0) {
		if (loc + len >= pt->seq_len)
			return 0;
		else
			return __seedmatch_further_fw(pt->ref_string + loc, seed, pt->seed_len, len);
	} else {
		if (loc < len)
			return 0;
		else
			return __seedmatch_further_rc(pt->ref_string + loc, seed, pt->seed_len, len);
	}
}

// from MurmurHash3, public domain, http://code.google.com/p/smhasher/wiki/MurmurHash3
#define BIG_CONSTANT(x) (x##ULL)

static inline uint64_t __fmix64(uint64_t k)

{
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
	k ^= k >> 33;
	return k;
}

#if 0
static inline uint64_t __get_hash(const uint64_t *seed, int bytes) {
	// very simple hash by cdkim
	int nblock = bytes / 8, rem = bytes % 8;
	int i;
	uint64_t h = 0;

	for (i = 0; i < nblock; i++)
		h ^= seed[i];

	if (rem > 0)
		h ^= seed[nblock] >> ((8 - rem) * 8);

	return __fmix64(h);
}
#endif

#if 0
static inline int64_t __get_hash_idx_fw(perfect_table_t *pt, const uint8_t *seed, int len) {
	uint64_t s = 0, h = 0;
	int i;

	for (i = 0; i < len; i++) {
		s = (s << 2) | (seed[i] & 0x3);
		if (i % 32 == (32 - 1)) {
			h ^= s;
			s = 0;
		}
	}

	if (len % 32 > 0) h ^= s;
	
	return (int64_t) (__fmix64(h) % pt->num_seed_entry); 
}
#endif

static inline int64_t __get_hash_idx_fw(perfect_table_t *pt, const uint8_t *seed, int len) {
	uint64_t s = 0, h = 0;
	int i;

	while (len >= 32) {
		s = (uint64_t) (seed[0] & 0x3) << 62
			| (uint64_t) (seed[1] & 0x3) << 60
			| (uint64_t) (seed[2] & 0x3) << 58
			| (uint64_t) (seed[3] & 0x3) << 56
			| (uint64_t) (seed[4] & 0x3) << 54
			| (uint64_t) (seed[5] & 0x3) << 52
			| (uint64_t) (seed[6] & 0x3) << 50
			| (uint64_t) (seed[7] & 0x3) << 48
			| (uint64_t) (seed[8] & 0x3) << 46
			| (uint64_t) (seed[9] & 0x3) << 44
			| (uint64_t) (seed[10] & 0x3) << 42
			| (uint64_t) (seed[11] & 0x3) << 40
			| (uint64_t) (seed[12] & 0x3) << 38
			| (uint64_t) (seed[13] & 0x3) << 36
			| (uint64_t) (seed[14] & 0x3) << 34
			| (uint64_t) (seed[15] & 0x3) << 32
			| (uint64_t) (seed[16] & 0x3) << 30
			| (uint64_t) (seed[17] & 0x3) << 28
			| (uint64_t) (seed[18] & 0x3) << 26
			| (uint64_t) (seed[19] & 0x3) << 24
			| (uint64_t) (seed[20] & 0x3) << 22
			| (uint64_t) (seed[21] & 0x3) << 20
			| (uint64_t) (seed[22] & 0x3) << 18
			| (uint64_t) (seed[23] & 0x3) << 16
			| (uint64_t) (seed[24] & 0x3) << 14
			| (uint64_t) (seed[25] & 0x3) << 12
			| (uint64_t) (seed[26] & 0x3) << 10
			| (uint64_t) (seed[27] & 0x3) << 8
			| (uint64_t) (seed[28] & 0x3) << 6
			| (uint64_t) (seed[29] & 0x3) << 4
			| (uint64_t) (seed[30] & 0x3) << 2
			| (uint64_t) (seed[31] & 0x3);

		h ^= s;

		seed += 32;
		len -= 32;
	}

	if (len == 0) 
		goto out;

	s = 0;
	while (len >= 8) {
		s = (s << 16) 
			| (uint64_t) (seed[0] & 0x3) << 14
			| (uint64_t) (seed[1] & 0x3) << 12
			| (uint64_t) (seed[2] & 0x3) << 10
			| (uint64_t) (seed[3] & 0x3) << 8
			| (uint64_t) (seed[4] & 0x3) << 6
			| (uint64_t) (seed[5] & 0x3) << 4
			| (uint64_t) (seed[6] & 0x3) << 2
			| (uint64_t) (seed[7] & 0x3);
		seed += 8;
		len -= 8;
	}

	while (len > 0) {
		s = (s << 2) | (uint64_t) ((*seed) & 0x3);
		seed++;
		len--;
	}

	h ^= s;

out:	
	return (int64_t) (__fmix64(h) % pt->num_seed_entry); 
}

#if 0
static inline int64_t __get_hash_idx_rc(perfect_table_t *pt, const uint8_t *seed, int len) {
	uint64_t s = 0, h = 0;
	int i;

	for (i = 0; i < len; i++) {
		s = (s << 2) | ((3 - seed[len - i - 1]) & 0x3);
		if (i % 32 == (32 - 1)) {
			h ^= s;
			s = 0;
		}
	}

	if (len % 32 > 0) h ^= s;
	
	return (int64_t) (__fmix64(h) % pt->num_seed_entry); 
}
#endif

static inline int64_t __get_hash_idx_rc(perfect_table_t *pt, const uint8_t *__seed, int len) {
	uint64_t s = 0, h = 0;
	const uint8_t *seed = __seed + len;
	int i;

	while (len >= 32) {
		seed -= 32;
		s = (uint64_t) (3 - (seed[0] & 0x3))
			| (uint64_t) (3 - (seed[1] & 0x3)) << 2
			| (uint64_t) (3 - (seed[2] & 0x3)) << 4
			| (uint64_t) (3 - (seed[3] & 0x3)) << 6
			| (uint64_t) (3 - (seed[4] & 0x3)) << 8
			| (uint64_t) (3 - (seed[5] & 0x3)) << 10
			| (uint64_t) (3 - (seed[6] & 0x3)) << 12
			| (uint64_t) (3 - (seed[7] & 0x3)) << 14
			| (uint64_t) (3 - (seed[8] & 0x3)) << 16
			| (uint64_t) (3 - (seed[9] & 0x3)) << 18
			| (uint64_t) (3 - (seed[10] & 0x3)) << 20
			| (uint64_t) (3 - (seed[11] & 0x3)) << 22
			| (uint64_t) (3 - (seed[12] & 0x3)) << 24
			| (uint64_t) (3 - (seed[13] & 0x3)) << 26
			| (uint64_t) (3 - (seed[14] & 0x3)) << 28
			| (uint64_t) (3 - (seed[15] & 0x3)) << 30
			| (uint64_t) (3 - (seed[16] & 0x3)) << 32
			| (uint64_t) (3 - (seed[17] & 0x3)) << 34
			| (uint64_t) (3 - (seed[18] & 0x3)) << 36
			| (uint64_t) (3 - (seed[19] & 0x3)) << 38
			| (uint64_t) (3 - (seed[20] & 0x3)) << 40
			| (uint64_t) (3 - (seed[21] & 0x3)) << 42
			| (uint64_t) (3 - (seed[22] & 0x3)) << 44
			| (uint64_t) (3 - (seed[23] & 0x3)) << 46
			| (uint64_t) (3 - (seed[24] & 0x3)) << 48
			| (uint64_t) (3 - (seed[25] & 0x3)) << 50
			| (uint64_t) (3 - (seed[26] & 0x3)) << 52
			| (uint64_t) (3 - (seed[27] & 0x3)) << 54
			| (uint64_t) (3 - (seed[28] & 0x3)) << 56
			| (uint64_t) (3 - (seed[29] & 0x3)) << 58
			| (uint64_t) (3 - (seed[30] & 0x3)) << 60
			| (uint64_t) (3 - (seed[31] & 0x3)) << 62;

		h ^= s;

		len -= 32;
	}

	if (len == 0) 
		goto out;

	s = 0;
	while (len >= 8) {
		seed -= 8;
		s = (s << 16) 
			| (uint64_t) (3 - (seed[0] & 0x3))
			| (uint64_t) (3 - (seed[1] & 0x3)) << 2
			| (uint64_t) (3 - (seed[2] & 0x3)) << 4
			| (uint64_t) (3 - (seed[3] & 0x3)) << 6
			| (uint64_t) (3 - (seed[4] & 0x3)) << 8
			| (uint64_t) (3 - (seed[5] & 0x3)) << 10
			| (uint64_t) (3 - (seed[6] & 0x3)) << 12
			| (uint64_t) (3 - (seed[7] & 0x3)) << 14;
		len -= 8;
	}

	while (len > 0) {
		seed--;
		s = (s << 2) | (uint64_t) (3 - ((*seed) & 0x3));
		len--;
	}

	h ^= s;

out:	
	return (int64_t) (__fmix64(h) % pt->num_seed_entry); 
}

static inline int __compare_fw_rc(const uint8_t *seed, int len);

static inline int64_t __get_hash_idx_seed(perfect_table_t *pt, uint8_t *seed, int fw_less) {
	if (fw_less)
		return __get_hash_idx_fw(pt, seed, pt->seed_len);
	else
		return __get_hash_idx_rc(pt, seed, pt->seed_len);
}

static inline int64_t get_hash_idx_seed(perfect_table_t *pt, uint8_t *seed) {
	if (__compare_fw_rc(seed, pt->seed_len))
		return __get_hash_idx_fw(pt, seed, pt->seed_len);
	else
		return __get_hash_idx_rc(pt, seed, pt->seed_len);
}

/* index-building only use ref_string */
static inline int64_t get_hash_idx_loc(perfect_table_t *pt, uint32_t loc) {
	return get_hash_idx_seed(pt, pt->ref_string + loc);
}

static inline int64_t get_hash_idx_ent(perfect_table_t *pt, seed_entry_t *ent) {
	return __get_hash_idx_seed(pt, pt->ref_string + ent->location, is_fw_less_entry(ent));
}

//#define pt_reference(pt) ((pt)->ref_string)

typedef struct {
	int64_t loc;
	int64_t pos;
	int rid;
	int flag;
	uint32_t is_rev:1, is_alt:1;
	//uint32_t mapq = 255, NM:22 = 0;
	//int n_cigar = 1;
	//uint32_t *cigar => [(l_seq << 4 | 0x0)];
	//char *XA;
	//int score = l_seq * opt->a;
	int sub;
} mem_aln_perfect_t;

typedef struct { size_t n, m; mem_aln_perfect_t *a; } mem_aln_perfect_v;

static inline seed_entry_t *get_seed_entry(perfect_table_t *pt, int64_t key) {
#ifdef MEMSCALE
	return key < pt->num_seed_load ? &pt->seed_table[key] : NULL;
#else
	return key == NO_ENTRY ? NULL : &pt->seed_table[key];
#endif
}

//int64_t get_hash_idx(perfect_table_t *pt, const uint64_t *seed);
void show_seed_entry(perfect_table_t *pt, uint32_t key);
void show_perfect_table_related(perfect_table_t *pt, uint32_t start);

extern perfect_table_t *perfect_table;
extern int perfect_table_seed_len;

void free_perfect_table();

/***************************************/
/* Load Perfect Table helper functions */
/***************************************/
static inline void __lpt_set_table_ptr(perfect_table_t *pt, perfect_table_t *chunk) {
	uint8_t *ptr = (uint8_t *) chunk;
	ptr += __aligned_size(sizeof(perfect_table_t), 64);
	pt->loc_table = (uint32_t *) ptr;
	ptr += __aligned_size(sizeof(uint32_t) * pt->num_loc_entry, 64);
	pt->seed_table = (seed_entry_t *) ptr;
}
static inline void __lpt_link_shm_to_pt(perfect_table_t *pt, perfect_table_t *shm) {
	memcpy(pt, shm, sizeof(perfect_table_t));
	__lpt_set_table_ptr(pt, shm);
}

static inline void __lpt_load_head(perfect_table_t *pt, FILE *fp) {
	assert(sizeof(perfect_table_t) % 64 == 0);
	err_fread_noeof(pt, sizeof(perfect_table_t), 1, fp);
#ifdef MEMSCALE
	pt->num_seed_load = pt->num_seed_entry;
#endif
}

#ifdef MEMSCALE
static inline void __lpt_set_num_seed_load(perfect_table_t *pt, uint32_t num_seed_load) {
	if (num_seed_load > 0 && num_seed_load < pt->num_seed_entry)
		pt->num_seed_load = num_seed_load;
	else
		pt->num_seed_load = pt->num_seed_entry;
}
#endif

#define ____lpt_shm_size(num_loc, num_seed) \
			( __aligned_size(sizeof(perfect_table_t), 64) \
			+ __aligned_size(sizeof(uint32_t) * (num_loc), 64) \
			+ __aligned_size(sizeof(seed_entry_t) * (num_seed), 64))

static inline size_t __lpt_shm_size(perfect_table_t *pt) {
#ifdef MEMSCALE
	return ____lpt_shm_size(pt->num_loc_entry, pt->num_seed_load);
#else
	return ____lpt_shm_size(pt->num_loc_entry, pt->num_seed_entry);
#endif
}

static inline size_t ____lpt_file_size(uint32_t num_loc_entry, uint32_t num_seed_entry) {
	return sizeof(perfect_table_t) 
					+ num_loc_entry * sizeof(uint32_t) 
					+ num_seed_entry * sizeof(seed_entry_t);
}

static inline size_t __lpt_file_size(perfect_table_t *pt) {
	return ____lpt_file_size(pt->num_loc_entry, pt->num_seed_entry);
}

static inline void __lpt_show_info(perfect_table_t *pt) {
	fprintf(stderr, "Reading perfect table size: %.2fGB seed_len: %u seq_len: %u #seed: %u #used: %u (%.2f%%) #key: %u (%.2f%%) #loc: %u\n",
			(float) __lpt_file_size(pt) / (1024 * 1024 * 1024),
			pt->seed_len, pt->seq_len, pt->num_seed_entry,
			pt->num_seed_used, (float) pt->num_seed_used * 100 / (float) pt->num_seed_entry,
			pt->num_seed_key, (float) pt->num_seed_key * 100 / (float) pt->num_seed_entry,
			pt->num_loc_entry);
	fflush(stderr);
}

static inline void __lpt_load_loc_table(perfect_table_t *pt, FILE *fp) {
	/* offset of fp should be the start of loc_table */
	uint32_t num_loaded = 0;
	int pct;

	for (pct = 0; pct < 10; pct++) {
		uint32_t chunk = pt->num_loc_entry / 10 + (pct < pt->num_loc_entry % 10 ? 1 : 0);
		err_fread_noeof(pt->loc_table + num_loaded, sizeof(uint32_t), chunk, fp);
		num_loaded += chunk;
		fprintf(stderr, "[Reading Location] %3u%% (%u/%u)\n", (pct + 1) * 10, num_loaded, pt->num_loc_entry);
		fflush(stderr);
	}
}

static inline void ____lpt_load_seed_table(perfect_table_t *pt, FILE *fp, 
												uint32_t beg __maybe_unused, 
												uint32_t end __maybe_unused) 
{
	uint32_t num_loaded = 0;
#ifdef MEMSCALE
	uint32_t num_to_load;
#else
#define num_to_load (pt->num_seed_entry)
#endif
	int pct;
	seed_entry_t *ptr = pt->seed_table;

#ifdef MEMSCALE
	if (end == 0)
		beg = 0, end = pt->num_seed_entry;
	else if (end > pt->num_seed_entry)
		end = pt->num_seed_entry;
	num_to_load = end - beg;
	fseek(fp, ____lpt_file_size(pt->num_loc_entry, beg), SEEK_SET);
	ptr += beg;
	
	if (beg > 0 || end < pt->num_seed_entry) 
		fprintf(stderr, "[Reading Table] part: %u ~ %u\n", beg, end);
#endif

	for (pct = 0; pct < 100; pct++) {
		uint32_t chunk = num_to_load / 100 + (pct < num_to_load % 100 ? 1 : 0);
		err_fread_noeof(ptr + num_loaded, sizeof(seed_entry_t), chunk, fp);
		num_loaded += chunk;
		fprintf(stderr, "[Reading Table] %3u%% (%u/%u)\n", pct + 1, num_loaded, num_to_load);
		fflush(stderr);
	}

#ifndef MEMSCALE
#undef num_to_load
#endif
}

static inline void __lpt_load_seed_table(perfect_table_t *pt, FILE *fp) {
#ifdef MEMSCALE
	____lpt_load_seed_table(pt, fp, 0, pt->num_seed_load);
#else
	____lpt_load_seed_table(pt, fp, 0, pt->num_seed_entry);
#endif
}


//#define PERFECT_PROFILE
/* perfect_map.cpp and bntseq.cpp also has some PERFECT_PROFILE related code */
#ifdef PERFECT_PROFILE
void set_perfect_profile_rid(int n_seq);
#endif

#define FIND_PERFECT_NO_TABLE 0
#define FIND_PERFECT_WITH_N 1
#define FIND_PERFECT_NOT_MATCHED 2
#define FIND_PERFECT_FW_MATCHED 3
#define FIND_PERFECT_RC_MATCHED 4
#define FIND_PERFECT_SEED_ONLY_MATCHED 5

#define PT_SEED_LEN_NO_TABLE (INT_MAX)
#define PT_SEED_LEN_AUTO_TABLE (INT_MAX - 1)

#endif /* PERFECT_MATCH */
#endif /* BWT_PERFECT_H */
