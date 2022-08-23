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

#include <stdio.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "sais.h"
#include "FMI_search.h"
#include "memcpy_bwamem.h"
#include "profiling.h"
#include "bwa_shm.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "safe_str_lib.h"
#ifdef __cplusplus
}
#endif

void *__load_file(const char *prefix, const char *postfix, void *buf, size_t *size);

#ifdef SMEM_ACCEL
static int building_smem_table = 0;

static inline int __seq_next(uint8_t *list, int len) {
	int i = len - 1;
	while (i >= 0) {
		 if (list[i] < 3) {
			 list[i]++;
			 for (i = i + 1; i < len; ++i)
				 list[i] = 0;
			 return 1;
		 }
		 i--;
	}

	return 0;
}

void FMI_search::__build_all_smem_table(uint8_t *seq, int len, all_smem_t *ent) {
	int i;
	uint8_t a;

	memset(ent, 0, sizeof(all_smem_t));

	a = seq[0];
	SMEM smem;
	smem.rid = 0;
	smem.m = 0;
	smem.n = 0;
	smem.k = count[a];
	smem.l = count[3 - a];
    smem.s = count[a + 1] - count[a];

	for (i = 1; i < len; ++i) {
		a = seq[i];
		SMEM smem_ = smem;

		// Forward extension is backward extension with the BWT of reverse complement
		smem_.k = smem.l;
		smem_.l = smem.k;
		SMEM newSmem_ = backwardExt(smem_, 3 - a);
		SMEM newSmem = newSmem_;
		newSmem.k = newSmem_.l;
		newSmem.l = newSmem_.k;
		newSmem.n = i;

		assert(newSmem.l >= count[3 - a]);
		ent->list[i - 1].l32 = (uint32_t) (newSmem.l - count[3 - a]);
		assert(newSmem.k >= smem.k);
		ent->list[i - 1].k32 = (uint32_t) (newSmem.k - smem.k);
		assert(newSmem.s <= UINT32_MAX);
		ent->list[i - 1].s32 = newSmem.s;
		
		if (newSmem.s > 0)
			ent->last_avail = i;
		else
			break;

		smem = newSmem;
	}
}

all_smem_t *FMI_search::build_all_smem_table(int len) {
	all_smem_t *table;
	int64_t num_entry, num_step;
	uint8_t seq[len];
	uint64_t i;

	num_entry = __num_smem_table_entry(len);
	num_step = num_entry / 2;
	table = (all_smem_t *) _mm_malloc(num_entry * sizeof(all_smem_t), 64);
	if (!table) {
		printf("ERROR: cannot allocate memory for all smem table\n");
		return NULL;
	}

	for (i = 0; i < len; ++i)
		seq[i] = 0;

	i = 0;
	__build_all_smem_table(seq, len, &table[i++]);
	while (__seq_next(seq, len)) {
		__build_all_smem_table(seq, len, &table[i++]);
		if (i % num_step == 0) {
			printf("%s: progress %ld/%ld (%.2f%%)\n",
					__func__, i, num_entry, 
					(double) i * 100 / (double) num_entry);
		}
	}
	
	assert(i == num_entry);

	return table;
}

void FMI_search::__build_last_smem_table(uint8_t *seq, int len, last_smem_t *ent) {
	int i = 0;
	uint8_t a;

	a = seq[0];
	SMEM smem;
	smem.rid = 0;
	smem.m = 0;
	smem.n = 0;
	smem.k = count[a];
	smem.l = count[3 - a];
    smem.s = count[a + 1] - count[a];

	for (i = 1; i < len; ++i) {
		a = seq[i];
		SMEM smem_ = smem;

		// Forward extension is backward extension with the BWT of reverse complement
		smem_.k = smem.l;
		smem_.l = smem.k;
		SMEM newSmem_ = backwardExt(smem_, 3 - a);
		SMEM newSmem = newSmem_;
		newSmem.k = newSmem_.l;
		newSmem.l = newSmem_.k;
		newSmem.n = i;

		if (newSmem.s == 0)
			break;

		smem = newSmem;
	}

	ent->bp = i;
	ent->kms = (int8_t) (smem.k >> 32);
	ent->kls = (uint32_t) (smem.k & 0xffffffff);
	ent->lms = (int8_t) (smem.l >> 32);
	ent->lls = (uint32_t) (smem.l & 0xffffffff);
	ent->sms = (int8_t) (smem.s >> 32);
	ent->sls = (uint32_t) (smem.s & 0xffffffff);
}

last_smem_t *FMI_search::build_last_smem_table(int len) {
	last_smem_t *table;
	int64_t num_entry, num_step;
	uint8_t seq[len];
	uint64_t i;

	num_entry = __num_smem_table_entry(len);
	num_step = num_entry / 100;
	table = (last_smem_t *)_mm_malloc(num_entry * sizeof(last_smem_t), 64);
	if (!table) {
		printf("ERROR: cannot allocate memory for last smem table\n");
		return NULL;
	}

	for (i = 0; i < len; ++i)
		seq[i] = 0;

	i = 0;
	__build_last_smem_table(seq, len, &table[i++]);
	while (__seq_next(seq, len)) {
		__build_last_smem_table(seq, len, &table[i++]);
		if (i % num_step == 0) {
			printf("%s: progress %ld/%ld (%.2f%%)\n",
					__func__, i, num_entry, 
					(double) i * 100 / (double) num_entry);
		}
	}
	
	assert(i == num_entry);

	return table;
}

int build_smem_tables(char *prefix) {
	char all_smem_fn[PATH_MAX];
	char last_smem_fn[PATH_MAX];
	FILE *fp;
	all_smem_t *all_smem_table;
	last_smem_t *last_smem_table;
	FMI_search *fmi;

	fmi = new FMI_search(prefix);
	building_smem_table = 1;
	fmi->load_index();

	/* all smem table */
	snprintf_s_si(all_smem_fn, PATH_MAX, "%s.all_smem.%d", prefix, ALL_SMEM_MAX_BP);
	printf("Build all smem table\n");
	all_smem_table = fmi->build_all_smem_table(ALL_SMEM_MAX_BP);
	if (all_smem_table == NULL) {
		printf("ERROR: failed to build all smem table\n");
		return -1;
	}
	printf("Write all smem table to %s\n", all_smem_fn);
	fp = xopen(all_smem_fn, "wb");
	err_fwrite(all_smem_table, sizeof(all_smem_t),
				__num_smem_table_entry(ALL_SMEM_MAX_BP), fp);
	err_fflush(fp);
	err_fclose(fp);
	_mm_free(all_smem_table);
	all_smem_table = NULL;

	/* last smem table */
	snprintf_s_si(last_smem_fn, PATH_MAX, "%s.last_smem.%d", prefix, LAST_SMEM_MAX_BP);
	printf("Build last smem table\n");
	last_smem_table = fmi->build_last_smem_table(LAST_SMEM_MAX_BP);
	if (last_smem_table == NULL) {
		printf("ERROR: failed to build last smem table\n");
		return -1;
	}
	printf("Write last smem table to %s\n", last_smem_fn);
	fp = xopen(last_smem_fn, "wb");
	err_fwrite(last_smem_table, sizeof(last_smem_t),
				__num_smem_table_entry(LAST_SMEM_MAX_BP), fp);
	err_fflush(fp);
	err_fclose(fp);
	_mm_free(last_smem_table);
	last_smem_table = NULL;
	
	delete(fmi);
	return 0;
}

static int __load_smem_table_from_file(const char *prefix, 
									all_smem_t **__all_smem_table, 
									last_smem_t **__last_smem_table)
{
	FILE *fp;
	all_smem_t *all_smem_table;
	last_smem_t *last_smem_table;

	/* all smem table */
	if (__all_smem_table) {
		char all_smem_fn[PATH_MAX];
		
		all_smem_table = *__all_smem_table;
		
		/* memory allocation if required */
		if (all_smem_table == NULL) 
			all_smem_table = (all_smem_t *)_mm_malloc(ALL_SMEM_TABLE_SIZE, 64);
		
		if (all_smem_table == NULL) {
			fprintf(stderr, "ERROR: cannot allocate memory for all smem table\n");
			exit(EXIT_FAILURE);
		}

		fprintf(stderr, "INFO: load all smem table from file (len: %d)\n", ALL_SMEM_MAX_BP);
		snprintf(all_smem_fn, PATH_MAX, "%s.all_smem.%d", prefix, ALL_SMEM_MAX_BP);
		fp = xopen(all_smem_fn, "rb");
		err_fread_noeof(all_smem_table, sizeof(all_smem_t),
					__num_smem_table_entry(ALL_SMEM_MAX_BP), fp);
		err_fclose(fp);
	
		/* return */
		*__all_smem_table = all_smem_table;
	}

	/* last smem table */
	if (__last_smem_table) {
		char last_smem_fn[PATH_MAX];
		
		last_smem_table = *__last_smem_table;
		
		/* memory allocation if required */
		if (last_smem_table == NULL) 
			last_smem_table = (last_smem_t *)_mm_malloc(LAST_SMEM_TABLE_SIZE, 64);
		
		if (last_smem_table == NULL) {
			fprintf(stderr, "ERROR: cannot lastocate memory for last smem table\n");
			exit(EXIT_FAILURE);
		}

		fprintf(stderr, "INFO: load last smem table from file (len: %d)\n", LAST_SMEM_MAX_BP);
		snprintf(last_smem_fn, PATH_MAX, "%s.last_smem.%d", prefix, LAST_SMEM_MAX_BP);
		fp = xopen(last_smem_fn, "rb");
		err_fread_noeof(last_smem_table, sizeof(last_smem_t),
					__num_smem_table_entry(LAST_SMEM_MAX_BP), fp);
		err_fclose(fp);

		/* return */
		*__last_smem_table = last_smem_table;
	}
	
	return 0;
}

#ifdef USE_SHM
static inline void __load_smem_table_set_shm_ptr(all_smem_t **__all,
												last_smem_t **__last,
												uint8_t *shm)
{
	if (__all) {
		*__all = (all_smem_t *) shm;
		shm += ALL_SMEM_TABLE_SIZE;
	}

	if (__last) {
		*__last = (last_smem_t *) shm;
	}
}

static int __load_smem_table_on_shm(const char *prefix,
									all_smem_t **__all_smem_table, 
									last_smem_t **__last_smem_table, 
									uint8_t *shm)
{
	__load_smem_table_set_shm_ptr(__all_smem_table, __last_smem_table, shm);
	return __load_smem_table_from_file(prefix, __all_smem_table, __last_smem_table);
}

static int __load_smem_table_from_shm(all_smem_t **__all_smem_table, 
									  last_smem_t **__last_smem_table, 
									  uint8_t *shm)
{
	__load_smem_table_set_shm_ptr(__all_smem_table, __last_smem_table, shm);
	return 0;
}

int _load_smem_table(const char *prefix, 
					all_smem_t **__all_smem_table, 
					last_smem_t **__last_smem_table) 
{
	fprintf(stderr, "INFO: load smem table (all_smem_len: %d last_smem_len: %d)\n",
						ALL_SMEM_MAX_BP, LAST_SMEM_MAX_BP);
	
	if (__all_smem_table &&
			__bwa_shm_load_file(prefix, STR_AND_VAL(".all_smem.", ALL_SMEM_MAX_BP),
								BWA_SHM_SALL, (void **) __all_smem_table) != 0) 
		return -1;
	
	if (__last_smem_table &&
			__bwa_shm_load_file(prefix, STR_AND_VAL(".last_smem.", LAST_SMEM_MAX_BP),
								BWA_SHM_SLAST, (void **) __last_smem_table) != 0) 
		return -1;

	return 0;
}

void FMI_search::load_smem_table() {
#ifdef MEMSCALE
	_load_smem_table(file_name,
						bwa_shm_info->smem_all_on ? &all_smem_table : NULL,
						bwa_shm_info->smem_last_on ? &last_smem_table : NULL);
#else
	//fprintf(stderr, "[DEBUG] %s all: %p last: %p\n", __func__, all_smem_table, last_smem_table);
	_load_smem_table(file_name, &all_smem_table, &last_smem_table);
#endif
}
#else
void FMI_search::load_smem_table() {
	all_smem_table = (all_smem_t *) __load_file(file_name, STR_AND_VAL(".all_smem.", ALL_SMEM_MAX_BP),
									NULL, NULL);
	last_smem_table = (last_smem_t *) __load_file(file_name, STR_AND_VAL(".last_smem.", LAST_SMEM_MAX_BP),
									NULL, NULL);
}
#endif

#endif

//#define SMEM_ACCEL_DEBUG
#if !defined(PERFECT_MATCH) && defined(SMEM_ACCEL_DEBUG)
void debug_smem_input(const char *header, int32_t rid, uint8_t *enc, int x, int readlength) {
	char base[readlength + 2];
	int i, e;

	e = 0;
	i = 0;
	while (e < readlength) {
		if (e == x)
			base[i++] = '!';
		base[i] = enc[e] < 4 ? "ACGT"[enc[e]] : 'N';
		i++; e++;
	}
	base[i] = '\0';

	fprintf(stderr, "[SMEM_INPUT]%s%s rid: %6d x: %3d len: %3d base: %s\n",
				header ? " " : "", header ? header : "", rid, x, readlength, base);

}

void debug_smem_output(const char *header, SMEM *array, int num_smem, int32_t *data) {
	if (num_smem == 1) { 
		SMEM *smem = &array[0];
		fprintf(stderr, "[SMEM_OUTPUT]%s%s rid: %6d m: %3d n: %3d k: %10ld l: %10ld s: %10ld",
					header ? " " : "", header ? header : "", 
											smem->rid, smem->m, smem->n, smem->k, smem->l, smem->s);
		if (data)
			fprintf(stderr, " data: %d\n", *data);
		else
			fprintf(stderr, "\n");
	} else {
		int i;
		fprintf(stderr, "[SMEM_OUTPUT]%s%s num_smem: %d",
					header ? " " : "", header ? header : "", num_smem);
		if (data)
			fprintf(stderr, " data: %d\n", *data);
		else
			fprintf(stderr, "\n");
			
		for (i = 0; i < num_smem; ++i) {
			SMEM *smem = &array[i];
			fprintf(stderr, "             [%03d] rid: %6d m: %3d n: %3d k: %10ld l: %10ld s: %10ld\n",
											i, smem->rid, smem->m, smem->n, smem->k, smem->l, smem->s);
		}
	}
}
#else
#define debug_smem_input(a, b, c, d, e) do{}while(0)
#define debug_smem_output(a, b, c, d) do{}while(0)
#endif

void FMI_search::init(const char *fname) 
{
    fprintf(stderr, "* Entering FMI_search\n");
    //strcpy(file_name, fname);
    strcpy_s(file_name, PATH_MAX, fname);
    reference_seq_len = 0;
    sentinel_index = 0;
#ifdef PERFECT_MATCH
	perfect_table = NULL;
#endif
    index_alloc = 0;
    sa_ls_word = NULL;
    sa_ms_byte = NULL;
    cp_occ = NULL;
#ifdef SMEM_ACCEL
	all_smem_table = NULL;
	last_smem_table = NULL;
#endif
	useErt = 0;
	kmer_offsets = NULL;
	mlt_table = NULL;

}

FMI_search::FMI_search(const char *fname)
{
	init(fname);
}

FMI_search::~FMI_search()
{
	//fprintf(stderr, "[DEBUG] %s all: %p last: %p\n", __func__, all_smem_table, last_smem_table);
#define _mm_free_safe(ptr) do { if (ptr) _mm_free(ptr); } while (0)
   	if (useErt) {
#ifdef USE_SHM 
		if (bwa_shm_unmap(BWA_SHM_KMER))
			_mm_free_safe(kmer_offsets);
		if (bwa_shm_unmap(BWA_SHM_MLT))
			_mm_free_safe(mlt_table);
#else
        _mm_free_safe(kmer_offsets);
        _mm_free_safe(mlt_table);
#endif
	} else {
#ifdef USE_SHM 
		if (bwa_shm_unmap(BWA_SHM_BWT)) {
			_mm_free_safe(sa_ms_byte);
			_mm_free_safe(sa_ls_word);
			_mm_free_safe(cp_occ);
		}
#ifdef SMEM_ACCEL
		if (bwa_shm_unmap(BWA_SHM_SALL))
			_mm_free_safe(all_smem_table);
		if (bwa_shm_unmap(BWA_SHM_SLAST))
			_mm_free_safe(last_smem_table);
#endif
#else
		_mm_free_safe(sa_ms_byte);
		_mm_free_safe(sa_ls_word);
		_mm_free_safe(cp_occ);
#ifdef SMEM_ACCEL
		_mm_free_safe(all_smem_table);
		_mm_free_safe(last_smem_table);
#endif
#endif
	}
#undef _mm_free_safe
}

int64_t FMI_search::pac_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	uint8_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	assert(c >= 0 && c <= 255);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

void FMI_search::pac2nt(const char *fn_pac, std::string &reference_seq)
{
	uint8_t *buf2;
	int64_t i, pac_size, seq_len;
	FILE *fp;

	// initialization
	seq_len = pac_seq_len(fn_pac);
    assert(seq_len > 0);
    assert(seq_len <= 0x7fffffffffL);
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (seq_len>>2) + ((seq_len&3) == 0? 0 : 1);
	buf2 = (uint8_t*)calloc(pac_size, 1);
    assert(buf2 != NULL);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	for (i = 0; i < seq_len; ++i) {
		int nt = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
        switch(nt)
        {
            case 0:
                reference_seq += "A";
            break;
            case 1: reference_seq += "C"; break;
            case 2:
                reference_seq += "G";
            break;
            case 3:
                reference_seq += "T";
            break;
            default:
                fprintf(stderr, "ERROR! Value of nt is not in 0,1,2,3!");
                exit(EXIT_FAILURE);
        }
	}
    for(i = seq_len - 1; i >= 0; i--)
    {
        char c = reference_seq[i];
        switch(c)
        {
            case 'A':
                reference_seq += "T";
            break;
            case 'C':
                reference_seq += "G";
            break;
            case 'G':
                reference_seq += "C";
            break;
            case 'T':
                reference_seq += "A";
            break;
	    default: /* cannot reach here */
	    break;
        }
    }
	free(buf2);
}

int FMI_search::build_fm_index(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count) {
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[PATH_MAX];

    strcpy_s(outname, PATH_MAX, ref_file_name);
    strcat_s(outname, PATH_MAX, CP_FILENAME_SUFFIX);
    //sprintf(outname, "%s.bwt.2bit.%d", ref_file_name, CP_BLOCK_SIZE);

    std::fstream outstream (outname, std::ios::out | std::ios::binary);
    outstream.seekg(0);	

    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char*)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE - 1) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;
    int64_t size = ref_seq_len_aligned * sizeof(uint8_t);
    bwt = (uint8_t *)_mm_malloc(size, 64);
    assert_not_null(bwt, size, index_alloc);

    int64_t sentinel_index = -1;
    for(i=0; i< ref_seq_len; i++)
    {
        if(sa_bwt[i] == 0)
        {
            bwt[i] = 4;
            printf("BWT[%ld] = 4\n", i);
            sentinel_index = i;
        }
        else
        {
            char c = binary_seq[sa_bwt[i]-1];
            switch(c)
            {
                case 0: bwt[i] = 0;
                          break;
                case 1: bwt[i] = 1;
                          break;
                case 2: bwt[i] = 2;
                          break;
                case 3: bwt[i] = 3;
                          break;
                default:
                        fprintf(stderr, "ERROR! i = %ld, c = %c\n", i, c);
                        exit(EXIT_FAILURE);
            }
        }
    }
    for(i = ref_seq_len; i < ref_seq_len_aligned; i++)
        bwt[i] = DUMMY_CHAR;


    printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT, CP_MASK);
    printf("sizeof CP_OCC = %ld\n", sizeof(CP_OCC));
    fflush(stdout);
    // create checkpointed occ
    int64_t cp_occ_size = (ref_seq_len >> CP_SHIFT) + 1;
    CP_OCC *cp_occ = NULL;

    size = cp_occ_size * sizeof(CP_OCC);
    cp_occ = (CP_OCC *)_mm_malloc(size, 64);
    assert_not_null(cp_occ, size, index_alloc);
    memset_s(cp_occ, cp_occ_size * sizeof(CP_OCC), 0);
    int64_t cp_count[16];

    memset_s(cp_count, 16 * sizeof(int64_t), 0);
    for(i = 0; i < ref_seq_len; i++)
    {
        if((i & CP_MASK) == 0)
        {
            CP_OCC cpo;
            cpo.cp_count[0] = cp_count[0];
            cpo.cp_count[1] = cp_count[1];
            cpo.cp_count[2] = cp_count[2];
            cpo.cp_count[3] = cp_count[3];

			int32_t j;
            cpo.one_hot_bwt_str[0] = 0;
            cpo.one_hot_bwt_str[1] = 0;
            cpo.one_hot_bwt_str[2] = 0;
            cpo.one_hot_bwt_str[3] = 0;

			for(j = 0; j < CP_BLOCK_SIZE; j++)
			{
                cpo.one_hot_bwt_str[0] = cpo.one_hot_bwt_str[0] << 1;
                cpo.one_hot_bwt_str[1] = cpo.one_hot_bwt_str[1] << 1;
                cpo.one_hot_bwt_str[2] = cpo.one_hot_bwt_str[2] << 1;
                cpo.one_hot_bwt_str[3] = cpo.one_hot_bwt_str[3] << 1;
				uint8_t c = bwt[i + j];
                //printf("c = %d\n", c);
                if(c < 4)
                {
                    cpo.one_hot_bwt_str[c] += 1;
                }
			}

            cp_occ[i >> CP_SHIFT] = cpo;
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC));
    _mm_free(cp_occ);
    _mm_free(bwt);

    #if SA_COMPRESSION  

    size = ((ref_seq_len >> SA_COMPX)+ 1)  * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    int64_t pos = 0;
    for(i = 0; i < ref_seq_len; i++)
    {
        if ((i & SA_COMPX_MASK) == 0)
        {
            sa_ls_word[pos] = sa_bwt[i] & 0xffffffff;
            sa_ms_byte[pos] = (sa_bwt[i] >> 32) & 0xff;
            pos++;
        }
    }
    fprintf(stderr, "pos: %ld, ref_seq_len__: %ld\n", pos, ref_seq_len >> SA_COMPX);
    outstream.write((char*)sa_ms_byte, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(uint32_t));
    
    #else
    
    size = ref_seq_len * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ref_seq_len * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    for(i = 0; i < ref_seq_len; i++)
    {
        sa_ls_word[i] = sa_bwt[i] & 0xffffffff;
        sa_ms_byte[i] = (sa_bwt[i] >> 32) & 0xff;
    }
    outstream.write((char*)sa_ms_byte, ref_seq_len * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ref_seq_len * sizeof(uint32_t));
    
    #endif

    outstream.write((char *)(&sentinel_index), 1 * sizeof(int64_t));
    outstream.close();
    printf("max_occ_ind = %ld\n", i >> CP_SHIFT);    
    fflush(stdout);

    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}

int FMI_search::build_index() {

    char *prefix = file_name;
    unsigned long long startTick;
    startTick = __rdtsc();
    index_alloc = 0;

    std::string reference_seq;
    char pac_file_name[PATH_MAX];
    strcpy_s(pac_file_name, PATH_MAX, prefix);
    strcat_s(pac_file_name, PATH_MAX, ".pac");
    //sprintf(pac_file_name, "%s.pac", prefix);
    pac2nt(pac_file_name, reference_seq);
	int64_t pac_len = reference_seq.length();
    int status;
    int64_t size = pac_len * sizeof(char);
    char *binary_ref_seq = (char *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(binary_ref_seq, size, index_alloc);
    char binary_ref_name[PATH_MAX];
    strcpy_s(binary_ref_name, PATH_MAX, prefix);
    strcat_s(binary_ref_name, PATH_MAX, ".0123");
    //sprintf(binary_ref_name, "%s.0123", prefix);
    std::fstream binary_ref_stream (binary_ref_name, std::ios::out | std::ios::binary);
    binary_ref_stream.seekg(0);
    fprintf(stderr, "init ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();
    int64_t i, count[16];
	memset(count, 0, sizeof(int64_t) * 16);
    for(i = 0; i < pac_len; i++)
    {
        switch(reference_seq[i])
        {
            case 'A':
            binary_ref_seq[i] = 0, ++count[0];
            break;
            case 'C':
            binary_ref_seq[i] = 1, ++count[1];
            break;
            case 'G':
            binary_ref_seq[i] = 2, ++count[2];
            break;
            case 'T':
            binary_ref_seq[i] = 3, ++count[3];
            break;
            default:
            binary_ref_seq[i] = 4;

        }
    }
    count[4]=count[0]+count[1]+count[2]+count[3];
    count[3]=count[0]+count[1]+count[2];
    count[2]=count[0]+count[1];
    count[1]=count[0];
    count[0]=0;
    fprintf(stderr, "ref seq len = %ld\n", pac_len);
    binary_ref_stream.write(binary_ref_seq, pac_len * sizeof(char));
    fprintf(stderr, "binary seq ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    size = (pac_len + 2) * sizeof(int64_t);
    int64_t *suffix_array=(int64_t *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(suffix_array, size, index_alloc);
    startTick = __rdtsc();
	//status = saisxx<const char *, int64_t *, int64_t>(reference_seq.c_str(), suffix_array + 1, pac_len, 4);
	status = saisxx(reference_seq.c_str(), suffix_array + 1, pac_len);
	suffix_array[0] = pac_len;
    fprintf(stderr, "build suffix-array ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

	build_fm_index(prefix, binary_ref_seq, pac_len, suffix_array, count);
    fprintf(stderr, "build fm-index ticks = %llu\n", __rdtsc() - startTick);
    _mm_free(binary_ref_seq);
    _mm_free(suffix_array);
    return 0;
}

void FMI_search::load_index_other_elements(int which) {
    char *ref_file_name = file_name;
    fprintf(stderr, "Reading other elements of the index from files %s\n", ref_file_name);
    bwa_idx_load_ele(ref_file_name, which);
}

#ifdef USE_SHM
int __load_BWT_from_file(char *cp_file_name, int64_t reference_seq_len, 
						int64_t *count, 
						CP_OCC *cp_occ, int64_t cp_occ_size,
						int8_t *sa_ms_byte, uint32_t *sa_ls_word,
						int64_t *_sentinel_index)
{
	FILE *cpstream = NULL;
	int64_t xx;
	
	cpstream = fopen(cp_file_name, "rb");
    
	if (cpstream == NULL) {
        fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
        exit(EXIT_FAILURE);
    } else
        fprintf(stderr, "* Index file found. Loading index from %s\n", cp_file_name);

    err_fread_noeof(&xx, sizeof(int64_t), 1, cpstream);
    assert(xx == reference_seq_len);
	
	err_fread_noeof(count, sizeof(int64_t), 5, cpstream);
    int64_t ii = 0;
    for(ii = 0; ii < 5; ii++)// update read count structure
    {
        count[ii] = count[ii] + 1;
    }

	// create checkpointed occ
	err_fread_noeof(cp_occ, sizeof(CP_OCC), cp_occ_size, cpstream);

    
	#if SA_COMPRESSION

    int64_t reference_seq_len_ = (reference_seq_len >> SA_COMPX) + 1;
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len_, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len_, cpstream);
    
    #else
    
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len, cpstream);

    #endif

    int64_t sentinel_index = -1;
    #if SA_COMPRESSION
    err_fread_noeof(&sentinel_index, sizeof(int64_t), 1, cpstream);
    fprintf(stderr, "* sentinel-index: %ld\n", sentinel_index);
    #endif
    fclose(cpstream);

    int64_t x;
    #if !SA_COMPRESSION
    for(x = 0; x < reference_seq_len; x++)
    {
        // fprintf(stderr, "x: %ld\n", x);
        #if SA_COMPRESSION
        if(get_sa_entry_compressed(x) == 0) {
            sentinel_index = x;
            break;
        }
        #else
        if(get_sa_entry(x) == 0) {
            sentinel_index = x;
            break;
        }
        #endif
    }
    fprintf(stderr, "\nsentinel_index: %ld\n", x);    
    #endif
	*_sentinel_index = sentinel_index;

    fprintf(stderr, "* Count:\n");
    for(x = 0; x < 5; x++)
    {
        fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)count[x]);
    }
    fprintf(stderr, "\n");  
	return 0;
}

static int64_t read_rlen_bwt(const char *bwt_file_name) {
	size_t n;
	int64_t ret;
	FILE *f = fopen(bwt_file_name, "rb");
	
	if (f == NULL) return -1;

	err_fread_noeof(&ret, sizeof(int64_t), 1, f);
	fclose(f);

	if (ret <= 0 || ret > (0xffffffffU * (int64_t)CP_BLOCK_SIZE))
		return -1;
	else
		return ret;
}

int64_t __load_BWT_rlen(const char *bwt_file_name) {
	int64_t ret;
	FILE *fp = xopen(bwt_file_name, "rb");
	if (fp == NULL) return -1;

	err_fread_noeof(&ret, sizeof(int64_t), 1, fp);
	err_fclose(fp);

	if (ret <= 0 || ret > (0xffffffffU * (int64_t)CP_BLOCK_SIZE))
		return -1;
	else
		return ret;
}

int __load_BWT_without_shm(char *ref_file_name, int64_t *_reference_seq_len, 
						int64_t *_count, 
						CP_OCC **__cp_occ, 
						int8_t **__sa_ms_byte, uint32_t **__sa_ls_word, 
						int64_t *_sentinel_index)
{
	char cp_file_name[PATH_MAX];
    strcpy_s(cp_file_name, PATH_MAX, ref_file_name);
    strcat_s(cp_file_name, PATH_MAX, CP_FILENAME_SUFFIX);
	
	int64_t reference_seq_len = __load_BWT_rlen(cp_file_name);
	*_reference_seq_len = reference_seq_len;

	int64_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
	CP_OCC *cp_occ = NULL;

    if ((cp_occ = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64)) == NULL) {
        fprintf(stderr, "ERROR! unable to allocated cp_occ memory\n");
        exit(EXIT_FAILURE);
    }

	int8_t *sa_ms_byte;
	uint32_t *sa_ls_word;
	#if SA_COMPRESSION
    int64_t reference_seq_len_ = (reference_seq_len >> SA_COMPX) + 1;
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len_ * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len_ * sizeof(uint32_t), 64);
    #else
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len * sizeof(uint32_t), 64);
    #endif
	
	if (sa_ms_byte == NULL || sa_ls_word == NULL) {
        fprintf(stderr, "ERROR! unable to allocated sa memory\n");
        exit(EXIT_FAILURE);
	}


	__load_BWT_from_file(cp_file_name, reference_seq_len, _count,
						 cp_occ, cp_occ_size,
						 sa_ms_byte, sa_ls_word,
						 _sentinel_index);
	
	*__cp_occ = cp_occ;	
	*__sa_ms_byte = sa_ms_byte;
	*__sa_ls_word = sa_ls_word;
	return 0;
}

int __load_BWT_on_shm(const char *ref_file_name, int64_t *_reference_seq_len, 
						int64_t *_count, 
						CP_OCC **__cp_occ, 
						int8_t **__sa_ms_byte, uint32_t **__sa_ls_word, 
						int64_t *_sentinel_index)
{
	shm_bwt_header_t *header;
	size_t shm_size;
	uint8_t *ptr = NULL;
	int i;
	int fd;
	char cp_file_name[PATH_MAX];
    strcpy_s(cp_file_name, PATH_MAX, ref_file_name);
    strcat_s(cp_file_name, PATH_MAX, CP_FILENAME_SUFFIX);
	

	int64_t reference_seq_len = bwa_shm_rlen();
	if (_reference_seq_len) *_reference_seq_len = reference_seq_len;

	shm_size = bwa_shm_size_bwt(reference_seq_len);

	fprintf(stderr, "INFO: shm_create for BWT index. hugetlb_flag: %x\n", bwa_shm_hugetlb_flags());
	fd = bwa_shm_create(BWA_SHM_BWT, shm_size);
	
	fprintf(stderr, "INFO: BWT shm_size: %ld fd: %d\n", shm_size, fd);
	if (fd >= 0)
		ptr = (uint8_t *) bwa_shm_map(BWA_SHM_BWT);
	else 
		fprintf(stderr, "[bwa_shm] failed to create BWA_SHM_BWT\n");

	if (!ptr)
		return -1;
   
   	header = (shm_bwt_header_t *) ptr;
	ptr += bwa_shm_size_bwt_header();
	header->reference_len = reference_seq_len;

	CP_OCC *cp_occ = (CP_OCC *) ptr;
	ptr += bwa_shm_size_bwt_cp_occ(reference_seq_len);
	int64_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
	
	int8_t *sa_ms_byte = (int8_t *) ptr;
	ptr += bwa_shm_size_bwt_sa_ms_byte(reference_seq_len);
	uint32_t *sa_ls_word = (uint32_t *) ptr;
	
	__load_BWT_from_file(cp_file_name, reference_seq_len,
						header->count,
						cp_occ, cp_occ_size,
						sa_ms_byte, sa_ls_word,
						&header->sentinel_index);

	if (__cp_occ) *__cp_occ = cp_occ;
	if (__sa_ms_byte) *__sa_ms_byte = sa_ms_byte;
	if (__sa_ls_word) *__sa_ls_word = sa_ls_word;
	if (_count) {
		for (i = 0; i < 5; ++i)
			_count[i] = header->count[i];
	}
	if (_sentinel_index) *_sentinel_index = header->sentinel_index;

	return 0;
}

int __load_BWT_from_shm(int64_t *_reference_seq_len,
						int64_t *_count, 
						CP_OCC **__cp_occ, 
						int8_t **__sa_ms_byte, uint32_t **__sa_ls_word, 
						int64_t *_sentinel_index)
{
	int64_t cp_occ_size;
	int64_t x;
	uint8_t *ptr;
	shm_bwt_header_t *header;

	ptr = (uint8_t *) bwa_shm_map(BWA_SHM_BWT);
	if (!ptr) return -1;

	header = (shm_bwt_header_t *) ptr;
	ptr += bwa_shm_size_bwt_header();

	*_reference_seq_len = header->reference_len;
	for (x = 0; x < 5; ++x)
		_count[x] = header->count[x];
	*_sentinel_index = header->sentinel_index;
    
	fprintf(stderr, "* sentinel-index: %ld\n", header->sentinel_index);
    
	fprintf(stderr, "* Count:\n");
    for(x = 0; x < 5; x++)
    {
        fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)_count[x]);
    }
    fprintf(stderr, "\n");  

	*__cp_occ = (CP_OCC *)ptr;
	ptr += bwa_shm_size_bwt_cp_occ(header->reference_len);
	*__sa_ms_byte = (int8_t *) ptr;
	ptr += bwa_shm_size_bwt_sa_ms_byte(header->reference_len);
	*__sa_ls_word = (uint32_t *) ptr;
	return 0;
}

static void load_BWT(char *ref_file_name, int64_t *_reference_seq_len, 
					int64_t *_count, 
					CP_OCC **__cp_occ, 
					int8_t **__sa_ms_byte, uint32_t **__sa_ls_word, 
					int64_t *_sentinel_index) 
{
	if (bwa_shm_mode == BWA_SHM_MATCHED) {
		if (__load_BWT_from_shm(_reference_seq_len, _count,
								__cp_occ, __sa_ms_byte, __sa_ls_word,
								_sentinel_index) == 0)
			return;
	} 
	
	if (bwa_shm_mode == BWA_SHM_RENEWAL) {
		if (__load_BWT_on_shm(ref_file_name, _reference_seq_len,
							_count, __cp_occ, __sa_ms_byte, __sa_ls_word,
							_sentinel_index) == 0)
			return;
	}

	__load_BWT_without_shm(ref_file_name, _reference_seq_len,
							_count, __cp_occ, __sa_ms_byte, __sa_ls_word,
							_sentinel_index);
}
#endif

#ifdef USE_SHM
size_t ____size_mlt(const char *prefix, const char *ref_file_name) {
	char path[PATH_MAX];
	FILE *fp = NULL;
	long ret = 0;
	
	if (prefix) {
		strcpy_s(path, PATH_MAX, prefix);
	} else if (ref_file_name) {
		size_t len;
		strcpy_s(path, PATH_MAX, ref_file_name);
		len = strnlen_s(path, PATH_MAX);
		if (len < 5) /* while initialization */
			return 0;
		// path = prefix + ".0123"
		path[len - 5] = '\0';
	} else 
		return 0;
	
	strcat_s(path, PATH_MAX, ".mlt_table");
	
	fp = fopen(path, "rb");
	if (!fp) {
		fprintf(stderr, "%s: failed to fopen %s errno: %d\n", __func__, path, errno);
		return 0L;
	}
	
	if (fseek(fp, 0L, SEEK_END) != 0) {
		fprintf(stderr, "%s: failed to fseek %s errno: %d\n", __func__, path, errno);
		return 0L;
	}
	
	ret = ftell(fp);
	if (ret < 0) {
		fprintf(stderr, "%s: failed to ftell %s errno: %d\n", __func__, path, errno);
		return 0L;
	}
	
	fclose(fp);
	return ret;
}

#define __size_mlt(prefix) ____size_mlt(prefix, NULL)

int _load_ert_index(const char *prefix, uint64_t **__kmer_offsets, uint8_t **__mlt_table) {
	if (__bwa_shm_load_file(prefix, ".kmer_table", BWA_SHM_KMER, (void **) __kmer_offsets))
		return -1;

   	if (__bwa_shm_load_file(prefix, ".mlt_table", BWA_SHM_MLT, (void **) __mlt_table))
		return -1;

	return 0;
}

int _load_kmer_table(const char *prefix, uint64_t **__kmer) {
	return __bwa_shm_load_file(prefix, ".kmer_table", BWA_SHM_KMER, (void **) __kmer);
}

int _load_mlt_table(const char *prefix, uint8_t **__mlt) {
   	return __bwa_shm_load_file(prefix, ".mlt_table", BWA_SHM_MLT, (void **) __mlt);
}


void FMI_search::load_ert_index() {
	int64_t allocMem = 0;
	size_t mlt_size = 0;
    double ctime, rtime;
    ctime = cputime(); rtime = realtime();

    fprintf(stderr, "[M::%s::ERT] Reading kmer index to memory\n", __func__);
	
	allocMem = numKmers * sizeof(uint64_t) + __size_mlt(file_name); 
	if (_load_kmer_table(file_name, &kmer_offsets))
		exit(EXIT_FAILURE);
	
	if (_load_mlt_table(file_name, &mlt_table))
		exit(EXIT_FAILURE);

    fprintf(stderr, "[M::%s::ERT] Index tables (%0.4lfGB) loaded in %.3f CPU sec, %.3f real sec...\n", 
					__func__, allocMem/1e9, cputime() - ctime, realtime() - rtime);

	useErt = 1;
}
#else
void FMI_search::load_ert_index() {
	int64_t allocMem = 0;
	size_t mlt_size = 0;
    double ctime, rtime;
    ctime = cputime(); rtime = realtime();
    
	fprintf(stderr, "[M::%s::ERT] Reading kmer index to memory\n", __func__);

	allocMem = numKmers * sizeof(uint64_t); 
	kmer_offsets = (uint64_t *) __load_file(file_name, ".kmer_table", NULL, NULL);

   	mlt_table = (uint8_t *) __load_file(file_name, ".mlt_table", NULL, &mlt_size);
	allocMem += mlt_size;

    fprintf(stderr, "[M::%s::ERT] Index tables (%0.4lfGB) loaded in %.3f CPU sec, %.3f real sec...\n", 
					__func__, allocMem/1e9, cputime() - ctime, realtime() - rtime);

	useErt = 1;
}
#endif

void FMI_search::load_index()
{
    one_hot_mask_array = (uint64_t *)_mm_malloc(64 * sizeof(uint64_t), 64);
    one_hot_mask_array[0] = 0;
    uint64_t base = 0x8000000000000000L;
    one_hot_mask_array[1] = base;
    int64_t i = 0;
    for(i = 2; i < 64; i++)
    {
        one_hot_mask_array[i] = (one_hot_mask_array[i - 1] >> 1) | base;
    }

    char *ref_file_name = file_name;
#ifdef USE_SHM
	load_BWT(ref_file_name, &reference_seq_len, count, 
			&cp_occ, &sa_ms_byte, &sa_ls_word, &sentinel_index);
#else
    //beCalls = 0;
    char cp_file_name[PATH_MAX];
    strcpy_s(cp_file_name, PATH_MAX, ref_file_name);
    strcat_s(cp_file_name, PATH_MAX, CP_FILENAME_SUFFIX);

    // Read the BWT and FM index of the reference sequence
    FILE *cpstream = NULL;
    cpstream = fopen(cp_file_name,"rb");
    if (cpstream == NULL)
    {
        fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stderr, "* Index file found. Loading index from %s\n", cp_file_name);
    }

    err_fread_noeof(&reference_seq_len, sizeof(int64_t), 1, cpstream);
    assert(reference_seq_len > 0);
    assert(reference_seq_len <= 0x7fffffffffL);

    fprintf(stderr, "* Reference seq len for bi-index = %ld\n", reference_seq_len);

    // create checkpointed occ
    int64_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
    cp_occ = NULL;

    err_fread_noeof(&count[0], sizeof(int64_t), 5, cpstream);
    if ((cp_occ = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64)) == NULL) {
        fprintf(stderr, "ERROR! unable to allocated cp_occ memory\n");
        exit(EXIT_FAILURE);
    }

    err_fread_noeof(cp_occ, sizeof(CP_OCC), cp_occ_size, cpstream);
    int64_t ii = 0;
    for(ii = 0; ii < 5; ii++)// update read count structure
    {
        count[ii] = count[ii] + 1;
    }

    #if SA_COMPRESSION

    int64_t reference_seq_len_ = (reference_seq_len >> SA_COMPX) + 1;
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len_ * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len_ * sizeof(uint32_t), 64);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len_, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len_, cpstream);
    
    #else
    
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len * sizeof(uint32_t), 64);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len, cpstream);

    #endif

    sentinel_index = -1;
    #if SA_COMPRESSION
    err_fread_noeof(&sentinel_index, sizeof(int64_t), 1, cpstream);
    fprintf(stderr, "* sentinel-index: %ld\n", sentinel_index);
    #endif
    fclose(cpstream);

    int64_t x;
    #if !SA_COMPRESSION
    for(x = 0; x < reference_seq_len; x++)
    {
        if(get_sa_entry(x) == 0) {
            sentinel_index = x;
            break;
        }
    }
    fprintf(stderr, "\nsentinel_index: %ld\n", x);    
    #endif

    fprintf(stderr, "* Count:\n");
    for(x = 0; x < 5; x++)
    {
        fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)count[x]);
    }
    fprintf(stderr, "\n");  
#endif /* !USE_SHM */

#ifdef SMEM_ACCEL
	if (building_smem_table == 0) {
		fprintf(stderr, "* Reading data for smem acceleration\n");
		load_smem_table();
	} else {
		all_smem_table = NULL;
		last_smem_table = NULL;
	}
#endif

    fprintf(stderr, "* Reading other elements of the index from files %s\n",
            ref_file_name);
    bwa_idx_load_ele(ref_file_name, BWA_IDX_ALL);

	useErt = 0;
    
	fprintf(stderr, "* Done reading Index!!\n");
}

void FMI_search::getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                         int16_t *query_pos_array,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    int64_t numTotalSmem = *__numTotalSmem;
    SMEM prevArray[max_readlength];

    uint32_t i;
    // Perform SMEM for original reads
    for(i = 0; i < numReads; i++)
    {
        int x = query_pos_array[i];
        int32_t rid = rid_array[i];
        int next_x = x + 1;

        int readlength = seq_[rid].l_seq;
        int offset = query_cum_len_ar[rid];
        // uint8_t a = enc_qdb[rid * readlength + x];
        uint8_t a = enc_qdb[offset + x];

        if(a < 4)
        {
			debug_smem_input("OnePosOne", rid, &enc_qdb[offset], x, readlength);
            SMEM smem;
            smem.rid = rid;
            smem.m = x;
            smem.n = x;
            smem.k = count[a];
            smem.l = count[3 - a];
            smem.s = count[a+1] - count[a];
            int numPrev = 0;
            
			int j;
#ifdef SMEM_ACCEL
#ifdef MEMSCALE
			if (all_smem_table && readlength - x >= ALL_SMEM_MAX_BP) 
#else
			if (readlength - x >= ALL_SMEM_MAX_BP) 
#endif
			{
				uint64_t all_smem_idx = 0;
				int k, last_idx, with_N = 0;
				uint8_t *enc = &enc_qdb[offset + x];
				for (k = 0; k < ALL_SMEM_MAX_BP; ++k, ++enc) {
					if ((*enc) >= 4) break;
					all_smem_idx = all_smem_idx | ((*enc) << (((ALL_SMEM_MAX_BP - 1) - k) * 2));
				}

				all_smem_t *ent = &all_smem_table[all_smem_idx];
				with_N = k < ALL_SMEM_MAX_BP ? 1 : 0;
				last_idx = (k > ent->last_avail ? ent->last_avail : k) - 1;

				for (j = x + 1, k = 0; k < last_idx; ++j, ++k) {
					a = enc_qdb[offset + j];
					next_x = j + 1;

					SMEM newSmem = smem;
					newSmem.k = smem.k + ent->list[k].k32;
					newSmem.l = count[3 - a] + ent->list[k].l32;
					newSmem.s = ent->list[k].s32;
					newSmem.n = j;

					int32_t s_neq_mask = newSmem.s != smem.s;

					prevArray[numPrev] = smem;
					debug_smem_output("OnePosOne-fwd-prev", &smem, 1, &s_neq_mask);
					numPrev += s_neq_mask;
					if (newSmem.s < min_intv_array[i])
					{
						next_x = j;
						j = readlength; // to skip the below loop
						break;
					}
					smem = newSmem;
				}

				if (with_N) {
					next_x = j + 1;
					j = readlength;
				}
			} else {
				j = x + 1;
			}

            for( ; j < readlength; j++)
#else
            for(j = x + 1; j < readlength; j++)
#endif
            {
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];
                next_x = j + 1;
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;
                    SMEM newSmem_ = backwardExt(smem_, 3 - a);
                    //SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

                    int32_t s_neq_mask = newSmem.s != smem.s;

                    prevArray[numPrev] = smem;
					debug_smem_output("OnePosOne-fwd-prev", &smem, 1, &s_neq_mask);
                    numPrev += s_neq_mask;
                    if(newSmem.s < min_intv_array[i])
                    {
                        next_x = j;
                        break;
                    }
                    smem = newSmem;
#ifdef ENABLE_PREFETCH
                    _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                }
                else
                {
                    break;
                }
            }
            if(smem.s >= min_intv_array[i])
            {

                prevArray[numPrev] = smem;
				debug_smem_output("OnePosOne-fwd-prev2", &smem, 1, &min_intv_array[i]);
                numPrev++;
            }

            SMEM *prev;
            prev = prevArray;

            int p;
            for(p = 0; p < (numPrev/2); p++)
            {
                SMEM temp = prev[p];
                prev[p] = prev[numPrev - p - 1];
                prev[numPrev - p - 1] = temp;
            }

            // Backward search
            int cur_j = readlength;
            for(j = x - 1; j >= 0; j--)
            {
                int numCurr = 0;
                int curr_s = -1;
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];

                if(a > 3)
                {
                    break;
                }
                for(p = 0; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;

                    if((newSmem.s < min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                    {
                        cur_j = j;

                        matchArray[numTotalSmem++] = smem;
						debug_smem_output("OnePosOne-match", &smem, 1, &minSeedLen);
                        break;
                    }
                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                        break;
                    }
                }
                p++;
                for(; p < numPrev; p++)
                {
                    SMEM smem = prev[p];

                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;


                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
						debug_smem_output("OnePosOne-BWD", &newSmem, 1, &min_intv_array[i]);
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                    }
                }
                numPrev = numCurr;
                if(numCurr == 0)
                {
                    break;
                }
            }
            if(numPrev != 0)
            {
                SMEM smem = prev[0];
                if(((smem.n - smem.m + 1) >= minSeedLen))
                {

                    matchArray[numTotalSmem++] = smem;
					debug_smem_output("OnePosOne-match", &smem, 1, &minSeedLen);
                }
                numPrev = 0;
            }
        }
        query_pos_array[i] = next_x;
    }
    (*__numTotalSmem) = numTotalSmem;
}

void FMI_search::getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    int16_t *query_pos_array = (int16_t *)_mm_malloc(numReads * sizeof(int16_t), 64);
    
    int32_t i;
    for(i = 0; i < numReads; i++)
        query_pos_array[i] = 0;

    int32_t numActive = numReads;
    (*__numTotalSmem) = 0;

    do
    {
        int32_t head = 0;
        int32_t tail = 0;
        for(head = 0; head < numActive; head++)
        {
            int readlength = seq_[rid_array[head]].l_seq;
            if(query_pos_array[head] < readlength)
            {
                rid_array[tail] = rid_array[head];
                query_pos_array[tail] = query_pos_array[head];
                min_intv_array[tail] = min_intv_array[head];
                tail++;             
            }               
        }
        getSMEMsOnePosOneThread(enc_qdb,
                                query_pos_array,
                                min_intv_array,
                                rid_array,
                                tail,
                                batch_size,
                                seq_,
                                query_cum_len_ar,
                                max_readlength,
                                minSeedLen,
                                matchArray,
                                __numTotalSmem);
        numActive = tail;
    } while(numActive > 0);

    _mm_free(query_pos_array);
}

int64_t FMI_search::bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                                   int32_t *max_intv_array,
                                                   int32_t numReads,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray)
{
#if defined(PERFECT_MATCH) && !defined(DO_NORMAL)
	int32_t pos = -1; /* for max_intv_array */
#endif
    int32_t i;

    int64_t numTotalSeed = 0;

    for(i = 0; i < numReads; i++)
    {
#if defined(PERFECT_MATCH) && !defined(DO_NORMAL)
		if (seq_[i].perfect.exist)
			continue;
		pos++;
#endif
        int readlength = seq_[i].l_seq;
        int16_t x = 0;
        while(x < readlength)
        {
            int next_x = x + 1;

            // Forward search
            SMEM smem;
            smem.rid = i;
            smem.m = x;
            smem.n = x;
            
            int offset = query_cum_len_ar[i];
            uint8_t a = enc_qdb[offset + x];
            // uint8_t a = enc_qdb[i * readlength + x];

            if(a < 4)
            {
				debug_smem_input("BWTSeed", i, &enc_qdb[offset], x, readlength);
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];

                int j;
#ifdef SMEM_ACCEL
#ifdef MEMSCALE
				if (last_smem_table && readlength - x >= LAST_SMEM_MAX_BP) 
#else
				if (readlength - x >= LAST_SMEM_MAX_BP) 
#endif
				{
					// TODO: use last_smem_table
					uint64_t last_smem_idx = 0;
					int with_N = 0;
					int k = 0;
					uint8_t *enc = &enc_qdb[offset + x];

					for (k = 0; k < LAST_SMEM_MAX_BP; ++k, ++enc) {
						last_smem_idx = last_smem_idx | ((*enc) << (((LAST_SMEM_MAX_BP - 1) - k) * 2));
						with_N += ((*enc) >> 2); // check *enc >= 4
					}

					if (with_N == 0) {
						last_smem_t *ent = &last_smem_table[last_smem_idx];
						
						j = x + ent->bp;
						next_x = j; /* if seed_end == readlength, next_x = seed_end,
									   otherwise, next_x will be set at the below */

						smem.k = __combine_ms_ls(ent->kms, ent->kls);
						smem.l = __combine_ms_ls(ent->lms, ent->lls);
						smem.s = __combine_ms_ls(ent->sms, ent->sls);
						smem.n = j - 1;

#if defined(PERFECT_MATCH) && !defined(DO_NORMAL)
                        if((smem.s < max_intv_array[pos]) && ((smem.n - smem.m + 1) >= minSeedLen))
#else
                        if((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
#endif
                        {

                            if(smem.s > 0) {
                                matchArray[numTotalSeed++] = smem;
								//debug_smem_output("BWTSeed-MATCH1", &smem, 1, &minSeedLen);
								debug_smem_output("BWTSeed-MATCH", &smem, 1, &minSeedLen);
							}
                        }

					} else {
						/* ignore the entry */
						j = x + 1;
					}
				} else {
					j = x + 1;
				}

                for( ; j < readlength; j++)
#else
                for(j = x + 1; j < readlength; j++)
#endif
                {
                    next_x = j + 1;
                    // a = enc_qdb[i * readlength + j];
                    a = enc_qdb[offset + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        //SMEM smem = backwardExt(smem, 3 - a);
                        //smem.n = j;
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;
                        smem = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif

#if defined(PERFECT_MATCH) && !defined(DO_NORMAL)
                        if((smem.s < max_intv_array[pos]) && ((smem.n - smem.m + 1) >= minSeedLen))
#else
                        if((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
#endif
                        {

                            if(smem.s > 0)
                            {
                                matchArray[numTotalSeed++] = smem;
								//debug_smem_output("BWTSeed-MATCH2", &smem, 1, &minSeedLen);
								debug_smem_output("BWTSeed-MATCH", &smem, 1, &minSeedLen);
                            }
                            break;
                        }
                    }
                    else
                    {

                        break;
                    }
                }

            }
            x = next_x;
        }
    }
    return numTotalSeed;
}


void FMI_search::getSMEMs(uint8_t *enc_qdb,
        int32_t numReads,
        int32_t batch_size,
        int32_t readlength,
        int32_t minSeedLen,
        int32_t nthreads,
        SMEM *matchArray,
        int64_t *numTotalSmem)
{
    SMEM *prevArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);
    SMEM *currArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);


// #pragma omp parallel num_threads(nthreads)
    {
        int tid = 0; //omp_get_thread_num();   // removed omp
        numTotalSmem[tid] = 0;
        SMEM *myPrevArray = prevArray + tid * readlength;
        SMEM *myCurrArray = prevArray + tid * readlength;

        int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
        int32_t first = tid * perThreadQuota;
        int32_t last  = (tid + 1) * perThreadQuota;
        if(last > numReads) last = numReads;
        SMEM *myMatchArray = matchArray + first * readlength;

        uint32_t i;
        // Perform SMEM for original reads
        for(i = first; i < last; i++)
        {
            int x = readlength - 1;
            int numPrev = 0;
            int numSmem = 0;

            while (x >= 0)
            {
                // Forward search
                SMEM smem;
                smem.rid = i;
                smem.m = x;
                smem.n = x;
                uint8_t a = enc_qdb[i * readlength + x];

                if(a > 3)
                {
                    x--;
                    continue;
                }
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];

                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    a = enc_qdb[i * readlength + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;

                        if(newSmem.s != smem.s)
                        {
                            myPrevArray[numPrev] = smem;
                            numPrev++;
                        }
                        smem = newSmem;
                        if(newSmem.s == 0)
                        {
                            break;
                        }
                    }
                    else
                    {
                        myPrevArray[numPrev] = smem;
                        numPrev++;
                        break;
                    }
                }
                if(smem.s != 0)
                {
                    myPrevArray[numPrev++] = smem;
                }

                SMEM *curr, *prev;
                prev = myPrevArray;
                curr = myCurrArray;

                int p;
                for(p = 0; p < (numPrev/2); p++)
                {
                    SMEM temp = prev[p];
                    prev[p] = prev[numPrev - p - 1];
                    prev[numPrev - p - 1] = temp;
                }

                int next_x = x - 1;

                // Backward search
                int cur_j = readlength;
                for(j = x - 1; j >= 0; j--)
                {
                    int numCurr = 0;
                    int curr_s = -1;
                    a = enc_qdb[i * readlength + j];
                    //printf("a = %d\n", a);
                    if(a > 3)
                    {
                        next_x = j - 1;
                        break;
                    }
                    for(p = 0; p < numPrev; p++)
                    {
                        SMEM smem = prev[p];
                        SMEM newSmem = backwardExt(smem, a);
                        newSmem.m = j;

                        if(newSmem.s == 0)
                        {
                            if((numCurr == 0) && (j < cur_j))
                            {
                                cur_j = j;
                                if((smem.n - smem.m + 1) >= minSeedLen)
                                    myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                            }
                        }
                        if((newSmem.s != 0) && (newSmem.s != curr_s))
                        {
                            curr_s = newSmem.s;
                            curr[numCurr++] = newSmem;
                        }
                    }
                    SMEM *temp = prev;
                    prev = curr;
                    curr = temp;
                    numPrev = numCurr;
                    if(numCurr == 0)
                    {
                        next_x = j;
                        break;
                    }
                    else
                    {
                        next_x = j - 1;
                    }
                }
                if(numPrev != 0)
                {
                    SMEM smem = prev[0];
                    if((smem.n - smem.m + 1) >= minSeedLen)
                        myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                    numPrev = 0;
                }
                x = next_x;
            }
            numTotalSmem[tid] += numSmem;
        }
    }

    _mm_free(prevArray);
    _mm_free(currArray);
}


int compare_smem(const void *a, const void *b)
{
    SMEM *pa = (SMEM *)a;
    SMEM *pb = (SMEM *)b;

    if(pa->rid < pb->rid)
        return -1;
    if(pa->rid > pb->rid)
        return 1;

    if(pa->m < pb->m)
        return -1;
    if(pa->m > pb->m)
        return 1;
    if(pa->n > pb->n)
        return -1;
    if(pa->n < pb->n)
        return 1;
    return 0;
}

void FMI_search::sortSMEMs(SMEM *matchArray,
        int64_t numTotalSmem[],
        int32_t numReads,
        int32_t readlength,
        int nthreads)
{
    int tid;
    int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
    for(tid = 0; tid < nthreads; tid++)
    {
        int32_t first = tid * perThreadQuota;
        SMEM *myMatchArray = matchArray + first * readlength;
        qsort(myMatchArray, numTotalSmem[tid], sizeof(SMEM), compare_smem);
    }
}


SMEM FMI_search::backwardExt(SMEM smem, uint8_t a)
{
    //beCalls++;
    uint8_t b;

    int64_t k[4], l[4], s[4];
    for(b = 0; b < 4; b++)
    {
        int64_t sp = (int64_t)(smem.k);
        int64_t ep = (int64_t)(smem.k) + (int64_t)(smem.s);
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        GET_OCC(ep, b, occ_id_ep, y_ep, occ_ep, one_hot_bwt_str_c_ep, match_mask_ep);
        k[b] = count[b] + occ_sp;
        s[b] = occ_ep - occ_sp;
    }

    int64_t sentinel_offset = 0;
    if((smem.k <= sentinel_index) && ((smem.k + smem.s) > sentinel_index)) sentinel_offset = 1;
    l[3] = smem.l + sentinel_offset;
    l[2] = l[3] + s[3];
    l[1] = l[2] + s[2];
    l[0] = l[1] + s[1];

    smem.k = k[a];
    smem.l = l[a];
    smem.s = s[a];
    return smem;
}

int64_t FMI_search::get_sa_entry(int64_t pos)
{
    int64_t sa_entry = sa_ms_byte[pos];
    sa_entry = sa_entry << 32;
    sa_entry = sa_entry + sa_ls_word[pos];
    return sa_entry;
}

void FMI_search::get_sa_entries(int64_t *posArray, int64_t *coordArray, uint32_t count, int32_t nthreads)
{
    uint32_t i;
// #pragma omp parallel for num_threads(nthreads)
    for(i = 0; i < count; i++)
    {
        int64_t pos = posArray[i];
        int64_t sa_entry = sa_ms_byte[pos];
        sa_entry = sa_entry << 32;
        sa_entry = sa_entry + sa_ls_word[pos];
        //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD), _MM_HINT_T0);
        coordArray[i] = sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ)
{
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = sa_ms_byte[pos];
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[pos];
            //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD * step), _MM_HINT_T0);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        coordCountArray[i] = c;
        totalCoordCount += c;
    }
}

// sa_compression
int64_t FMI_search::get_sa_entry_compressed(int64_t pos, int tid)
{
    if ((pos & SA_COMPX_MASK) == 0) {
        
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[pos >> SA_COMPX];
        #else
        int64_t sa_entry = sa_ms_byte[pos];     // simulation
        #endif
        
        sa_entry = sa_entry << 32;
        
        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];
        #else
        sa_entry = sa_entry + sa_ls_word[pos];   // simulation
        #endif
        
        return sa_entry;        
    }
    else {
        // tprof[MEM_CHAIN][tid] ++;
        int64_t offset = 0; 
        int64_t sp = pos;
        while(true)
        {
            int64_t occ_id_pp_ = sp >> CP_SHIFT;
            int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
            uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
            uint8_t b;

            if((one_hot_bwt_str[0] >> y_pp_) & 1)
                b = 0;
            else if((one_hot_bwt_str[1] >> y_pp_) & 1)
                b = 1;
            else if((one_hot_bwt_str[2] >> y_pp_) & 1)
                b = 2;
            else if((one_hot_bwt_str[3] >> y_pp_) & 1)
                b = 3;
            else
                b = 4;

            if (b == 4) {
                return offset;
            }

            GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);

            sp = count[b] + occ_sp;
            
            offset ++;
            // tprof[ALIGN1][tid] ++;
            if ((sp & SA_COMPX_MASK) == 0) break;
        }
        // assert((reference_seq_len >> SA_COMPX) - 1 >= (sp >> SA_COMPX));
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[sp >> SA_COMPX];
        #else
        int64_t sa_entry = sa_ms_byte[sp];      // simultion
        #endif
        
        sa_entry = sa_entry << 32;

        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
        #else
        sa_entry = sa_entry + sa_ls_word[sp];      // simulation
        #endif
        
        sa_entry += offset;
        return sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ, int tid)
{
    
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        // coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
}

// SA_COPMRESSION w/ PREFETCH
int64_t FMI_search::call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset)
{
    if ((pos & SA_COMPX_MASK) == 0) {        
        sa_entry = sa_ms_byte[pos >> SA_COMPX];        
        sa_entry = sa_entry << 32;        
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];        
        // return sa_entry;
        return 1;
    }
    else {
        // int64_t offset = 0; 
        int64_t sp = pos;

        int64_t occ_id_pp_ = sp >> CP_SHIFT;
        int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
        uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
        uint8_t b;

        if((one_hot_bwt_str[0] >> y_pp_) & 1)
            b = 0;
        else if((one_hot_bwt_str[1] >> y_pp_) & 1)
            b = 1;
        else if((one_hot_bwt_str[2] >> y_pp_) & 1)
            b = 2;
        else if((one_hot_bwt_str[3] >> y_pp_) & 1)
            b = 3;
        else
            b = 4;
        if (b == 4) {
            sa_entry = 0;
            return 1;
        }
        
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        
        sp = count[b] + occ_sp;
        
        offset ++;
        if ((sp & SA_COMPX_MASK) == 0) {
    
            sa_entry = sa_ms_byte[sp >> SA_COMPX];        
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
            
            sa_entry += offset;
            // return sa_entry;
            return 1;
        }
        else {
            sa_entry = sp;
            return 0;
        }
    } // else
}

void FMI_search::get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                         int64_t *coordCountArray, int64_t count,
                                         const int32_t max_occ, int tid, int64_t &id_)
{
    
    // uint32_t i;
    int32_t totalCoordCount = 0;
    int32_t mem_lim = 0, id = 0;
    
    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        mem_lim += smem.s;
    }

    int64_t *pos_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);
    int64_t *map_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);

    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
             pos_ar[id]  = pos;
             map_ar[id++] = totalCoordCount + c;
            // int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            // coordArray[totalCoordCount + c] = sa_entry;
        }
        //coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
    
    id_ += id;
    
    const int32_t sa_batch_size = 20;
    int64_t working_set[sa_batch_size], map_pos[sa_batch_size];;
    int64_t offset[sa_batch_size] = {-1};
    
    int i = 0, j = 0;    
    while(i<id && j<sa_batch_size)
    {
        int64_t pos =  pos_ar[i];
        working_set[j] = pos;
        map_pos[j] = map_ar[i];
        offset[j] = 0;
        
        if (pos & SA_COMPX_MASK == 0) {
            _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
            _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
        }
        else {
            int64_t occ_id_pp_ = pos >> CP_SHIFT;
            _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
        }
        i++;
        j++;
    }
        
    int lim = j, all_quit = 0;
    while (all_quit < id)
    {
        
        for (int k=0; k<lim; k++)
        {
            int64_t sp = 0, pos = 0;
            bool quit;
            if (offset[k] >= 0) {
                quit = call_one_step(working_set[k], sp, offset[k]);
            }
            else
                continue;
            
            if (quit) {
                coordArray[map_pos[k]] = sp;
                all_quit ++;
                
                if (i < id)
                {
                    pos = pos_ar[i];
                    working_set[k] = pos;
                    map_pos[k] = map_ar[i++];
                    offset[k] = 0;
                    
                    if (pos & SA_COMPX_MASK == 0) {
                        _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
                        _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
                    }
                    else {
                        int64_t occ_id_pp_ = pos >> CP_SHIFT;
                        _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                    }
                }
                else
                    offset[k] = -1;
            }
            else {
                working_set[k] = sp;
                if (sp & SA_COMPX_MASK == 0) {
                    _mm_prefetch(&sa_ms_byte[sp >> SA_COMPX], _MM_HINT_T0);
                    _mm_prefetch(&sa_ls_word[sp >> SA_COMPX], _MM_HINT_T0);
                }
                else {
                    int64_t occ_id_pp_ = sp >> CP_SHIFT;
                    _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                }                
            }
        }
    }
    
    _mm_free(pos_ar);
    _mm_free(map_ar);
}
