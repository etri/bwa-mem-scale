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

#ifdef USE_SHM
#include "macro.h"
#include "bwa_shm.h"
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/mount.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <mntent.h>
#include <errno.h>
#include "FMI_search.h"
#include "fastmap.h"
#include "safe_lib.h"
#include <fcntl.h>
#ifdef PERFECT_MATCH
#include "perfect.h"
#endif
#include <errno.h>

#define KB_UNIT_STR_1GB "1048576kB"
#define KB_UNIT_STR_2MB "2048kB"

/* TODO: default huge page size from the system */
#define KB_UNIT_STR_DEFAULT KB_UNIT_STR_2MB
#define DEFAULT_HUGETLB_PAGESIZE (2L << 20) 

#define B2GB(x) ((x + (1 << 30) - 1) >> 30)
#define B2GB_DOUBLE(x) (((double) x) / (1LL << 30))

bwa_shm_info_t *loading_info = NULL;
int opt_bwa_shm_map_touch = 0;

static const char *mmap_prefix = NULL;

void *__load_file(const char *prefix, const char *postfix, void *buf, size_t *size);
int __bwa_shm_load(const char *prefix, enum hugetlb_mode huge_mode, int huge_force,
			int pt_seed_len, int pt_mmap, size_t gb_limit);

int use_mmap(int m) {
#ifdef PERFECT_MATCH
	if (m == BWA_SHM_PERFECT) {
		if (loading_info)
			return loading_info->pt_mmap;
		else if (bwa_shm_info)
			return bwa_shm_info->pt_mmap;
		else
			return 0;
	}
#endif
	return 0;
}

inline void set_mmap_prefix(const char *prefix) {
	mmap_prefix = prefix;
}

char *bwa_shm_mmap_filename(int m, char *buf) {
	if (m == BWA_SHM_BWT) {
		strcpy_s(buf, PATH_MAX, mmap_prefix);
		strcat_s(buf, PATH_MAX, CP_FILENAME_SUFFIX);
		return buf;
	} else if (m == BWA_SHM_REF) {
		strcpy_s(buf, PATH_MAX, mmap_prefix);
		strcat_s(buf, PATH_MAX, ".0123");
		return buf;
	} else if (m == BWA_SHM_PAC) {
		strcpy_s(buf, PATH_MAX, mmap_prefix);
		strcat_s(buf, PATH_MAX, ".pac");
		return buf;
	} else if (m == BWA_SHM_KMER) {
		strcpy_s(buf, PATH_MAX, mmap_prefix);
		strcat_s(buf, PATH_MAX, ".kmer_table");
		return buf;
	} else if (m == BWA_SHM_MLT) {
		strcpy_s(buf, PATH_MAX, mmap_prefix);
		strcat_s(buf, PATH_MAX, ".mlt_table");
		return buf;
#ifdef PERFECT_MATCH
	} else if (m == BWA_SHM_PERFECT) {
		if (perfect_table_seed_len <= 0)
			return NULL;
		snprintf(buf, PATH_MAX, "%s.perfect.%d", mmap_prefix, perfect_table_seed_len);
		return buf;
#endif
#ifdef SMEM_ACCEL
	} else if (m == BWA_SHM_SALL) {
		snprintf(buf, PATH_MAX, "%s.all_smem.%d", mmap_prefix, ALL_SMEM_MAX_BP);
		return buf;
	} else if (m == BWA_SHM_SLAST) {
		snprintf(buf, PATH_MAX, "%s.all_smem.%d", mmap_prefix, LAST_SMEM_MAX_BP);
		return buf;
#endif
	} else {
		return NULL;
	}
}

static const char *bwa_shm_type_str[] = {
	"INFO",
	"BWT",
	"REF",
	"PAC",
	"KMER",
	"MLT",
#ifdef PERFECT_MATCH
	"PERFECT",
#endif
#ifdef SMEM_ACCEL
	"SMEM_ALL",
	"SMEM_LAST",
#endif
	"others"
};

static const char *bwa_shm_name_str[] = {
	"bwa_mem_large_shm",
	"bwa_mem_large_bwt",
	"bwa_mem_large_ref",
	"bwa_mem_large_pac",
	"bwa_mem_large_kmer",
	"bwa_mem_large_mlt",
#ifdef PERFECT_MATCH
	"bwa_mem_large_perfect",
#endif
#ifdef SMEM_ACCEL
	"bwa_mem_large_smem_all",
	"bwa_mem_large_smem_last",
#endif
};

#define BWA_SHM_HUGE_DIR "/bwa_shm_huge"
static const char *bwa_shm_huge_name_str[] = {
	BWA_SHM_HUGE_DIR "/bwa_mem_large_shm", /* actually, not used */
	BWA_SHM_HUGE_DIR "/bwa_mem_large_bwt",
	BWA_SHM_HUGE_DIR "/bwa_mem_large_ref",
	BWA_SHM_HUGE_DIR "/bwa_mem_large_pac",
	BWA_SHM_HUGE_DIR "/bwa_mem_large_kmer",
	BWA_SHM_HUGE_DIR "/bwa_mem_large_mlt",
#ifdef PERFECT_MATCH
	BWA_SHM_HUGE_DIR "/bwa_mem_large_perfect",
#endif
#ifdef SMEM_ACCEL
	BWA_SHM_HUGE_DIR "/bwa_mem_large_smem_all",
	BWA_SHM_HUGE_DIR "/bwa_mem_large_smem_last",
#endif
};

enum bwa_shm_mode bwa_shm_mode = BWA_SHM_DISABLE;

int shm_fd[NUM_BWA_SHM];
void *shm_ptr[NUM_BWA_SHM];
#define BWA_SHM_INFO_FD  shm_fd[BWA_SHM_INFO]
#define BWA_SHM_INFO_PTR shm_ptr[BWA_SHM_INFO]
bwa_shm_info_t *bwa_shm_info;

static inline int lock_bwa_shm_info() {
	int ret;	
	if (!bwa_shm_info) return 0;
	__sync_synchronize();
	if (bwa_shm_info->state == BWA_SHM_STATE_NOT_INIT)
		return 0;
retry:
	ret = __sync_val_compare_and_swap(&bwa_shm_info->lock, 0, 1);
	if (ret == 0) {
		__sync_synchronize();
		return 1;
	} else
		goto retry;
}

static inline int unlock_bwa_shm_info() {
	int ret;	
	if (!bwa_shm_info) return 0;
	__sync_synchronize();
	ret = __sync_val_compare_and_swap(&bwa_shm_info->lock, 1, 0);
	__sync_synchronize();
	if (ret == 1)
		return 1;
	else
		return 0;
}

void show_bwa_shm_info(bwa_shm_info_t *info, const char *name) {
	fprintf(stderr, "[BWA_SHM_INFO]%s%s\n", name ? " name: " : "", name ? name : "");
	fprintf(stderr, "[BWA_SHM_INFO] state: %d num_read: %d num_manager: %d hugetlb_flags: %x useErt: %d\n",
					info->state, info->num_map_read, info->num_map_manager,
					info->hugetlb_flags, info->useErt);
#ifdef MEMSCALE
	fprintf(stderr, "[BWA_SHM_INFO] [memscale] bwt: %d pac: %d ref: %d kmer: %d mlt: %d perfect: %d smem_all: %d smem_last: %d\n",
					info->bwt_on, info->pac_on, info->ref_on, 
					info->kmer_on, info->mlt_on,
					info->perfect_on,
					info->smem_all_on, info->smem_last_on);
#endif
#ifdef PERFECT_MATCH
	fprintf(stderr, "[BWA_SHM_INFO] perfect_mmap: %d perfect_seed_len: %d perfect_num_loc: %u perfect_num_seed: %u\n",
					info->pt_mmap,
					info->pt_seed_len,
					info->pt_num_loc_entry,
					info->pt_num_seed_entry);
#endif
#ifdef MEMSCALE
	fprintf(stderr, "[BWA_SHM_INFO] [memscale] perfect_num_seed_load: %u\n",
					info->pt_num_seed_entry_loaded);
#endif
	fprintf(stderr, "[BWA_SHM_INFO] reference_len: %ld ref_file_name(%d): %s\n",
					info->reference_len, info->ref_file_name_len, info->ref_file_name);
}


/*********************************************************/
/* hugeTLB helper functions                              */
/*********************************************************/
static enum hugetlb_mode parse_hugetlb_mode(const char *arg) {
	if (strcmp(arg, "huge") == 0)
		return BWA_SHM_HUGE_PAGE;
	else if (strcmp(arg, "2mb") == 0)
		return BWA_SHM_HUGE_2MB;
	else if (strcmp(arg, "1gb") == 0)
		return BWA_SHM_HUGE_1GB;
	else if (strcmp(arg, "normal") == 0)
		return BWA_SHM_NORMAL_PAGE;
	else {
		fprintf(stderr, "ERROR: hugetlb_option %s is not supported. the normal page will be used.\n", arg);
		return BWA_SHM_NORMAL_PAGE;
	}
}

static inline int get_hugetlb_flag(enum hugetlb_mode mode) {
	switch(mode) {
	case BWA_SHM_NORMAL_PAGE: return 0;
	case BWA_SHM_HUGE_PAGE:   return MAP_HUGETLB;
	case BWA_SHM_HUGE_2MB:    return MAP_HUGETLB | MAP_HUGE_2MB;
	case BWA_SHM_HUGE_1GB:    return MAP_HUGETLB | MAP_HUGE_1GB;
	default:                  return 0;
	}
}

static inline size_t get_hugetlb_unit(enum hugetlb_mode mode) {
	switch(mode) {
	case BWA_SHM_NORMAL_PAGE: return (4L * 1024);
	case BWA_SHM_HUGE_PAGE:   return DEFAULT_HUGETLB_PAGESIZE;
	case BWA_SHM_HUGE_2MB:    return (2L << 20);
	case BWA_SHM_HUGE_1GB:    return (1L << 30);
	default:                  return 0;
	}
}

static inline size_t page_aligned_size(size_t size) {
	return __aligned_size(size, 4L << 10);
}

static inline size_t hugetlb_aligned_size(size_t size) {
	int flags = bwa_shm_hugetlb_flags();
	if (flags == 0)
		return __aligned_size(size, 4L << 10);
	else if (flags == MAP_HUGETLB)
		return __aligned_size(size, DEFAULT_HUGETLB_PAGESIZE);
	else
		return __aligned_size(size, 1L << (flags >> MAP_HUGE_SHIFT));
}

#ifdef MEMSCALE
static size_t get_bwa_shm_size(int m) {
	bwa_shm_info_t *info = loading_info ? loading_info : bwa_shm_info;
	size_t size = 0;
	
	if (!info) {
		if (m != BWA_SHM_INFO)
			return 0;
		else if (BWA_SHM_INFO_FD >= 0)
			return __get_shm_size(BWA_SHM_INFO_FD);
		else
			return page_aligned_size(sizeof(bwa_shm_info_t));
	}

	switch(m) {
	case BWA_SHM_INFO:	size = bwa_shm_size_info(info->ref_file_name_len);
						return page_aligned_size(size);
						break;

	case BWA_SHM_BWT:	if (info->bwt_on)
							size = bwa_shm_size_bwt(info->reference_len);
						break;

	case BWA_SHM_PAC:	if (info->pac_on)
							size = bwa_shm_size_pac(info->reference_len);
						break;

	case BWA_SHM_REF:	if (info->ref_on)
							size = bwa_shm_size_ref(info->reference_len);
						break;

	case BWA_SHM_KMER:  if (info->kmer_on)
							size = bwa_shm_size_kmer();
						break;

	case BWA_SHM_MLT:   if (info->mlt_on)
							size = bwa_shm_size_mlt(info);
						break;

#ifdef PERFECT_MATCH
	case BWA_SHM_PERFECT: if (info->perfect_on)
							size = ____lpt_shm_size(info->pt_num_loc_entry, 
												  info->pt_num_seed_entry_loaded);
						break;
#endif
#ifdef SMEM_ACCEL
	case BWA_SHM_SALL:	if (info->smem_all_on)
							size = bwa_shm_size_sall();
						break;
	case BWA_SHM_SLAST: if (info->smem_last_on)
							size = bwa_shm_size_slast();
						break;
#endif
	default:
					return 0;
	}
	
	if (size > 0)
		return hugetlb_aligned_size(size);
	else
		return 0;
}
#else
static size_t get_bwa_shm_size(int m) {
	bwa_shm_info_t *info = loading_info ? loading_info : bwa_shm_info;
	size_t size = 0;

	if (!info) {
		if (m != BWA_SHM_INFO)
			return 0;
		else if (BWA_SHM_INFO_FD >= 0)
			return __get_shm_size(BWA_SHM_INFO_FD);
		else
			return page_aligned_size(sizeof(bwa_shm_info_t));
	}

	switch(m) {
	case BWA_SHM_INFO: size = bwa_shm_size_info(info->ref_file_name_len);
						return page_aligned_size(size);
						break;
	case BWA_SHM_BWT: size = bwa_shm_size_bwt(info->reference_len);
						break;
	case BWA_SHM_PAC: size = bwa_shm_size_pac(info->reference_len);
						break;
	case BWA_SHM_REF: size = bwa_shm_size_ref(info->reference_len);
						break;
	case BWA_SHM_KMER: size = bwa_shm_size_kmer();
						break;
	case BWA_SHM_MLT: size = bwa_shm_size_mlt(info);
						break;
#ifdef PERFECT_MATCH
	case BWA_SHM_PERFECT: size = ____lpt_shm_size(info->pt_num_loc_entry, 
												  info->pt_num_seed_entry);
						break;
#endif
#ifdef SMEM_ACCEL
	case BWA_SHM_SALL: size = bwa_shm_size_sall();
						break;
	case BWA_SHM_SLAST: size = bwa_shm_size_slast();
						break;
#endif
	default:
					return 0;
	}

	return hugetlb_aligned_size(size);
}
#endif

static inline const char *bwa_shm_filename(int m) {
	return use_hugetlb(m) ? bwa_shm_huge_name_str[m]
						  : bwa_shm_name_str[m];
}

int bwa_shm_create(int m, size_t size) {
	int fd = shm_fd[m];
	if (fd >= 0) {
		if (close(shm_fd[m]))
			return -1;
		else
			shm_fd[m] = -1;
	}

	if (use_mmap(m)) {
		fd = bwa_shm_open(m);
		if (fd < 0) {
			char fn[PATH_MAX];
			fprintf(stderr, "[bwa_shm] %s: failed to open %s. errno: %d\n",
								__func__, bwa_shm_mmap_filename(m, fn), errno);
		}
		return fd;
	} else if (use_hugetlb(m)) {
		fd = open(bwa_shm_filename(m), 
						bwa_shm_create_flags, 
						bwa_shm_create_mode);
	} else {
		fd = shm_open(bwa_shm_filename(m), 
						bwa_shm_create_flags, 
						bwa_shm_create_mode);
	}

	if (fd < 0) {
		fprintf(stderr, "[bwa_shm] %s: failed to open %s. errno: %d\n", 
					__func__, bwa_shm_filename(m), errno);
		return fd;
	}

	if ((!use_mmap(m)) && (!use_hugetlb(m))) {
		if (ftruncate(fd, page_aligned_size(size))) {
			fprintf(stderr, "[bwa_shm] %s: failed to truncate %s to 0x%lx. errno: %d\n", 
					__func__, bwa_shm_filename(m), page_aligned_size(size), errno);
			return -1;
		}
	}
	
	shm_fd[m] = fd;
	return fd;
}

int bwa_shm_open(int m) {
	if (shm_fd[m] >= 0)
		return shm_fd[m];
	if (use_mmap(m)) {
		char buf[PATH_MAX];
		char *fn = bwa_shm_mmap_filename(m, buf);
		/* when fn == NULL, open() will return -1 and set errno. */
		shm_fd[m] = open(fn, O_RDONLY | O_DIRECT | O_SYNC);
	} else if (use_hugetlb(m))
		shm_fd[m] = open(bwa_shm_filename(m), O_RDWR, 0666);
	else
		shm_fd[m] = shm_open(bwa_shm_filename(m), O_RDWR, 0666);

	return shm_fd[m];
}

void *bwa_shm_map(int m) {
	int fd;
	void *ptr;
	size_t size;
	if (m < 0 || m >= NUM_BWA_SHM)
		return NULL;

	fd = bwa_shm_open(m);
	if (fd < 0)
		return NULL;

	ptr = shm_ptr[m];
	if (ptr)
		return ptr;
	
	size = get_bwa_shm_size(m);
	if (use_mmap(m) == 0)
		ptr = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED | bwa_shm_hugetlb_flags(),
					fd, 0);
	else
		ptr = mmap(NULL, size, PROT_READ, MAP_SHARED | bwa_shm_hugetlb_flags(),
					fd, 0);
		
	if (ptr == MAP_FAILED) {
		fprintf(stderr, "%s: mmap failed. type: %d fd: %d size: %ld hugetlb: %x errno: %d\n",
					__func__, m, fd, size, bwa_shm_hugetlb_flags(), errno);
		return NULL;
	}
	if (opt_bwa_shm_map_touch == 1) {
		uint8_t x = 0, *p = (uint8_t *)ptr;
		size_t i;
		for (i = 0; i < size; ++i)
			x = x ^ p[i];
		fprintf(stderr, "INFO: shm_map_touch type: %d size: 0x%lx result: 0x%x\n", m, size, x);
	}
	fprintf(stderr, "INFO: shm_map. type: %d fd: %d size: 0x%lx\n", m, fd, size);
	shm_ptr[m] = ptr;
	if (m == BWA_SHM_INFO)
		bwa_shm_info = (bwa_shm_info_t *) ptr;
	
	return ptr;
}

int bwa_shm_close(int m) {
	int ret = close(shm_fd[m]);
	if (ret == 0)
		shm_fd[m] = -1;
	return ret;
}

int bwa_shm_unmap(int m) {
	size_t size = get_bwa_shm_size(m);
	if (!shm_ptr[m]) return -1;
	munmap(shm_ptr[m], size);
	fprintf(stderr, "INFO: shm_unmap type: %d size: 0x%lx\n", m, size);
	shm_ptr[m] = NULL;
	if (m == BWA_SHM_INFO)
		bwa_shm_info = NULL;
	return bwa_shm_close(m);
}

static inline void __bwa_shm_unmap_all() {
	int m;
	for (m = NUM_BWA_SHM - 1; m >= 0; m--)
		bwa_shm_unmap(m);
}

/* return a negative number if error occurs. */
int __bwa_shm_remove(int m) {
	int fd;
	bwa_shm_unmap(m);
	if (use_mmap(m)) {
		fprintf(stderr, "remove_shm: type: %d DO_NOT_REMOVE_MMAPED_REGION\n", m);
		return 0;
	}

	fd = bwa_shm_open(m);
	if (fd < 0) {
		fprintf(stderr, "remove_shm: type: %d name: %s NOT_EXIST\n", m, bwa_shm_filename(m));
		return 0;
	}

	if (bwa_shm_close(m)) {
		fprintf(stderr, "remove_shm: type: %d name: %s FAILED (close) errno: %d\n", m, bwa_shm_filename(m), errno);
		return -1;
	}

	if (use_hugetlb(m) ? unlink(bwa_shm_filename(m))
					   : shm_unlink(bwa_shm_filename(m))) {
		fprintf(stderr, "remove_shm: type: %d name: %s FAILED (unlink) errno: %d\n", m, bwa_shm_filename(m), errno);
		return -1;
	}
	
	fprintf(stderr, "remove_shm: type: %d name: %s SUCCEED\n", m, bwa_shm_filename(m));
	return 0;
}

static int __check_mount_hugetlbfs(const char *path, enum hugetlb_mode mode, int print_error);

int __bwa_shm_remove_all() {
	int ret = 0;
	int m;
	bwa_shm_info_t info;
	bwa_shm_info_t *bak;

	memset(&info, 0, sizeof(bwa_shm_info_t));
	info.hugetlb_flags = MAP_HUGETLB;
	
	bak = loading_info;
	loading_info = &info;
	for (m = NUM_BWA_SHM - 1; m > BWA_SHM_INFO; m--) {
		if (__bwa_shm_remove(m) < 0)
			ret = -1;
	}

	info.hugetlb_flags = 0;
	for (m = NUM_BWA_SHM - 1; m >= BWA_SHM_INFO; m--) {
		if (__bwa_shm_remove(m) < 0)
			ret = -1;
	}
	loading_info = bak;

	return ret;
}

static int __bwa_shm_remove_hugetlb() {
	int ret = 0, r;
	r = __check_mount_hugetlbfs(BWA_SHM_HUGE_DIR, BWA_SHM_HUGE_PAGE, 0);

	if (r == 0 || r == EINVAL) {
		ret = umount(BWA_SHM_HUGE_DIR);
		ret = ret == 0 ? 0 : (errno == EINVAL ? 0 : ret); /* if BWA_SHM_HUGE_DIR is not a mount point, try to remove it */
		fprintf(stderr, "remove_shm: unmount %s: %s\n", BWA_SHM_HUGE_DIR, ret == 0 ? "SUCCEED" : "FAILED");
		if (ret == 0) {
			ret = rmdir(BWA_SHM_HUGE_DIR);
			if (ret < 0)
				fprintf(stderr, "remove_shm remove %s: FAILED (erro: %d)\n", BWA_SHM_HUGE_DIR, errno);
			else
				fprintf(stderr, "remove_shm: remove %s: SUCCEED\n", BWA_SHM_HUGE_DIR);
		} 
	} else {
		fprintf(stderr, "remove_shm: unmount %s: %s\n", BWA_SHM_HUGE_DIR, r == ENOENT ? "NOT_EXIST" : "FAILED");
		ret -1;
	}

	return ret;
}

static void __bwa_shm_init_data(const char *prefix) {
	int m;

	/* initialize data structure */
	bwa_shm_info = NULL;
	for (m = 0; m < NUM_BWA_SHM; ++m) {
		shm_fd[m] = -1;
		shm_ptr[m] = NULL;
	}
	if (prefix)
		set_mmap_prefix(prefix);
}

void bwa_shm_init(const char *prefix, int *useErt, int pt_seed_len,
							enum bwa_shm_init_mode mode) 
{
	struct stat ref_stat;		
	struct timespec mtim_ref;
	char ref_file_name[PATH_MAX];
	char *abs_path;
	size_t abs_path_len, size;
	bwa_shm_info_t *info;
	int64_t rlen;
	int m;
	int fd;
	int state = BWA_SHM_STATE_NOT_INIT, num_second;

	/* initialize data structure */
	__bwa_shm_init_data(ref_file_name);

    strcpy_s(ref_file_name, PATH_MAX, prefix);
    strcat_s(ref_file_name, PATH_MAX, ".0123");
	
	fprintf(stderr, "ref_file: %s\n", ref_file_name);

	if (stat(ref_file_name, &ref_stat) != 0) {
		fprintf(stderr, "ERROR! Unable to stat the file: %s\n", ref_file_name);
		exit(EXIT_FAILURE);
	}

	mtim_ref = ref_stat.st_mtim;
	rlen = ref_stat.st_size + 1;

	abs_path = realpath(ref_file_name, NULL);
	abs_path_len = strlen(abs_path);

	fd = bwa_shm_open(BWA_SHM_INFO);
	
	if (fd < 0) {
		if (mode != BWA_SHM_INIT_NEW) 
			fprintf(stderr, "[bwa_shm] the previous info does not exist\n");
		goto renewal;
	} else if (mode == BWA_SHM_INIT_NEW) {
		fprintf(stderr, "[bwa_shm] the previous info still exist.\n");
		goto renewal;
	}
		
	info = (bwa_shm_info_t *) bwa_shm_map(BWA_SHM_INFO);
	if (info == NULL) {
		fprintf(stderr, "[bwa_shm] failed to map INFO\n");
		goto renewal;
	}

	num_second = 0;
	state = BWA_SHM_STATE_NOT_INIT;
	while (state != BWA_SHM_STATE_AVAIL) {
		if (lock_bwa_shm_info() == 1) {
			state = bwa_shm_info->state;
			if (state != BWA_SHM_STATE_AVAIL) {
				unlock_bwa_shm_info();
				fprintf(stderr, "[bwa_shm] an index manager is modifying the indexes...(%d)\n", num_second);
				sleep(1);
				num_second++;
			} else 
				break; /* DON'T unlock */
		} else {
			fprintf(stderr, "[bwa_shm] an index manager is initializing the indexes...(%d)\n", num_second);
			sleep(1);
			num_second++;
		}
	}

	/* locked and state == AVAIL */
	if (mode == BWA_SHM_INIT_READ) {
		bwa_shm_info->num_map_read++;
	} else {
		bwa_shm_info->num_map_manager++;
		if (bwa_shm_info->num_map_read > 0)
			bwa_shm_info->state = BWA_SHM_STATE_WAIT;
		else
			bwa_shm_info->state = BWA_SHM_STATE_MODIFY;
	}
	state = bwa_shm_info->state;
	unlock_bwa_shm_info();

	if (mode == BWA_SHM_INIT_NEW)
		goto renewal;

	if (info->mtim_ref.tv_sec != mtim_ref.tv_sec
			|| info->mtim_ref.tv_nsec != mtim_ref.tv_nsec) {
		fprintf(stderr, "[bwa_shm] the last modified time of reference file is changed.\n");
		goto renewal;
	}

	if (bwa_shm_info->ref_file_name_len != abs_path_len) {
		fprintf(stderr, "[bwa_shm] you use a different reference file from before.\n");
		goto renewal;
	}

	if (strcmp(bwa_shm_info->ref_file_name, abs_path) != 0) {
		fprintf(stderr, "[bwa_shm] you use a different reference file from before.\n");
		goto renewal;
	}

	if (*useErt >= 0 && bwa_shm_info->useErt != *useErt) {
		fprintf(stderr, "[bwa_shm] you previously %suse ERT, but now you %suse ERT.\n",
				bwa_shm_info->useErt == 1 ? "" : "don't ",
				*useErt == 1 ? "" : "don't ");
		goto renewal;
	} else {
		*useErt = bwa_shm_info->useErt;
	}
	
	/* get other shm_fds */
	for (m = BWA_SHM_INFO + 1; m < NUM_BWA_SHM; ++m) {
		if (use_mmap(m)) /* for mmap'ed region, don't check shm_fd here */
			continue;

		if (bwa_shm_info->useErt) {
			if (m == BWA_SHM_BWT)
				continue;
#ifdef SMEM_ACCEL
			if (m == BWA_SHM_SALL || m == BWA_SHM_SLAST)
				continue;
#endif
		} else {
			if (m == BWA_SHM_KMER || m == BWA_SHM_MLT)
				continue;
		}

		fd = bwa_shm_open(m);
		if (fd < 0) {
#ifdef MEMSCALE
			if (m == BWA_SHM_PERFECT || m == BWA_SHM_SALL || m == BWA_SHM_SLAST)
				continue;
#endif
			fprintf(stderr, "[bwa_shm] failed to get shared memory of %s\n",
						bwa_shm_type_str[m]);
			goto renewal;
		}
	}
	
	for (m = BWA_SHM_INFO + 1; m < NUM_BWA_SHM; ++m) {
		if (use_mmap(m)) /* for mmap'ed region, don't check mapping here */
			continue;

		if (bwa_shm_info->useErt) {
			if (m == BWA_SHM_BWT)
				continue;
#ifdef SMEM_ACCEL
			if (m == BWA_SHM_SALL || m == BWA_SHM_SLAST)
				continue;
#endif
		} else {
			if (m == BWA_SHM_KMER || m == BWA_SHM_MLT)
				continue;
		}
		
		if (bwa_shm_map(m) == NULL) 
		{
#ifdef MEMSCALE
			if (m == BWA_SHM_PERFECT || m == BWA_SHM_SALL || m == BWA_SHM_SLAST)
				continue;
#endif
			fprintf(stderr, "[bwa_shm] failed to map shared memory of %s\n",
						bwa_shm_type_str[m]);
			goto disable;
		}
	}

	bwa_shm_mode = BWA_SHM_MATCHED;
	fprintf(stderr, "BWA_SHM_MODE: MATCHED\n");
	free(abs_path);
	return;

renewal:
	/* unmap current bwa_shm_info */
	if (bwa_shm_info) {
		/* if bwa_shm_info exists,
		   wait for other mapping processes, and re-init the info */

		/* wait for other processes */
		if (state == BWA_SHM_STATE_WAIT) {
			int num_map_read = 1;
			num_second = 0;
			
			while (num_map_read > 0) {
				if (lock_bwa_shm_info() == 0) {
					num_map_read = bwa_shm_info->num_map_read;
					unlock_bwa_shm_info();
					if (num_map_read > 0) {
						fprintf(stderr, "[bwa_shm] %d mappers are using the indexes...(%d)\n", 
														num_map_read, num_second);
						sleep(1);
						num_second++;
					} else 
						break;
				} else 
					break;
			}
		}

		if (ftruncate(BWA_SHM_INFO_FD, bwa_shm_size_info(abs_path_len))) {
			if (__get_shm_size(BWA_SHM_INFO_FD) < bwa_shm_size_info(abs_path_len)) {
				fprintf(stderr, "[bwa_shm] failed to increase the size of shared memory for information\n");
				goto disable;
			}
		}
		/* The following memset() unlockes the info,
		   and info->state is now "NOT_INIT", so other process will not use this info. */
		memset(bwa_shm_info, 0, bwa_shm_size_info(abs_path_len));
	}

	/* remove shms */
	if (__bwa_shm_remove_all()) {
		fprintf(stderr, "[bwa_shm] failed to remove the previous shared memories\n");
		goto disable; /* failed to remove some shms */
	}
	
	/* re-create bwa_shm_info and fill it except for refernece_len.
	   reference_len will be filled at load_BWT_and_FMI() */
	fd = bwa_shm_create(BWA_SHM_INFO, bwa_shm_size_info(abs_path_len));
	if (fd < 0) {
		fprintf(stderr, "[bwa_shm] failed to create BWA_SHM_INFO. errno: %d\n", errno);
		goto disable;
	}

	info = (bwa_shm_info_t *) bwa_shm_map(BWA_SHM_INFO);
	if (!info) {
		fprintf(stderr, "[bwa_shm] failed to map BWA_SHM_INFO\n");
		goto disable;
	}
	
	info->num_map_read = 0;
	info->num_map_manager = 1; /* me */
	info->hugetlb_flags = 0;
	if ((*useErt) < 0) {
		info->useErt = DEFAULT_USE_ERT;
		*useErt = DEFAULT_USE_ERT;
	} else {
		info->useErt = *useErt;
	}
#ifdef MEMSCALE
	info->bwt_on = 0;
	info->pac_on = 0;
	info->ref_on = 0;
	info->kmer_on = 0;
	info->mlt_on = 0;
	info->perfect_on = 0;
	info->smem_all_on = 0;
	info->smem_last_on = 0;
	info->pt_num_seed_entry_loaded = 0;
#endif
#ifdef PERFECT_MATCH
	info->pt_num_loc_entry = 0;
	info->pt_num_seed_entry = 0;
	if (pt_seed_len > 0
			&& pt_seed_len != PT_SEED_LEN_NO_TABLE
			&& pt_seed_len != PT_SEED_LEN_AUTO_TABLE)
		info->pt_seed_len = pt_seed_len;
	else
		info->pt_seed_len = 0;
	info->pt_mmap = DEFAULT_MMAP_PERFECT;
#endif

	info->reference_len = rlen;
	info->mtim_ref = mtim_ref;
	info->ref_file_name_len = abs_path_len;
	strncpy(bwa_shm_info->ref_file_name, abs_path, abs_path_len);

	bwa_shm_info->lock = 0;
	
	__sync_synchronize();
	bwa_shm_info->state = BWA_SHM_STATE_MODIFY;
	
	bwa_shm_mode = BWA_SHM_RENEWAL;
	free(abs_path);
	fprintf(stderr, "BWA_SHM_MODE: RENEWAL\n");

#ifdef MEMSCALE
	if (mode == BWA_SHM_INIT_READ) {
		__bwa_shm_load(prefix, BWA_SHM_NORMAL_PAGE, 0,
				info->pt_seed_len, info->pt_mmap, 0);
		bwa_shm_mode = BWA_SHM_MATCHED;
	}
#endif
	return;

disable:
	__bwa_shm_unmap_all();
	bwa_shm_mode = BWA_SHM_DISABLE;
	free(abs_path);
	fprintf(stderr, "BWA_SHM_MODE: DISABLE\n");
	return;
}

/* called after index initialization of 'mem' command */
void bwa_shm_complete(enum bwa_shm_init_mode mode) {
	if (mode == BWA_SHM_INIT_READ && bwa_shm_mode == BWA_SHM_RENEWAL) {
		int locked = lock_bwa_shm_info();
		if (locked) {
			bwa_shm_info->num_map_manager--;
			bwa_shm_info->num_map_read++;
			bwa_shm_info->state = BWA_SHM_STATE_AVAIL;
			bwa_shm_mode = BWA_SHM_MATCHED;
			show_bwa_shm_info(bwa_shm_info, "complete");
			unlock_bwa_shm_info();
		}
	}
}

void bwa_shm_final(enum bwa_shm_init_mode mode) {
	int locked;
	
	if (bwa_shm_mode == BWA_SHM_DISABLE)
		return;

	locked = lock_bwa_shm_info();
	if (locked) {
		if (mode == BWA_SHM_INIT_READ) {
			/* when bwa_shm_mode == BWA_SHM_RENEWAL, 
			   bwa_shm_complete() decrement num_map_manager counter */
			bwa_shm_info->num_map_read--;
		} else { /* BWA_SHM_INIT_NEW or BWA_SHM_INIT_MODIFY */
			bwa_shm_info->num_map_manager--;
			bwa_shm_info->state = BWA_SHM_STATE_AVAIL;
			show_bwa_shm_info(bwa_shm_info, "final");
		}
		unlock_bwa_shm_info();
	}
	__bwa_shm_unmap_all();
}

static inline size_t __get_hugetlbfs_pagesize(const char *opts) {
	char *buf = strdup(opts);
	char *ptr;
	size_t ret = 0;


	if (buf == NULL)
		return SIZE_MAX;
	
	ptr = strstr(buf, "pagesize=");
	if (ptr == NULL) {
		ret = DEFAULT_HUGETLB_PAGESIZE;
		goto out;
	}

	ptr = ptr + 9; /* skip "pagesize=" */

	while ((*ptr != ',') && (*ptr != '\0')) {
		if (*ptr >= '0' && *ptr <= '9') {
			ret = ret * 10 + (*ptr - '0');
		} else if (*ptr == 'k' || *ptr == 'K') {
			ret = ret * (1 << 10);
			break;
		} else if (*ptr == 'm' || *ptr == 'M') {
			ret = ret * (1 << 20);
			break;
		} else if (*ptr == 'g' || *ptr == 'G') {
			ret = ret * (1 << 30);
			break;
		} else {
			ret = SIZE_MAX;
			break;
		}
		ptr++;
	}

out:
	free(buf);
	return ret;
}

static inline const char *__hugetlbfs_opts(enum hugetlb_mode mode) {
	switch(mode) {
	case BWA_SHM_HUGE_PAGE:   
	case BWA_SHM_HUGE_2MB:    return "pagesize=2M";
	case BWA_SHM_HUGE_1GB:    return "pagesize=1024M";
	case BWA_SHM_NORMAL_PAGE: 
	default:                  return "NEVER_A_REAL_OPTION_STRING";
	}
}
/* path is a mount point. check the mounted fs is hugetlbfs matching with @mode 
 * return <0 if error occurs
 * return 0 if path is a hugetlbfs mount point and the option is matched.
 * return EINVAL if path is a hugetlbfs mount point but the option is not matched.
 * return EAGAIN if path is a mount point but its filesystem is not hugetlbfs.
 * return EISDIR if path is a directory but it is not a mount point.
 * return ENOENT if path is invalid.
 */
static int __check_mount_hugetlbfs(const char *path, enum hugetlb_mode mode, int print_error) {
	struct mntent *ent = NULL, *ent_buf = NULL;
	char *str_buf = NULL;
	static size_t str_buf_size = 2048 * sizeof(char);
	FILE *f = NULL;
	DIR *dir = NULL;
	int ret = 0;

	dir = opendir(path);
	if (dir == NULL) {
		if (errno == ENOENT)
			return ENOENT;
		else if (errno == ENOTDIR) {
			if (print_error)
				fprintf(stderr, "ERROR: %s already exists and not a directory.\n", path);
			return -ENOTDIR;
		} else { 
			ret = -errno;
			if (print_error)
				fprintf(stderr, "ERROR: failed to open %s. errno: %d\n", path, errno);
			return ret;
		}
	} else
		closedir(dir);
	
	f = setmntent("/etc/mtab", "r");
	if (!f) {
		fprintf(stderr, "ERROR: failed to open /etc/mtab to check the filesystem of %s\n", path);
		ret = -ENXIO;
		goto out;
	}
	
	ent_buf = (struct mntent *) malloc(sizeof(struct mntent));
	str_buf = (char *) malloc(str_buf_size);

	if (ent_buf == NULL || str_buf == NULL) {
		fprintf(stderr, "ERROR: failed to allocate memory to check the filesystem of %s\n", path);
		ret = -ENOMEM;
		goto out;
	}

	while ((ent = getmntent_r(f, ent_buf, str_buf, str_buf_size)) != NULL) {
		if (strcmp(ent->mnt_dir, path) == 0)
			break;
	}

	if (ent == NULL) {
		fprintf(stderr, "ERROR: %s is not a mount point\n", path);
		ret = EISDIR;
		goto out;
	}

	fprintf(stderr, "INFO: [mount] path: %s type: %s fs: %s opts: %s\n", 
					path, ent->mnt_type, ent->mnt_fsname, ent->mnt_opts);

	if (strcmp(ent->mnt_type, "hugetlbfs") == 0
			&& strcmp(ent->mnt_fsname, "nodev") == 0) {
		if (__get_hugetlbfs_pagesize(ent->mnt_opts)
					== get_hugetlb_unit(mode))
			ret = 0;
		else
			ret = EINVAL;
	} else
		ret = EAGAIN;

out:
	if (ent_buf) free(ent_buf);
	if (str_buf) free(str_buf);
	if (f) endmntent(f);

	return ret;
}

static inline int __mount_hugetlbfs(const char *path, enum hugetlb_mode mode) {
	int ret; 

	ret = __check_mount_hugetlbfs(path, mode, 1);

	if (ret <= 0)
		return ret; /* already mounted as we want, or error */
	
	if (ret == ENOENT) {
		if (mkdir(BWA_SHM_HUGE_DIR, 0777)) {
			ret = errno;
			fprintf(stderr, "ERROR: failed to mkdir %s. errno: %d\n",
						path, ret);
			return -ret;
		} else
			fprintf(stderr, "INFO: %s is created\n", path);
		
	}

	if (ret == EINVAL || ret == EAGAIN) { /* umount first */
		if (umount(path)) {
			ret = -errno;
			fprintf(stderr, "ERROR: failed umount %s to re-mount. errno: %d\n",
								path, errno);
			return ret;
		}
	}

	ret = mount("nodev", path, "hugetlbfs", 0, __hugetlbfs_opts(mode));
	if (ret < 0)
		fprintf(stderr, "ERROR: failed to mount %s. errno: %d\n",
					path, errno);
	else
		fprintf(stderr, "INFO: mount hugetlbfs at %s with %s\n",
					path, __hugetlbfs_opts(mode));
	return ret;
}

static int __check_hugetlb(const char *unit) {
	char fn[PATH_MAX];
	FILE *f;

    strcpy_s(fn, PATH_MAX, "/sys/kernel/mm/hugepages/hugepages-");
    strcat_s(fn, PATH_MAX, unit);
    strcat_s(fn, PATH_MAX, "/nr_hugepages");
	f = fopen(fn, "r+");
	if (f == NULL) {
		fprintf(stderr, "ERROR: failed to open a file for %s hugepage. errno: %d\n",
							unit, errno);
		return -1;
	}
	fclose(f);
	return 0;
}

static int check_hugetlb(enum hugetlb_mode mode) {
	int ret;

	if (mode == BWA_SHM_NORMAL_PAGE)
		return 0;

	if (getuid() != 0) {
		fprintf(stderr, "ERROR: hugeTLB requires the root permission.\n");
		return -EPERM;
	}

	ret = __mount_hugetlbfs(BWA_SHM_HUGE_DIR, mode);
	if (ret)
		return ret;
		
	if (mode == BWA_SHM_HUGE_PAGE) {
		ret = __check_hugetlb(KB_UNIT_STR_DEFAULT);
		if (ret) 
			fprintf(stderr, "ERROR: hugetlb is not available\n");
		return ret;
	} else if (mode == BWA_SHM_HUGE_2MB) {
		ret = __check_hugetlb(KB_UNIT_STR_2MB);
		if (ret) 
			fprintf(stderr, "ERROR: huge_2mb is not available\n");
		return ret;
	} else { // if (mode == BWA_SHM_HUGE_1GB)
		ret = __check_hugetlb(KB_UNIT_STR_1GB);
		if (ret)
			fprintf(stderr, "ERROR: huge_1gb is not available.\n");
		return ret;
	}
	
	return 0;
}

static int __get_hugepages(const char *unit, size_t num_pages) {
	char fn[PATH_MAX];
	FILE *f;
	int ret;
	size_t num_curr, num_free, num_final;

    strcpy_s(fn, PATH_MAX, "/sys/kernel/mm/hugepages/hugepages-");
    strcat_s(fn, PATH_MAX, unit);
    strcat_s(fn, PATH_MAX, "/free_hugepages");
	f = fopen(fn, "r");
	if (f == NULL) {
		fprintf(stderr, "ERROR: failed to open %s. errno: %d\n", fn, errno);
		return -1;
	}
	ret = fscanf(f, "%ld", &num_free);
	fclose(f);

	if (ret != 1) {
		fprintf(stderr, "ERROR: failed to read from %s. errno: %d\n", fn, errno);
		return -1;
	}
	
    strcpy_s(fn, PATH_MAX, "/sys/kernel/mm/hugepages/hugepages-");
    strcat_s(fn, PATH_MAX, unit);
    strcat_s(fn, PATH_MAX, "/nr_hugepages");
	f = fopen(fn, "r");
	if (f == NULL) {
		fprintf(stderr, "ERROR: failed to open %s. errno: %d\n", fn, errno);
		return -1;
	}
	ret = fscanf(f, "%ld", &num_curr);
	fclose(f);

	if (ret != 1) {
		fprintf(stderr, "ERROR: failed to read from %s. errno: %d\n", fn, errno);
		return -1;
	}

	num_final = num_curr + (num_pages - num_free);

	if (num_final <= num_curr)
		return 0;

	f = fopen(fn, "w");
	if (f == NULL) {
		fprintf(stderr, "ERROR: failed to open %s for write. errno: %d\n", fn, errno);
		return -1;
	}
	fprintf(stderr, "INFO: nr_hugepages(%s) will increase from %ld to %ld\n", unit, num_curr, num_final);
	ret = fprintf(f, "%ld", num_final);
	fclose(f);

	if (ret < 0) {
		fprintf(stderr, "ERROR: failed to write on %s. errno: %d\n", fn, errno);
		return -1;
	}
	
	f = fopen(fn, "r");
	if (f == NULL) {
		fprintf(stderr, "ERROR: failed to open %s. errno: %d\n", fn, errno);
		return -1;
	}
	ret = fscanf(f, "%ld", &num_curr);
	fclose(f);

	if (ret != 1) {
		fprintf(stderr, "ERROR: failed to read from %s. errno: %d\n", fn, errno);
		return -1;
	}

	if (num_curr != num_final) {
		fprintf(stderr, "ERROR: failed to set %s.\n", fn);
		return -1;
	}

	return 0;
}

static int get_hugepages(enum hugetlb_mode *mode_in, size_t shm_size, int opt_force) {
	size_t num_pages;
	enum hugetlb_mode mode = *mode_in;
	int ret = 0;
	
	switch(mode) {
	case BWA_SHM_HUGE_1GB:
							num_pages = (shm_size + (1 << 30) - 1) / (1 << 30);
							fprintf(stderr, "INFO: required 1GB pages: %ld\n", num_pages);
							if (__get_hugepages(KB_UNIT_STR_1GB, num_pages) == 0) {
								fprintf(stderr, "INFO: succeed to get %ld pages\n", num_pages);
								break;
							} else {
								fprintf(stderr, "INFO: failed to get %ld 1GB pages. Use 2MB pages.\n", num_pages);
								ret = -1;
								mode = BWA_SHM_HUGE_2MB;
							}

							if (opt_force)
								break;
							/* fall-through */
						
	case BWA_SHM_HUGE_2MB:
							num_pages = (shm_size + (1 << 21) - 1) / (1 << 21);
							fprintf(stderr, "INFO: required 2MB pages: %ld\n", num_pages);
							if (__get_hugepages(KB_UNIT_STR_2MB, num_pages) == 0) {
								fprintf(stderr, "INFO: succeed to get %ld pages\n", num_pages);
								break;
							} else {
								fprintf(stderr, "INFO: failed to get %ld 2MB pages. Use normal pages.\n", num_pages);
								ret = -1;
								mode = BWA_SHM_NORMAL_PAGE;
							}

							if (opt_force)
								break;
							/* fall-through */
	
	case BWA_SHM_HUGE_PAGE:
							num_pages = (shm_size + DEFAULT_HUGETLB_PAGESIZE - 1) / DEFAULT_HUGETLB_PAGESIZE;
							fprintf(stderr, "INFO: required huge pages: %ld\n", num_pages);
							if (__get_hugepages(KB_UNIT_STR_DEFAULT, num_pages) == 0) {
								fprintf(stderr, "INFO: succeed to get %ld pages\n", num_pages);
								break;
							} else {
								fprintf(stderr, "INFO: failed to get %ld huge pages. Use normal pages.\n", num_pages);
								ret = -1;
								mode = BWA_SHM_NORMAL_PAGE;
							}

							if (opt_force)
								break;
							/* fall-through */

	case BWA_SHM_NORMAL_PAGE:
							break;
	default:
							ret = -1;
							break;
	}

	*mode_in = mode;

	return ret;
}

#ifdef PERFECT_MATCH
static inline void get_perfect_table_filename(char *buf, const char *prefix, int seed_len) {
	snprintf(buf, PATH_MAX, "%s.perfect.%d", prefix, seed_len);
	return;
}

static size_t get_perfect_table_size(const char *prefix, int seed_len, 
				size_t *_size_head, size_t *_size_loc, size_t *_size_seed,
				uint32_t *_num_loc, uint32_t *_num_seed) {
	char file_name[PATH_MAX];
	FILE *fp;
	uint32_t dummy;
	size_t size_head, size_loc, size_seed;
	perfect_table_t pt;

	get_perfect_table_filename(file_name, prefix, seed_len);
	
	fp = fopen(file_name, "rb");
	if (fp == NULL) {
		fprintf(stderr, "ERROR: failed to open %s\n", file_name);
		return 0;
	}

	__lpt_load_head(&pt, fp);

	fclose(fp);

	size_head = __aligned_size(sizeof(perfect_table_t), 64);
	size_loc = __aligned_size(sizeof(uint32_t) * pt.num_loc_entry, 64);
	size_seed = __aligned_size(sizeof(seed_entry_t) * pt.num_seed_entry, 64);
	
	if (_size_head) *_size_head = size_head;
	if (_size_loc) *_size_loc = size_loc;
	if (_size_seed) *_size_seed = size_seed;
	if (_num_loc) *_num_loc = pt.num_loc_entry;
	if (_num_seed) *_num_seed = pt.num_seed_entry;

	return size_head + size_loc + size_seed;
}
#endif

static void usage_bwa_shm_load(void)
{
    fprintf(stderr, "Usage: bwa-mem2 load-shm [options] <idxbase>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -Z 0 or 1                Use ERT(Enumarated Radix Tree) index [%d]\n"
					"                             Show better performance but use 60GB more memory.\n", DEFAULT_USE_ERT);
    fprintf(stderr, "    -f                       Force using hugetlb. Exit with failure if setting hugh TLB fails.\n"
					"                             Default: fallback to normal pages.\n");
    fprintf(stderr, "    -H normal,huge,2mb,1gb   huge TLB options [normal]\n");
#ifdef MEMSCALE
#define HG38_RLEN (3209286105LL * 2 + 1)
    fprintf(stderr, "    -m                       Modify the loaded index\n");
    fprintf(stderr, "    -g                       The number of gigabytes of memory for index. [0]\n"
					"                             0 for unlimited. Range for hg38: %lld ~ %lld\n",
					B2GB(bwa_shm_size_bwt(HG38_RLEN) 
						+ bwa_shm_size_ref(HG38_RLEN)
						+ bwa_shm_size_pac(HG38_RLEN)),
					B2GB(bwa_shm_size_kmer() 
						+ 55067950773 // mlt size (hard coded)
						+ bwa_shm_size_ref(HG38_RLEN)
						+ bwa_shm_size_pac(HG38_RLEN)
						+ __aligned_size(sizeof(perfect_table_t), 64) 
						+ (((HG38_RLEN - 1) / 2) * 11 / 10) * sizeof(seed_entry_t) 
						+ (((HG38_RLEN - 1) / 2) / 100) * sizeof(uint32_t) * 2)
					);
#endif
#ifdef PERFECT_MATCH
	fprintf(stderr, "    -l INT                   load perfect hash table with the specified seed length\n");
#endif
}

int __bwa_shm_load_file(const char *prefix, const char *postfix, int m, void **ret_ptr) {
	void *ptr;
	int fd;

	if (bwa_shm_mode == BWA_SHM_RENEWAL) {
		size_t shm_size = get_bwa_shm_size(m);
		fprintf(stderr, "INFO: shm_create for %s. size: %ld hugetlb_flag: 0x%x\n",
					bwa_shm_type_str[m], shm_size, bwa_shm_hugetlb_flags());
		fd = bwa_shm_create(m, shm_size);

		if (fd >= 0) {
			/* if bwa_shm_map() fails, it returns NULL,
				and __load_file() runs in the allocation mode */
			ptr = __load_file(prefix, postfix, bwa_shm_map(m), &shm_size);
			goto out;
		} else {
			fprintf(stderr, "[bwa_shm] failed to create BWA_SHM_%s\n", bwa_shm_type_str[m]);
		}
	} else if (bwa_shm_mode == BWA_SHM_MATCHED) {
		ptr = bwa_shm_map(m);
		if (ptr) goto out;
	}

	ptr = __load_file(prefix, postfix, NULL, NULL); /* ERROR or BWA_SHM_DISABLE */

out:
	if (ptr == NULL)
		return -1;
	if (ret_ptr)
		*ret_ptr = ptr;
	return 0;
}

static int __bwa_shm_load_BWT(const char *prefix, size_t size) {
	int ret = 0;
	
	fprintf(stderr, "INFO: load BWT index (size: %ld)\n", size);

	int __load_BWT_on_shm(const char *ref_file_name, int64_t *_reference_seq_len, 
							int64_t *_count, 
							CP_OCC **__cp_occ, 
							int8_t **__sa_ms_byte, uint32_t **__sa_ls_word, 
							int64_t *_sentinel_index);
	
	if (__load_BWT_on_shm(prefix, NULL, NULL, NULL, NULL, NULL, NULL) != 0) {
		fprintf(stderr, "ERROR: failed to load shm for BWT index\n");
		ret = -1;
	}
	
	return ret;
}

static int __bwa_shm_load_pac(const char *prefix, size_t size) {
	fprintf(stderr, "INFO: load pac (size: %ld)\n", bwa_shm_size_pac(bwa_shm_rlen()));
	
	uint8_t *load_pac_file(const char *prefix, int64_t l_pac);
	if (load_pac_file(prefix, (bwa_shm_rlen() - 1) / 2) == NULL) {
		fprintf(stderr, "ERROR: failed to load shm for PAC\n");
		return -1;
	}

	return 0;
}

static int __bwa_shm_load_ref(const char *prefix, size_t size) {
	int _load_ref_string(const char *prefix, uint8_t **ret_ptr);
	fprintf(stderr, "INFO: load reference string (size: %ld)\n", size);

	if (_load_ref_string(prefix, NULL)) {
		fprintf(stderr, "ERROR: failed to load shm for reference string\n");
		return -1;
	} 

	return 0;
}

static int __bwa_shm_load_kmer(const char *prefix, size_t size) {
	int _load_kmer_table(const char *prefix, uint64_t **ret_ptr);
	fprintf(stderr, "INFO: load kmer table (size: %ld)\n", size);

	if (_load_kmer_table(prefix, NULL)) {
		fprintf(stderr, "ERROR: failed to load shm for kmer table\n");
		return -1;
	} 

	return 0;
}

static int __bwa_shm_load_mlt(const char *prefix, size_t size) {
	int _load_mlt_table(const char *prefix, uint8_t **ret_ptr);
	fprintf(stderr, "INFO: load mlt table (size: %ld)\n", size);

	if (_load_mlt_table(prefix, NULL)) {
		fprintf(stderr, "ERROR: failed to load shm for mlt table\n");
		return -1;
	} 

	return 0;
}

#ifdef PERFECT_MATCH
static int __bwa_shm_load_perfect(const char *prefix, int pt_seed_len, 
									uint32_t num_seed_load) 
{
	int ____load_perfect_table_on_shm(char *file_name, int len, uint32_t num_seed_load, perfect_table_t **ret_ptr);
	
	char filename[PATH_MAX];

	/* to directly call ____load_perfect_table_on_shm(), it is required to set perfect_table_seed_len. */
	perfect_table_seed_len = pt_seed_len; 
	get_perfect_table_filename(filename, prefix, pt_seed_len);
	if (____load_perfect_table_on_shm(filename, pt_seed_len, num_seed_load, NULL)) {
		fprintf(stderr, "ERROR: failed to load shm for perfect hash table with seedlen=%d\n", pt_seed_len);
		return -1;
	}
	
	return 0;
}
#endif

#ifdef MEMSCALE
static int __bwa_shm_resize_perfect(const char *prefix, int pt_seed_len,
									size_t size_head, size_t size_loc,
									uint32_t old_num, uint32_t new_num)
{
	/* resize the shm */
	size_t size = size_head + size_loc + sizeof(seed_entry_t) * new_num;

	/* we cannot resize a shared memory region if hugetlb enabled. */
	assert(!use_hugetlb(BWA_SHM_PERFECT));

	if (ftruncate(shm_fd[BWA_SHM_PERFECT], size)) {
		fprintf(stderr, "[bwa_shm] %s: failed to truncate memory for perfect table to 0x%lx. errno: %d\n", 
			__func__, size, errno);
		return -1;
	} 
	
	if (new_num > old_num) {
		/* load if needed */
		char filename[PATH_MAX];
		FILE *fp;
		perfect_table_t pt;
		get_perfect_table_filename(filename, prefix, pt_seed_len);
	
		fp = xopen(filename, "rb");
		if (!fp) {
			fprintf(stderr, "%s: failed to open %s\n", __func__, filename);
			return -1;
		}
		
		__lpt_link_shm_to_pt(&pt, (perfect_table_t *) shm_ptr[BWA_SHM_PERFECT]);
		____lpt_load_seed_table(&pt, fp, old_num, new_num);
	}
	return 0;
}
#endif

#ifdef SMEM_ACCEL
static int __bwa_shm_load_accel(const char *prefix, 
				const int smem_all_on, const int smem_last_on) 
{
	int _load_smem_table(const char *prefix, 
						all_smem_t **__all_smem_table, 
						last_smem_t **__last_smem_table);
	
	all_smem_t *all_smem_table = NULL;
	last_smem_t *last_smem_table = NULL;

	if (_load_smem_table(prefix,
						smem_all_on ? &all_smem_table : NULL,
						smem_last_on ? &last_smem_table : NULL)) {
		fprintf(stderr, "ERROR: failed to load shm for smem accel index\n");
		return -1;
	} else
		return 0;
}
#endif

int __bwa_shm_load(const char *prefix, 
						enum hugetlb_mode huge_mode, int huge_force, 
						int pt_seed_len __maybe_unused, int pt_mmap __maybe_unused,
						size_t gb_limit __maybe_unused) 
{
	int ret = 0;
	bwa_shm_info_t *old_info;
	bwa_shm_info_t *new_info;
	int got_huge = 0;
	size_t huge_unit;
	size_t size_total = 0;
	size_t size_bwt, size_pac, size_ref;
	size_t size_kmer, size_mlt;
#ifdef PERFECT_MATCH
	size_t size_pt, size_pt_head, size_pt_loc, size_pt_seed;
	uint32_t num_pt_loc, num_pt_seed;
#endif
#ifdef SMEM_ACCEL
	size_t size_all_smem, size_last_smem;
#endif
#ifdef MEMSCALE
	size_t size_load = 0;
	size_t limit = 0, limit_min = 0, limit_max = 0;
	ssize_t rem;
#else
#define size_load size_total
#endif
	
	lock_bwa_shm_info();
	old_info = (bwa_shm_info_t *) malloc(bwa_shm_size_info(bwa_shm_info->ref_file_name_len));
	new_info = (bwa_shm_info_t *) malloc(bwa_shm_size_info(bwa_shm_info->ref_file_name_len));
	if (new_info)	
		memcpy(new_info, bwa_shm_info, 
					bwa_shm_size_info(bwa_shm_info->ref_file_name_len));
	unlock_bwa_shm_info();

	if (old_info == NULL || new_info == NULL) {
		fprintf(stderr, "ERROR: cannot allocate memory to buffer for bwa_shm_info\n");
		ret = -1;
		goto out;
	}

reset_newinfo:
	/* set new info */
	new_info->hugetlb_flags = get_hugetlb_flag(huge_mode);

	huge_unit = get_hugetlb_unit(huge_mode);

	size_bwt = __aligned_size(bwa_shm_size_bwt(bwa_shm_rlen()), huge_unit);
	size_pac = __aligned_size(bwa_shm_size_pac(bwa_shm_rlen()), huge_unit);
	size_ref = __aligned_size(bwa_shm_size_ref(bwa_shm_rlen()), huge_unit);
	size_kmer = __aligned_size(bwa_shm_size_kmer(), huge_unit);
	size_mlt = __aligned_size(bwa_shm_size_mlt(new_info), huge_unit);
	if (!new_info->useErt)
		size_total = size_bwt + size_pac + size_ref;
	else
		size_total = size_pac + size_ref + size_kmer + size_mlt;

#ifdef PERFECT_MATCH
	if (pt_seed_len > 0) {
		size_pt = get_perfect_table_size(prefix, pt_seed_len,
									&size_pt_head, &size_pt_loc, &size_pt_seed,
									&num_pt_loc, &num_pt_seed);
		size_pt = __aligned_size(size_pt, huge_unit);
		size_total += size_pt;
		new_info->pt_num_loc_entry = num_pt_loc;  
		new_info->pt_num_seed_entry = num_pt_seed;  
	} else {
		fprintf(stderr, "[memscale] Read length is not given. Perfect match table will not be used.\n");
		size_pt = 0;
		size_pt_head = 0;
		size_pt_loc = 0;
		size_pt_seed = 0;
	}
	new_info->pt_mmap = pt_mmap;
	new_info->pt_seed_len = pt_seed_len;
	new_info->pt_num_loc_entry = num_pt_loc;  
	new_info->pt_num_seed_entry = num_pt_seed;  
#endif
#ifdef SMEM_ACCEL
	size_all_smem = __aligned_size(ALL_SMEM_TABLE_SIZE, huge_unit);
	size_last_smem = __aligned_size(LAST_SMEM_TABLE_SIZE, huge_unit);
	size_total += size_all_smem + size_last_smem;
#endif

#ifdef MEMSCALE
	limit = gb_limit << 30;
	limit_min = size_bwt + size_pac + size_ref;
	limit_max = size_pac + size_ref + size_pt + size_kmer + size_mlt;
	if (limit == 0) {
		fprintf(stderr, "[memscale] gb_limit is set to the max value (%.1f)\n",
						B2GB_DOUBLE(limit_max));
		limit = limit_max;
	} else if (limit < limit_min) {
		fprintf(stderr, "[memscale] gb_limit (%ld) is less than the min value (%.1f). "
						"Set it to the min value.\n", gb_limit,  B2GB_DOUBLE(limit_min));
		limit = limit_min;
	} else if (limit >= limit_max) {
		fprintf(stderr, "[memscale] gb_limit (%ld) is larger than the max value (%.1f). "
						"Set it to the max value.\n", gb_limit, B2GB_DOUBLE(limit_max));
		limit = limit_max;
	}

	if ((limit & ((1L << 30) - 1)) == 0L)
		fprintf(stderr, "[memscale] gb_limit set: %ldGB\n", B2GB(limit));
	else
		fprintf(stderr, "[memscale] gb_limit set: %.1fGB\n", B2GB_DOUBLE(limit));

	rem = limit;
	new_info->bwt_on = 1;
	rem -= size_bwt;
	new_info->pac_on = 1;
	rem -= size_pac;
	new_info->ref_on = 1;
	rem -= size_ref;
	size_load = size_bwt + size_pac + size_ref;

	/* among the optional indices, all_smem and last_smem have the best capacity-performance ratio. */
	if (rem >= size_all_smem) {
		new_info->smem_all_on = 1;
		rem -= size_all_smem;
		size_load += size_all_smem;
	} else {
		new_info->smem_all_on = 0;
	}
	
	if (rem >= size_last_smem) {
		new_info->smem_last_on = 1;
		rem -= size_last_smem;
		size_load += size_last_smem;
	} else {
		new_info->smem_last_on = 0;
	}

	/* the second best is the perfect matching. */
	if (pt_seed_len > 0
		&& rem >= __aligned_size(size_pt_head + size_pt_loc
					+ __aligned_size(sizeof(seed_entry_t), 64), 
					huge_unit)) {
		size_t size_load_pt;
		size_t num_seed_load;
		size_t rem_pt = (rem / huge_unit) * huge_unit;
		new_info->perfect_on = 1;
		rem_pt -= size_pt_head + size_pt_loc;
		rem_pt -= (rem_pt % 64); /* aligned range for seed entries, in fact, unnecessary now. */
		num_seed_load = rem_pt / sizeof(seed_entry_t);
		if (num_seed_load > new_info->pt_num_seed_entry)
			num_seed_load = new_info->pt_num_seed_entry;
		size_load_pt = __aligned_size(size_pt_head + size_pt_loc
										+ __aligned_size(num_seed_load * sizeof(seed_entry_t), 64), 
									huge_unit);
		rem -= size_load_pt;
		size_load += size_load_pt;
		new_info->pt_num_seed_entry_loaded = num_seed_load;
	} else {
		new_info->perfect_on = 0;
		new_info->pt_num_seed_entry_loaded = 0;
	}
	
	/* check whether loading ERT tables is possible */
	if (size_kmer + size_mlt <= rem + size_bwt
								+ (new_info->smem_all_on ? size_all_smem : 0)
								+ (new_info->smem_last_on ? size_last_smem : 0)) {
		new_info->bwt_on = 0;
		rem += size_bwt;
		size_load -= size_bwt;
		if (new_info->smem_all_on) {
			new_info->smem_all_on = 0;
			rem += size_all_smem;
			size_load -= size_all_smem;
		}
		if (new_info->smem_last_on) {
			new_info->smem_last_on = 0;
			rem += size_last_smem;
			size_load -= size_last_smem;
		}
		rem -= size_kmer + size_mlt;
		size_load += size_kmer + size_mlt;
		new_info->kmer_on = 1;
		new_info->mlt_on = 1;

		new_info->useErt = 1;
	} else {
		new_info->useErt = 0;
	}
#endif
	
	if (got_huge == 0 &&
				get_hugepages(&huge_mode, size_load, huge_force)) 
	{
		fprintf(stderr, "ERROR: failed to get hugepages\n");
		
		if (huge_force) {
			ret = -1;
			goto out;
		} else {
			got_huge = 1;
			goto reset_newinfo;
		}
	}

	/* remove if needed */
#ifdef MEMSCALE
	lock_bwa_shm_info();
	if (new_info->hugetlb_flags != bwa_shm_info->hugetlb_flags) {
		fprintf(stderr, "[memscale] hugetlb_flags changes. Reload all.\n");
		if (bwa_shm_info->bwt_on == 1) {
			__bwa_shm_remove(BWA_SHM_BWT);
			bwa_shm_info->bwt_on = 0;
		}

		if (bwa_shm_info->pac_on == 1) {
			__bwa_shm_remove(BWA_SHM_PAC);
			bwa_shm_info->pac_on = 0;
		}

		if (bwa_shm_info->ref_on == 1) {
			__bwa_shm_remove(BWA_SHM_REF);
			bwa_shm_info->ref_on = 0;
		}
		
		if (bwa_shm_info->kmer_on == 1) {
			__bwa_shm_remove(BWA_SHM_KMER);
			bwa_shm_info->kmer_on = 0;
		}
		
		if (bwa_shm_info->mlt_on == 1) {
			__bwa_shm_remove(BWA_SHM_MLT);
			bwa_shm_info->mlt_on = 0;
		}

		if (bwa_shm_info->perfect_on == 1) {
			__bwa_shm_remove(BWA_SHM_PERFECT);
			bwa_shm_info->perfect_on = 0;
		}

		if (bwa_shm_info->smem_all_on == 1) {
			__bwa_shm_remove(BWA_SHM_SALL);
			bwa_shm_info->smem_all_on = 0;
		}

		if (bwa_shm_info->smem_last_on == 1) {
			__bwa_shm_remove(BWA_SHM_SLAST);
			bwa_shm_info->smem_last_on = 0;
		}
		bwa_shm_info->hugetlb_flags = new_info->hugetlb_flags;
	}

	if (bwa_shm_info->perfect_on == 1
			&& new_info->pt_seed_len != bwa_shm_info->pt_seed_len) {
		fprintf(stderr, "[memscale] pt_seed_len is changed. Reload perfect_table.\n");
		__bwa_shm_remove(BWA_SHM_PERFECT);
		bwa_shm_info->perfect_on = 0;
	}
	
	if (bwa_shm_info->perfect_on == 1
			&& new_info->pt_mmap != bwa_shm_info->pt_mmap) {
		fprintf(stderr, "[memscale] pt_mmap is changed. Reload perfect_table.\n");
		__bwa_shm_remove(BWA_SHM_PERFECT);
		bwa_shm_info->perfect_on = 0;
	}
	
	if (bwa_shm_info->perfect_on == 1 && bwa_shm_info->hugetlb_flags != 0
			&& new_info->pt_num_seed_entry_loaded != bwa_shm_info->pt_num_seed_entry_loaded) {
		fprintf(stderr, "[memscale] with hugetlb, truncate shared memory size is not supported. Since perfect_num_seed_load changes, reload perfect_table.\n");
		__bwa_shm_remove(BWA_SHM_PERFECT);
		bwa_shm_info->perfect_on = 0;
	}

	if (bwa_shm_info->bwt_on == 1 && new_info->bwt_on == 0) {
		__bwa_shm_remove(BWA_SHM_BWT);
		bwa_shm_info->bwt_on = 0;
	}
	
	if (bwa_shm_info->kmer_on == 1 && new_info->kmer_on == 0) {
		__bwa_shm_remove(BWA_SHM_KMER);
		bwa_shm_info->kmer_on = 0;
	}
	
	if (bwa_shm_info->mlt_on == 1 && new_info->mlt_on == 0) {
		__bwa_shm_remove(BWA_SHM_MLT);
		bwa_shm_info->mlt_on = 0;
	}
	
	if (bwa_shm_info->perfect_on == 1 && new_info->perfect_on == 0) {
		__bwa_shm_remove(BWA_SHM_PERFECT);
		bwa_shm_info->perfect_on = 0;
	}
	
	if (bwa_shm_info->smem_all_on == 1 && new_info->smem_all_on == 0) {
		__bwa_shm_remove(BWA_SHM_SALL);
		bwa_shm_info->smem_all_on = 0;
	}
	
	if (bwa_shm_info->smem_last_on == 1 && new_info->smem_last_on == 0) {
		__bwa_shm_remove(BWA_SHM_SLAST);
		bwa_shm_info->smem_last_on = 0;
	}

	memcpy(old_info, bwa_shm_info, 
				bwa_shm_size_info(bwa_shm_info->ref_file_name_len));
	unlock_bwa_shm_info();
#else
	__bwa_shm_remove(BWA_SHM_BWT);
	__bwa_shm_remove(BWA_SHM_PAC);
	__bwa_shm_remove(BWA_SHM_REF);
	__bwa_shm_remove(BWA_SHM_KMER);
	__bwa_shm_remove(BWA_SHM_MLT);
#ifdef PERFECT_MATCH
	__bwa_shm_remove(BWA_SHM_PERFECT);
#endif
#ifdef SMEM_ACCEL
	__bwa_shm_remove(BWA_SHM_SALL);
	__bwa_shm_remove(BWA_SHM_SLAST);
#endif
#endif
	
#ifndef MEMSCALE
#undef size_load
#endif

#ifdef MEMSCALE
	show_bwa_shm_info(old_info, "old_info");
#endif
	show_bwa_shm_info(new_info, "new_info");

	loading_info = new_info;
#ifdef MEMSCALE
	if (old_info->bwt_on == 0 && new_info->bwt_on == 1
			&& __bwa_shm_load_BWT(prefix, size_bwt)) {
		ret = -1;
		goto out;
	}

	if (old_info->pac_on == 0 && __bwa_shm_load_pac(prefix, size_pac)) {
		ret = -1;
		goto out;
	}

	if (old_info->ref_on == 0 && __bwa_shm_load_ref(prefix, size_ref)) {
		ret = -1;
		goto out;
	}

	if (old_info->kmer_on == 0 && new_info->kmer_on == 1
			&& __bwa_shm_load_kmer(prefix, size_kmer)) {
		ret = -1;
		goto out;
	}
	
	if (old_info->mlt_on == 0 && new_info->mlt_on == 1
			&& __bwa_shm_load_mlt(prefix, size_mlt)) {
		ret = -1;
		goto out;
	}

	if (new_info->perfect_on) {
		if (old_info->perfect_on == 0) {
			if (__bwa_shm_load_perfect(prefix, pt_seed_len, 
						new_info->pt_num_seed_entry_loaded)) {
					ret = -1;
					goto out;
				}
		} else {
			if (__bwa_shm_resize_perfect(prefix, pt_seed_len,
						size_pt_head, size_pt_loc,
						old_info->pt_num_seed_entry_loaded,
						new_info->pt_num_seed_entry_loaded)) {
					ret = -1;
					goto out;
				}
		}
	}

	if (new_info->smem_all_on == 1 || new_info->smem_last_on == 1) {
		if (__bwa_shm_load_accel(prefix, new_info->smem_all_on, new_info->smem_last_on)) {
			ret = -1;
			goto out;
		}
	}

#else /* !MEMSCALE */
	if (__bwa_shm_load_pac(prefix, size_pac)) {
		ret = -1;
		goto out;
	}

	if (__bwa_shm_load_ref(prefix, size_ref)) {
		ret = -1;
		goto out;
	}
	
	if (!new_info->useErt) {
		if (__bwa_shm_load_BWT(prefix, size_bwt)) {
			ret = -1;
			goto out;
		}

#ifdef SMEM_ACCEL
		if (__bwa_shm_load_accel(prefix, 1, 1)) {
			ret = -1;
			goto out;
		}
#endif
	} else {
		if (__bwa_shm_load_kmer(prefix, size_kmer)) {
			ret = -1;
			goto out;
		}

		if (__bwa_shm_load_mlt(prefix, size_mlt)) {
			ret = -1;
			goto out;
		}
	}
#ifdef PERFECT_MATCH
	if (__bwa_shm_load_perfect(prefix, pt_seed_len, 0)) {
		ret = -1;
		goto out;
	}
#endif
#endif /* !MEMSCALE */
	
	lock_bwa_shm_info();
#define copy_struct_var(dst, src, var) (dst)->var = (src)->var
	copy_struct_var(bwa_shm_info, new_info, hugetlb_flags);
	copy_struct_var(bwa_shm_info, new_info, useErt);
#ifdef MEMSCALE
	copy_struct_var(bwa_shm_info, new_info, bwt_on);
	copy_struct_var(bwa_shm_info, new_info, pac_on);
	copy_struct_var(bwa_shm_info, new_info, ref_on);
	copy_struct_var(bwa_shm_info, new_info, kmer_on);
	copy_struct_var(bwa_shm_info, new_info, mlt_on);
	copy_struct_var(bwa_shm_info, new_info, perfect_on);
	copy_struct_var(bwa_shm_info, new_info, smem_all_on);
	copy_struct_var(bwa_shm_info, new_info, smem_last_on);
	copy_struct_var(bwa_shm_info, new_info, pt_num_seed_entry_loaded);
#endif
#ifdef PERFECT_MATCH
	copy_struct_var(bwa_shm_info, new_info, pt_num_loc_entry);
	copy_struct_var(bwa_shm_info, new_info, pt_num_seed_entry);
	copy_struct_var(bwa_shm_info, new_info, pt_seed_len);
	copy_struct_var(bwa_shm_info, new_info, pt_mmap);
#endif
#undef copy_struct_var
	unlock_bwa_shm_info();
out:
	loading_info = NULL;
	if (new_info) free(new_info);
	if (old_info) free(old_info);
	return ret;
}

int bwa_shm_load(int argc, char *argv[]) {
	int ret = -1;	
	enum hugetlb_mode hugetlb_mode;
	int opt_force = 0;
	int opt_modify = 0;
	int useErt = DEFAULT_USE_ERT;
#ifdef MEMSCALE
	enum bwa_shm_init_mode init_mode = BWA_SHM_INIT_NEW;
	int opt_gb = 0;
#else
	const enum bwa_shm_init_mode init_mode = BWA_SHM_INIT_NEW;
	const int opt_gb = 0;
#endif
	int i;
	char c;
	char *prefix;
#ifdef PERFECT_MATCH
	int pt_seed_len = 0;
	int pt_mmap = DEFAULT_MMAP_PERFECT;
#else
	const int pt_seed_len = 0;
	const int pt_mmap = 0;
#endif

	hugetlb_mode = BWA_SHM_NORMAL_PAGE;

    /* Parse input arguments */
    while ((c = getopt(argc, argv, "fH:mg:l:p:Z:")) >= 0)
    {
		if (c == 'f') opt_force = 1;
        else if (c == 'H') hugetlb_mode = parse_hugetlb_mode(optarg);
		else if (c == 'Z') useErt = atoi(optarg) ? 1 : 0;
#ifdef MEMSCALE
        else if (c == 'm') {
			opt_modify = 1;
			init_mode = BWA_SHM_INIT_MODIFY;
        } else if (c == 'g') opt_gb = atoi(optarg);
#endif
		else if (c == 'l') {
#ifdef PERFECT_MATCH
			int seed_len = atoi(optarg);
			if (seed_len <= 0) {
				fprintf(stderr, "ERROR: The hash seed length for perfect match should be larger than 0.\n");
				exit(EXIT_FAILURE);
			}
			pt_seed_len = seed_len;
		} else if (c == 'p') {
			pt_mmap = atoi(optarg);
			if (pt_mmap)
				fprintf(stderr, "INFO: The perfect table will be mmap()ed directly from the filesystem.\n");
			else
				fprintf(stderr, "INFO: The perfect table will be loaded on memory.\n");
#endif
		} else {
			fprintf(stderr, "ERROR: Unknown option: %c\n", c);
			exit(EXIT_FAILURE);
		}
	}

	if (optind != argc - 1) {
		if (optind == argc) {
			fprintf(stderr, "ERROR: you may not give index base.\n");
		} else if (optind > argc) {
			fprintf(stderr, "ERROR: you may add options after the index base.\n");
		}
		usage_bwa_shm_load();
		exit(EXIT_FAILURE);
	}
	
	prefix = argv[optind];

	/* check the hugetlb availability and modify the modes if required */
	if (check_hugetlb(hugetlb_mode)) {
		if (opt_force) {
			fprintf(stderr, "ERROR: hugetlb is not available.\n");
			exit(EXIT_FAILURE);
		} else {
			fprintf(stderr, "WARN: hugetlb is not available. use normal pages.\n");
			hugetlb_mode = BWA_SHM_NORMAL_PAGE;
		}
	}
	
	bwa_shm_init(prefix, &useErt, pt_seed_len, init_mode);
	
	if (bwa_shm_mode == BWA_SHM_DISABLE) {
		fprintf(stderr, "ERROR: failed to init shm\n");
		ret = -1;
		goto out;
	}

	if (opt_modify == 0 && bwa_shm_mode != BWA_SHM_RENEWAL) {
		fprintf(stderr, "ERROR: failed to remove the previous data\n");
		ret = -1;
		goto out;
	} 
	
	if (opt_modify == 1 && bwa_shm_mode != BWA_SHM_MATCHED) {
		fprintf(stderr, "ERROR: failed to load the previous data to modify\n");
		ret = -1;
		goto out;
	}
	fprintf(stderr, "========BWA_SHM_LOAD_BEGIN==========================================\n");
	ret = __bwa_shm_load(prefix, hugetlb_mode, opt_force, 
					pt_seed_len, pt_mmap, opt_gb);
	fprintf(stderr, "========BWA_SHM_LOAD_END============================================\n");

out:
	bwa_shm_final(init_mode);

	if (ret == 0)
		return 0;

	__bwa_shm_remove_all();
	__bwa_shm_remove_hugetlb();
	return ret;
}

int bwa_shm_remove() {
	__bwa_shm_init_data(NULL);
	__bwa_shm_remove_all();
	__bwa_shm_remove_hugetlb();
	return 0;
}
#endif
