/* The MIT License

   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>

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

   Modified Copyright (C) 2020 Intel Corporation, Heng Li.
   Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
   Heng Li <hli@jimmy.harvard.edu> 
*/


#ifndef AC_KSEQ_H
#define AC_KSEQ_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "memcpy_bwamem.h"
#include <pthread.h>
#include <semaphore.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifdef OPT_RW
#define KSEQ_BUF_SIZE (1 << 24)

#define KSTREAM_NUM_BUF 8

typedef struct {
	unsigned char *buf;
	int begin, end, is_eof;
	sem_t sem_r, sem_w;
} __kstream_buf_t;

#define __KS_TYPE(type_t)						\
	typedef struct __kstream_t {				\
		__kstream_buf_t buf[KSTREAM_NUM_BUF];	\
		int curr; /* buf id to read */			\
		int done; /* done. stop the thread */	\
		pthread_t thread;						\
		type_t f;								\
	} kstream_t;

#define ks_eof(ks) ((ks)->buf[(ks)->curr].is_eof \
				&& (ks)->buf[(ks)->curr].begin >= (ks)->buf[(ks)->curr].end)

#define __KS_BASIC(type_t, __read, __bufsize)						\
	static void *ks_read_thread(void *arg) {						\
		kstream_t *ks = (kstream_t *)arg;							\
		int curr = 0; /* write */									\
		int is_eof = 0;												\
		__kstream_buf_t *buf;										\
		/* small first read for early initialization */				\
		buf = &ks->buf[curr];										\
		sem_wait(&buf->sem_w);										\
		buf->begin = 0;												\
		buf->end = __read(ks->f, buf->buf,							\
						__bufsize < 16384 ? __bufsize : 16384);		\
		if (buf->end > 0) {											\
			buf->is_eof = 0;										\
		} else {													\
			buf->is_eof = 1;										\
			is_eof = 1;												\
		}															\
		sem_post(&buf->sem_r);										\
		curr = (curr + 1) % (KSTREAM_NUM_BUF);						\
		while (ks->done == 0 && is_eof == 0) {						\
			buf = &ks->buf[curr];									\
			sem_wait(&buf->sem_w);									\
			buf->begin = 0;											\
			buf->end = __read(ks->f, buf->buf, __bufsize);			\
			if (buf->end > 0) {										\
				buf->is_eof = 0;									\
			} else {												\
				buf->is_eof = 1;									\
				is_eof = 1;											\
			}														\
			sem_post(&buf->sem_r);									\
			curr = (curr + 1) % (KSTREAM_NUM_BUF);					\
		}															\
		return NULL;												\
	}																\
	static inline kstream_t *ks_init(type_t f)						\
	{																\
		int i;														\
		kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));	\
        assert(ks != NULL);                                         \
		ks->f = f;													\
		ks->done = 0;												\
		for (i = 0; i < KSTREAM_NUM_BUF; ++i) {						\
			ks->buf[i].begin = 0;									\
			ks->buf[i].end = 0;										\
			ks->buf[i].is_eof = 0;									\
			ks->buf[i].buf = (unsigned char*)malloc(__bufsize);		\
        	assert(ks->buf[i].buf != NULL);                         \
			sem_init(&(ks->buf[i].sem_r), 0, 0);					\
			sem_init(&(ks->buf[i].sem_w), 0, 1);					\
		}															\
		pthread_create(&ks->thread, NULL, ks_read_thread, ks);		\
		ks->curr = 0;												\
		sem_wait(&ks->buf[0].sem_r);								\
		return ks;													\
	}																\
	static inline void ks_destroy(kstream_t *ks)					\
	{																\
		if (ks) {													\
			int i;													\
			for (i = 0; i < KSTREAM_NUM_BUF; ++i) 					\
				sem_post(&(ks->buf[i].sem_w));						\
			(ks)->done = 1;											\
			pthread_join((ks)->thread, NULL);						\
			for (i = 0; i < KSTREAM_NUM_BUF; ++i) {					\
				free(ks->buf[i].buf);								\
				sem_destroy(&(ks->buf[i].sem_r));					\
				sem_destroy(&(ks->buf[i].sem_w));					\
			}														\
			free(ks);												\
		}															\
	}

#define __KS_GETC()												\
	static inline int ks_getc(kstream_t *ks)					\
	{															\
		__kstream_buf_t *buf = &ks->buf[ks->curr];				\
		if (buf->begin < buf->end)								\
			return (int)buf->buf[buf->begin++];					\
		else if (buf->is_eof == 0) {							\
			sem_post(&buf->sem_w);								\
			ks->curr = (ks->curr + 1) % KSTREAM_NUM_BUF;		\
			buf = &ks->buf[ks->curr];							\
			sem_wait(&buf->sem_r);								\
			if (buf->end > 0) 									\
				return (int)buf->buf[buf->begin++];				\
			else { assert(buf->is_eof == 1); return -1; }		\
		} else /* buf->is_eof == 1 */							\
			return -1;											\
	}

#define __KS_GETUNTIL_LINE()											\
	static int ks_getuntil_line2(kstream_t *ks, kstring_t *str, int *dret, int append) \
	{																	\
		int gotany = 0;													\
		if (dret) *dret = 0;											\
		str->l = append? str->l : 0;									\
		for (;;) {														\
			int i; 														\
			__kstream_buf_t *buf = &ks->buf[ks->curr];					\
			if (buf->begin >= buf->end) {								\
				if (buf->is_eof == 1)									\
					break;												\
				sem_post(&buf->sem_w);									\
				ks->curr = (ks->curr + 1) % KSTREAM_NUM_BUF;			\
				buf = &ks->buf[ks->curr];								\
				sem_wait(&buf->sem_r);									\
				if (buf->end == 0) {									\
					assert(buf->is_eof == 1);							\
					break;												\
				}														\
			}															\
			/*for (i = buf->begin; i < buf->end; ++i) */				\
				/*if (buf->buf[i] == '\n') break; */ 					\
			for (i = buf->begin; i < buf->end - 8; i += 8) 				\
				if (buf->buf[i] == '\n' ||								\
						buf->buf[i+1] == '\n' ||						\
						buf->buf[i+2] == '\n' ||						\
						buf->buf[i+3] == '\n' ||						\
						buf->buf[i+4] == '\n' ||						\
						buf->buf[i+5] == '\n' ||						\
						buf->buf[i+6] == '\n' ||						\
						buf->buf[i+7] == '\n')							\
					break;												\
			for ( ; i < buf->end; ++i)									\
				if (buf->buf[i] == '\n') break; 						\
			if (str->m - str->l < (size_t)(i - buf->begin + 1)) {		\
				str->m = str->l + (i - buf->begin) + 1;					\
				kroundup32(str->m);										\
				str->s = (char*)realloc(str->s, str->m);				\
			}															\
			gotany = 1;													\
			/*memcpy_bwamem(str->s + str->l, str->m - str->l, buf->buf + buf->begin, i - buf->begin, __FILE__, __LINE__);*/ \
			memcpy(str->s + str->l, buf->buf + buf->begin, i - buf->begin); \
			str->l = str->l + (i - buf->begin);							\
			buf->begin = i + 1;											\
			if (i < buf->end) {											\
				if (dret) *dret = buf->buf[i];							\
				break;													\
			}															\
		}																\
		if (!gotany && ks_eof(ks)) return -1;							\
		if (str->s == 0) {												\
			str->m = 1;													\
			str->s = (char*)calloc(1, 1);								\
            assert(str->s != NULL);                                     \
		} else if (str->l > 1 && str->s[str->l-1] == '\r') --str->l; 	\
		str->s[str->l] = '\0';											\
		return str->l;													\
	}																	\
	static inline int ks_getuntil_line(kstream_t *ks, kstring_t *str, int *dret) \
	{ return ks_getuntil_line2(ks, str, dret, 0); }

#define __KS_GETUNTIL_SPACE()											\
	static int ks_getuntil_space2(kstream_t *ks, kstring_t *str, int *dret, int append) \
	{																	\
		int gotany = 0;													\
		if (dret) *dret = 0;											\
		str->l = append? str->l : 0;									\
		for (;;) {														\
			int i;														\
			__kstream_buf_t *buf = &ks->buf[ks->curr];					\
			if (buf->begin >= buf->end) {								\
				if (buf->is_eof == 1)									\
					break;												\
				sem_post(&buf->sem_w);									\
				ks->curr = (ks->curr + 1) % KSTREAM_NUM_BUF;			\
				buf = &ks->buf[ks->curr];								\
				sem_wait(&buf->sem_r);									\
				if (buf->end == 0) {									\
					assert(buf->is_eof == 1);							\
					break;												\
				}														\
			}															\
			/*for (i = buf->begin; i < buf->end; ++i)*/					\
				/* if (isspace(buf->buf[i])) break;	*/					\
			for (i = buf->begin; i < buf->end - 8; i += 8) 				\
				if (isspace(buf->buf[i]) ||								\
						isspace(buf->buf[i+1]) ||						\
						isspace(buf->buf[i+2]) ||						\
						isspace(buf->buf[i+3]) ||						\
						isspace(buf->buf[i+4]) ||						\
						isspace(buf->buf[i+5]) ||						\
						isspace(buf->buf[i+6]) ||						\
						isspace(buf->buf[i+7]))							\
					break;												\
			for ( ; i < buf->end; ++i)									\
				if (isspace(buf->buf[i])) break; 						\
			if (str->m - str->l < (size_t)(i - buf->begin + 1)) {		\
				str->m = str->l + (i - buf->begin) + 1;					\
				kroundup32(str->m);										\
				str->s = (char*)realloc(str->s, str->m);				\
			}															\
			gotany = 1;													\
			/*memcpy_bwamem(str->s + str->l, str->m - str->l, buf->buf + buf->begin, i - buf->begin, __FILE__, __LINE__);*/ \
			memcpy(str->s + str->l, buf->buf + buf->begin, i - buf->begin); \
			str->l = str->l + (i - buf->begin);							\
			buf->begin = i + 1;											\
			if (i < buf->end) {											\
				if (dret) *dret = buf->buf[i];							\
				break;													\
			}															\
		}																\
		if (!gotany && ks_eof(ks)) return -1;							\
		if (str->s == 0) {												\
			str->m = 1;													\
			str->s = (char*)calloc(1, 1);								\
            assert(str->s != NULL);                                     \
		} 																\
		str->s[str->l] = '\0';											\
		return str->l;													\
	}																	\
	static inline int ks_getuntil_space(kstream_t *ks, kstring_t *str, int *dret) \
	{ return ks_getuntil_space2(ks, str, dret, 0); }

#define KSTREAM_INIT(type_t, __read, __bufsize) \
	__KS_TYPE(type_t)							\
	__KS_BASIC(type_t, __read, __bufsize)		\
	__KS_GETC()									\
	__KS_GETUNTIL_LINE()						\
	__KS_GETUNTIL_SPACE()	

#define __ks_rewind(ks, func_rewind) do {					\
	int i; 													\
	for (i = 0; i < KSTREAM_NUM_BUF; ++i) 					\
		sem_post(&(ks->buf[i].sem_w));						\
	(ks)->done = 1;											\
	pthread_join((ks)->thread, NULL);						\
	func_rewind((ks)->f);									\
	for (i = 0; i < KSTREAM_NUM_BUF; ++i) {					\
		sem_destroy(&(ks->buf[i].sem_r));					\
		sem_destroy(&(ks->buf[i].sem_w));					\
		ks->buf[i].begin = 0;								\
		ks->buf[i].end = 0;									\
		ks->buf[i].is_eof = 0;								\
		assert(ks->buf[i].buf != NULL);                     \
		sem_init(&(ks->buf[i].sem_r), 0, 0);				\
		sem_init(&(ks->buf[i].sem_w), 0, 1);				\
	}														\
	(ks)->curr = 0;											\
	(ks)->done = 0;											\
	pthread_create(&ks->thread, NULL, ks_read_thread, ks);	\
} while (0)

#define kseq_rewind(ks, func_rewind) do {	\
	__ks_rewind(ks->f); 					\
	(ks)->last_char = 0;					\
} while (0)

#define __KSEQ_BASIC(SCOPE, type_t)										\
	SCOPE kseq_t *kseq_init(type_t fd)									\
	{																	\
		kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));					\
        assert(s != NULL);                                              \
		s->f = ks_init(fd);												\
		return s;														\
	}																	\
	SCOPE void kseq_destroy(kseq_t *ks)									\
	{																	\
		if (!ks) return;												\
		free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s); \
		ks_destroy(ks->f);												\
		free(ks);														\
	}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
#define __KSEQ_READ(SCOPE) \
	SCOPE int64_t kseq_read(kseq_t *seq) \
	{ \
		int c; \
		kstream_t *ks = seq->f; \
		if (seq->last_char == 0) { /* then jump to the next header line */ \
			while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@'); \
			if (c == -1) return -1; /* end of file */ \
			seq->last_char = c; \
		} /* else: the first header char has been read in the previous call */ \
		seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */ \
		if (ks_getuntil_space(ks, &seq->name, &c) < 0) return -1; /* normal exit: EOF */ \
		if (c != '\n') ks_getuntil_line(ks, &seq->comment, 0); /* read FASTA/Q comment */ \
		if (seq->seq.s == 0) { /* we can do this in the loop below, but that is slower */ \
			seq->seq.m = 256; \
			seq->seq.s = (char*)malloc(seq->seq.m); \
            assert(seq->seq.s != NULL);             \
		} \
		while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') { \
			if (c == '\n') continue; /* skip empty lines */ \
			seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ \
			ks_getuntil_line2(ks, &seq->seq, 0, 1); /* read the rest of the line */ \
		} \
		if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */	\
		if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ \
			/*seq->seq.m = seq->seq.l + 2;*/							\
			seq->seq.m = seq->seq.l + 16;								\
			kroundup32(seq->seq.m); /* rounded to the next closest 2^k */ \
			seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); \
            assert(seq->seq.s != NULL); \
		} \
		seq->seq.s[seq->seq.l] = 0;	/* null terminated string */ \
		if (c != '+') return seq->seq.l; /* FASTA */ \
		if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */ \
			seq->qual.m = seq->seq.m; \
			seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); \
		} \
		while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */ \
		if (c == -1) return -2; /* error: no quality string */ \
		while (ks_getuntil_line2(ks, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l); \
		seq->last_char = 0;	/* we have not come to the next header line */ \
		if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */ \
		return seq->seq.l; \
	}

#else /* !OPT_RW */
#define KSEQ_BUF_SIZE (16384)

#define __KS_TYPE(type_t)						\
	typedef struct __kstream_t {				\
		unsigned char *buf;						\
		int begin, end, is_eof;					\
		type_t f;								\
	} kstream_t;

#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)


#define __KS_BASIC(type_t, __bufsize)								\
	static inline kstream_t *ks_init(type_t f)						\
	{																\
		kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));	\
        assert(ks != NULL);                                         \
		ks->f = f;													\
		ks->buf = (unsigned char*)malloc(__bufsize);				\
        assert(ks->buf != NULL);                                    \
		return ks;													\
	}																\
	static inline void ks_destroy(kstream_t *ks)					\
	{																\
		if (ks) {													\
			free(ks->buf);											\
			free(ks);												\
		}															\
	}

#define __KS_GETC(__read, __bufsize)						\
	static inline int ks_getc(kstream_t *ks)				\
	{														\
		if (ks->is_eof && ks->begin >= ks->end) return -1;	\
		if (ks->begin >= ks->end) {							\
			ks->begin = 0;									\
			ks->end = __read(ks->f, ks->buf, __bufsize);	\
			if (ks->end == 0) { ks->is_eof = 1; return -1;}	\
		}													\
		return (int)ks->buf[ks->begin++];					\
	}

#define __KS_GETUNTIL(__read, __bufsize)								\
	static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append) \
	{																	\
		int gotany = 0;													\
		if (dret) *dret = 0;											\
		str->l = append? str->l : 0;									\
		for (;;) {														\
			int i;														\
			if (ks->begin >= ks->end) {									\
				if (!ks->is_eof) {										\
					ks->begin = 0;										\
					ks->end = __read(ks->f, ks->buf, __bufsize);		\
					if (ks->end == 0) { ks->is_eof = 1; break; }		\
				} else break;											\
			}															\
			if (delimiter == KS_SEP_LINE) { \
				for (i = ks->begin; i < ks->end; ++i) \
					if (ks->buf[i] == '\n') break; \
			} else if (delimiter > KS_SEP_MAX) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (ks->buf[i] == delimiter) break;					\
			} else if (delimiter == KS_SEP_SPACE) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (isspace(ks->buf[i])) break;						\
			} else if (delimiter == KS_SEP_TAB) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break; \
			} else i = 0; /* never come to here! */						\
			if (str->m - str->l < (size_t)(i - ks->begin + 1)) {		\
				str->m = str->l + (i - ks->begin) + 1;					\
				kroundup32(str->m);										\
				str->s = (char*)realloc(str->s, str->m);				\
			}															\
			gotany = 1;													\
			memcpy_bwamem(str->s + str->l, str->m - str->l, ks->buf + ks->begin, i - ks->begin, __FILE__, __LINE__); \
			str->l = str->l + (i - ks->begin);							\
			ks->begin = i + 1;											\
			if (i < ks->end) {											\
				if (dret) *dret = ks->buf[i];							\
				break;													\
			}															\
		}																\
		if (!gotany && ks_eof(ks)) return -1;							\
		if (str->s == 0) {												\
			str->m = 1;													\
			str->s = (char*)calloc(1, 1);								\
            assert(str->s != NULL);                                     \
		} else if (delimiter == KS_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l; \
		str->s[str->l] = '\0';											\
		return str->l;													\
	}																	\
	static inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret) \
	{ return ks_getuntil2(ks, delimiter, str, dret, 0); }

#define KSTREAM_INIT(type_t, __read, __bufsize) \
	__KS_TYPE(type_t)							\
	__KS_BASIC(type_t, __bufsize)				\
	__KS_GETC(__read, __bufsize)				\
	__KS_GETUNTIL(__read, __bufsize)

#define kseq_rewind(ks) ((ks)->last_char = (ks)->f->is_eof = (ks)->f->begin = (ks)->f->end = 0)

#define __KSEQ_BASIC(SCOPE, type_t)										\
	SCOPE kseq_t *kseq_init(type_t fd)									\
	{																	\
		kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));					\
        assert(s != NULL);                                              \
		s->f = ks_init(fd);												\
		return s;														\
	}																	\
	SCOPE void kseq_destroy(kseq_t *ks)									\
	{																	\
		if (!ks) return;												\
		free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s); \
		ks_destroy(ks->f);												\
		free(ks);														\
	}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
#define __KSEQ_READ(SCOPE) \
	SCOPE int64_t kseq_read(kseq_t *seq) \
	{ \
		int c; \
		kstream_t *ks = seq->f; \
		if (seq->last_char == 0) { /* then jump to the next header line */ \
			while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@'); \
			if (c == -1) return -1; /* end of file */ \
			seq->last_char = c; \
		} /* else: the first header char has been read in the previous call */ \
		seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */ \
		if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1; /* normal exit: EOF */ \
		if (c != '\n') ks_getuntil(ks, KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */ \
		if (seq->seq.s == 0) { /* we can do this in the loop below, but that is slower */ \
			seq->seq.m = 256; \
			seq->seq.s = (char*)malloc(seq->seq.m); \
            assert(seq->seq.s != NULL);             \
		} \
		while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') { \
			if (c == '\n') continue; /* skip empty lines */ \
			seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ \
			ks_getuntil2(ks, KS_SEP_LINE, &seq->seq, 0, 1); /* read the rest of the line */ \
		} \
		if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */	\
		if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ \
			/*seq->seq.m = seq->seq.l + 2;*/							\
			seq->seq.m = seq->seq.l + 16;								\
			kroundup32(seq->seq.m); /* rounded to the next closest 2^k */ \
			seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); \
            assert(seq->seq.s != NULL); \
		} \
		seq->seq.s[seq->seq.l] = 0;	/* null terminated string */ \
		if (c != '+') return seq->seq.l; /* FASTA */ \
		if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */ \
			seq->qual.m = seq->seq.m; \
			seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); \
		} \
		while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */ \
		if (c == -1) return -2; /* error: no quality string */ \
		while (ks_getuntil2(ks, KS_SEP_LINE, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l); \
		seq->last_char = 0;	/* we have not come to the next header line */ \
		if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */ \
		return seq->seq.l; \
	}
#endif /* !OPT_RW */

#define __KSEQ_TYPE(type_t)						\
	typedef struct {							\
		kstring_t name, comment, seq, qual;		\
		int last_char;							\
		kstream_t *f;							\
	} kseq_t;

#define KSEQ_INIT2(SCOPE, type_t, __read)		\
	KSTREAM_INIT(type_t, __read, KSEQ_BUF_SIZE)			\
	__KSEQ_TYPE(type_t)							\
	__KSEQ_BASIC(SCOPE, type_t)					\
	__KSEQ_READ(SCOPE)

#define KSEQ_INIT(type_t, __read) KSEQ_INIT2(static, type_t, __read)

#define KSEQ_DECLARE(type_t)											\
	__KS_TYPE(type_t)													\
		__KSEQ_TYPE(type_t)												\
		extern kseq_t *kseq_init(type_t fd);							\
	void kseq_destroy(kseq_t *ks);										\
	int64_t kseq_read(kseq_t *seq);											\
	static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append);
	
#endif
