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

#include "read_index_ele.h"
#include "bwa_shm.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "safe_mem_lib.h"
#include "safe_str_lib.h"
#ifdef __cplusplus
}
#endif

indexEle::indexEle()
{
    idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
    assert(idx != NULL);
}

indexEle::~indexEle()
{
    if (idx == 0) return;
    if (idx->mem == 0)
    {
        if (idx->bns) bns_destroy(idx->bns);


#ifdef USE_SHM
        if (idx->pac && bwa_shm_unmap(BWA_SHM_PAC))
			free(idx->pac);
#else /* !USE_SHM */
        if (idx->pac)
			free(idx->pac);
#endif /* !USE_SHM */
    } else {
        free(idx->bns->anns); free(idx->bns);
        if (!idx->is_shm) free(idx->mem);
    }
    free(idx);
}

#ifdef USE_SHM
uint8_t *load_pac_file(const char *prefix, int64_t l_pac __maybe_unused) {
	uint8_t *pac = NULL;
	if (__bwa_shm_load_file(prefix, ".pac", BWA_SHM_PAC, (void **) &pac) == 0)
		return pac;
	else {
		printf("ERROR: [%s] failed to load packed reference\n", __func__);
		exit(EXIT_FAILURE);
	}
}
#else
uint8_t *load_pac_file(const char *prefix, int64_t l_pac) {
	char pac_filename[PATH_MAX];
	FILE *fp;
	uint8_t *pac;
	strcpy_s(pac_filename, PATH_MAX, prefix);
	strcat_s(pac_filename, PATH_MAX, ".pac");
	fp = xopen(pac_filename, "rb");

	pac = (uint8_t*) calloc(l_pac/4+1, 1);
	assert(pac != NULL);
	err_fread_noeof(pac, 1, l_pac/4+1, fp); // concatenated 2-bit encoded sequence
	err_fclose(fp);
	return pac;
}
#endif

void indexEle::bwa_idx_load_ele(const char *hint, int which)
{
    char *prefix;
    int l_hint = strlen(hint);
    prefix = (char *) malloc(l_hint + 3 + 4 + 1);
    assert(prefix != NULL);
    strcpy_s(prefix, l_hint + 3 + 4 + 1, hint);

    fprintf(stderr, "* Index prefix: %s\n", prefix);
    
    // idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
    if (which & BWA_IDX_BNS) {
        int i, c;
        idx->bns = bns_restore(prefix);
        if (idx->bns == 0) {
            printf("Error!! : [%s] bns is NULL!!\n", __func__);
            exit(EXIT_FAILURE);
        }
        for (i = c = 0; i < idx->bns->n_seqs; ++i)
            if (idx->bns->anns[i].is_alt) ++c;
        
        fprintf(stderr, "* Read %d ALT contigs\n", c);
        
        if (which & BWA_IDX_PAC)
            idx->pac = load_pac_file(prefix, idx->bns->l_pac); // concatenated 2-bit encoded sequence
    }
    free(prefix);
}

#include <sys/file.h>
char* indexEle::bwa_idx_infer_prefix(const char *hint)
{
    char *prefix;
    int l_hint;
    FILE *fp;
    l_hint = strlen(hint);
    prefix = (char *) malloc(l_hint + 3 + 4 + 1);
    assert(prefix != NULL);
    strcpy_s(prefix, l_hint + 3 + 4 + 1, hint);
    strcpy_s(prefix + l_hint, 8, ".64.bwt");
    if ((fp = fopen(prefix, "rb")) != 0)
    {
        fclose(fp);
        prefix[l_hint + 3] = 0;
        return prefix;
    } else {
        strcpy_s(prefix + l_hint, 8, ".bwt");
        if ((fp = fopen(prefix, "rb")) == 0)
        {
            free(prefix);
            return 0;
        } else {
            //flock(fileno(fp), 1);
            //flock(fileno(fp), 1);  // Unlock the file
            fclose(fp);
            prefix[l_hint] = 0;
            return prefix;
        }
    }
}

#if TEST
//int main(int argc, char* argv[])
//{
//  printf("Testing read_index_ele...\n");
//  indexEle *bwaEle = new indexEle();
//  
//  bwaEle->bwa_idx_load_ele("/projects/PCL-GBB/wasim/read_and_ref_data_1/hgaa.fa",
//                          BWA_IDX_ALL);
//
//  delete bwaEle;
//}
#endif
