#include "structure.h"

static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
{
	const unsigned int bufsize = 0x1000000; // 16M block
	bwtint_t offset = 0;
	while (size) {
		bwtint_t x = bufsize < size ? bufsize : size;
		x = fread(((char*)a + offset), 1, x, fp);
		size -= x; offset += x;
	}
	return offset;
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = fopen(fn, "rb");
	fread(&primary, sizeof(bwtint_t), 1, fp);
	//xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	fread(skipped, sizeof(bwtint_t), 4, fp); // skip
	fread(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	fread(&primary, sizeof(bwtint_t), 1, fp);
	//xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	bwt->sa[0] = -1;

	fread_fix(fp, sizeof(bwtint_t) * (bwt->n_sa - 1), bwt->sa + 1);
	fclose(fp);
}

bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[10240];
	FILE *fp;
	const char *fname;
	bntseq_t *bns;
	long long xx;
	int i;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = fopen(fname = ann_filename, "r");
		(void)fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			(void)fscanf(fp, "%u%s", &p->gi, str);
			p->name = strdup(str);
			// read fasta comments 
			while (str - q < (int)(sizeof(str) - 1) && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			while (c != '\n' && c != EOF) c = fgetc(fp);

			*q = 0;
			if (q - str > 1) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			(void)fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			p->offset = xx;
		}
		fclose(fp);
	}
	{ // read .amb
		int32_t n_seqs;
		fp = fopen(fname = amb_filename, "r");
		(void)fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		//xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = bns->n_holes? (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t)) : 0;
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			(void)fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			p->offset = xx;
			p->amb = str[0];
		}
		fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = fopen(pac_filename, "rb");
	}
	return bns;
}

void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = fopen(fn, "rb");
	fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
	fseek(fp, 0, SEEK_SET);
	fread(&bwt->primary, sizeof(bwtint_t), 1, fp);
	fread(bwt->L2 + 1, sizeof(bwtint_t), 4, fp);
	fread_fix(fp, bwt->bwt_size << 2, bwt->bwt);
	bwt->seq_len = bwt->L2[4];
	fclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}

bwt_t *bwa_idx_load_bwt(const char *hint)
{
	char *tmp;
	bwt_t *bwt;

	tmp = (char*)calloc(strlen(hint) + 5, 1);
	strcat(strcpy(tmp, hint), ".bwt"); // FM-index
	bwt = bwt_restore_bwt(tmp);
	strcat(strcpy(tmp, hint), ".sa");  // partial suffix array (SA)
	bwt_restore_sa(tmp, bwt);
	free(tmp); tmp = NULL;

	return bwt;
}

bntseq_t *bns_restore(const char *prefix)
{  
	char ann_filename[256], amb_filename[256], pac_filename[256];
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");
	return bns_restore_core(ann_filename, amb_filename, pac_filename);
}

bwaidx_t *bwa_idx_load(const char *hint)
{
	bwaidx_t *idx;

	//fprintf(stderr, "\tLoad the reference index files...");
	idx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));
	idx->bwt = bwa_idx_load_bwt(hint);
	idx->bns = bns_restore(hint);
	idx->pac = (uint8_t*)calloc(idx->bns->l_pac/4+1, 1);
	//fprintf(stderr, "Done!\n");

	return idx;
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}

void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		if (bns->fp_pac) fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}

void bwa_idx_destroy(bwaidx_t *idx)
{
	if (idx == 0) return;
	if (idx->bwt) bwt_destroy(idx->bwt);
	if (idx->bns) bns_destroy(idx->bns);
	if (idx->pac) free(idx->pac);
	free(idx);
}

void *IdvLoadReferenceSequences(void *arg)
{
	int base, *my_id;
	int64_t fPos, rPos;

	my_id = (int*)arg;
	for (fPos = *my_id, rPos = TwoGenomeSize - fPos - 1 ; fPos < GenomeSize; fPos += iThreadNum, rPos-= iThreadNum)
	{
		base = RefIdx->pac[fPos >> 2] >> ((~fPos & 3) << 1) & 3;
		switch (base)
		{
		case 0: RefSequence[fPos] = 'A'; RefSequence[rPos] = 'T'; break;
		case 1: RefSequence[fPos] = 'C'; RefSequence[rPos] = 'G'; break;
		case 2: RefSequence[fPos] = 'G'; RefSequence[rPos] = 'C'; break;
		case 3: RefSequence[fPos] = 'T'; RefSequence[rPos] = 'A'; break;
		default:RefSequence[fPos] = RefSequence[rPos] = 'N';
		}
	}
	return (void*)(1);
}

void RestoreReferenceSequences()
{
	int i, *JobIDArr = new int[iThreadNum];

	pthread_t *ThreadArr = new pthread_t[iThreadNum];
	for (i = 0; i < iThreadNum; i++)
	{
		JobIDArr[i] = i;
		pthread_create(&ThreadArr[i], NULL, IdvLoadReferenceSequences, JobIDArr + i);
	}
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

	delete[] ThreadArr; delete[] JobIDArr;
}

void RestoreReferenceInfo()
{
	int i;
	int64_t iTotalLength = 0;

	GenomeSize = RefIdx->bns->l_pac; TwoGenomeSize = (GenomeSize << 1);
	ChromosomeVec.resize((iChromsomeNum = RefIdx->bns->n_seqs));

	fprintf(stderr, "\tLoad the reference sequences (%d %s)\n", RefIdx->bns->n_seqs, RefIdx->bns->n_seqs > 1 ? "chromosomes" : "chromosome");
	fseek(RefIdx->bns->fp_pac, 0, SEEK_SET);
	fread(RefIdx->pac, 1, GenomeSize / 4 + 1, RefIdx->bns->fp_pac);

	for (i = 0; i < iChromsomeNum; i++)
	{
		ChromosomeVec[i].len = RefIdx->bns->anns[i].len;
		ChromosomeVec[i].name = RefIdx->bns->anns[i].name;
		//fprintf(stderr, "%d --> %s\n", i, ChromosomeVec[i].name);

		ChromosomeVec[i].FowardLocation = iTotalLength; iTotalLength += ChromosomeVec[i].len;
		ChromosomeVec[i].ReverseLocation = TwoGenomeSize - iTotalLength;
		//printf("%s: %lld %lld\n", ChromosomeVec[i].name, ChromosomeVec[i].FowardLocation, ChromosomeVec[i].ReverseLocation);

		ChrLocMap.insert(make_pair(ChromosomeVec[i].FowardLocation + ChromosomeVec[i].len - 1, i));
		ChrLocMap.insert(make_pair(ChromosomeVec[i].ReverseLocation + ChromosomeVec[i].len - 1, i));
	}
	RefSequence = new char[TwoGenomeSize + 1]; RefSequence[TwoGenomeSize] = '\0';
	RestoreReferenceSequences();
	//fprintf(stderr, "\t(Current memory consumption: %d MB)\n\n", CheckMemoryUsage());

	//if (bDebugMode)
	//{
	//	printf("Reference sequences:\n");
	//	for (map<int64_t, int>::iterator iter = ChrLocMap.begin(); iter != ChrLocMap.end(); iter++)
	//		printf("\t%s [%lld -- %lld]\n", ChromosomeVec[iter->second].name, iter->first - ChromosomeVec[iter->second].len + 1, iter->first);
	//}
}
