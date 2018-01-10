#include <cmath>
#include "structure.h"

#define WindowSize 80
#define MaximumGaps 300
#define SeedExplorationChunk 10000

int QueryChrIdx;
FILE *AlnFile = stdout;
vector<FragPair_t> SeedVec;
static pthread_mutex_t Lock;
vector<AlnBlock_t> AlnBlockVec;
int64_t LocalAlignmentNum = 0, SNP_num = 0, IND_num = 0, SVS_num = 0;

bool CompByQueryPos(const FragPair_t& p1, const FragPair_t& p2)
{
	if (p1.qPos == p2.qPos) return p1.rPos < p2.rPos;
	else return p1.qPos < p2.qPos;
}

void ShowAlnBlockDistance(int idx1, int idx2)
{
	int n, qPos1, qPos2;
	int64_t rPos1, rPos2;

	n = (int)AlnBlockVec[idx1].FragPairVec.size(); 
	qPos1 = AlnBlockVec[idx1].FragPairVec[n - 1].qPos + AlnBlockVec[idx1].FragPairVec[n - 1].qLen - 1;
	rPos1 = AlnBlockVec[idx1].FragPairVec[n - 1].rPos + AlnBlockVec[idx1].FragPairVec[n - 1].rLen - 1;

	qPos2 = AlnBlockVec[idx2].FragPairVec[0].qPos;
	rPos2 = AlnBlockVec[idx2].FragPairVec[0].rPos;

	printf("qDistance = %d, rDistance = %lld\n", qPos2 - qPos1, rPos2 - rPos1);
}

void ShowFragSeqs(FragPair_t& FragPair)
{
	char *frag1, *frag2;

	frag1 = new char[FragPair.qLen + 1]; strncpy(frag1, QueryChrVec[QueryChrIdx].seq.c_str() + FragPair.qPos, FragPair.qLen); frag1[FragPair.qLen] = '\0';
	frag2 = new char[FragPair.rLen + 1]; strncpy(frag2, RefSequence + FragPair.rPos, FragPair.rLen); frag2[FragPair.rLen] = '\0';

	printf("q[%d-%d]=%d r[%lld-%lld]=%d\n%s\n%s\n\n", FragPair.qPos, FragPair.qPos + FragPair.qLen - 1, FragPair.qLen, FragPair.rPos, FragPair.rPos + FragPair.rLen - 1, FragPair.rLen, frag1, frag2);
}

void ShowFragPairVec(vector<FragPair_t>& FragPairVec)
{
	printf("FragPairVec (N=%d)\n", (int)FragPairVec.size());
	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if(iter->bSeed) printf("\t\tq[%d-%d] r[%lld-%lld] len=%d\n", iter->qPos, iter->qPos + iter->qLen - 1, iter->rPos, iter->rPos + iter->rLen - 1, iter->qLen);
		else printf("\t\tq[%d-%d]=%d r[%lld-%lld]=%d\n", iter->qPos, iter->qPos + iter->qLen - 1, iter->qLen, iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen);
	}
}

void ShowAlnBlock(int idx)
{
	int qBeg, qEnd;
	int64_t rBeg, rEnd;

	AlnBlock_t& AlnBlock = AlnBlockVec[idx];

	int i, num = (int)AlnBlock.FragPairVec.size();

	printf("AlnBlock: score = %d\n", AlnBlock.score);
	if (idx > 0) ShowAlnBlockDistance(idx - 1, idx);

	qBeg = AlnBlock.FragPairVec[0].qPos; qEnd = AlnBlock.FragPairVec[num - 1].qPos + AlnBlock.FragPairVec[num - 1].qLen - 1;
	rBeg = AlnBlock.FragPairVec[0].rPos; rEnd = AlnBlock.FragPairVec[num - 1].rPos + AlnBlock.FragPairVec[num - 1].rLen - 1;

	printf("q[%d-%d]=%d r[%lld-%lld]=%lld\n", qBeg, qEnd, qEnd - qBeg + 1, rBeg, rEnd, rEnd - rBeg + 1);
	
	//ShowFragPairVec(AlnBlock.FragPairVec);

	printf("\n\n\n");
}

void *IdentifyLocalMEM(void *arg)
{
	FragPair_t seed;
	vector<FragPair_t> vec;
	bwtSearchResult_t bwtSearchResult;
	int i, pos, start, stop, num, len, *my_id = (int*)arg;

	string& seq = QueryChrVec[QueryChrIdx].seq; len = (int)seq.length(); seed.bSeed = true;
	for (pos = (*my_id*SeedExplorationChunk); pos < len; pos += (iThreadNum * SeedExplorationChunk))
	{
		start = pos; if((stop = start + SeedExplorationChunk) > len) stop = len;
		while (start < stop)
		{
			if (nst_nt4_table[(int)seq[start]] > 3) start++;
			else
			{
				bwtSearchResult = BWT_Search(seq, start, stop);
				if (bwtSearchResult.freq > 0)
				{
					seed.rLen = seed.qLen = bwtSearchResult.len; seed.qPos = start;
					for (i = 0; i != bwtSearchResult.freq; i++)
					{
						seed.rPos = bwtSearchResult.LocArr[i];
						seed.PosDiff = seed.rPos - seed.qPos;
						vec.push_back(seed);
					}
					delete[] bwtSearchResult.LocArr;
				}
				start += (bwtSearchResult.len + 1);
			}
		}
		if (*my_id == 0) fprintf(stderr, "\r\t\tSeed exploration: %d / %d (%d%%)...", stop, len, (int)(100 * ((1.0*stop / len))));
	}
	if (*my_id == 0) fprintf(stderr, "\r\t\tSeed exploration: %d / %d (100%%)...done!\n", len, len);
	std::sort(vec.begin(), vec.end(), CompByQueryPos);

	if (iThreadNum == 1) SeedVec.swap(vec);
	else
	{
		pthread_mutex_lock(&Lock);

		num = (int)SeedVec.size();
		copy(vec.begin(), vec.end(), back_inserter(SeedVec));
		inplace_merge(SeedVec.begin(), SeedVec.begin() + num, SeedVec.end(), CompByQueryPos);

		pthread_mutex_unlock(&Lock);
	}

	return (void*)(1);
}

int CalAlnBlockScore(vector<FragPair_t>& FragPairVec)
{
	int score = 0;
	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++) score += iter->qLen;

	return score;
}

void AlignmentBlockClustering()
{
	bool *VisitArr;
	AlnBlock_t AlnBlock;
	int i, j, curIdx, num;

	num = (int)SeedVec.size(); SeedVec.resize(num + 1); SeedVec[num].qLen = 0; SeedVec[num].PosDiff = 0; num += 1;

	VisitArr = new bool[num]();
	for (i = 0; i < num; i++)
	{
		if (VisitArr[i] == false)
		{
			//printf("ClusterHead: q[%d-%d] r[%lld-%lld] len=%d PD=%lld\n", SeedVec[i].qPos, SeedVec[i].qPos + SeedVec[i].qLen - 1, SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].qLen, SeedVec[i].PosDiff);
			VisitArr[i] = true; AlnBlock.FragPairVec.clear(); AlnBlock.FragPairVec.push_back(SeedVec[i]);
			for (curIdx = i, j = i + 1; j < num; j++)
			{
				if (VisitArr[j] == false)
				{
					if ((abs(SeedVec[j].PosDiff - SeedVec[curIdx].PosDiff)) < MaximumGaps && (SeedVec[j].qPos - (SeedVec[curIdx].qPos + SeedVec[curIdx].qLen)) < MaximumGaps)
					{
						VisitArr[(curIdx = j)] = true;
						AlnBlock.FragPairVec.push_back(SeedVec[j]);
					}
				}
				if (j - curIdx > MaxSeedFreq) break;
			}
			if((AlnBlock.score = CalAlnBlockScore(AlnBlock.FragPairVec)) > 100) AlnBlockVec.push_back(AlnBlock);
		}
	}
	delete[] VisitArr;
}

void RemoveRedundantAlnBlocks()
{
	int i, j, n, thr, num, iBeg, iEnd, jBeg, jEnd;

	num = (int)AlnBlockVec.size();

	for (i = 0; i < num; i++)
	{
		if (AlnBlockVec[i].score == 0) continue;

		n = (int)AlnBlockVec[i].FragPairVec.size();

		iBeg = AlnBlockVec[i].FragPairVec[0].qPos; iEnd = AlnBlockVec[i].FragPairVec[n - 1].qPos + AlnBlockVec[i].FragPairVec[n - 1].qLen - 1;

		for (j = i + 1; j < num; j++)
		{
			if (AlnBlockVec[i].score == 0) continue;

			n = (int)AlnBlockVec[j].FragPairVec.size(); thr = (AlnBlockVec[i].score >> 1);

			jBeg = AlnBlockVec[j].FragPairVec[0].qPos; jEnd = AlnBlockVec[j].FragPairVec[n - 1].qPos + AlnBlockVec[j].FragPairVec[n - 1].qLen - 1;

			if (jBeg >= iBeg && jEnd <= iEnd) // inclusive
			{
				if(AlnBlockVec[j].score < thr) AlnBlockVec[j].score = 0;
			}
			else break;
		}
	}
	vector<AlnBlock_t> TmpVec;
	for (vector<AlnBlock_t>::iterator iter = AlnBlockVec.begin(); iter != AlnBlockVec.end(); iter++)
	{
		if (iter->score > 0) TmpVec.push_back(*iter);
	}
	AlnBlockVec.swap(TmpVec);
}

void RemoveSmallAlnBlocks()
{
	vector<AlnBlock_t> TmpVec;
	vector<AlnBlock_t>::iterator ABiter;
	vector<FragPair_t>::iterator FragPairIter;

	for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		ABiter->score = 0;
		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++) ABiter->score += FragPairIter->qLen;
		if (ABiter->score < 50) ABiter->score = 0;
		else TmpVec.push_back(*ABiter);
	}
	AlnBlockVec.swap(TmpVec);
}


bool RemoveOverlaps(vector<FragPair_t>& FragPairVec)
{
	bool bNullPair;
	int64_t irBeg, irEnd, jrBeg, jrEnd;
	int iqBeg, iqEnd, iLen, jqBeg, jqEnd, jLen;
	int i, j, overlap_size, num = (int)FragPairVec.size();

	for (i = 0; i < num; i++)
	{
		if (FragPairVec[i].qLen == 0) continue;

		iLen = FragPairVec[i].qLen;
		iqBeg = FragPairVec[i].qPos; iqEnd = iqBeg + iLen - 1;
		irBeg = FragPairVec[i].rPos; irEnd = irBeg + iLen - 1;

		for (j = i + 1; j < num; j++)
		{
			if (FragPairVec[j].qLen == 0) continue;

			jLen = FragPairVec[j].qLen;
			jqBeg = FragPairVec[j].qPos; jqEnd = jqBeg + jLen - 1;
			jrBeg = FragPairVec[j].rPos; jrEnd = jrBeg + jLen - 1;

			if (jqBeg <= iqEnd)
			{
				if ((overlap_size = iqEnd - jqBeg + 1) > jLen) overlap_size = jLen;
				//printf("overlap@query\n\ti[%d-%d, %lld-%lld]=%d\n\tj[%d-%d, %lld-%lld]=%d\n\tOverlap=%d\n", iqBeg, iqEnd, irBeg, irEnd, iLen, jqBeg, jqEnd, jrBeg, jrEnd, jLen, overlap_size);

				if (iLen <= jLen) // shrink block i
				{
					if ((FragPairVec[i].qLen -= overlap_size) < 10) FragPairVec[i].qLen = 0;
					FragPairVec[i].rLen = FragPairVec[i].qLen;
					//printf("update:i[%d-%d, %lld-%lld]=%d\n", iqBeg, iqBeg + FragPairVec[i].qLen - 1, irBeg, irBeg + FragPairVec[i].rLen - 1, FragPairVec[i].qLen);
				}
				else // shrink block j
				{
					FragPairVec[j].qPos += overlap_size; FragPairVec[j].rPos += overlap_size;
					if ((FragPairVec[j].qLen -= overlap_size) < 10) FragPairVec[j].qLen = 0;
					FragPairVec[j].rLen = FragPairVec[j].qLen;
					//printf("update:j[%d-%d, %lld-%lld]=%d\n", FragPairVec[j].qPos, FragPairVec[j].qPos + FragPairVec[j].qLen - 1, FragPairVec[j].rPos, FragPairVec[j].rPos + FragPairVec[j].rLen - 1, FragPairVec[j].qLen);
				}
			}
			else if (jrBeg <= irEnd)
			{
				if ((overlap_size = irEnd - jrBeg + 1) > jLen) overlap_size = jLen;
				//printf("overlap@ref\n\ti[%d-%d, %lld-%lld]=%d\n\tj[%d-%d, %lld-%lld]=%d\n\tOverlap=%d\n", iqBeg, iqEnd, irBeg, irEnd, iLen, jqBeg, jqEnd, jrBeg, jrEnd, jLen, overlap_size);
				if (iLen <= jLen) // shrink block i
				{
					if ((FragPairVec[i].qLen -= overlap_size) < 10) FragPairVec[i].qLen = 0;
					FragPairVec[i].rLen = FragPairVec[i].qLen;
					//printf("update:i[%d-%d, %lld-%lld]=%d\n", iqBeg, iqBeg + FragPairVec[i].qLen - 1, irBeg, irBeg + FragPairVec[i].rLen - 1, FragPairVec[i].qLen);
				}
				else // shrink block j
				{
					FragPairVec[j].qPos += overlap_size; FragPairVec[j].rPos += overlap_size;
					if ((FragPairVec[j].qLen -= overlap_size) < 10) FragPairVec[j].qLen = 0;
					FragPairVec[j].rLen = FragPairVec[j].qLen;
					//printf("update:j[%d-%d, %lld-%lld]=%d\n", FragPairVec[j].qPos, FragPairVec[j].qPos + FragPairVec[j].qLen - 1, FragPairVec[j].rPos, FragPairVec[j].rPos + FragPairVec[j].rLen - 1, FragPairVec[j].qLen);
				}
			}
			else break;
		}
	}
	for (bNullPair = false, i = 0; i < num; i++)
	{
		if (FragPairVec[i].qLen == 0)
		{
			bNullPair = true;
			break;
		}
	}
	return bNullPair;
}

void RemoveNullFragPairs(vector<FragPair_t>& FragPairVec)
{
	vector<FragPair_t> vec;

	vec.reserve((int)FragPairVec.size());
	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (iter->qLen > 0) vec.push_back(*iter);
	}
	FragPairVec.swap(vec);
}

void IdentifyNormalPairs(vector<FragPair_t>& FragPairVec)
{
	FragPair_t FragPair;
	int i, j, qGaps, rGaps, num = (int)FragPairVec.size();

	FragPair.bSeed = false;
	for (i = 0, j = 1; j < num; i++, j++)
	{
		if ((qGaps = FragPairVec[j].qPos - (FragPairVec[i].qPos + FragPairVec[i].qLen)) < 0) qGaps = 0;
		if ((rGaps = FragPairVec[j].rPos - (FragPairVec[i].rPos + FragPairVec[i].rLen)) < 0) rGaps = 0;

		if (qGaps > 0 || rGaps > 0)
		{
			FragPair.qPos = FragPairVec[i].qPos + FragPairVec[i].qLen;
			FragPair.rPos = FragPairVec[i].rPos + FragPairVec[i].rLen;
			FragPair.PosDiff = FragPair.rPos - FragPair.qPos;
			FragPair.qLen = qGaps; FragPair.rLen = rGaps;
			FragPairVec.push_back(FragPair);
			//if (qGaps > MinSeedLength || rGaps > MinSeedLength) ShowFragSeqs(FragPair);
			//printf("insert a normal pair: r[%d-%d] g[%lld-%lld] and r[%d-%d] g[%lld-%lld]: r[%d-%d] g[%lld-%lld]\n", FragPairVec[i].rPos, FragPairVec[i].rPos + FragPairVec[i].rLen - 1, FragPairVec[i].gPos, FragPairVec[i].gPos + FragPairVec[i].gLen - 1, FragPairVec[j].rPos, FragPairVec[j].rPos + FragPairVec[j].rLen - 1, FragPairVec[j].gPos, FragPairVec[j].gPos + FragPairVec[j].gLen - 1, FragPair.rPos, FragPair.rPos + FragPair.rLen - 1, FragPair.gPos, FragPair.gPos + FragPair.gLen - 1);
		}
	}
	if ((int)FragPairVec.size() > num) inplace_merge(FragPairVec.begin(), FragPairVec.begin() + num, FragPairVec.end(), CompByQueryPos);
}

void CheckFragPairContinuity(vector<FragPair_t>& FragPairVec)
{
	int64_t rPos;
	bool bChecked = true;
	int i, qPos, num = (int)FragPairVec.size();

	if (num == 0) return;

	qPos = FragPairVec[0].qPos + FragPairVec[0].qLen;
	rPos = FragPairVec[0].rPos + FragPairVec[0].rLen;

	for (i = 1; i < num; i++)
	{
		if (FragPairVec[i].qPos != qPos || FragPairVec[i].rPos != rPos)
		{
			printf("Stop!\n"); ShowFragSeqs(FragPairVec[i]);
			bChecked = false;
			break;
		}
		else
		{
			qPos = FragPairVec[i].qPos + FragPairVec[i].qLen;
			rPos = FragPairVec[i].rPos + FragPairVec[i].rLen;
		}
	}
}

int CheckFragPairMismatch(FragPair_t* FragPair)
{
	int i, mismatch = 0;
	char *TemplateSeq = RefSequence + FragPair->rPos;
	char *QuerySeq = (char*)QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos;

	for (i = 0; i < FragPair->qLen; i++)
	{
		if (QuerySeq[i] != TemplateSeq[i] && QuerySeq[i] != 'N' && TemplateSeq[i] != 'N') mismatch++;
	}
	return mismatch;
}

void *GenerateFragAlignment(void *arg)
{
	FragPair_t *FragPair;
	int i, j, AlnBlockNum, FragPairNum, TailIdx, mismatch, *my_id = (int*)arg;

	AlnBlockNum = (int)AlnBlockVec.size();

	for (i = 0; i < AlnBlockNum; i++)
	{
		FragPairNum = (int)AlnBlockVec[i].FragPairVec.size(); TailIdx = FragPairNum - 1;
		for (j = *my_id; j < FragPairNum; j+=iThreadNum)
		{
			if (AlnBlockVec[i].FragPairVec[j].bSeed) continue;

			FragPair = &AlnBlockVec[i].FragPairVec[j];
			if (FragPair->qLen == 0)
			{
				FragPair->aln1.resize(FragPair->rLen); strncpy((char*)FragPair->aln1.c_str(), RefSequence + FragPair->rPos, FragPair->rLen);
				FragPair->aln2.assign(FragPair->rLen, '-');
			}
			else if (FragPair->rLen == 0)
			{
				FragPair->aln1.assign(FragPair->qLen, '-');
				FragPair->aln2.resize(FragPair->qLen); strncpy((char*)FragPair->aln2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos, FragPair->qLen);
			}
			else if (FragPair->qLen == FragPair->rLen && (mismatch = CheckFragPairMismatch(FragPair)) <= 5)
			{
				FragPair->aln1.resize(FragPair->rLen); strncpy((char*)FragPair->aln1.c_str(), RefSequence + FragPair->rPos, FragPair->rLen);
				FragPair->aln2.resize(FragPair->qLen); strncpy((char*)FragPair->aln2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos, FragPair->qLen);
			}
			else
			{
				if (bDebugMode) printf("GenAln: %d vs %d\n", FragPair->qLen, FragPair->rLen), fflush(stdout);
				FragPair->aln1.resize(FragPair->rLen); strncpy((char*)FragPair->aln1.c_str(), RefSequence + FragPair->rPos, FragPair->rLen);
				FragPair->aln2.resize(FragPair->qLen); strncpy((char*)FragPair->aln2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos, FragPair->qLen);
				
				nw_alignment(FragPair->rLen, FragPair->aln1, FragPair->qLen, FragPair->aln2);
			}

			if (j == 0 && !AlnBlockVec[i].FragPairVec[j].bSeed)
			{
				int k, rGaps = 0, qGaps = 0;
				for (k = 0; k < (int)FragPair->aln1.length(); k++)
				{
					if (FragPair->aln1[k] == '-') rGaps++;
					else if (FragPair->aln2[k] == '-') qGaps++;
					else break;
				}
				if(k > 0)
				{
					FragPair->aln1 = FragPair->aln1.substr(k);
					FragPair->aln2 = FragPair->aln2.substr(k);
					FragPair->rLen -= rGaps; FragPair->rPos += rGaps;
					FragPair->qLen -= qGaps; FragPair->qPos += qGaps;
				}
			}
			else if (j == TailIdx && !AlnBlockVec[i].FragPairVec[j].bSeed)
			{
				int k, rGaps = 0, qGaps = 0;
				for (k = FragPair->aln1.length() - 1; k >= 0; k--)
				{
					if (FragPair->aln1[k] == '-') rGaps++;
					else if (FragPair->aln2[k] == '-') qGaps++;
					else break;
				}
				if (++k < FragPair->aln1.length())
				{
					FragPair->aln1.resize(k);
					FragPair->aln2.resize(k);
					FragPair->rLen -= rGaps;
					FragPair->qLen -= qGaps;
				}
			}
		}
	}
	return (void*)(1);
}

Coordinate_t GenCoordinateInfo(int64_t rPos)
{
	Coordinate_t coordinate;
	map<int64_t, int>::iterator iter;

	if (rPos < GenomeSize)
	{
		coordinate.bDir = true;

		iter = ChrLocMap.lower_bound(rPos);
		coordinate.ChromosomeIdx = iter->second;
		coordinate.gPos = rPos + 1 - ChromosomeVec[iter->second].FowardLocation;
	}
	else
	{
		coordinate.bDir = false;

		iter = ChrLocMap.lower_bound(rPos);
		coordinate.ChromosomeIdx = iter->second;
		coordinate.gPos = iter->first - rPos + 1;
	}
	return coordinate;
}

int CountBaseNum(string& frag)
{
	int n = (int)frag.length();

	for (string::iterator iter = frag.begin(); iter != frag.end(); iter++) if (*iter == '-') n--;

	return n;
}

int CountIdenticalPairs(string& aln1, string& aln2)
{
	int i, n, len = (int)aln1.length();

	for (n = len, i = 0; i < len; i++)
	{
		if (aln1[i] != aln2[i]) n--;
	}
	return n;
}

void OutputAlignment()
{
	char* frag;
	FILE *outFile;
	int64_t RefPos;
	int i, p, q, n, RefIdx, QueryPos, idy;
	vector<FragPair_t>::iterator FragPairIter;
	string QueryChrName, RefChrName, aln1, aln2, frag1, frag2;

	outFile = fopen(alnFileName, "a"); 
	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		aln1.clear(); aln2.clear(); idy = 0;
		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (FragPairIter->bSeed)
			{
				idy += FragPairIter->qLen; 
				frag = new char[FragPairIter->qLen + 1]; frag[FragPairIter->qLen] = '\0';
				strncpy(frag, QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos, FragPairIter->qLen);
				aln1 += frag; aln2 += frag; delete[] frag;
			}
			else
			{
				idy += CountIdenticalPairs(FragPairIter->aln1, FragPairIter->aln2);
				aln1 += FragPairIter->aln1;
				aln2 += FragPairIter->aln2;
			}
		}
		if ((n = (int)aln1.length()) > 200)
		{
			LocalAlignmentNum++;
			RefIdx = ABiter->coor.ChromosomeIdx;
			QueryChrName = QueryChrVec[QueryChrIdx].name; RefChrName = ChromosomeVec[RefIdx].name;
			if (QueryChrName.length() > RefChrName.length()) RefChrName += string().assign((QueryChrName.length() - RefChrName.length()), ' ');
			else QueryChrName += string().assign((RefChrName.length() - QueryChrName.length()), ' ');

			fprintf(outFile, "#Identity = %d / %d (%.2f%%) Orientation = %s\n\n", idy, n, (int)(10000 * (1.0*idy / n)) / 100.0, ABiter->coor.bDir? "Forward":"Reverse");
			//ShowFragPairVec(ABiter->FragPairVec); printf("\n\n");
			i = 0; QueryPos = ABiter->FragPairVec[0].qPos + 1; RefPos = ABiter->coor.gPos;
			while (i < n)
			{
				frag1 = aln1.substr(i, WindowSize); frag2 = aln2.substr(i, WindowSize);
				p = CountBaseNum(frag1); q = CountBaseNum(frag2);

				fprintf(outFile, "%s\t%12d\t%s\n%s\t%12d\t%s\n\n", RefChrName.c_str(), RefPos, frag1.c_str(), QueryChrName.c_str(), QueryPos, frag2.c_str());
				i += WindowSize; RefPos += (ABiter->coor.bDir ? p : 0 - p); QueryPos += q;
			}
			fprintf(outFile, "%s\n", string().assign(100, '*').c_str());
		}
	}
	std::fclose(outFile);
	//printf("Total alignment length = %lld\n", TotalAlnLen);
}

void OutputMAF()
{
	char* frag;
	FILE *outFile;
	int i, p, q, n, RefIdx, idy, ns;
	vector<FragPair_t>::iterator FragPairIter;
	string QueryChrName, RefChrName, aln1, aln2, frag1, frag2;

	outFile = fopen(mafFileName, "a");
	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		aln1.clear(); aln2.clear(); idy = 0;
		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (FragPairIter->bSeed)
			{
				idy += FragPairIter->qLen;
				frag = new char[FragPairIter->qLen + 1]; frag[FragPairIter->qLen] = '\0';
				strncpy(frag, QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos, FragPairIter->qLen);
				aln1 += frag; aln2 += frag; delete[] frag;
			}
			else
			{
				idy += CountIdenticalPairs(FragPairIter->aln1, FragPairIter->aln2);
				aln1 += FragPairIter->aln1;
				aln2 += FragPairIter->aln2;
			}
		}
		if ((n = (int)aln1.length()) > 200)
		{
			LocalAlignmentNum++;
			RefIdx = ABiter->coor.ChromosomeIdx;
			QueryChrName = QueryChrVec[QueryChrIdx].name; RefChrName = ChromosomeVec[RefIdx].name;
			if (QueryChrName.length() > RefChrName.length()) RefChrName += string().assign((QueryChrName.length() - RefChrName.length()), ' ');
			else QueryChrName += string().assign((RefChrName.length() - QueryChrName.length()), ' ');

			fprintf(outFile, "a score=%d\n", idy);
			if (ABiter->coor.bDir)
			{
				fprintf(outFile, "s %s %lld %lld + %d %s\n", ChromosomeVec[RefIdx].name, ABiter->coor.gPos, (long long)aln1.length(), ChromosomeVec[RefIdx].len, (char*)aln1.c_str());
				fprintf(outFile, "s %s %lld %lld + %d %s\n\n", (char*)QueryChrName.c_str(), ABiter->FragPairVec[0].qPos, (long long)aln2.length(), (int)QueryChrVec[QueryChrIdx].seq.length(), (char*)aln2.c_str());
			}
			else
			{
				i = (int)ABiter->FragPairVec.size() - 1;
				int64_t rPos = ABiter->FragPairVec[i].rPos + ABiter->FragPairVec[i].rLen - 1;
				SelfComplementarySeq((int)aln1.length(), (char*)aln1.c_str());
				SelfComplementarySeq((int)aln2.length(), (char*)aln2.c_str());
				fprintf(outFile, "s %s %lld %lld + %d %s\n", ChromosomeVec[RefIdx].name, GenCoordinateInfo(rPos).gPos, (long long)aln1.length(), ChromosomeVec[RefIdx].len, (char*)aln1.c_str());
				fprintf(outFile, "s %s %lld %lld - %d %s\n\n", (char*)QueryChrName.c_str(), ABiter->FragPairVec[0].qPos, (long long)aln2.length(), (int)QueryChrVec[QueryChrIdx].seq.length(), (char*)aln2.c_str());
			}
		}
	}
	std::fclose(outFile);
}

void OutputVariantCallingFile()
{
	int64_t rPos;
	FILE *outFile;
	string RefChrName, frag1, frag2;
	int i, j, qPos, aln_len, ind_len;
	vector<FragPair_t>::iterator FragPairIter;

	outFile = fopen(vcfFileName, "a");
	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (!ABiter->coor.bDir) continue;
		RefChrName = ChromosomeVec[ABiter->coor.ChromosomeIdx].name;

		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (!FragPairIter->bSeed)
			{
				if (FragPairIter->qLen == 0) // delete
				{
					frag1.resize(FragPairIter->rLen + 1); strncpy((char*)frag1.c_str(), RefSequence + FragPairIter->rPos - 1, FragPairIter->rLen + 1);
					fprintf(outFile, "%s\t%d\t.\t%s\t%c\t100\tPASS\tmt=DELETE\n", RefChrName.c_str(), GenCoordinateInfo(FragPairIter->rPos - 1).gPos, (char*)frag1.c_str(), QueryChrVec[QueryChrIdx].seq[FragPairIter->qPos - 1]);
				}
				else if (FragPairIter->rLen == 0) // insert
				{
					frag2.resize(FragPairIter->qLen + 1); strncpy((char*)frag2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos - 1, FragPairIter->qLen + 1);
					fprintf(outFile, "%s\t%d\t.\t%c\t%s\t100\tPASS\tmt=INSERT\n", RefChrName.c_str(), GenCoordinateInfo(FragPairIter->rPos - 1).gPos, RefSequence[FragPairIter->rPos - 1], (char*)frag2.c_str());
				}
				else if (FragPairIter->qLen == 1 && FragPairIter->rLen == 1) // substitution
				{
					fprintf(outFile, "%s\t%d\t.\t%c\t%c\t100\tPASS\tmt=SUBSTITUTE\n", RefChrName.c_str(), GenCoordinateInfo(FragPairIter->rPos).gPos, FragPairIter->aln2[0], FragPairIter->aln1[0]);
				}
				else
				{
					//fprintf(stdout, "ref=%s\nqry=%s\n", FragPairIter->aln2.c_str(), FragPairIter->aln1.c_str());
					for (aln_len = (int)FragPairIter->aln1.length(), rPos = FragPairIter->rPos, qPos = FragPairIter->qPos, i = 0; i < aln_len; i++)
					{
						if (FragPairIter->aln1[i] == FragPairIter->aln2[i])
						{
							rPos++; qPos++;
						}
						else if (FragPairIter->aln1[i] == '-') // insert
						{
							ind_len = 1; while (FragPairIter->aln1[i + ind_len] == '-') ind_len++;
							frag2 = QueryChrVec[QueryChrIdx].seq.substr(qPos - 1, ind_len + 1);
							fprintf(outFile, "%s\t%d\t.\t%c\t%s\t100\tPASS\tmt=INSERT\n", RefChrName.c_str(), GenCoordinateInfo(rPos - 1).gPos, frag2[0], (char*)frag2.c_str());
							qPos += ind_len; i += ind_len - 1;
						}
						else if (FragPairIter->aln2[i] == '-') // delete
						{
							ind_len = 1; while (FragPairIter->aln2[i + ind_len] == '-') ind_len++;
							frag1.resize(ind_len + 2); frag1[ind_len + 1] = '\0';
							strncpy((char*)frag1.c_str(), RefSequence + rPos - 1, ind_len + 1);
							fprintf(outFile, "%s\t%d\t.\t%s\t%c\t100\tPASS\tmt=DELETE\n", RefChrName.c_str(), GenCoordinateInfo(rPos - 1).gPos, (char*)frag1.c_str(), frag1[0]);
							rPos += ind_len; i += ind_len - 1;
						}
						else // substitute
						{
							fprintf(outFile, "%s\t%d\t.\t%c\t%c\t100\tPASS\tmt=SUBSTITUTE\n", RefChrName.c_str(), GenCoordinateInfo(rPos).gPos, FragPairIter->aln1[i], FragPairIter->aln2[i]);
							rPos++; qPos++;
						}
					}
				}
			}
		}
	}
	std::fclose(outFile);
}

void OutputSNPs()
{
	char* frag;
	FILE *outFile;
	int64_t RefPos;
	int i, snpPos, RefIdx, QueryPos;
	vector<FragPair_t>::iterator FragPairIter;
	string QueryChrName, RefChrName, frag1, frag2, aln;

	outFile = fopen(snpFileName, "a");
	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		RefIdx = ABiter->coor.ChromosomeIdx;
		//QueryPos = ABiter->FragPairVec[0].qPos; RefPos = ABiter->coor.gPos;

		QueryChrName = QueryChrVec[QueryChrIdx].name; RefChrName = ChromosomeVec[RefIdx].name;

		//if (QueryChrName.length() > RefChrName.length()) RefChrName += string().assign((QueryChrName.length() - RefChrName.length()), ' ');
		//else QueryChrName += string().assign((RefChrName.length() - QueryChrName.length()), ' ');

		frag1.resize(21); frag2.resize(21);
		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (!FragPairIter->bSeed && FragPairIter->qLen == FragPairIter->rLen)
			{
				if (FragPairIter->qLen == 1)
				{
					QueryPos = FragPairIter->qPos - 10; RefPos = FragPairIter->rPos - 10;
					strncpy((char*)frag1.c_str(), RefSequence + RefPos, 21);
					strncpy((char*)frag2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + QueryPos, 21);
					fprintf(outFile, "SNP#%d\t%s%s:%d vs %s:%d\n", ++SNP_num, RefChrName.c_str(), (ABiter->coor.bDir ? "" : "_rev"), GenCoordinateInfo(RefPos).gPos, QueryChrName.c_str(), QueryPos + 1);
					fprintf(outFile, "R: %s\nQ: %s\n             ^             \n\n", (char*)frag1.c_str(), (char*)frag2.c_str());
				}
				else if(FragPairIter->qLen == (int)FragPairIter->aln2.size())
				{
					aln.assign(FragPairIter->qLen, ' '); snpPos = -1;
					for (i = 0; i < FragPairIter->qLen; i++)
					{
						if (FragPairIter->aln1[i] != FragPairIter->aln2[i])
						{
							aln[i] = '^';
							if (snpPos == -1) snpPos = i;
						}
					}

					QueryPos = FragPairIter->qPos - 10; RefPos = FragPairIter->rPos - 10;
					fprintf(outFile, "SNP#%d\t%s%s:%d vs %s:%d\n", ++SNP_num, RefChrName.c_str(), (ABiter->coor.bDir ? "" : "_rev"), GenCoordinateInfo(RefPos).gPos, QueryChrName.c_str(), QueryPos + 1);

					strncpy((char*)frag1.c_str(), RefSequence + RefPos, 10); frag1[10] = '\0';
					strncpy((char*)frag2.c_str(), RefSequence + FragPairIter->rPos + FragPairIter->rLen, 10); frag2[10] = '\0';
					fprintf(outFile, "R: %s%s%s\n", (char*)frag1.c_str(), (char*)FragPairIter->aln2.c_str(), (char*)frag2.c_str());

					strncpy((char*)frag1.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + QueryPos, 10); frag1[10] = '\0';
					strncpy((char*)frag2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos + FragPairIter->qLen, 10); frag2[10] = '\0';
					fprintf(outFile, "Q: %s%s%s\n", (char*)frag1.c_str(), (char*)FragPairIter->aln1.c_str(), (char*)frag2.c_str());

					fprintf(outFile, "             %s          \n\n", (char*)aln.c_str());
				}
			}
		}
	}
	std::fclose(outFile);
}

void OutputIndeles()
{
	char* frag;
	FILE *outFile;
	int64_t RefPos, TotalAlnLen = 0;
	int i, p, q, n, RefIdx, QueryPos;
	vector<FragPair_t>::iterator FragPairIter;
	string QueryChrName, RefChrName, head, tail, frag1, frag2;

	outFile = fopen(indFileName, "a"); head.resize(10); tail.resize(10);
	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		RefIdx = ABiter->coor.ChromosomeIdx;
		QueryPos = ABiter->FragPairVec[0].qPos; RefPos = ABiter->coor.gPos;

		QueryChrName = QueryChrVec[QueryChrIdx].name; RefChrName = ChromosomeVec[RefIdx].name;
		if (QueryChrName.length() > RefChrName.length()) RefChrName += string().assign((QueryChrName.length() - RefChrName.length()), ' ');
		else QueryChrName += string().assign((RefChrName.length() - QueryChrName.length()), ' ');

		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (!FragPairIter->bSeed)
			{
				if (FragPairIter->qLen == 0 || FragPairIter->rLen == 0 || FragPairIter->qLen != FragPairIter->rLen)
				{
					QueryPos = FragPairIter->qPos - 10; RefPos = FragPairIter->rPos - 10;

					strncpy((char*)head.c_str(), RefSequence + RefPos, 10);
					strncpy((char*)tail.c_str(), RefSequence + FragPairIter->rPos + FragPairIter->rLen, 10);
					frag1 = head + FragPairIter->aln2 + tail;

					strncpy((char*)head.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + QueryPos, 10);
					strncpy((char*)tail.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos + FragPairIter->qLen, 10);
					frag2 = head + FragPairIter->aln1 + tail;

					fprintf(outFile, "Ind#%d\t%s%s:%d vs %s:%d\n", ++IND_num, RefChrName.c_str(), (ABiter->coor.bDir ? "" : "_rev"), GenCoordinateInfo(RefPos).gPos, QueryChrName.c_str(), QueryPos + 1);
					fprintf(outFile, "R: %s\nQ: %s\n\n", (char*)frag1.c_str(), (char*)frag2.c_str());
				}
			}
		}
	}
	std::fclose(outFile);
}

void OutputDotplot()
{
	int64_t last_ref_end;
	vector<int> ChrScoreVec;
	vector<AlnBlock_t>::iterator ABiter;
	FILE *dataFile1, *dataFile2, *outFile;
	string QueryChrName, RefChrName, cmd, DataFileName;
	vector<vector<vector<AlnBlock_t>::iterator> > ChrClusterVec;
	int i, j, last_query_end, FragNum, num, iCluster, ChrIdx, thr;

	ChrScoreVec.resize(iChromsomeNum); ChrClusterVec.resize(iChromsomeNum);
	for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) ChrScoreVec[ABiter->coor.ChromosomeIdx] += ABiter->score;
	for (thr = i = 0; i < iChromsomeNum; i++) if (thr < ChrScoreVec[i]) thr = ChrScoreVec[i]; thr /= 4;
	for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (ABiter->score >= 100 && ChrScoreVec[ABiter->coor.ChromosomeIdx] > thr)
			ChrClusterVec[ABiter->coor.ChromosomeIdx].push_back(ABiter);
	}
	for (iCluster = i = 0; i < iChromsomeNum; i++) if (ChrScoreVec[i] > thr) iCluster++;

	outFile = fopen(gpFileName, "w"); QueryChrName = QueryChrVec[QueryChrIdx].name;
	fprintf(outFile, "set terminal postscript color solid 'Courier' 8\nset output '%s-%s.ps'\n", OutputPrefix, QueryChrVec[QueryChrIdx].name.c_str());
	fprintf(outFile, "set multiplot layout %d,1\nset grid\nset border 1\n", iCluster);
	fprintf(outFile, "set style line 1 lw 4 pt 0 ps 0.5 lc rgb 'red'\nset style line 2 lw 4 pt 0 ps 0.5 lc rgb 'blue'\n");

	for (i = 0; i < iChromsomeNum; i++)
	{
		if ((num = (int)ChrClusterVec[i].size()) > 0)
		{
			DataFileName = (string)OutputPrefix + "." + QueryChrName + "-" + ChromosomeVec[i].name;
			fprintf(outFile, "set xrange[1:*]\nset yrange[1:*]\nset xlabel 'Query (%s)'\nset ylabel 'Ref (%s)'\n", (char*)QueryChrName.c_str(), ChromosomeVec[i].name);
			fprintf(outFile, "plot '%s' title 'Forward' with lp ls 1, '%s' title 'Reverse' with lp ls 2\n\n", (char*)(DataFileName + ".1").c_str(), (char*)(DataFileName + ".2").c_str());

			dataFile1 = fopen((char*)(DataFileName + ".1").c_str(), "w"); dataFile2 = fopen((char*)(DataFileName + ".2").c_str(), "w");
			fprintf(dataFile1, "0 0\n0 0\n\n"); fprintf(dataFile2, "0 0\n0 0\n\n");

			for (j = 0; j < num; j++)
			{
				FragNum = (int)ChrClusterVec[i][j]->FragPairVec.size() - 1;
				last_query_end = ChrClusterVec[i][j]->FragPairVec[FragNum].qPos + ChrClusterVec[i][j]->FragPairVec[FragNum].qLen - 1;
				last_ref_end = ChrClusterVec[i][j]->FragPairVec[FragNum].rPos + ChrClusterVec[i][j]->FragPairVec[FragNum].rLen - 1;
				if (ChrClusterVec[i][j]->coor.bDir) fprintf(dataFile1, "%d %d\n%d %d\n\n", ChrClusterVec[i][j]->FragPairVec[0].qPos + 1, GenCoordinateInfo(ChrClusterVec[i][j]->FragPairVec[0].rPos).gPos, last_query_end + 1, GenCoordinateInfo(last_ref_end).gPos);
				else fprintf(dataFile2, "%d %d\n%d %d\n\n", ChrClusterVec[i][j]->FragPairVec[0].qPos + 1, GenCoordinateInfo(ChrClusterVec[i][j]->FragPairVec[0].rPos).gPos, last_query_end + 1, GenCoordinateInfo(last_ref_end).gPos);
			}
			std::fclose(dataFile1); std::fclose(dataFile2);
		}
	}
	fprintf(outFile, "#unset multiplot\n"); std::fclose(outFile);
	cmd = (string)GnuPlotPath + " " + (string)gpFileName; system((char*)cmd.c_str());
	DataFileName = (string)OutputPrefix + "." + QueryChrName + "-*.*";
	cmd = "rm " + DataFileName; system(cmd.c_str());
}

void GenomeComparison()
{
	int i, n;
	bool* CoverageArr;
	vector<AlnBlock_t>::iterator ABiter;
	int64_t iTotalLength = 0, iCoverage = 0;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];
	
	vector<int> vec(iThreadNum); for (i = 0; i < iThreadNum; i++) vec[i] = i;

	if (bDebugMode) iThreadNum = 1;

	fprintf(stderr, "Step2. Sequence analysis for all query chromosomes\n");
	for (QueryChrIdx = 0; QueryChrIdx != iQueryChrNum; QueryChrIdx++)
	{
		fprintf(stderr, "\tProcess query chromsomoe: %s...\n", QueryChrVec[QueryChrIdx].name.c_str());
		CoverageArr = new bool[QueryChrVec[QueryChrIdx].seq.length()]();

		SeedVec.clear(); AlnBlockVec.clear();
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, IdentifyLocalMEM, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

		//for (vector<FragPair_t>::iterator iter = SeedVec.begin(); iter != SeedVec.end(); iter++)
		//	printf("q[%d-%d] r[%lld-%lld] len=%d PD=%lld\n", iter->qPos, iter->qPos + iter->qLen - 1, iter->rPos, iter->rPos + iter->rLen - 1, iter->qLen, iter->PosDiff);

		AlignmentBlockClustering();
		RemoveRedundantAlnBlocks();

		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
		{
			if (RemoveOverlaps(ABiter->FragPairVec)) RemoveNullFragPairs(ABiter->FragPairVec);
		}
		RemoveSmallAlnBlocks();

		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) IdentifyNormalPairs(ABiter->FragPairVec);
		//for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) ShowFragPairVec(ABiter->FragPairVec);
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, GenerateFragAlignment, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
		{
			n = (int)ABiter->FragPairVec.size() - 1;
			memset(CoverageArr + ABiter->FragPairVec[0].qPos, true, ABiter->FragPairVec[n].qPos + ABiter->FragPairVec[n].qLen - ABiter->FragPairVec[0].qPos);
			ABiter->coor = GenCoordinateInfo(ABiter->FragPairVec[0].rPos);
		}
		//fprintf(stderr, "\tOutput the alignment for query chromosome %s in the file: %s...\n", QueryChrVec[QueryChrIdx].name.c_str(), alnFileName);	OutputAlignment();
		fprintf(stderr, "\tOutput the MAF for query chromosome %s in the file: %s...\n", QueryChrVec[QueryChrIdx].name.c_str(), mafFileName);	OutputMAF();
		fprintf(stderr, "\tOutput the variants for query chromosome %s in the file: %s...\n", QueryChrVec[QueryChrIdx].name.c_str(), vcfFileName); OutputVariantCallingFile();
		//if (bShowSubstitution) fprintf(stderr, "\tOutput the SNPs for query chromosome %s in the file: %s...\n", QueryChrVec[QueryChrIdx].name.c_str(), snpFileName), OutputSNPs();
		//if (bShowIndel) fprintf(stderr, "\tOutput the indels for query chromosome %s in the file: %s...\n", QueryChrVec[QueryChrIdx].name.c_str(), indFileName), OutputIndeles();
		if (GnuPlotPath != NULL) fprintf(stderr, "\tGenerate the dotplot for query chromosome %s in the file: %s-%s.ps...\n", QueryChrVec[QueryChrIdx].name.c_str(), OutputPrefix, QueryChrVec[QueryChrIdx].name.c_str()), OutputDotplot();

		iTotalLength += (n = (int)QueryChrVec[QueryChrIdx].seq.length());
		for (i = 0; i < n; i++) if (CoverageArr[i]) iCoverage++;
		delete[] CoverageArr;
	}
	fprintf(stderr, "\nAlignment # = %lld, coverage = %.4f\n", (long long)LocalAlignmentNum, 1.0*iCoverage / iTotalLength);

	delete[] ThreadArr;
}
