#include "structure.h"

static pthread_mutex_t Lock;

bool CompByPosDiff(const FragPair_t& p1, const FragPair_t& p2)
{
	if (p1.PosDiff == p2.PosDiff) return p1.qPos < p2.qPos;
	else return p1.PosDiff < p2.PosDiff;
}

bool CompByQueryPos(const FragPair_t& p1, const FragPair_t& p2)
{
	if (p1.qPos == p2.qPos) return p1.rPos < p2.rPos;
	else return p1.qPos < p2.qPos;
}

bool CompByRemoval(const FragPair_t& p1, const FragPair_t& p2)
{
	if (p1.bSeed && p2.bSeed) return p1.qPos < p2.qPos;
	else return (p1.bSeed > p2.bSeed);
}

bool CompByAlnBlockScore(const AlnBlock_t& p1, const AlnBlock_t& p2)
{
	return p1.score > p2.score;
}

int CalAlnBlockScore(vector<FragPair_t>& FragPairVec)
{
	if (FragPairVec.size() == 0) return 0;
	else if ((FragPairVec.rbegin()->qPos + FragPairVec.rbegin()->qLen - FragPairVec.begin()->qPos) < MinAlnLength) return 0;
	else
	{
		int score = 0;
		for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++) score += iter->qLen;
		return score;
	}
}

int CountIdenticalPairs(string& aln1, string& aln2)
{
	int i, n, len = (int)aln1.length();

	for (n = 0, i = 0; i < len; i++)
	{
		if (nst_nt4_table[(int)aln1[i]] == nst_nt4_table[(int)aln2[i]]) n++;
	}
	return n;
}

int CheckFragPairMismatch(FragPair_t* FragPair)
{
	int i, mismatch = 0;
	char *TemplateSeq = RefSequence + FragPair->rPos;
	char *QuerySeq = (char*)QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos;

	for (i = 0; i < FragPair->qLen; i++)
	{
		if (nst_nt4_table[QuerySeq[i]] != nst_nt4_table[TemplateSeq[i]]) mismatch++;
	}
	return mismatch;
}

void RemoveBadSeeds(vector<FragPair_t>& FragPairVec)
{
	int num = (int)FragPairVec.size();
	sort(FragPairVec.begin(), FragPairVec.end(), CompByRemoval);
	while (!FragPairVec[num - 1].bSeed) num--;
	//printf("size: %d --> %d\n", (int)FragPairVec.size(), num);
	FragPairVec.resize(num);
}

void RemoveBadAlnBlocks()
{
	int num;
	num = AlnBlockVec.size(); sort(AlnBlockVec.begin(), AlnBlockVec.end(), CompByAlnBlockScore);
	while (AlnBlockVec[num -1].score == 0) num--;
	//fprintf(stderr, "\t\tAlnBlockNum = %d --> %d\n", AlnBlockVec.size(), num);
	AlnBlockVec.resize(num);
}

void CheckGapsBetweenSeeds(AlnBlock_t& AlnBlock)
{
	bool bSim;
	AlnBlock_t SubAlnBlock;
	vector<int> BreakPointVec;
	vector<int>::iterator iter;
	int p, i, j, num, qGap, rGap;

	num = (int)AlnBlock.FragPairVec.size();
	for (i = 0, j = 1; j < num; i++, j++)
	{
		qGap = AlnBlock.FragPairVec[j].qPos - AlnBlock.FragPairVec[i].qPos - AlnBlock.FragPairVec[i].qLen;
		rGap = AlnBlock.FragPairVec[j].rPos - AlnBlock.FragPairVec[i].rPos - AlnBlock.FragPairVec[i].rLen;
		if ((qGap > 300 || rGap > 300))
		{
			if (abs(qGap - rGap) > 100 || qGap > MaxSeedGap || rGap > MaxSeedGap || CalGapSimilarity(AlnBlock.FragPairVec[i].qPos + AlnBlock.FragPairVec[i].qLen, AlnBlock.FragPairVec[j].qPos, AlnBlock.FragPairVec[i].rPos + AlnBlock.FragPairVec[i].rLen, AlnBlock.FragPairVec[j].rPos) == false) 
				BreakPointVec.push_back(j);
			//printf("gap_size = %d\n", gap); //ShowFragPair(AlnBlock.FragPairVec[i - 1]); ShowFragPair(AlnBlock.FragPairVec[i]);
		}
	}
	if ((num = (int)BreakPointVec.size()) > 0) // split this alignment block
	{
		pthread_mutex_lock(&Lock);
		AlnBlock.score = 0;
		for (i = 0, iter = BreakPointVec.begin(); iter != BreakPointVec.end();iter++)
		{
			j = *iter;
			//printf("split [%d-%d]\n", i, j); ShowFragPair(AlnBlock.FragPairVec[j - 1]); ShowFragPair(AlnBlock.FragPairVec[j]);
			SubAlnBlock.FragPairVec.clear(); copy(AlnBlock.FragPairVec.begin() + i, AlnBlock.FragPairVec.begin() + j, back_inserter(SubAlnBlock.FragPairVec));
			if ((SubAlnBlock.score = CalAlnBlockScore(AlnBlock.FragPairVec)) > MinAlnBlockScore) AlnBlockVec.push_back(SubAlnBlock);
			i = j;
		}
		//printf("split [%d-%d]\n", i, (int)AlnBlock.FragPairVec.size()); ShowFragPair(AlnBlock.FragPairVec[i - 1]); ShowFragPair(AlnBlock.FragPairVec[i]);
		SubAlnBlock.FragPairVec.clear(); copy(AlnBlock.FragPairVec.begin() + i, AlnBlock.FragPairVec.end(), back_inserter(SubAlnBlock.FragPairVec));
		if ((SubAlnBlock.score = CalAlnBlockScore(AlnBlock.FragPairVec)) > MinAlnBlockScore) AlnBlockVec.push_back(SubAlnBlock);
		pthread_mutex_unlock(&Lock);
	}
}

void *CheckAlnBlockLargeGaps(void *arg)
{
	int i, *my_id = (int*)arg;
	for (i = *my_id; i < AlnBlockNum; i += iThreadNum) CheckGapsBetweenSeeds(AlnBlockVec[i]);

	return (void*)(1);
}

//void CheckOverlaps(vector<FragPair_t>& FragPairVec)
//{
//	int i, j, num;
//	bool bOverlap;
//
//	num = (int)FragPairVec.size();
//	for (i = 0, j = 1; j < num; i++, j++)
//	{
//		bOverlap = false;
//		if (FragPairVec[j].qPos >= FragPairVec[i].qPos && FragPairVec[j].qPos < (FragPairVec[i].qPos + FragPairVec[i].qLen)) bOverlap = true;
//		if (FragPairVec[j].rPos >= FragPairVec[i].rPos && FragPairVec[j].rPos < (FragPairVec[i].rPos + FragPairVec[i].rLen)) bOverlap = true;
//		if (bOverlap) printf("Overlap!\n"), ShowFragPair(FragPairVec[i]), ShowFragPair(FragPairVec[j]);
//	}
//}

void RemoveOverlaps(vector<FragPair_t>& FragPairVec)
{
	bool bModified;
	int i, j, overlap_size, num;

	while (true)
	{
		bModified = false; num = (int)FragPairVec.size();
		for (i = 0, j = 1; j < num; i++, j++)
		{
			if (FragPairVec[j].rPos <= FragPairVec[i].rPos)
			{
				//printf("case1\n"); ShowFragPair(FragPairVec[i]); ShowFragPair(FragPairVec[j]);
				bModified = true; FragPairVec[i].bSeed = false;
				continue;
			}
			if ((overlap_size = FragPairVec[i].rPos + FragPairVec[i].rLen - FragPairVec[j].rPos) > 0)
			{
				//printf("case2\n"); ShowFragPair(FragPairVec[i]); ShowFragPair(FragPairVec[j]);
				FragPairVec[i].qLen -= overlap_size; FragPairVec[i].rLen -= overlap_size;
				if (FragPairVec[i].qLen <= 0 || FragPairVec[i].rLen <= 0)
				{
					bModified = true; FragPairVec[i].bSeed = false;
					continue;
				}
				//printf("update:\n"); ShowFragPair(FragPairVec[i]); ShowFragPair(FragPairVec[j]);
			}
			if ((overlap_size = FragPairVec[i].qPos + FragPairVec[i].qLen - FragPairVec[j].qPos) > 0)
			{
				//printf("case3\n"); ShowFragPair(FragPairVec[i]); ShowFragPair(FragPairVec[j]);
				FragPairVec[i].qLen -= overlap_size; FragPairVec[i].rLen -= overlap_size;
				if (FragPairVec[i].qLen <= 0 || FragPairVec[i].rLen <= 0)
				{
					bModified = true; FragPairVec[i].bSeed = false;
					continue;
				}
				//printf("update:\n"); ShowFragPair(FragPairVec[i]); ShowFragPair(FragPairVec[j]);
			}
		}
		if (bModified) RemoveBadSeeds(FragPairVec);
		else break;
	}
}
void *CheckAlnBlockOverlaps(void *arg)
{
	int i, *my_id = (int*)arg;

	for (i = *my_id; i < AlnBlockNum; i += iThreadNum) RemoveOverlaps(AlnBlockVec[i].FragPairVec);

	return (void*)(1);
}

void IdentifyNormalPairs(vector<FragPair_t>& FragPairVec)
{
	bool bBreak = false;
	FragPair_t FragPair;
	int i, j, qGaps, rGaps, num;

	if ((num = (int)FragPairVec.size()) == 1) return;

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
			//if(qGaps > MaxSeedGap || rGaps > MaxSeedGap) printf("Gaps=%d, %d\n", qGaps, rGaps), ShowFragPair(FragPairVec[i]), ShowFragPair(FragPairVec[j]);
		}
	}
	if ((int)FragPairVec.size() > num) inplace_merge(FragPairVec.begin(), FragPairVec.begin() + num, FragPairVec.end(), CompByQueryPos);
}

void *FillAlnBlockGaps(void *arg)
{
	int i, *my_id = (int*)arg;

	for (i = *my_id; i < AlnBlockNum; i += iThreadNum)
	{
		if (AlnBlockVec[i].score > 0) IdentifyNormalPairs(AlnBlockVec[i].FragPairVec);
	}
	return (void*)(1);
}

void *GenerateFragAlignment(void *arg)
{
	FragPair_t *FragPair;
	int i, j, aln_len, score, FragPairNum, mismatch, *my_id = (int*)arg;

	for (i = 0; i < AlnBlockNum; i++)
	{
		aln_len = score = 0;
		FragPairNum = (int)AlnBlockVec[i].FragPairVec.size();
		//for (j = 0; j < FragPairNum;j++)
		for (j = *my_id; j < FragPairNum; j+= iThreadNum)
		{
			if (AlnBlockVec[i].FragPairVec[j].bSeed)
			{
				aln_len += AlnBlockVec[i].FragPairVec[j].qLen;
				score += AlnBlockVec[i].FragPairVec[j].qLen;
			}
			else
			{
				FragPair = &AlnBlockVec[i].FragPairVec[j];
				if (FragPair->qLen == 0)
				{
					aln_len += FragPair->rLen;
					FragPair->aln1.resize(FragPair->rLen); strncpy((char*)FragPair->aln1.c_str(), RefSequence + FragPair->rPos, FragPair->rLen);
					FragPair->aln2.assign(FragPair->rLen, '-');
				}
				else if (FragPair->rLen == 0)
				{
					aln_len += FragPair->qLen;
					FragPair->aln1.assign(FragPair->qLen, '-');
					FragPair->aln2.resize(FragPair->qLen); strncpy((char*)FragPair->aln2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos, FragPair->qLen);
				}
				//else if (FragPair->qLen == FragPair->rLen && ((mismatch = CheckFragPairMismatch(FragPair)) <= 5 || 1.0* mismatch / FragPair->qLen < 0.1))
				else if (FragPair->qLen == FragPair->rLen && (mismatch = CheckFragPairMismatch(FragPair)) <= 5)
				{
					FragPair->aln1.resize(FragPair->rLen); strncpy((char*)FragPair->aln1.c_str(), RefSequence + FragPair->rPos, FragPair->rLen);
					FragPair->aln2.resize(FragPair->qLen); strncpy((char*)FragPair->aln2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos, FragPair->qLen);
					aln_len += FragPair->qLen;
					score += (FragPair->qLen - mismatch);
					//if (FragPair->qPos == 91210692) printf("case1, mismatch=%d\n", mismatch);
				}
				else
				{
					//if (bDebugMode) printf("GenAln: %d vs %d\n", FragPair->qLen, FragPair->rLen), fflush(stdout);
					//ShowFragPair(*FragPair); fflush(stdout);
					FragPair->aln1.resize(FragPair->rLen); strncpy((char*)FragPair->aln1.c_str(), RefSequence + FragPair->rPos, FragPair->rLen);
					FragPair->aln2.resize(FragPair->qLen); strncpy((char*)FragPair->aln2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos, FragPair->qLen);
					ksw2_alignment(FragPair->rLen, FragPair->aln1, FragPair->qLen, FragPair->aln2);
					//if (FragPair->qPos == 1231410) ShowFragPair(*FragPair), printf("%s\n%s\n", FragPair->aln1.c_str(), FragPair->aln2.c_str());
					aln_len += FragPair->aln1.length();
					score += CountIdenticalPairs(FragPair->aln1, FragPair->aln2);
					//if (FragPair->qPos == 91210692) printf("case2\n");
				}
			}
		}
		pthread_mutex_lock(&Lock);
		AlnBlockVec[i].aln_len += aln_len;
		AlnBlockVec[i].score += score;
		pthread_mutex_unlock(&Lock);
	}
	return (void*)(1);
}

void CheckAlnBlockCompleteness(vector<FragPair_t>& FragPairVec)
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
			printf("Gaps\n");
			ShowFragPair(FragPairVec[i-1]); ShowFragPair(FragPairVec[i]);
			break;
		}
		else
		{
			qPos = FragPairVec[i].qPos + FragPairVec[i].qLen;
			rPos = FragPairVec[i].rPos + FragPairVec[i].rLen;
		}
	}
}
