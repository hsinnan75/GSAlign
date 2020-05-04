#include <queue>
#include "structure.h"

#define Min_PD_Freq 3
#define SeedExplorationChunk 10000

int ObrPos = -1;
pthread_mutex_t Lock;
int64_t *RefChrScoreArr;
vector<FragPair_t> SeedVec;
vector<AlnBlock_t> AlnBlockVec;
uint32_t QrySeqPos, QryChrLength;
vector<pair<int, int> > SeedGroupVec;
int SeedNum, SeedGroupNum, GroupID, AlnBlockNum, DupAlnNum = 0;
int64_t TotalAlignmentLength = 0, TotalAlignmentMatches = 0, LocalAlignmentNum = 0, SNP_num = 0, IND_num = 0, SVS_num = 0;

bool CompByAlnBlockQueryPos(const AlnBlock_t& p1, const AlnBlock_t& p2)
{
	if (p1.FragPairVec.begin()->qPos == p2.FragPairVec.begin()->qPos) return p1.score > p2.score;
	else return p1.FragPairVec.begin()->qPos < p2.FragPairVec.begin()->qPos;
}

bool CompByAlnBlockRefPos(const AlnBlock_t& p1, const AlnBlock_t& p2)
{
	if (p1.FragPairVec.begin()->rPos == p2.FragPairVec.begin()->rPos) return p1.score > p2.score;
	else return p1.FragPairVec.begin()->rPos < p2.FragPairVec.begin()->rPos;
}

void AddAlnBlock(int i, int j)
{
	int region_size;
	AlnBlock_t AlnBlock;

	AlnBlock.score = 0; AlnBlock.bDup = false;
	copy(SeedVec.begin() + i, SeedVec.begin() + j, back_inserter(AlnBlock.FragPairVec));
	for (vector<FragPair_t>::iterator iter = AlnBlock.FragPairVec.begin(); iter != AlnBlock.FragPairVec.end(); iter++) AlnBlock.score += iter->qLen;
	region_size = (AlnBlock.FragPairVec.rbegin()->qPos + AlnBlock.FragPairVec.rbegin()->qLen) - AlnBlock.FragPairVec.begin()->qPos;
	if (AlnBlock.score < MinAlnBlockScore || region_size < MinAlnLength || (AlnBlock.score < 1000 && AlnBlock.score < region_size*0.05))
	{
		//printf("discard!\n"); ShowAlnBlockBoundary(AlnBlock.score, AlnBlock.FragPairVec);
	}
	else
	{
		//printf("add\n"); ShowAlnBlockBoundary(AlnBlock.score, AlnBlock.FragPairVec);
		pthread_mutex_lock(&Lock);
		AlnBlockVec.push_back(AlnBlock);
		pthread_mutex_unlock(&Lock);
	}
}

void *IdentifyLocalMEM(void *arg)
{
	int i;
	FragPair_t seed;
	uint32_t start, stop;
	vector<FragPair_t> seed_vec;
	bwtSearchResult_t bwtSearchResult;

	string& seq = QueryChrVec[QueryChrIdx].seq; seed.bSeed = true;

	while (true)
	{
		pthread_mutex_lock(&Lock);
		start = QrySeqPos; stop = (QrySeqPos += SeedExplorationChunk);
		if (stop > QryChrLength) stop = QryChrLength;
		fprintf(stderr, "\r\t\tSeed exploration: %d / %d (%d%%)...", stop, QryChrLength, (int)(100 * ((1.0*stop / QryChrLength))));
		pthread_mutex_unlock(&Lock);

		if (start >= QryChrLength) break;
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
						seed_vec.push_back(seed);
					}
					delete[] bwtSearchResult.LocArr;
					if (bSensitive) start += 5;
					//else start += bwtSearchResult.len;
					else start += (bwtSearchResult.len + 1);
				}
				else start++;
			}
		}
	}
	sort(seed_vec.begin(), seed_vec.end(), CompByPosDiff);

	if (iThreadNum == 1) SeedVec.swap(seed_vec);
	else if (seed_vec.size() > 0)
	{
		pthread_mutex_lock(&Lock);
		int num = (int)SeedVec.size();
		copy(seed_vec.begin(), seed_vec.end(), back_inserter(SeedVec));
		inplace_merge(SeedVec.begin(), SeedVec.begin() + num, SeedVec.end(), CompByPosDiff);
		pthread_mutex_unlock(&Lock);
	}
	return (void*)(1);
}

void FindSpecificLocalMEM(int start, int stop, int64_t rPos1, int64_t rPos2)
{
	FragPair_t seed;

	string& seq = QueryChrVec[QueryChrIdx].seq; seed.bSeed = true;
	while (start < stop)
	{
		if (nst_nt4_table[(int)seq[start]] > 3) start++;
		else
		{
			seed = Specific_BWT_Search(seq, start, stop, rPos1, rPos2);
			if (seed.qLen > 0) start += seed.qLen;
			else start++;
		}
	}
}

int SeedGrouping()
{
	int i, j, p;

	for (p = i = 0, j = 1; j < SeedNum; i++, j++)
	{
		//printf("check: "); ShowFragPair(SeedVec[i]);
		if ((SeedVec[j].PosDiff - SeedVec[i].PosDiff) > MaxIndelSize)
		{
			//printf("Break!!\n"); ShowFragPair(SeedVec[i]); ShowFragPair(SeedVec[j]);
			SeedGroupVec.push_back(make_pair(p, j));
			p = j;
		}
	}
	if (p < j) SeedGroupVec.push_back(make_pair(p, j));

	return (int)SeedGroupVec.size();
}

bool Check_PD_Frequency(int PD, map<int, int>& PDFmap)
{
	map<int, int>::iterator iter, prev_iter, next_iter;

	iter = PDFmap.find(PD); if (iter->second >= Min_PD_Freq) return true;
	//else if ((prev_iter = iter)-- != PDFmap.end() && iter->first - prev_iter->first == 1 && prev_iter->second >= Min_PD_Freq) return true;
	//else if ((next_iter = iter)++ != PDFmap.end() && next_iter->first - iter->first == 1 && next_iter->second >= Min_PD_Freq) return true;
	else return false;
}

int FindPivot(int i, int j, int BegIdx, int EndIdx, bool *UniqueArr)
{
	int p1, p2, pivot = i;

	p1 = i - 1; p2 = j;
	while (p1 >= BegIdx && (!UniqueArr[p1 - BegIdx] || !SeedVec[p1].bSeed)) p1--;
	while (p2 < EndIdx && (!UniqueArr[p2 - BegIdx] || !SeedVec[p2].bSeed)) p2++;

	if (p1 >= BegIdx || p2 < EndIdx)
	{
		//if (p1 >= BegIdx) printf("pivot1: "), ShowFragPair(SeedVec[p1]);
		//if (p2 < EndIdx)  printf("pivot2: "), ShowFragPair(SeedVec[p2]);
		if (p1 < BegIdx) return pivot = p2;
		if (p2 == EndIdx) return pivot = p1;
		else
		{
			if (SeedVec[i].qPos - SeedVec[p1].qPos < SeedVec[p2].qPos - SeedVec[i].qPos) pivot = p1;
			else pivot = p2;
		}
	}
	return pivot;
}

int64_t FindNeighboringPosDiffAvg(int i, int j, int BegIdx, int EndIdx, bool *UniqueArr)
{
	int64_t sum1, sum2;
	int p1, p2, n1, n2;

	n1 = 0; sum1 = 0;
	//printf("Repeat[%d-%d]\n", i, j-1);
	for (p1 = i - 1; p1 >= BegIdx; p1--)
	{
		if (UniqueArr[p1 - BegIdx] && SeedVec[p1].bSeed)
		{
			//printf("L-%d:", p1); ShowFragPair(SeedVec[p1]);
			n1++; sum1 += SeedVec[p1].PosDiff;
			if (n1 == 5) break;
		}
	}
	n2 = 0; sum2 = 0;
	for (p2 = j; p2 < EndIdx && p2 > BegIdx; p2++)
	{
		if (UniqueArr[p2 - BegIdx] && SeedVec[p2].bSeed)
		{
			//printf("R-%d:", p2); ShowFragPair(SeedVec[p2]);
			n2++; sum2 += SeedVec[p2].PosDiff;
			if (n2 == 5) break;
		}
	}
	if (n1 > 0 || n2 > 0) return (sum1 + sum2) / (n1 + n2);
	else return SeedVec[i].PosDiff;
}

void RemoveRedundantSeeds(int i, int j, int64_t AvgPosDiff)
{
	int k, idx = -1;
	int64_t diff, min_diff = GenomeSize;

	//printf("Avg: %lld\n", AvgPosDiff);
	for (k = i; k < j; k++)
	{
		//printf("\t"); ShowFragPair(SeedVec[k]);
		if ((diff = abs(SeedVec[k].PosDiff - AvgPosDiff)) < MaxIndelSize && diff < min_diff)
		{
			min_diff = diff;
			idx = k;
		}
	}
	//if (idx != -1) printf("\tKeep "), ShowFragPair(SeedVec[idx]);
	for (k = i; k < j; k++) if (k != idx) SeedVec[k].bSeed = false;
}

bool CheckAvgReliability(int64_t avg, int BegIdx, int EndIdx, bool *UniqueArr)
{
	int i, n = 0;
	int64_t diff_sum = 0;

	for (i = BegIdx; i < EndIdx; i++)
	{
		if (UniqueArr[i - BegIdx])
		{
			n++;
			diff_sum += abs(avg - SeedVec[i].PosDiff);
		}
	}
	//printf("avg=%lld, diff_avg=%lld\n", avg, diff_sum / n);
	if ((diff_sum / n) > MaxIndelSize) return false;
	else return true;
}

void RefinePDFmap(int BegIdx, int EndIdx, map<int, int>& PDFmap)
{
	pair<int, int> MaxFreqPair;
	map<int, int>::iterator iter;

	MaxFreqPair = make_pair(0, 0);
	for (iter = PDFmap.begin(); iter != PDFmap.end(); iter++) if (iter->second > MaxFreqPair.second) MaxFreqPair = *iter;

	//Reset frequency for bad PDs
	for (iter = PDFmap.begin(); iter != PDFmap.end(); iter++)
	{
		if (abs(iter->first - MaxFreqPair.first) >= 3) iter->second = 0;
	}
}

void RemoveOutlierSeeds(int BegIdx, int EndIdx, bool *UniqueArr)
{
	int n, i, pd;
	int64_t sum, avg;
	map<int, int> PDFmap; //PosDiff frequency

	n = 0; sum = 0;
	for (i = BegIdx; i < EndIdx; i++)
	{
		if (UniqueArr[i - BegIdx]) PDFmap[((int)(SeedVec[i].PosDiff >> 4))]++;
	}
	RefinePDFmap(BegIdx, EndIdx, PDFmap); 
	for (i = BegIdx; i < EndIdx; i++)
	{
		if (UniqueArr[i - BegIdx])
		{
			if (PDFmap[((int)(SeedVec[i].PosDiff >> 4))] > 0)
			{
				sum += SeedVec[i].PosDiff; n++;
			}
		}
	}
	avg = n > 0 ? sum / n : GenomeSize;

	//printf("avg = %lld (N=%d)\n", avg, n);
	for (i = BegIdx; i < EndIdx; i++)
	{
		if (UniqueArr[i - BegIdx])
		{
			pd = (int)(SeedVec[i].PosDiff >> 4);
			if (abs(avg - SeedVec[i].PosDiff) > MaxIndelSize && !Check_PD_Frequency(pd, PDFmap)) SeedVec[i].bSeed = false;
			//if (abs(avg - SeedVec[i].PosDiff) > MaxIndelSize && !Check_PD_Frequency(pd, PDFmap)) printf("remove %d: ", i), ShowFragPair(SeedVec[i]);
			//else printf("%d: ", i), ShowFragPair(SeedVec[i]);
		}
	}
	//printf("break\n\n");
}

int FindSeedGroupScore(int BegIdx, int EndIdx)
{
	int score = 0;
	for (int i = BegIdx; i < EndIdx; i++) score += SeedVec[i].qLen;
	return score;
}

void SeedGroupAnalysis(int BegIdx, int EndIdx)
{
	bool *UniqueArr;
	int n, p, i, j, k;
	AlnBlock_t AlnBlock;
	vector<FragPair_t> SeedCandidateVec;

	sort(SeedVec.begin() + BegIdx, SeedVec.begin() + EndIdx, CompByQueryPos);
	//if (EndIdx - BegIdx < 10000) return; printf("Group[%d-%d]\n", BegIdx, EndIdx - 1);
	//for (i = BegIdx; i < EndIdx; i++) printf("%d:", i),ShowFragPair(SeedVec[i]);
	UniqueArr = new bool[(EndIdx - BegIdx)]();
	for (i = BegIdx, j = i + 1; i < EndIdx; i++, j++)
	{
		if (j < EndIdx && SeedVec[i].qPos == SeedVec[j].qPos)
		{
			while (++j < EndIdx && SeedVec[i].qPos == SeedVec[j].qPos);
			//for (; i < j; i++) printf("\t%d: ", i), ShowFragPair(SeedVec[i]);
			i = j - 1;
		}
		else UniqueArr[(i - BegIdx)] = true;
	}
	for (n = UniqueArr[0] ? 1 : 0, i = BegIdx, j = BegIdx + 1; j < EndIdx; j++)
	{
		if (UniqueArr[(j - BegIdx)])
		{
			if (SeedVec[j].PosDiff == SeedVec[j - 1].PosDiff) n++;
			else if (++n >= 30 && SeedVec[j].qPos - SeedVec[i].qPos > 3000) //interval
			{
				RemoveOutlierSeeds(i, j, UniqueArr + (i - BegIdx));
				i = j; n = 0;
			}
		}
	}
	RemoveOutlierSeeds(i, EndIdx, UniqueArr + (i - BegIdx));

	//Remove outliers in repetitive regions
	for (i = BegIdx, j = i + 1; i < EndIdx; i++, j++)
	{
		if (j < EndIdx && SeedVec[i].qPos == SeedVec[j].qPos)
		{
			while (++j < EndIdx && SeedVec[i].qPos == SeedVec[j].qPos);
			//RemoveRedundantSeeds(i, j, SeedVec[FindPivot(i, j, BegIdx, EndIdx, UniqueArr)].PosDiff);
			RemoveRedundantSeeds(i, j, FindNeighboringPosDiffAvg(i, j, BegIdx, EndIdx, UniqueArr));
			i = j - 1;
		}
	}
	delete[] UniqueArr;

	sort(SeedVec.begin() + BegIdx, SeedVec.begin() + EndIdx, CompByRemoval); while (!SeedVec[EndIdx - 1].bSeed) EndIdx--;
	//printf("Group\n"); for (i = BegIdx; i < EndIdx; i++) ShowFragPair(SeedVec[i]); printf("End\n\n");
	for (i = BegIdx, j = i + 1, k = j + 1; k < EndIdx; i++, j++, k++)
	{
		if (abs(SeedVec[j].PosDiff - SeedVec[i].PosDiff) > 5 && abs(SeedVec[j].PosDiff - SeedVec[k].PosDiff) > 5)
		{
			//printf("noise\n"),ShowFragPair(SeedVec[i]), ShowFragPair(SeedVec[j]), ShowFragPair(SeedVec[k]);
			SeedVec[j].bSeed = false;
		}
	}
	sort(SeedVec.begin() + BegIdx, SeedVec.begin() + EndIdx, CompByRemoval); while (!SeedVec[EndIdx - 1].bSeed) EndIdx--;
	for (p = i = BegIdx, j = i + 1; j < EndIdx; i++, j++)
	{
		//ShowFragPair(SeedVec[i]);
		if (SeedVec[j].qPos - SeedVec[i].qPos - SeedVec[i].qLen > MaxSeedGap || abs(SeedVec[i].PosDiff - SeedVec[j].PosDiff) > 100)
		{
			//printf("break!\n");
			AddAlnBlock(p, j);
			p = j;
		}
	}
	AddAlnBlock(p, j);
}

void *GenerateAlignmentBlocks(void *arg)
{
	int i;
	while (true)
	{
		pthread_mutex_lock(&Lock);
		i = GroupID++;
		pthread_mutex_unlock(&Lock);

		if (i >= SeedGroupNum) break;
		else if (FindSeedGroupScore(SeedGroupVec[i].first, SeedGroupVec[i].second) < MinAlnBlockScore) continue;
		else if (SeedGroupVec[i].first < SeedGroupVec[i].second) SeedGroupAnalysis(SeedGroupVec[i].first, SeedGroupVec[i].second);
	}
	return (void*)(1);
}

void EstChromosomeSimilarity()
{
	int i, chr_idx;
	map<int64_t, int>::iterator iter;

	AlnBlockNum = (int)AlnBlockVec.size();
	memset(RefChrScoreArr, 0, sizeof(int64_t)*iChromsomeNum);
	for (i = 0; i < AlnBlockNum; i++)
	{
		iter = ChrLocMap.lower_bound(AlnBlockVec[i].FragPairVec.begin()->rPos);
		chr_idx = iter->second;
		RefChrScoreArr[chr_idx] += AlnBlockVec[i].score;
	}
	//printf("Query=%s\n", QueryChrVec[QueryChrIdx].name.c_str()); for (i = 0; i < iChromsomeNum; i++) printf("%s: score=%lld\n", ChromosomeVec[i].name, RefChrScoreArr[i]);
}

bool CheckDuplicatedChrScore(int score1, int score2)
{
	if (score1 > score2 && score1 >= score2 * 2) return true;
	else return false;
}

void RemoveRedundantAlnBlocks(int type) //type1:query, type2:ref
{
	float f1, f2;
	int i, j, chr_idx1, chr_idx2;
	int64_t HeadPos1, HeadPos2, TailPos1, TailPos2, overlap;

	AlnBlockNum = (int)AlnBlockVec.size(); 
	if (type == 1) sort(AlnBlockVec.begin(), AlnBlockVec.end(), CompByAlnBlockQueryPos);
	else sort(AlnBlockVec.begin(), AlnBlockVec.end(), CompByAlnBlockRefPos);

	for (i = 0; i < AlnBlockNum; i++)
	{
		///AlnBlockVec[i].bDup = false;
		if (AlnBlockVec[i].score == 0) continue;
		HeadPos1 = (type == 1 ? AlnBlockVec[i].FragPairVec.begin()->qPos : AlnBlockVec[i].FragPairVec.begin()->rPos);
		TailPos1 = (type == 1 ? AlnBlockVec[i].FragPairVec.rbegin()->qPos + AlnBlockVec[i].FragPairVec.rbegin()->qLen - 1 : AlnBlockVec[i].FragPairVec.rbegin()->rPos + AlnBlockVec[i].FragPairVec.rbegin()->rLen - 1);
		chr_idx1 = ChrLocMap.lower_bound(AlnBlockVec[i].FragPairVec.begin()->rPos)->second;
		if (type == 2 && HeadPos1 >= GenomeSize) ReverseRefCoordinate(HeadPos1, TailPos1);

		for (j = i + 1; j < AlnBlockNum; j++)
		{
			if (AlnBlockVec[j].score == 0) continue;
			HeadPos2 = (type == 1 ? AlnBlockVec[j].FragPairVec.begin()->qPos : AlnBlockVec[j].FragPairVec.begin()->rPos);
			TailPos2 = (type == 1 ? AlnBlockVec[j].FragPairVec.rbegin()->qPos + AlnBlockVec[j].FragPairVec.rbegin()->qLen - 1 : AlnBlockVec[j].FragPairVec.rbegin()->rPos + AlnBlockVec[j].FragPairVec.rbegin()->rLen - 1);
			
			if (type == 1 && HeadPos1 == HeadPos2 && TailPos1 == TailPos2)
			{
				//printf("\n%s: h1=%d t1=%d, h2=%d t2=%d\n", ChromosomeVec[chr_idx1].name, HeadPos1, TailPos1, HeadPos2, TailPos2);
				AlnBlockVec[i].bDup = true;
				AlnBlockVec[j].score = 0;
				continue;
			}

			chr_idx2 = ChrLocMap.lower_bound(AlnBlockVec[j].FragPairVec.begin()->rPos)->second;
			if (type == 2 && HeadPos2 >= GenomeSize) ReverseRefCoordinate(HeadPos2, TailPos2);

			if (HeadPos2 < TailPos1) // overlap
			{
				overlap = (TailPos2 > TailPos1 ? TailPos1 - HeadPos2 : TailPos2 - HeadPos2);
				f1 = 1.*overlap / (TailPos1 - HeadPos1); f2 = 1.*overlap / (TailPos2 - HeadPos2);
				//printf("\nType=%d Overlap=%d, f1=%.4f f2=%.4f\n", type, overlap, f1, f2); ShowAlnBlockBoundary(AlnBlockVec[i].score, AlnBlockVec[i].FragPairVec); ShowAlnBlockBoundary(AlnBlockVec[j].score, AlnBlockVec[j].FragPairVec);
				//if (((HeadPos1 < ObrPos && TailPos1 > ObrPos) || (HeadPos2 < ObrPos && TailPos2 > ObrPos))) printf("!!!\n");
				if ((f1 > f2 && f1 >= 0.9) || (OneOnOneMode && CheckDuplicatedChrScore(RefChrScoreArr[chr_idx2], RefChrScoreArr[chr_idx1])))
				{
					AlnBlockVec[i].score = 0; //printf("Remove first\n");
					break;
				}
				if ((f2 > f1 && f2 >= 0.9) || (OneOnOneMode && CheckDuplicatedChrScore(RefChrScoreArr[chr_idx1], RefChrScoreArr[chr_idx2])))
				{
					AlnBlockVec[j].score = 0; //printf("Remove second\n");
				}
			}
			else break;
		}
	}
	RemoveBadAlnBlocks();
}

void GenomeComparison()
{
	int i;
	vector<AlnBlock_t>::iterator ABiter;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];

	vector<int> vec(iThreadNum); for (i = 0; i < iThreadNum; i++) vec[i] = i;

	RefChrScoreArr = new int64_t[iChromsomeNum]; pthread_mutex_init(&Lock, NULL);
	fprintf(stderr, "Step2. Sequence analysis for all query chromosomes\n");
	for (QueryChrIdx = 0; QueryChrIdx != iQueryChrNum; QueryChrIdx++)
	{
		//if (QueryChrVec[QueryChrIdx].name != "") continue;
		//printf("Query=%s\n", QueryChrVec[QueryChrIdx].name.c_str());
		fprintf(stderr, "\tProcess query chromsomoe: %s...\n", QueryChrVec[QueryChrIdx].name.c_str());
		QrySeqPos = 0; QryChrLength = (uint32_t)QueryChrVec[QueryChrIdx].seq.length();

		SeedVec.clear(); SeedGroupVec.clear(); AlnBlockVec.clear(); //reset data structures

		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, IdentifyLocalMEM, NULL);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL); fprintf(stderr, "\n");

		SeedNum = (int)SeedVec.size(); GroupID = 0; SeedGroupNum = SeedGrouping();
		fprintf(stderr, "\t\tSeed clustering and chaining..."); //iThreadNum = 1;
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, GenerateAlignmentBlocks, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL); fprintf(stderr, "\n");
		fprintf(stderr, "\t\tFix overlapping seeds and close gaps between gaps...");
		AlnBlockNum = (int)AlnBlockVec.size();// iThreadNum = 1;
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, CheckAlnBlockOverlaps, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);
		AlnBlockNum = (int)AlnBlockVec.size();
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, CheckAlnBlockLargeGaps, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL); RemoveBadAlnBlocks();
		AlnBlockNum = (int)AlnBlockVec.size();
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, CheckAlnBlockSpanMultiSeqs, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL); RemoveBadAlnBlocks();

		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) ABiter->bDup = false;
		EstChromosomeSimilarity(); RemoveRedundantAlnBlocks(1); RemoveRedundantAlnBlocks(2);
		AlnBlockNum = (int)AlnBlockVec.size(); //iThreadNum = 1;
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, FillAlnBlockGaps, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL); fprintf(stderr, "\n");
		//sort(AlnBlockVec.begin(), AlnBlockVec.end(), CompByAlnBlockRefPos);
		//for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) ShowAlnBlockBoundary(ABiter->score, ABiter->FragPairVec);

		//sort(AlnBlockVec.begin(), AlnBlockVec.end(), CompByAlnBlockRefPos);
		//for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) ShowAlnBlockBoundary(ABiter->score, ABiter->FragPairVec);
		//for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) ShowFragPairVec(ABiter->FragPairVec);
		//continue;

		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) ABiter->aln_len = ABiter->score = 0;
		fprintf(stderr, "\t\tGenreate sequence alignment..."); //iThreadNum = 1;
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, GenerateFragAlignment, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL); fprintf(stderr, "\n");

		int n = 0; int64_t aln_score = 0, aln_len = 0;
		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
		{
			if ((int)(100 * (1.0*ABiter->score / ABiter->aln_len)) < MinSeqIdy) ABiter->score = 0;
			else
			{
				if (ABiter->bDup) DupAlnNum++;
				n++; aln_len += ABiter->aln_len; aln_score += ABiter->score;
				LocalAlignmentNum++; TotalAlignmentLength += ABiter->aln_len; TotalAlignmentMatches += ABiter->score;
				ABiter->coor = GenCoordinateInfo(ABiter->FragPairVec[0].rPos);
			}
		}
		RemoveBadAlnBlocks();
		if (n == 0) continue;
		fprintf(stderr, "\t\tProduce %d local alignments (length = %lld), ANI=%.2f%%\n", n, (long long)aln_len, (100*(1.0*aln_score/aln_len)));
		if (OutputFormat == 1) fprintf(stderr, "\t\tOutput alignments for query sequence (%s)\n", mafFileName), OutputMAF();
		if (OutputFormat == 2) fprintf(stderr, "\t\tOutput alignments for query sequence (%s)\n", alnFileName), OutputAlignment();
		if (bVCF) fprintf(stderr, "\t\tIdentify sequence variants for query sequence...\n"), VariantIdentification();
		if (bShowPlot && GnuPlotPath != NULL) fprintf(stderr, "\t\tGenerate dotplot for query sequence (%s-%s.ps)\n", OutputPrefix, QueryChrVec[QueryChrIdx].name.c_str()), OutputDotplot();
		fprintf(stderr, "\n");
	}
	if (LocalAlignmentNum > 0) fprintf(stderr, "\tAlignment#=%d (total alignment length=%lld) ANI=%.2f%%, unique alignment#=%d\n", (int)LocalAlignmentNum, (long long)TotalAlignmentLength, (100 * (1.0*TotalAlignmentMatches / TotalAlignmentLength)), int(LocalAlignmentNum - DupAlnNum));
	fprintf(stderr, "\tIt took %lld seconds for genome sequence alignment.\n", (long long)(time(NULL) - StartProcessTime));
	delete[] RefChrScoreArr; delete[] ThreadArr;
}
