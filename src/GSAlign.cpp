#include <cmath>
#include "structure.h"

#define MaxGapSize 300
#define SeedExplorationChunk 10000

int ObrPos = -1;
int64_t *RefChrScoreArr;
vector<Variant_t> VarVec;
static pthread_mutex_t Lock;
static vector<FragPair_t> SeedVec;
static int QrySeqPos, QryChrLength;
int64_t TotalAlignmentLength = 0, LocalAlignmentNum = 0, SNP_num = 0, IND_num = 0, SVS_num = 0;
string LineColorArr[10] = {"red", "blue", "web-green", "dark-magenta", "orange", "yellow", "turquoise", "dark-yellow", "violet", "dark-grey"};

bool CompByPosDiff(const FragPair_t& p1, const FragPair_t& p2)
{
	if (p1.PosDiff == p2.PosDiff) return p1.qPos < p2.qPos;
	else return p1.PosDiff < p2.PosDiff;
}

static bool CompByQueryPos(const FragPair_t& p1, const FragPair_t& p2)
{
	if (p1.qPos == p2.qPos) return p1.rPos < p2.rPos;
	else return p1.qPos < p2.qPos;
}

//static bool CompByRefPos(const FragPair_t& p1, const FragPair_t& p2)
//{
//	if (p1.rPos == p2.rPos) return p1.qPos < p2.qPos;
//	else return p1.rPos < p2.rPos;
//}

bool CompByAlnBlockQueryPos(const AlnBlock_t& p1, const AlnBlock_t& p2)
{
	return p1.FragPairVec.begin()->qPos < p2.FragPairVec.begin()->qPos;
}

bool CompByAlnBlockRefPos(const AlnBlock_t& p1, const AlnBlock_t& p2)
{
	return p1.FragPairVec.begin()->rPos < p2.FragPairVec.begin()->rPos;
}

bool CompByVariantPos(const Variant_t& p1, const Variant_t& p2)
{
	if(p1.chr_idx == p2.chr_idx) return p1.pos < p2.pos;
	else return p1.chr_idx < p2.chr_idx;
}

bool CompByChrScore(const pair<int, int64_t>& p1, const pair<int, int64_t>& p2)
{
	return p1.second > p2.second;
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

void *IdentifyLocalMEM(void *arg)
{
	FragPair_t seed;
	vector<FragPair_t> vec;
	int i, start, stop, num;
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
						vec.push_back(seed);
					}
					delete[] bwtSearchResult.LocArr;
					if (bSensitive) start += 5;
					else start += (bwtSearchResult.len + 1);
				}
				else start++;
			}
		}
	}
	sort(vec.begin(), vec.end(), CompByPosDiff);

	if (iThreadNum == 1) SeedVec.swap(vec);
	else
	{
		pthread_mutex_lock(&Lock);
		num = (int)SeedVec.size();
		copy(vec.begin(), vec.end(), back_inserter(SeedVec));
		inplace_merge(SeedVec.begin(), SeedVec.begin() + num, SeedVec.end(), CompByPosDiff);
		pthread_mutex_unlock(&Lock);
	}
	return (void*)(1);
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

bool CalGapSimilarity(int qPos1, int qPos2, int64_t rPos1, int64_t rPos2)
{
	int q_len, r_len;
	bool bSimilar = false;
	string query_frag, ref_frag;
	vector<unsigned int> KmerVec1, KmerVec2, vec;

	q_len = qPos2 - qPos1; r_len = rPos2 - rPos1; 
	if (q_len < 1000 && r_len < 1000 && abs(abs(rPos1 - qPos1)  - abs(rPos2 - qPos2)) < 50)
	{
		query_frag = QueryChrVec[QueryChrIdx].seq.substr(qPos1, q_len);
		ref_frag.resize(r_len); strncpy((char*)ref_frag.c_str(), RefSequence + rPos1, r_len);
		//printf("PosDiff: %lld vs %lld\n", rPos1 - qPos1, rPos2 - qPos2);
		//printf("len = %d vs %d\n", q_len, r_len);
		//printf("CheckSimilarity\nqry:%s\nref:%s\n\n", query_frag.c_str(), ref_frag.c_str());

		KmerVec1 = CreateKmerVecFromReadSeq(q_len, (char*)query_frag.c_str());
		KmerVec2 = CreateKmerVecFromReadSeq(r_len, (char*)ref_frag.c_str());
		set_intersection(KmerVec1.begin(), KmerVec1.end(), KmerVec2.begin(), KmerVec2.end(), back_inserter(vec));

		//if ((bLowSimilarity && (int)vec.size() > (q_len + r_len)*0.1)
		//	|| (int)vec.size() > (q_len + r_len)*0.25) 
		//	bSimilar = true;
		if ((int)vec.size() > (q_len + r_len)*0.25) bSimilar = true;
	}
	return bSimilar;
}

void AlignmentBlockClustering()
{
	int64_t rTailPos;
	AlnBlock_t AlnBlock;
	vector<pair<int, int> > GroupVec;
	int headIdx, i, j, ci, num, qTailPos, diff;

	if ((num = (int)SeedVec.size()) == 0) return;
	
	for (headIdx = i = 0, j = 1; j < num; i++, j++)
	{
		if (SeedVec[i].qPos == SeedVec[j].qPos)
		{
			GroupVec.push_back(make_pair(headIdx, i));
			GroupVec.push_back(make_pair(i, j));
			headIdx = j;
		}
		if (abs(SeedVec[j].PosDiff - SeedVec[i].PosDiff) > MaxIndelSize)
		{
			GroupVec.push_back(make_pair(headIdx, j));
			headIdx = j;
			//printf("Break!! %lld vs %lld\n\n", SeedVec[i].PosDiff, SeedVec[j].PosDiff);
		}
	}
	GroupVec.push_back(make_pair(headIdx, j));

	bool *VisitSeedArr = new bool[num]();
	for (vector<pair<int, int> >::iterator iter = GroupVec.begin(); iter != GroupVec.end(); iter++)
	{
		sort(SeedVec.begin() + iter->first, SeedVec.begin() + iter->second, CompByQueryPos);
		//for (i = iter->first; i < iter->second;i++) ShowFragPair(SeedVec[i]); printf("GroupBreak!!\n\n");
		
		for (i = iter->first; i < iter->second; i++)
		{
			if (!VisitSeedArr[i])
			{
				//printf("HeadSeed:"); ShowFragPair(SeedVec[i]);
				VisitSeedArr[i] = true; AlnBlock.FragPairVec.clear(); AlnBlock.FragPairVec.push_back(SeedVec[i]);
				//chr_boundary = ChrLocMap.lower_bound(SeedVec[i].rPos)->first;
				qTailPos = SeedVec[i].qPos + SeedVec[i].qLen - 1; rTailPos = SeedVec[i].rPos + SeedVec[i].rLen - 1;
				for (ci = i, j = i + 1; j < iter->second; j++)
				{
					if (!VisitSeedArr[j])
					{
						//if (SeedVec[j].rPos > chr_boundary) break;
						if (SeedVec[j].qPos - qTailPos > MaxGapSize && SeedVec[j].rPos - rTailPos > MaxGapSize && CalGapSimilarity(qTailPos, SeedVec[j].qPos, rTailPos, SeedVec[j].rPos) == false)
						{
							//printf("Break:   "); ShowFragPair(SeedVec[j]); printf("qGaps=%d, rGaps=%d\n", SeedVec[j].qPos - qTailPos, SeedVec[j].rPos - rTailPos);
							break;
						}
						else if (abs(SeedVec[j].PosDiff - SeedVec[ci].PosDiff) < MaxIndelSize)
						{
							//printf("\t "); ShowFragPair(SeedVec[j]);
							AlnBlock.FragPairVec.push_back(SeedVec[j]); VisitSeedArr[j] = true;
							qTailPos = SeedVec[j].qPos + SeedVec[j].qLen - 1; rTailPos = SeedVec[j].rPos + SeedVec[j].rLen - 1;
							ci = j; 
						}
					}
				}
				i = ci;
				if ((AlnBlock.score = CalAlnBlockScore(AlnBlock.FragPairVec)) > MinClusterSize) AlnBlockVec.push_back(AlnBlock);
			}
		}
	}
	delete[] VisitSeedArr;
}

void EstChromosomeSimilarity()
{
	int i, chr_idx, num;
	map<int64_t, int>::iterator iter;

	memset(RefChrScoreArr, 0, sizeof(int64_t)*iChromsomeNum);

	num = (int)AlnBlockVec.size();
	for (i = 0; i < num; i++)
	{
		iter = ChrLocMap.lower_bound(AlnBlockVec[i].FragPairVec.begin()->rPos);
		chr_idx = iter->second;
		RefChrScoreArr[chr_idx] += AlnBlockVec[i].score;
	}
	//for (i = 0; i < iChromsomeNum; i++) printf("%s: score=%lld\n", ChromosomeVec[i].name, RefChrScoreArr[i]);
}

void RemoveRedundantAlnBlocksByQueryPos()
{
	int i, j, num, HeadPos1, HeadPos2, TailPos1, TailPos2, chr_idx1, chr_idx2;

	num = (int)AlnBlockVec.size(); sort(AlnBlockVec.begin(), AlnBlockVec.end(), CompByAlnBlockQueryPos);
	for (i = 0; i < num; i++)
	{
		if (AlnBlockVec[i].score == 0) continue;
		HeadPos1 = AlnBlockVec[i].FragPairVec.begin()->qPos;
		TailPos1 = AlnBlockVec[i].FragPairVec.rbegin()->qPos + AlnBlockVec[i].FragPairVec.rbegin()->qLen - 1;
		chr_idx1 = ChrLocMap.lower_bound(AlnBlockVec[i].FragPairVec.begin()->rPos)->second;

		for (j = i + 1; j < num; j++)
		{
			if (AlnBlockVec[j].score == 0) continue;
			HeadPos2 = AlnBlockVec[j].FragPairVec.begin()->qPos;
			TailPos2 = AlnBlockVec[j].FragPairVec.rbegin()->qPos + AlnBlockVec[j].FragPairVec.rbegin()->qLen - 1;
			chr_idx2 = ChrLocMap.lower_bound(AlnBlockVec[j].FragPairVec.begin()->rPos)->second;

			if (HeadPos2 > TailPos1) break;
			else
			{
				//printf("Overlap!\n"); ShowAlnBlockBoundary(AlnBlockVec[i].score, AlnBlockVec[i].FragPairVec); ShowAlnBlockBoundary(AlnBlockVec[j].score, AlnBlockVec[j].FragPairVec);
				if (TailPos2 < TailPos1) AlnBlockVec[j].score = 0;
				else if (HeadPos1 == HeadPos2)
				{
					if (TailPos1 == TailPos2)
					{
						if (chr_idx1 == chr_idx2) continue;
						else if (RefChrScoreArr[chr_idx1] > RefChrScoreArr[chr_idx2]) AlnBlockVec[j].score = 0;
						else
						{
							//printf("remove first\n");
							AlnBlockVec[i].score = 0;
							break;
						}
					}
					else if (TailPos1 > TailPos2) AlnBlockVec[j].score = 0;
					else
					{
						//printf("remove first\n");
						AlnBlockVec[i].score = 0;
						break;
					}
				}
				else if (TailPos1 == TailPos2) AlnBlockVec[j].score = 0;
				//if (AlnBlockVec[j].score == 0) printf("remove second\n");
				//else printf("Do nothing!\n");
			}
		}
	}
	vector<AlnBlock_t> NewVec; NewVec.reserve(AlnBlockVec.size());
	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) if (ABiter->score > 0) NewVec.push_back(*ABiter);
	AlnBlockVec.swap(NewVec);
}

void RemoveRedundantAlnBlocksByRefPos()
{
	int i, j, num, chr_idx1, chr_idx2;
	int64_t HeadPos1, HeadPos2, TailPos1, TailPos2;

	num = (int)AlnBlockVec.size(); sort(AlnBlockVec.begin(), AlnBlockVec.end(), CompByAlnBlockRefPos);
	for (i = 0; i < num; i++)
	{
		if (AlnBlockVec[i].score == 0) continue;
		HeadPos1 = AlnBlockVec[i].FragPairVec.begin()->rPos;
		TailPos1 = AlnBlockVec[i].FragPairVec.rbegin()->rPos + AlnBlockVec[i].FragPairVec.rbegin()->rLen - 1;
		chr_idx1 = ChrLocMap.lower_bound(AlnBlockVec[i].FragPairVec.begin()->rPos)->second;

		for (j = i + 1; j < num; j++)
		{
			if (AlnBlockVec[j].score == 0) continue;
			HeadPos2 = AlnBlockVec[j].FragPairVec.begin()->rPos;
			TailPos2 = AlnBlockVec[j].FragPairVec.rbegin()->rPos + AlnBlockVec[j].FragPairVec.rbegin()->rLen - 1;
			chr_idx2 = ChrLocMap.lower_bound(AlnBlockVec[j].FragPairVec.begin()->rPos)->second;

			if (HeadPos2 > TailPos1) break;
			else
			{
				//printf("Overlap!\n"); ShowAlnBlockBoundary(AlnBlockVec[i].score, AlnBlockVec[i].FragPairVec); ShowAlnBlockBoundary(AlnBlockVec[j].score, AlnBlockVec[j].FragPairVec);
				if (TailPos2 < TailPos1) AlnBlockVec[j].score = 0;
				else if (HeadPos1 == HeadPos2)
				{
					if (TailPos1 == TailPos2)
					{
						if (chr_idx1 == chr_idx2) continue;
						else if (RefChrScoreArr[chr_idx1] > RefChrScoreArr[chr_idx2]) AlnBlockVec[j].score = 0;
						else
						{
							//printf("remove first\n");
							AlnBlockVec[i].score = 0;
							break;
						}
					}
					else if (TailPos1 > TailPos2) AlnBlockVec[j].score = 0;
					else
					{
						//printf("remove first\n");
						AlnBlockVec[i].score = 0;
						break;
					}
				}
				else if (TailPos1 == TailPos2) AlnBlockVec[j].score = 0;
				//if (AlnBlockVec[j].score == 0) printf("remove second\n");
				//else printf("Do nothing!\n");
			}
		}
	}
	vector<AlnBlock_t> NewVec;
	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) if (ABiter->score > 0) NewVec.push_back(*ABiter);
	AlnBlockVec.swap(NewVec);
}

void CheckOverlaps(vector<FragPair_t>& FragPairVec)
{
	int i, qPos;
	int64_t rPos;
	map<int, bool> qMap;
	bool bOverlap = false;
	map<int64_t, bool> rMap;

	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		for (qPos = iter->qPos, i = 0; i < iter->qLen; i++, qPos++)
		{
			if (qMap.find(qPos) == qMap.end()) qMap.insert(make_pair(qPos, false));
			else
			{
				bOverlap = true;
				//printf("q_error@%d, FragID=%d\n", qPos, iter - FragPairVec.begin());
				break;
			}
		}
		for (rPos = iter->rPos, i = 0; i < iter->rLen; i++, rPos++)
		{
			if (rMap.find(rPos) == rMap.end()) rMap.insert(make_pair(rPos, false));
			else
			{
				bOverlap = true;
				//printf("r_error@%lld, FragID=%d\n", rPos, iter - FragPairVec.begin());
				break;
			}
		}
	}
	if (bOverlap) ShowFragPairVec(FragPairVec);
}

void RemoveOverlaps(vector<FragPair_t>& FragPairVec)
{
	bool bNullPair;
	//bool bShow = false;
	bool bRepeat, bOverlap = false;
	int i, j, pivot, q_overlap_size, r_overlap_size, num;

	//if (ObrPos != -1 && FragPairVec.begin()->rPos < ObrPos && FragPairVec.rbegin()->rPos > ObrPos) bShow=true,ShowFragPairVec(FragPairVec);
	if ((num = (int)FragPairVec.size()) == 1) return;
	for (pivot = i = 0; i < num; i++)
	{
		if (FragPairVec[i].qLen == 0) continue;
		if (i > 0 && i < num - 1 && abs(FragPairVec[i].PosDiff - FragPairVec[i - 1].PosDiff) >= 5 && abs(FragPairVec[i].PosDiff - FragPairVec[i + 1].PosDiff) >= 5)
		{
			FragPairVec[i].qLen = FragPairVec[i].rLen = 0, bNullPair = true;
			continue;
		}
		bRepeat = false;
		for (j = i + 1; j < num; j++)
		{
			if (FragPairVec[j].qLen == 0) continue;
			if (FragPairVec[j].qPos == FragPairVec[i].qPos && FragPairVec[j].qLen == FragPairVec[i].qLen) // repeats
			{
				bRepeat = true;
				FragPairVec[j].qLen = FragPairVec[j].rLen = 0, bNullPair = true; continue;
				//if (abs(FragPairVec[j].PosDiff > abs(FragPairVec[i].PosDiff)))
				//{
				//	FragPairVec[j].qLen = FragPairVec[j].rLen = 0, bNullPair = true; continue;
				//}
				//else
				//{
				//	FragPairVec[i].qLen = FragPairVec[i].rLen = 0, bNullPair = true; break;
				//}
			}
			else if (FragPairVec[j].rPos <= FragPairVec[i].rPos)
			{
				FragPairVec[j].qLen = FragPairVec[j].rLen = 0, bNullPair = true; continue;
			}
			if ((r_overlap_size = FragPairVec[i].rPos + FragPairVec[i].rLen - FragPairVec[j].rPos) > 0)
			{
				//if (bMsg) printf("overlap_size=%d\n", r_overlap_size);
				if (r_overlap_size >= FragPairVec[j].rLen)
				{
					FragPairVec[j].qLen = FragPairVec[j].rLen = 0, bNullPair = true;
				}
				else if (r_overlap_size >= FragPairVec[i].rLen)
				{
					FragPairVec[i].qLen = FragPairVec[i].rLen = 0, bNullPair = true; break;
				}
				else
				{
					bOverlap = true; FragPairVec[i].rLen = (FragPairVec[i].qLen -= r_overlap_size); // shrink block i
				}
			}
			if ((q_overlap_size = FragPairVec[i].qPos + FragPairVec[i].qLen - FragPairVec[j].qPos) > 0)
			{
				if (q_overlap_size >= FragPairVec[j].qLen)
				{
					FragPairVec[j].qLen = FragPairVec[j].rLen = 0, bNullPair = true;
				}
				else if(q_overlap_size >= FragPairVec[i].qLen)
				{
					FragPairVec[i].qLen = FragPairVec[i].rLen = 0, bNullPair = true; break;
				}
				else
				{
					bOverlap = true; FragPairVec[i].rLen = (FragPairVec[i].qLen -= q_overlap_size); // shrink block i
				}
			}
			if (r_overlap_size <= 0 && q_overlap_size <= 0) break;
			//if (FragPairVec[i].PosDiff <= FragPairVec[j].PosDiff) break;
		}
		if (bRepeat) bNullPair = true, FragPairVec[i].qLen = FragPairVec[i].rLen = 0;
	}
	if (bNullPair)
	{
		vector<FragPair_t> vec;

		vec.reserve((int)FragPairVec.size()); num = 0;
		for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
		{
			if (iter->qLen == 0) continue;
			if (num > 0 && abs(iter->PosDiff - vec[num - 1].PosDiff) > MaxIndelSize) continue;
			else
			{
				num++; vec.push_back(*iter);
			}
		}
		FragPairVec.swap(vec);
	}
	//if (bShow) printf("After removing overlaps\n"), ShowFragPairVec(FragPairVec);
}

void IdentifyNormalPairs(vector<FragPair_t>& FragPairVec)
{
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
		if (QuerySeq[i] != TemplateSeq[i] && QuerySeq[i] != 'N') mismatch++;
	}
	return mismatch;
}

void *GenerateFragAlignment(void *arg)
{
	FragPair_t *FragPair;
	int i, j, AlnBlockNum, FragPairNum, TailIdx, mismatch, *my_id = (int*)arg;

	AlnBlockNum = (int)AlnBlockVec.size();

	for (i = *my_id; i < AlnBlockNum; i+= iThreadNum)
	{
		AlnBlockVec[i].score = AlnBlockVec[i].aln_len = 0;

		FragPairNum = (int)AlnBlockVec[i].FragPairVec.size(); TailIdx = FragPairNum - 1;
		for (j = 0; j < FragPairNum; j++)
		{
			if (AlnBlockVec[i].FragPairVec[j].bSeed)
			{
				AlnBlockVec[i].score += AlnBlockVec[i].FragPairVec[j].qLen;
				AlnBlockVec[i].aln_len += AlnBlockVec[i].FragPairVec[j].qLen;
				continue;
			}
			FragPair = &AlnBlockVec[i].FragPairVec[j];
			if (FragPair->qLen == 0)
			{
				AlnBlockVec[i].aln_len += FragPair->rLen;
				FragPair->aln1.resize(FragPair->rLen); strncpy((char*)FragPair->aln1.c_str(), RefSequence + FragPair->rPos, FragPair->rLen);
				FragPair->aln2.assign(FragPair->rLen, '-');
			}
			else if (FragPair->rLen == 0)
			{
				AlnBlockVec[i].aln_len += FragPair->qLen;
				FragPair->aln1.assign(FragPair->qLen, '-');
				FragPair->aln2.resize(FragPair->qLen); strncpy((char*)FragPair->aln2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos, FragPair->qLen);
			}
			else if (FragPair->qLen == FragPair->rLen && ((mismatch = CheckFragPairMismatch(FragPair)) <= 5 || 1.0* mismatch / FragPair->qLen < 0.15))
			{
				FragPair->aln1.resize(FragPair->rLen); strncpy((char*)FragPair->aln1.c_str(), RefSequence + FragPair->rPos, FragPair->rLen);
				FragPair->aln2.resize(FragPair->qLen); strncpy((char*)FragPair->aln2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos, FragPair->qLen);
				AlnBlockVec[i].score += CountIdenticalPairs(FragPair->aln1, FragPair->aln2);
				AlnBlockVec[i].aln_len += FragPair->qLen;
			}
			else
			{
				//if (bDebugMode) printf("GenAln: %d vs %d\n", FragPair->qLen, FragPair->rLen), fflush(stdout);
				FragPair->aln1.resize(FragPair->rLen); strncpy((char*)FragPair->aln1.c_str(), RefSequence + FragPair->rPos, FragPair->rLen);
				FragPair->aln2.resize(FragPair->qLen); strncpy((char*)FragPair->aln2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPair->qPos, FragPair->qLen);
				//if (FragPair->qPos >= 106500 && FragPair->qPos + FragPair->qLen <= 10800)
				//{
				//	printf("Alignment pair:\n%s\n%s\n", FragPair->aln1.c_str(), FragPair->aln2.c_str());
				//}
				nw_alignment(FragPair->rLen, FragPair->aln1, FragPair->qLen, FragPair->aln2);
				//if (FragPair->qPos >= 10600 && FragPair->qPos + FragPair->qLen <= 10800)
				//{
				//	printf("Alignment pair:\n%s\n%s\n", FragPair->aln1.c_str(), FragPair->aln2.c_str());
				//}

				AlnBlockVec[i].score += CountIdenticalPairs(FragPair->aln1, FragPair->aln2);
				AlnBlockVec[i].aln_len += FragPair->aln1.length();
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
					AlnBlockVec[i].aln_len -= (rGaps + qGaps);
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
				if (++k < (int)FragPair->aln1.length())
				{
					FragPair->aln1.resize(k);
					FragPair->aln2.resize(k);
					FragPair->rLen -= rGaps;
					FragPair->qLen -= qGaps;
					AlnBlockVec[i].aln_len -= (rGaps + qGaps);
				}
			}
		}
		if (AlnBlockVec[i].aln_len < MinAlnLength || (int)(100 * (1.0*AlnBlockVec[i].score / AlnBlockVec[i].aln_len)) < MinSeqIdy) AlnBlockVec[i].score = 0;
		else AlnBlockVec[i].coor = GenCoordinateInfo(AlnBlockVec[i].FragPairVec[0].rPos);

		//if (ObrPos != -1 && AlnBlockVec[i].FragPairVec.begin()->rPos < ObrPos && AlnBlockVec[i].FragPairVec.rbegin()->rPos> ObrPos) ShowFragPairVec(AlnBlockVec[i].FragPairVec);
	}
	return (void*)(1);
}

void VariantIdentification()
{
	int64_t rPos;
	Variant_t Variant;
	string frag1, frag2;
	int i, qPos, aln_len, ind_len;
	vector<FragPair_t>::iterator FragPairIter;

	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (!ABiter->coor.bDir || ABiter->score == 0) continue;

		Variant.chr_idx = ABiter->coor.ChromosomeIdx;
		//if (ObrPos != -1 && ABiter->coor.gPos < ObrPos && ABiter->coor.gPos + ABiter->aln_len > ObrPos) OutputDesiredAlignment(*ABiter);
		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (!FragPairIter->bSeed)
			{
				if (FragPairIter->qLen == 0 && FragPairIter->rLen == 0) continue;
				else if (FragPairIter->qLen == 0) // delete
				{
					frag1.resize(FragPairIter->rLen + 1); strncpy((char*)frag1.c_str(), RefSequence + FragPairIter->rPos - 1, FragPairIter->rLen + 1);
					Variant.type = 2;
					Variant.pos = GenCoordinateInfo(FragPairIter->rPos - 1).gPos;
					Variant.ref_frag = frag1;
					Variant.alt_frag.resize(1); Variant.alt_frag[0] = QueryChrVec[QueryChrIdx].seq[FragPairIter->qPos - 1];
					VarVec.push_back(Variant);
					//if (GenCoordinateInfo(rPos - 1).gPos == 1238929) OutputDesiredAlignment(*ABiter);
					//fprintf(outFile, "%s\t%d\t.\t%s\t%c\t100\tPASS\tmt=DELETE\n", RefChrName.c_str(), GenCoordinateInfo(FragPairIter->rPos - 1).gPos, (char*)frag1.c_str(), QueryChrVec[QueryChrIdx].seq[FragPairIter->qPos - 1]);
				}
				else if (FragPairIter->rLen == 0) // insert
				{
					frag2.resize(FragPairIter->qLen + 1); strncpy((char*)frag2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos - 1, FragPairIter->qLen + 1);
					Variant.type = 1;
					Variant.pos = GenCoordinateInfo(FragPairIter->rPos - 1).gPos;
					Variant.ref_frag.resize(1); Variant.ref_frag[0] = RefSequence[FragPairIter->rPos - 1];
					Variant.alt_frag = frag2;
					VarVec.push_back(Variant);
					//fprintf(outFile, "%s\t%d\t.\t%c\t%s\t100\tPASS\tmt=INSERT\n", RefChrName.c_str(), GenCoordinateInfo(FragPairIter->rPos - 1).gPos, RefSequence[FragPairIter->rPos - 1], (char*)frag2.c_str());
				}
				else if (FragPairIter->qLen == 1 && FragPairIter->rLen == 1) // substitution
				{
					if (nst_nt4_table[(int)FragPairIter->aln1[0]] != nst_nt4_table[(int)FragPairIter->aln2[0]] && nst_nt4_table[(int)FragPairIter->aln1[0]] != 4 && nst_nt4_table[(int)FragPairIter->aln2[0]] != 4)
					{
						Variant.type = 0;
						Variant.pos = GenCoordinateInfo(FragPairIter->rPos).gPos;
						Variant.ref_frag = FragPairIter->aln1;
						Variant.alt_frag = FragPairIter->aln2;
						VarVec.push_back(Variant);
						//fprintf(outFile, "%s\t%d\t.\t%c\t%c\t100\tPASS\tmt=SUBSTITUTE\n", RefChrName.c_str(), GenCoordinateInfo(FragPairIter->rPos).gPos, FragPairIter->aln1[0], FragPairIter->aln2[0]);
					}
				}
				else
				{
					//fprintf(stdout, "ref=%s\nqry=%s\n", FragPairIter->aln2.c_str(), FragPairIter->aln1.c_str());
					for (aln_len = (int)FragPairIter->aln1.length(), rPos = FragPairIter->rPos, qPos = FragPairIter->qPos, i = 0; i < aln_len; i++)
					{
						if (FragPairIter->aln1[i] == '-') // insert
						{
							ind_len = 1; while (FragPairIter->aln1[i + ind_len] == '-') ind_len++;
							frag2 = QueryChrVec[QueryChrIdx].seq.substr(qPos - 1, ind_len + 1);
							//fprintf(outFile, "%s\t%d\t.\t%c\t%s\t100\tPASS\tmt=INSERT\n", RefChrName.c_str(), GenCoordinateInfo(rPos - 1).gPos, frag2[0], (char*)frag2.c_str());
							Variant.type = 1;
							Variant.pos = GenCoordinateInfo(rPos - 1).gPos;
							Variant.ref_frag.resize(1); Variant.ref_frag[0] = frag2[0];
							Variant.alt_frag = frag2;
							VarVec.push_back(Variant);

							qPos += ind_len; i += ind_len - 1;
						}
						else if (FragPairIter->aln2[i] == '-') // delete
						{
							//if (GenCoordinateInfo(rPos - 1).gPos == 1238929) OutputDesiredAlignment(*ABiter);
							ind_len = 1; while (FragPairIter->aln2[i + ind_len] == '-') ind_len++;
							frag1.resize(ind_len + 2); frag1[ind_len + 1] = '\0';
							strncpy((char*)frag1.c_str(), RefSequence + rPos - 1, ind_len + 1);
							//fprintf(outFile, "%s\t%d\t.\t%s\t%c\t100\tPASS\tmt=DELETE\n", RefChrName.c_str(), GenCoordinateInfo(rPos - 1).gPos, (char*)frag1.c_str(), frag1[0]);
							Variant.type = 2;
							Variant.pos = GenCoordinateInfo(rPos - 1).gPos;
							Variant.ref_frag = frag1;
							Variant.alt_frag.resize(1); Variant.alt_frag[0] = frag1[0];
							VarVec.push_back(Variant);

							rPos += ind_len; i += ind_len - 1;
						}
						else if (nst_nt4_table[(int)FragPairIter->aln1[i]] != nst_nt4_table[(int)FragPairIter->aln2[i]])// substitute
						{
							if (nst_nt4_table[(int)FragPairIter->aln1[i]] != 4 && nst_nt4_table[(int)FragPairIter->aln2[i]] != 4)
							{
								//fprintf(outFile, "%s\t%d\t.\t%c\t%c\t100\tPASS\tmt=SUBSTITUTE\n", RefChrName.c_str(), GenCoordinateInfo(rPos).gPos, FragPairIter->aln1[i], FragPairIter->aln2[i]);
								Variant.type = 0;
								Variant.pos = GenCoordinateInfo(rPos).gPos;
								Variant.ref_frag.resize(1); Variant.ref_frag[0] = FragPairIter->aln1[i];
								Variant.alt_frag.resize(1); Variant.alt_frag[0] = FragPairIter->aln2[i];
								VarVec.push_back(Variant);
							}
							rPos++; qPos++;
						}
						else
						{
							rPos++; qPos++;
						}
					}
				}
			}
		}
	}
}

void OutputSequenceVariants()
{
	FILE *outFile;
	const char *MutType[3] = { "SUBSTITUTE", "INSERT", "DELETE" };
	
	sort(VarVec.begin(), VarVec.end(), CompByVariantPos);
	
	outFile = fopen(vcfFileName, "w");
	fprintf(outFile, "##fileformat=VCFv4.3\n");
	if(IndexFileName != NULL) fprintf(outFile, "##reference=%s\n", IndexFileName);
	else fprintf(outFile, "##reference=%s\n", RefSeqFileName);
	fprintf(outFile, "##source=GSAlign %s\n", VersionStr);
	fprintf(outFile, "##INFO=<ID=TYPE,Type=String,Description=\"The type of allele, either SUBSTITUTE, INSERT, DELETE, or BND.\">\n");
	for (vector<Variant_t>::iterator iter = VarVec.begin(); iter != VarVec.end(); iter++)
	{
		fprintf(outFile, "%s\t%d\t.\t%s\t%s\t100\tPASS\tTYPE=%s\n", ChromosomeVec[iter->chr_idx].name, iter->pos, (char*)iter->ref_frag.c_str(), (char*)iter->alt_frag.c_str(), MutType[iter->type]);
	}
	std::fclose(outFile);
}

void OutputDotplot()
{
	FILE *outFile;
	vector<int> vec;
	int64_t last_ref_end;
	char tmpFileName[256];
	string cmd, DataFileName;
	map<int, int> ChrColorMap;
	int i, last_query_end, thr;
	map<int, FILE*> ChrFileHandle;
	vector<AlnBlock_t>::iterator ABiter;
	vector<pair<int, int64_t> > ChrScoreVec;

	if (AlnBlockVec.size() == 0) return;

	vec.clear(); vec.resize(iChromsomeNum);
	outFile = fopen(gpFileName, "w"); DataFileName = (string)OutputPrefix + "." + QueryChrVec[QueryChrIdx].name;
	for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (ABiter->score > 0) vec[ABiter->coor.ChromosomeIdx] += ABiter->score;
	}
	for (i = 0; i < iChromsomeNum; i++)
	{
		if (vec[i] >= 1000) ChrScoreVec.push_back(make_pair(i, vec[i]));
	}
	if (ChrScoreVec.size() == 0) return;
	sort(ChrScoreVec.begin(), ChrScoreVec.end(), CompByChrScore);
	if (ChrScoreVec.size() > 5) ChrScoreVec.resize(5);
	thr = ChrScoreVec.rbegin()->second; if (thr < ChrScoreVec.begin()->second / 10) thr = ChrScoreVec.begin()->second / 10;

	for (i = 0; i < (int)ChrScoreVec.size(); i++)
	{
		ChrColorMap[ChrScoreVec[i].first] = i + 1;
		//printf("first=%d, name=%s\n", ChrScoreVec[i].first, ChromosomeVec[ChrScoreVec[i].first].name); fflush(stdout);
		sprintf(tmpFileName, "%svs%s", DataFileName.c_str(), ChromosomeVec[ChrScoreVec[i].first].name);
		ChrFileHandle[ChrScoreVec[i].first] = fopen(tmpFileName, "w");
		fprintf(ChrFileHandle[ChrScoreVec[i].first], "0 0\n0 0\n\n");
	}
	fprintf(outFile, "set terminal postscript color solid 'Courier' 15\nset output '%s-%s.ps'\nset grid\nset border 1\n", OutputPrefix, QueryChrVec[QueryChrIdx].name.c_str());
	for (i = 0; i < (int)ChrScoreVec.size(); i++) fprintf(outFile, "set style line %d lw 4 pt 0 ps 0.5 lc '%s'\n", i + 1, LineColorArr[i].c_str());
	fprintf(outFile, "set xrange[1:*]\nset yrange[1:*]\nset xlabel 'Query (%s)'\nset ylabel 'Ref'\n", (char*)QueryChrVec[QueryChrIdx].name.c_str());
	fprintf(outFile, "plot "); 
	for (i = 0; i < (int)ChrScoreVec.size(); i++)
	{
		sprintf(tmpFileName, "%svs%s", DataFileName.c_str(), ChromosomeVec[ChrScoreVec[i].first].name);
		fprintf(outFile, "'%s' title '%s' with lp ls %d%s", tmpFileName, ChromosomeVec[ChrScoreVec[i].first].name, ChrColorMap[ChrScoreVec[i].first], (i != (int)ChrScoreVec.size() - 1 ? ", " : "\n\n"));
	}
	for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (ABiter->score > 0 && ChrFileHandle.find(ABiter->coor.ChromosomeIdx) != ChrFileHandle.end())
		{
			//FragNum = (int)ABiter->FragPairVec.size() - 1;
			last_query_end = ABiter->FragPairVec.rbegin()->qPos + ABiter->FragPairVec.rbegin()->qLen - 1;
			last_ref_end = ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen - 1;
			if (ABiter->coor.bDir) fprintf(ChrFileHandle[ABiter->coor.ChromosomeIdx], "%d %d\n%d %d\n\n", ABiter->FragPairVec.begin()->qPos + 1, GenCoordinateInfo(ABiter->FragPairVec.begin()->rPos).gPos, last_query_end + 1, GenCoordinateInfo(last_ref_end).gPos);
			else fprintf(ChrFileHandle[ABiter->coor.ChromosomeIdx], "%d %d\n%d %d\n\n", ABiter->FragPairVec.begin()->qPos + 1, GenCoordinateInfo(ABiter->FragPairVec[0].rPos).gPos, last_query_end + 1, GenCoordinateInfo(last_ref_end).gPos);
		}
	}
	for (i = 0; i < (int)ChrScoreVec.size(); i++) fclose(ChrFileHandle[ChrScoreVec[i].first]); fclose(outFile);
	cmd = (string)GnuPlotPath + " " + (string)gpFileName; i = system((char*)cmd.c_str());
	cmd = "rm " + DataFileName + "*"; i = system(cmd.c_str());
}

void GenomeComparison()
{
	int i, n;
	//bool* CoverageArr;
	vector<AlnBlock_t>::iterator ABiter;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];
	int64_t iTotalQueryLength = 0, iCoverage = 0;
	
	vector<int> vec(iThreadNum); for (i = 0; i < iThreadNum; i++) vec[i] = i;

	//if (bDebugMode) iThreadNum = 1;
	if (OutputFormat == 0)
	{
		FILE *outFile;
		outFile = fopen(mafFileName, "w");
		fprintf(outFile, "##maf version=1\n");
		fclose(outFile);
	}
	RefChrScoreArr = new int64_t[iChromsomeNum];
	fprintf(stderr, "Step2. Sequence analysis for all query chromosomes\n");
	for (QueryChrIdx = 0; QueryChrIdx != iQueryChrNum; QueryChrIdx++)
	{
		fprintf(stderr, "\tProcess query chromsomoe: %s...\n", QueryChrVec[QueryChrIdx].name.c_str());
		QrySeqPos = 0; QryChrLength = QueryChrVec[QueryChrIdx].seq.length();

		SeedVec.clear(); AlnBlockVec.clear();
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, IdentifyLocalMEM, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);
		fprintf(stderr, "\n");
		//fprintf(stderr, "\t\tChaining...\n");
		AlignmentBlockClustering();
		if ((int)AlnBlockVec.size() == 0)
		{
			fprintf(stderr, "\t\tGSAlign did not find any similarity between the reference and query (%s).\n", QueryChrVec[QueryChrIdx].name.c_str());
			continue;
		}
		//fprintf(stderr, "\t\tIdentify candidate alignments...\n");
		EstChromosomeSimilarity();
		RemoveRedundantAlnBlocksByQueryPos(); RemoveRedundantAlnBlocksByRefPos();
		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) RemoveOverlaps(ABiter->FragPairVec);
		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) IdentifyNormalPairs(ABiter->FragPairVec);
		//for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) ShowFragPairVec(ABiter->FragPairVec); //ShowAlnBlockBoundary(ABiter->score, ABiter->FragPairVec); // ShowFragPairVec(ABiter->FragPairVec);
		//for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) CheckOverlaps(ABiter->FragPairVec);
		fprintf(stderr, "\t\tGenreate sequence alignment...\n");
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, GenerateFragAlignment, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

		//CoverageArr = new bool[QueryChrVec[QueryChrIdx].seq.length()]();
		//for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
		//{
		//	if (ABiter->score == 0) continue;
		//	memset((CoverageArr + ABiter->FragPairVec.begin()->qPos), true, (int)(ABiter->FragPairVec.rbegin()->qPos + ABiter->FragPairVec.rbegin()->qLen - ABiter->FragPairVec.begin()->qPos));
		//	//ABiter->coor = GenCoordinateInfo(ABiter->FragPairVec[0].rPos);
		//}
		//sort(AlnBlockVec.begin(), AlnBlockVec.end(), CompByAlnBlockCoordinate);
		if (OutputFormat == 0) fprintf(stderr, "\t\tOutput the MAF for query chromosome %s in the file: %s\n", QueryChrVec[QueryChrIdx].name.c_str(), mafFileName), OutputMAF();
		if (OutputFormat == 1) fprintf(stderr, "\t\tOutput the alignment for query chromosome %s in the file: %s\n", QueryChrVec[QueryChrIdx].name.c_str(), alnFileName), OutputAlignment();
		if (bVCF) fprintf(stderr, "\t\tIdentify sequence variants for query chromosome %s\n", QueryChrVec[QueryChrIdx].name.c_str()), VariantIdentification();
		if (bShowPlot && GnuPlotPath != NULL) fprintf(stderr, "\t\tGenerate the dotplot for query chromosome %s in the file: %s-%s.ps\n", QueryChrVec[QueryChrIdx].name.c_str(), OutputPrefix, QueryChrVec[QueryChrIdx].name.c_str()), OutputDotplot();

		iTotalQueryLength += (int)QueryChrVec[QueryChrIdx].seq.length();
		//iTotalQueryLength += (n = (int)QueryChrVec[QueryChrIdx].seq.length());
		//for (i = 0; i < n; i++) if (CoverageArr[i]) iCoverage++;
		//delete[] CoverageArr;
	}
	//if (LocalAlignmentNum > 0 && iTotalQueryLength > 0) fprintf(stderr, "\n\tAlignment # = %lld, Total alignment length = %lld (avgLen=%lld), coverage = %.2f%%\n", (long long)LocalAlignmentNum, (long long)TotalAlignmentLength, (long long)(TotalAlignmentLength/ LocalAlignmentNum), 100*(1.0*iCoverage / iTotalQueryLength));
	if (LocalAlignmentNum > 0 && iTotalQueryLength > 0) fprintf(stderr, "\n\tAlignment#=%lld (total length:%lld)\n", (long long)LocalAlignmentNum, (long long)TotalAlignmentLength);
	delete[] RefChrScoreArr; delete[] ThreadArr;
}
