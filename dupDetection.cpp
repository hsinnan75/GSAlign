#include "structure.h"

#define DupBlockSize 30
#define SeedExplorationChunk 10000

static int QryChrLength;
static pthread_mutex_t Lock;
unsigned short* RefCoverageArr;
static vector<SeedFragPair_t> SeedVec;
static vector<FragPair_t>  VarDupSeedVec;

static bool CompByQueryPos(const SeedFragPair_t& p1, const SeedFragPair_t& p2)
{
	if (p1.qPos == p2.qPos) return p1.rPos < p2.rPos;
	else return p1.qPos < p2.qPos;
}

static bool CompByRefPos(const SeedFragPair_t& p1, const SeedFragPair_t& p2)
{
	if (p1.rPos == p2.rPos) return p1.qPos < p2.qPos;
	else return p1.rPos < p2.rPos;
}

void *IdentifyRepetitiveMEM(void *arg)
{
	SeedFragPair_t seed;
	vector<SeedFragPair_t> Qvec, Rvec;
	bwtSearchResult_t bwtSearchResult;
	int i, pos, start, stop, num, *my_id = (int*)arg;

	string& seq = QueryChrVec[QueryChrIdx].seq;
	for (pos = (*my_id*SeedExplorationChunk); pos < QryChrLength; pos += (iThreadNum * SeedExplorationChunk))
	{
		start = pos; if ((stop = start + SeedExplorationChunk) > QryChrLength) stop = QryChrLength;
		while (start < stop)
		{
			if (nst_nt4_table[(int)seq[start]] > 3) start++;
			else
			{
				bwtSearchResult = BWT_Search(seq, start, (start + DupBlockSize));

				if (bwtSearchResult.freq > 0 && bwtSearchResult.len == DupBlockSize)
				{
					seed.qPos = start;
					for (i = 0; i != bwtSearchResult.freq; i++)
					{
						seed.rPos = bwtSearchResult.LocArr[i]; Qvec.push_back(seed); 
						//printf("q[%d] r[%lld]\n", seed.qPos, seed.rPos);
					}
					delete[] bwtSearchResult.LocArr;
				}
				start += DupBlockSize;
			}
		}
		if (*my_id == 0) fprintf(stderr, "\r\t\tSeed exploration: %d / %d (%d%%)...", stop, QryChrLength, (int)(100 * ((1.0*stop / QryChrLength))));
	}
	if (*my_id == 0) fprintf(stderr, "\r\t\tSeed exploration: %d / %d (100%%)...done!\n", QryChrLength, QryChrLength);

	sort(Qvec.begin(), Qvec.end(), CompByQueryPos);

	if (iThreadNum == 1) SeedVec.swap(Qvec);
	else
	{
		pthread_mutex_lock(&Lock);
		num = (int)SeedVec.size();
		copy(Qvec.begin(), Qvec.end(), back_inserter(SeedVec));
		inplace_merge(SeedVec.begin(), SeedVec.begin() + num, SeedVec.end(), CompByQueryPos);
		pthread_mutex_unlock(&Lock);
	}
	return (void*)(1);
}

void Identify_VarDup_Seeds()
{
	int i, j, num, m, n;
	FragPair_t FragPair;
	vector<SeedFragPair_t>::iterator iter;

	num = (int)SeedVec.size(); FragPair.bSeed = true;
	for (i = 0; i < num;)
	{
		for (m = 1, j = i + 1; j < num; j++)
		{
			if (SeedVec[j].qPos != SeedVec[i].qPos) break;
			else m++;
		}
		for (; i < j; i++)
		{
			if (m != RefCoverageArr[(SeedVec[i].rPos / DupBlockSize)])
			{
				FragPair.qPos = SeedVec[i].qPos;
				FragPair.rPos = SeedVec[i].rPos;
				FragPair.PosDiff = SeedVec[i].rPos - SeedVec[i].qPos;
				FragPair.qLen = FragPair.rLen = DupBlockSize;
				VarDupSeedVec.push_back(FragPair);
			}
		}
	}
	sort(VarDupSeedVec.begin(), VarDupSeedVec.end(), CompByPosDiff);
}

void ExtendLeftEnd(FragPair_t& FragPair)
{
	int qPos = FragPair.qPos;
	int64_t rPos = FragPair.rPos;
	string& seq = QueryChrVec[QueryChrIdx].seq;

	while (qPos > 0 && rPos > 0 && nst_nt4_table[seq[qPos - 1]] == nst_nt4_table[RefSequence[rPos - 1]])
	{
		qPos--; rPos--;
	}
	FragPair.qPos = qPos; FragPair.rPos = rPos;
}

void ExtendRightEnd(FragPair_t& FragPair)
{
	int qPos = FragPair.qPos, qStop;
	int64_t rPos = FragPair.rPos, rStop;
	string& seq = QueryChrVec[QueryChrIdx].seq;

	qStop = QryChrLength - 1; rStop = ChrLocMap.lower_bound(rPos)->first - 1;
	while (qPos < QryChrLength && rPos < rStop && nst_nt4_table[seq[qPos + 1]] == nst_nt4_table[RefSequence[rPos + 1]])
	{
		qPos++; rPos++;
	}
	FragPair.qPos = qPos; FragPair.rPos = rPos;
}

void Cluster_VarDup_Seeds()
{
	AlnBlock_t AlnBlock;
	int i, j, num, headIdx;

	num = (int)VarDupSeedVec.size();
	for (i = 0, j = 1; j < num; i++, j++)
	{
		AlnBlock.FragPairVec.push_back(VarDupSeedVec[i]);

		//printf("q[%d] r[%lld] PosDiff=%lld\n", VarDupSeedVec[i].qPos, VarDupSeedVec[i].rPos, VarDupSeedVec[i].PosDiff);
		if (VarDupSeedVec[j].PosDiff != VarDupSeedVec[i].PosDiff || VarDupSeedVec[j].qPos - VarDupSeedVec[i].qPos > 100)
		{
			ExtendLeftEnd(*AlnBlock.FragPairVec.begin()); ExtendRightEnd(*AlnBlock.FragPairVec.rbegin());
			//printf("Break!!\n\n");
			if ((AlnBlock.score = CalAlnBlockScore(AlnBlock.FragPairVec)) > MinClusterSize) AlnBlockVec.push_back(AlnBlock);
			//printf("AlnBlockScore=%d\n", AlnBlock.score); ShowFragPairVec(AlnBlock.FragPairVec);
			AlnBlock.FragPairVec.clear();
		}
	}
	if ((AlnBlock.score = CalAlnBlockScore(AlnBlock.FragPairVec)) > MinClusterSize) AlnBlockVec.push_back(AlnBlock);
}

void dupDetection()
{
	int i;
	vector<AlnBlock_t>::iterator ABiter;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];
	int64_t obs_pos, iTotalQueryLength = 0, iCoverage = 0;

	vector<int> vec(iThreadNum); for (i = 0; i < iThreadNum; i++) vec[i] = i;

	if (OutputFormat == 0)
	{
		FILE *outFile;
		outFile = fopen(mafFileName, "w");
		fprintf(outFile, "##maf version=1\n");
		fclose(outFile);
	}
	RefCoverageArr = new unsigned short[(TwoGenomeSize / DupBlockSize)]();

	fprintf(stderr, "Step2. Sequence analysis for all query chromosomes\n");
	for (QueryChrIdx = 0; QueryChrIdx != iQueryChrNum; QueryChrIdx++)
	{
		fprintf(stderr, "\tProcess query chromsomoe: %s...\n", QueryChrVec[QueryChrIdx].name.c_str());

		QryChrLength = QueryChrVec[QueryChrIdx].seq.length();

		SeedVec.clear(); VarDupSeedVec.clear();  AlnBlockVec.clear();
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, IdentifyRepetitiveMEM, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

		for (vector<SeedFragPair_t>::iterator iter = SeedVec.begin(); iter != SeedVec.end(); iter++) RefCoverageArr[(iter->rPos / DupBlockSize)]++;

		Identify_VarDup_Seeds(); Cluster_VarDup_Seeds();
		//for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) ShowFragPairVec(ABiter->FragPairVec);
		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) RemoveOverlaps(ABiter->FragPairVec);
		for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) IdentifyNormalPairs(ABiter->FragPairVec);
		//for (ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++) CheckOverlaps(ABiter->FragPairVec);
		fprintf(stderr, "\t\tGenreate sequence alignment...\n");
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, GenerateFragAlignment, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

		if (OutputFormat == 0) fprintf(stderr, "\tOutput the MAF for query chromosome %s in the file: %s\n", QueryChrVec[QueryChrIdx].name.c_str(), mafFileName), OutputMAF();
		if (OutputFormat == 1) fprintf(stderr, "\tOutput the alignment for query chromosome %s in the file: %s\n", QueryChrVec[QueryChrIdx].name.c_str(), alnFileName), OutputAlignment();
	}
	delete[] ThreadArr;
	delete[] RefCoverageArr;
}
