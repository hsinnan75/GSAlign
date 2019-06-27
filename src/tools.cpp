#include "structure.h" 

#define WindowSize 80

static const char ReverseMap[255] =
{
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
	'\0', '\0', '\0', '\0', '\0', '-', '\0', '\0', '\0', '\0', /*  40 -  49 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
	'\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
	'\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
	'\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', /*  90 -  99 */
	'\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
	'N',  '\0', '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', /* 110 - 119 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
	'\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
};

void SelfComplementarySeq(int len, char* seq)
{
	int i, j;
	char aa1, aa2;

	for (j = len - 1, i = 0; i<j; i++, j--)
	{
		aa1 = seq[i]; aa2 = seq[j];
		seq[i] = ReverseMap[(int)aa2]; seq[j] = ReverseMap[(int)aa1];
	}
	if (i == j) seq[i] = ReverseMap[(int)seq[i]];
}

int parseLine(char* line) {
	// This assumes that a digit will be found and the line ends in " Kb".
	int i = strlen(line);
	const char* p = line;
	while (*p <'0' || *p > '9') p++;
	line[i - 3] = '\0';
	i = atoi(p);
	return i;
}

// return value is in GB
int CheckMemoryUsage()
{
	FILE* file = fopen("/proc/self/status", "r");
	int iKB = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL) {
		if (strncmp(line, "VmRSS:", 6) == 0) {
			iKB = parseLine(line);
			break;
		}
	}
	fclose(file);

	if (iKB > 0) return (iKB >> 10);
	else return 0;
}

void ShowFragPair(FragPair_t& FragPair)
{
	printf("q[%d-%d]=%d r[%lld-%lld]=%d PD=%lld\n", FragPair.qPos, FragPair.qPos + FragPair.qLen - 1, FragPair.qLen, FragPair.rPos, FragPair.rPos + FragPair.rLen - 1, FragPair.rLen, FragPair.PosDiff);
}

void ShowAlnBlockBoundary(int score, vector<FragPair_t>& FragPairVec)
{
	int chr_idx = ChrLocMap.lower_bound(FragPairVec.begin()->rPos)->second;
	printf("AlnBlockBoundary Q[%d-%d] R[%lld-%lld] chr=%s score = %d size = %d\n", FragPairVec.begin()->qPos, FragPairVec.rbegin()->qPos + FragPairVec.rbegin()->qLen - 1, FragPairVec.begin()->rPos, FragPairVec.rbegin()->rPos + FragPairVec.rbegin()->rLen - 1, ChromosomeVec[chr_idx].name, score, FragPairVec.rbegin()->qPos + FragPairVec.rbegin()->qLen - FragPairVec.begin()->qPos);
}

void ShowFragPairVec(vector<FragPair_t>& FragPairVec)
{
	printf("FragPairVec (N=%d)\n", (int)FragPairVec.size());
	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (iter->bSeed)
		{
			printf("\t\tq[%d-%d] r[%lld-%lld] PosDiff=%lld, len=%d\n", iter->qPos, iter->qPos + iter->qLen - 1, iter->rPos, iter->rPos + iter->rLen - 1, iter->PosDiff, iter->qLen);
			//char *frag = new char[iter->qLen + 1];
			//strncpy(frag, QueryChrVec[QueryChrIdx].seq.c_str() + iter->qPos, iter->qLen); frag[iter->qLen] = '\0';
			//printf("\t\t%s\n", frag);
			//delete[] frag;
		}
		else
		{
			printf("\t\tq[%d-%d]=%d r[%lld-%lld]=%d\n", iter->qPos, iter->qPos + iter->qLen - 1, iter->qLen, iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen);
			printf("\t\t%s\n\t\t%s\n", iter->aln1.c_str(), iter->aln2.c_str());
		}
	}
	printf("\n\n");
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

void OutputMAF()
{
	char* frag;
	FILE *outFile;
	int i, RefIdx;
	vector<FragPair_t>::iterator FragPairIter;
	string QueryChrName, RefChrName, aln1, aln2, frag1, frag2;

	outFile = fopen(mafFileName, "a");
	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (ABiter->score == 0) continue;

		aln1.clear(); aln2.clear();
		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (FragPairIter->bSeed)
			{
				frag = new char[FragPairIter->qLen + 1]; frag[FragPairIter->qLen] = '\0';
				strncpy(frag, QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos, FragPairIter->qLen);
				aln1 += frag; aln2 += frag; delete[] frag;
			}
			else
			{
				aln1 += FragPairIter->aln1;
				aln2 += FragPairIter->aln2;
			}
		}
		LocalAlignmentNum++; TotalAlignmentLength += ABiter->aln_len;
		RefIdx = ABiter->coor.ChromosomeIdx;
		QueryChrName = QueryChrVec[QueryChrIdx].name; RefChrName = ChromosomeVec[RefIdx].name;
		if (QueryChrName.length() > RefChrName.length()) RefChrName += string().assign((QueryChrName.length() - RefChrName.length()), ' ');
		else QueryChrName += string().assign((RefChrName.length() - QueryChrName.length()), ' ');

		if (ABiter->coor.bDir)
		{
			fprintf(outFile, "a score=%d\n", ABiter->score);
			fprintf(outFile, "s %s %lld %lld + %d %s\n", ChromosomeVec[RefIdx].name, ABiter->coor.gPos, (long long)aln1.length(), ChromosomeVec[RefIdx].len, (char*)aln1.c_str());
			fprintf(outFile, "s %s %lld %lld + %d %s\n\n", (char*)QueryChrName.c_str(), ABiter->FragPairVec[0].qPos + 1, (long long)aln2.length(), (int)QueryChrVec[QueryChrIdx].seq.length(), (char*)aln2.c_str());
		}
		else
		{
			//ShowFragPairVec(ABiter->FragPairVec);
			i = (int)ABiter->FragPairVec.size() - 1;
			int64_t rPos = ABiter->FragPairVec[i].rPos + ABiter->FragPairVec[i].rLen - 1;
			SelfComplementarySeq((int)aln1.length(), (char*)aln1.c_str());
			SelfComplementarySeq((int)aln2.length(), (char*)aln2.c_str());
			fprintf(outFile, "a score=%d\n", ABiter->score);
			fprintf(outFile, "s %s %lld %lld + %d %s\n", ChromosomeVec[RefIdx].name, GenCoordinateInfo(rPos).gPos, (long long)aln1.length(), ChromosomeVec[RefIdx].len, (char*)aln1.c_str());
			fprintf(outFile, "s %s %lld %lld - %d %s\n\n", (char*)QueryChrName.c_str(), (int)QueryChrVec[QueryChrIdx].seq.length() - (ABiter->FragPairVec[i].qPos + ABiter->FragPairVec[i].qLen - 1), (long long)aln2.length(), (int)QueryChrVec[QueryChrIdx].seq.length(), (char*)aln2.c_str());
		}
	}
	std::fclose(outFile);
}

int CountBaseNum(string& frag)
{
	int n = (int)frag.length();

	for (string::iterator iter = frag.begin(); iter != frag.end(); iter++) if (*iter == '-') n--;

	return n;
}

void OutputAlignment()
{
	char* frag;
	FILE *outFile;
	int64_t RefPos;
	int i, p, q, RefIdx, QueryPos;
	vector<FragPair_t>::iterator FragPairIter;
	string QueryChrName, RefChrName, aln1, aln2, frag1, frag2;

	outFile = fopen(alnFileName, "a");
	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (ABiter->score == 0) continue;

		aln1.clear(); aln2.clear();
		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (FragPairIter->bSeed)
			{
				frag = new char[FragPairIter->qLen + 1]; frag[FragPairIter->qLen] = '\0';
				strncpy(frag, QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos, FragPairIter->qLen);
				aln1 += frag; aln2 += frag; delete[] frag;
			}
			else
			{
				aln1 += FragPairIter->aln1;
				aln2 += FragPairIter->aln2;
			}
		}
		LocalAlignmentNum++; TotalAlignmentLength += ABiter->aln_len;
		RefIdx = ABiter->coor.ChromosomeIdx;
		QueryChrName = QueryChrVec[QueryChrIdx].name; RefChrName = ChromosomeVec[RefIdx].name;
		if (QueryChrName.length() > RefChrName.length()) RefChrName += string().assign((QueryChrName.length() - RefChrName.length()), ' ');
		else QueryChrName += string().assign((RefChrName.length() - QueryChrName.length()), ' ');

		fprintf(outFile, "#Identity = %d / %d (%.2f%%) Orientation = %s\n\n", ABiter->score, ABiter->aln_len, (int)(10000 * (1.0*ABiter->score / ABiter->aln_len)) / 100.0, ABiter->coor.bDir ? "Forward" : "Reverse");
		//ShowFragPairVec(ABiter->FragPairVec); printf("\n\n");
		i = 0; QueryPos = ABiter->FragPairVec[0].qPos + 1; RefPos = ABiter->coor.gPos;
		while (i < ABiter->aln_len)
		{
			frag1 = aln1.substr(i, WindowSize); frag2 = aln2.substr(i, WindowSize);
			p = CountBaseNum(frag1); q = CountBaseNum(frag2);

			fprintf(outFile, "%s\t%12d\t%s\n%s\t%12d\t%s\n\n", RefChrName.c_str(), RefPos, frag1.c_str(), QueryChrName.c_str(), QueryPos, frag2.c_str());
			i += WindowSize; RefPos += (ABiter->coor.bDir ? p : 0 - p); QueryPos += q;
		}
		fprintf(outFile, "%s\n", string().assign(100, '*').c_str());
	}
	std::fclose(outFile);
	//printf("Total alignment length = %lld\n", TotalAlnLen);
}

void OutputDesiredAlignment(AlnBlock_t AlnBlock)
{
	char* frag;
	int64_t RefPos;
	int i, p, q, RefIdx, QueryPos;
	vector<FragPair_t>::iterator FragPairIter;
	string QueryChrName, RefChrName, aln1, aln2, frag1, frag2;

	for (FragPairIter = AlnBlock.FragPairVec.begin(); FragPairIter != AlnBlock.FragPairVec.end(); FragPairIter++)
	{
		if (FragPairIter->bSeed)
		{
			frag = new char[FragPairIter->qLen + 1]; frag[FragPairIter->qLen] = '\0';
			strncpy(frag, QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos, FragPairIter->qLen);
			aln1 += frag; aln2 += frag; delete[] frag;
		}
		else
		{
			aln1 += FragPairIter->aln1;
			aln2 += FragPairIter->aln2;
		}
	}
	RefIdx = AlnBlock.coor.ChromosomeIdx;
	QueryChrName = QueryChrVec[QueryChrIdx].name; RefChrName = ChromosomeVec[RefIdx].name;
	if (QueryChrName.length() > RefChrName.length()) RefChrName += string().assign((QueryChrName.length() - RefChrName.length()), ' ');
	else QueryChrName += string().assign((RefChrName.length() - QueryChrName.length()), ' ');

	fprintf(stdout, "#Identity = %d / %d (%.2f%%) Orientation = %s\n\n", AlnBlock.score, AlnBlock.aln_len, (int)(10000 * (1.0*AlnBlock.score / AlnBlock.aln_len)) / 100.0, AlnBlock.coor.bDir ? "Forward" : "Reverse");
	//ShowFragPairVec(ABiter->FragPairVec); printf("\n\n");
	i = 0; QueryPos = AlnBlock.FragPairVec[0].qPos + 1; RefPos = AlnBlock.coor.gPos;
	while (i < AlnBlock.aln_len)
	{
		frag1 = aln1.substr(i, WindowSize); frag2 = aln2.substr(i, WindowSize);
		p = CountBaseNum(frag1); q = CountBaseNum(frag2);

		fprintf(stdout, "%s\t%12d\t%s\n%s\t%12d\t%s\n\n", RefChrName.c_str(), RefPos, frag1.c_str(), QueryChrName.c_str(), QueryPos, frag2.c_str());
		i += WindowSize; RefPos += (AlnBlock.coor.bDir ? p : 0 - p); QueryPos += q;
	}
	fprintf(stdout, "%s\n", string().assign(100, '*').c_str());
}