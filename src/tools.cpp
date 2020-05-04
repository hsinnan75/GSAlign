#include "structure.h" 

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
	FILE* file;
	if ((file = fopen("/proc/self/status", "r")) != NULL)
	{
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
	else return 0;
}

void ShowFragPair(FragPair_t& FragPair)
{
	printf("q[%d-%d] r[%lld-%lld] L:%d D:%lld \n", FragPair.qPos, FragPair.qPos + FragPair.qLen - 1, (long long)FragPair.rPos, (long long)(FragPair.rPos + FragPair.rLen - 1), FragPair.qLen, (long long)FragPair.PosDiff);
	//if (FragPair.bSeed) printf("q[%d-%d] r[%lld-%lld] L:%d D:%lld \n", FragPair.qPos, FragPair.qPos + FragPair.qLen - 1, FragPair.rPos, FragPair.rPos + FragPair.rLen - 1, FragPair.qLen, FragPair.PosDiff);
	//else printf("q[%d-%d]=%d r[%lld-%lld]=%d\n", FragPair.qPos, FragPair.qPos + FragPair.qLen - 1, FragPair.qLen, FragPair.rPos, FragPair.rPos + FragPair.rLen - 1, FragPair.rLen);
}

void ShowAlnBlockBoundary(int score, vector<FragPair_t>& FragPairVec)
{
	int64_t r1, r2;
	int q1, q2, chr_idx = ChrLocMap.lower_bound(FragPairVec.begin()->rPos)->second;

	q1 = FragPairVec.begin()->qPos; q2 = FragPairVec.rbegin()->qPos + FragPairVec.rbegin()->qLen - 1;
	r1 = FragPairVec.begin()->rPos; r2 = FragPairVec.rbegin()->rPos + FragPairVec.rbegin()->rLen - 1;
	printf("AlnBlockBoundary Q[%d-%d] R[%lld-%lld] chr=%s score = %d size = %d\n", q1, q2, (long long)r1, (long long)r2, ChromosomeVec[chr_idx].name, score, FragPairVec.rbegin()->qPos + FragPairVec.rbegin()->qLen - FragPairVec.begin()->qPos);
	if (r2 - r1 < 100) ShowFragPairVec(FragPairVec);
}

void ShowFragPairVec(vector<FragPair_t>& FragPairVec)
{
	printf("FragPairVec (N=%d)\n", (int)FragPairVec.size());
	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		ShowFragPair(*iter);
		//if (iter->bSeed)
		//{
		//	//printf("\t\tq[%d-%d] r[%lld-%lld] len=%d\n", iter->qPos, iter->qPos + iter->qLen - 1, iter->rPos, iter->rPos + iter->rLen - 1, iter->qLen);
		//	//char *frag = new char[iter->qLen + 1];
		//	//strncpy(frag, QueryChrVec[QueryChrIdx].seq.c_str() + iter->qPos, iter->qLen); frag[iter->qLen] = '\0';
		//	//printf("\t\t%s\n", frag);
		//	//delete[] frag;
		//}
		//else
		//{
		//	printf("\t\tq[%d-%d]=%d r[%lld-%lld]=%d\n", iter->qPos, iter->qPos + iter->qLen - 1, iter->qLen, iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen);
		//	printf("\t\t%s\n\t\t%s\n", iter->aln1.c_str(), iter->aln2.c_str());
		//}
	}
	printf("End\n\n");
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

int CountGapNum(char *aln, int i, int stop)
{
	int n = 0;
	for (; i < stop; i++) if (aln[i] == '-') n++;
	return n;
}

void OutputMAF()
{
	FILE *outFile;
	uint32_t aln_len;
	char *aln1, *aln2;
	int RefIdx, iExtension;
	vector<FragPair_t>::iterator FragPairIter;
	string QueryChrName, RefChrName;

	if (QueryChrIdx == 0)
	{
		outFile = fopen(mafFileName, "w");
		fprintf(outFile, "##maf version=1\n");
	}
	else outFile = fopen(mafFileName, "a");

	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (bAllowDuplication == false && ABiter->bDup) continue;

		aln_len = 0; aln1 = new char[ABiter->aln_len + 1]; aln2 = new char[ABiter->aln_len + 1]; aln1[ABiter->aln_len] = aln2[ABiter->aln_len] = '\0';
		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (FragPairIter->bSeed)
			{
				strncpy(aln1 + aln_len, QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos, FragPairIter->qLen);
				strncpy(aln2 + aln_len, QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos, FragPairIter->qLen);
				aln_len += FragPairIter->qLen;
			}
			else
			{
				strcpy(aln1 + aln_len, FragPairIter->aln1.c_str());
				strcpy(aln2 + aln_len, FragPairIter->aln2.c_str());
				aln_len += (uint32_t)FragPairIter->aln1.length();
			}
		}

		RefIdx = ABiter->coor.ChromosomeIdx;
		QueryChrName = QueryChrVec[QueryChrIdx].name; RefChrName = ChromosomeVec[RefIdx].name;
		if (QueryChrName.length() > RefChrName.length()) RefChrName += string().assign((QueryChrName.length() - RefChrName.length()), ' ');
		else QueryChrName += string().assign((RefChrName.length() - QueryChrName.length()), ' ');

		//ShowAlnBlockBoundary(ABiter->score, ABiter->FragPairVec);
		iExtension = 0;
		if (ABiter->coor.bDir && ((ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen) > (ChromosomeVec[RefIdx].FowardLocation + ChromosomeVec[RefIdx].len))) iExtension = (int)((ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen) - (ChromosomeVec[RefIdx].FowardLocation + ChromosomeVec[RefIdx].len));
		else if (!ABiter->coor.bDir && ((ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen) > (ChromosomeVec[RefIdx].ReverseLocation + ChromosomeVec[RefIdx].len))) iExtension = (int)((ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen) - (ChromosomeVec[RefIdx].ReverseLocation + ChromosomeVec[RefIdx].len));
		if (iExtension > 0)
		{
			//printf("%s: %lld-%lld --> %lld (%d)\n", ChromosomeVec[RefIdx].name, (ABiter->coor.bDir ? ChromosomeVec[RefIdx].FowardLocation : ChromosomeVec[RefIdx].ReverseLocation), (ABiter->coor.bDir ? ChromosomeVec[RefIdx].FowardLocation : ChromosomeVec[RefIdx].ReverseLocation) + ChromosomeVec[RefIdx].len, (ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen), iExtension);
			ABiter->aln_len -= iExtension; ABiter->score -= iExtension;
			ABiter->FragPairVec.rbegin()->rLen -= iExtension;
			ABiter->FragPairVec.rbegin()->qLen -= iExtension;
			aln1[ABiter->aln_len] = aln2[ABiter->aln_len] = '\0';
		}
		if (ABiter->coor.bDir)
		{
			fprintf(outFile, "a score=%d\n", ABiter->bDup ? 1 : ABiter->score);
			fprintf(outFile, "s ref.%s %d %d + %d %s\n", ChromosomeVec[RefIdx].name, ABiter->coor.gPos - 1, ABiter->aln_len - CountGapNum(aln1, 0, ABiter->aln_len), ChromosomeVec[RefIdx].len, aln1);
			fprintf(outFile, "s qry.%s %d %d + %d %s\n\n", (char*)QueryChrName.c_str(), ABiter->FragPairVec[0].qPos, ABiter->aln_len - CountGapNum(aln2, 0, ABiter->aln_len), (uint32_t)QueryChrVec[QueryChrIdx].seq.length(), aln2);
		}
		else
		{
			int64_t rPos = ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen - 1;
			SelfComplementarySeq(ABiter->aln_len, aln1); SelfComplementarySeq(ABiter->aln_len, aln2);
			fprintf(outFile, "a score=%d\n", ABiter->bDup ? 1 : ABiter->score);
			fprintf(outFile, "s ref.%s %d %d + %d %s\n", ChromosomeVec[RefIdx].name, GenCoordinateInfo(rPos).gPos - 1, ABiter->aln_len - CountGapNum(aln1, 0, ABiter->aln_len), ChromosomeVec[RefIdx].len, aln1);
			fprintf(outFile, "s qry.%s %d %d - %d %s\n\n", (char*)QueryChrName.c_str(), (uint32_t)QueryChrVec[QueryChrIdx].seq.length() - (ABiter->FragPairVec.rbegin()->qPos + ABiter->FragPairVec.rbegin()->qLen), ABiter->aln_len - CountGapNum(aln2, 0, ABiter->aln_len), (uint32_t)QueryChrVec[QueryChrIdx].seq.length(), aln2);
		}
		delete[] aln1; delete[] aln2;
	}
	std::fclose(outFile);
}

void OutputAlignment()
{
	FILE *outFile;
	int64_t RefPos;
	char *aln1, *aln2;
	uint32_t pos, aln_len;
	string QueryChrName, RefChrName;
	int p, q, RefIdx, QueryPos, iExtension;
	vector<FragPair_t>::iterator FragPairIter;

	if (QueryChrIdx == 0) outFile = fopen(alnFileName, "w");
	else outFile = fopen(alnFileName, "a");

	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (bAllowDuplication == false && ABiter->bDup) continue;

		aln_len = 0; aln1 = new char[ABiter->aln_len + 1]; aln2 = new char[ABiter->aln_len + 1]; aln1[ABiter->aln_len] = aln2[ABiter->aln_len] = '\0';
		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (FragPairIter->bSeed)
			{
				strncpy(aln1 + aln_len, QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos, FragPairIter->qLen);
				strncpy(aln2 + aln_len, QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos, FragPairIter->qLen);
				aln_len += FragPairIter->qLen;
			}
			else
			{
				strcpy(aln1 + aln_len, FragPairIter->aln1.c_str());
				strcpy(aln2 + aln_len, FragPairIter->aln2.c_str());
				aln_len += (uint32_t)FragPairIter->aln1.length();
			}
		}
		RefIdx = ABiter->coor.ChromosomeIdx;
		QueryChrName = QueryChrVec[QueryChrIdx].name; RefChrName = ChromosomeVec[RefIdx].name;
		if (QueryChrName.length() > RefChrName.length()) RefChrName += string().assign((QueryChrName.length() - RefChrName.length()), ' ');
		else QueryChrName += string().assign((RefChrName.length() - QueryChrName.length()), ' ');

		iExtension = 0;
		if (ABiter->coor.bDir && ((ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen) > (ChromosomeVec[RefIdx].FowardLocation + ChromosomeVec[RefIdx].len))) iExtension = (int)((ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen) - (ChromosomeVec[RefIdx].FowardLocation + ChromosomeVec[RefIdx].len));
		else if (!ABiter->coor.bDir && ((ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen) > (ChromosomeVec[RefIdx].ReverseLocation + ChromosomeVec[RefIdx].len))) iExtension = (int)((ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen) - (ChromosomeVec[RefIdx].ReverseLocation + ChromosomeVec[RefIdx].len));
		if (iExtension > 0)
		{
			//printf("%s: %lld-%lld --> %lld (%d)\n", ChromosomeVec[RefIdx].name, (ABiter->coor.bDir ? ChromosomeVec[RefIdx].FowardLocation : ChromosomeVec[RefIdx].ReverseLocation), (ABiter->coor.bDir ? ChromosomeVec[RefIdx].FowardLocation : ChromosomeVec[RefIdx].ReverseLocation) + ChromosomeVec[RefIdx].len, (ABiter->FragPairVec.rbegin()->rPos + ABiter->FragPairVec.rbegin()->rLen), iExtension);
			ABiter->aln_len -= iExtension; ABiter->score -= iExtension;
			ABiter->FragPairVec.rbegin()->rLen -= iExtension;
			ABiter->FragPairVec.rbegin()->qLen -= iExtension;
			aln1[ABiter->aln_len] = aln2[ABiter->aln_len] = '\0';
		}

		fprintf(outFile, "#Identity = %d / %d (%.2f%%) Orientation = %s\n\n", ABiter->score, ABiter->aln_len, (int)(1000 * (1.0*ABiter->score / ABiter->aln_len)) / 10.0, ABiter->coor.bDir ? "Forward" : "Reverse");
		pos = 0; QueryPos = ABiter->FragPairVec[0].qPos + 1; RefPos = ABiter->coor.gPos;
		while (pos < aln_len)
		{
			p = 80 - CountGapNum(aln1, pos, (pos + 80 > aln_len ? aln_len : pos + 80));
			q = 80 - CountGapNum(aln2, pos, (pos + 80 > aln_len ? aln_len : pos + 80));

			fprintf(outFile, "ref.%s\t%12lld\t%.80s\nqry.%s\t%12d\t%.80s\n\n", RefChrName.c_str(), (long long)RefPos, aln1 + pos, QueryChrName.c_str(), QueryPos, aln2 + pos);
			pos += 80; RefPos += (ABiter->coor.bDir ? p : 0 - p); QueryPos += q;
		}
		delete[] aln1; delete[] aln2;
		fprintf(outFile, "%s\n", string().assign(100, '*').c_str());
	}
	std::fclose(outFile);
}

int64_t CalPosDiffAvg(vector<int64_t>& vec)
{
	int64_t sum = 0;
	int i, j, n, size;

	size = (int)vec.size(); 
	for (n = 0, i = 0, j = 1; j < size; i++, j++)
	{
		if (abs(vec[j] - vec[i]) <= MaxIndelSize)
		{
			n+=2; sum += (vec[i] + vec[j]);
		}
	}
	if (n > 0) return (sum / n);
	else return GenomeSize;
}

void ReverseRefCoordinate(int64_t &pos1, int64_t &pos2)
{
	int64_t tmp_pos;

	tmp_pos = pos1;
	pos1 = TwoGenomeSize - 1 - pos2;
	pos2 = TwoGenomeSize - 1 - tmp_pos;
}
