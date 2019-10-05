#include "structure.h"

string LineColorArr[10] = { "red", "blue", "web-green", "dark-magenta", "orange", "yellow", "turquoise", "dark-yellow", "violet", "dark-grey" };

bool CompByChrScore(const pair<int, int64_t>& p1, const pair<int, int64_t>& p2)
{
	return p1.second > p2.second;
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
