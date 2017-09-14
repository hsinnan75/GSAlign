#include "structure.h"
#include <sys/stat.h>

bwt_t *Refbwt;
bwaidx_t *RefIdx;
time_t StartProcessTime;
vector<QueryChr_t> QueryChrVec;
const char* VersionStr = "0.9.0";
int iThreadNum, iQueryChrNum, MinSeedLength;
bool bDebugMode, bShowSubstitution, bShowIndel;
char *RefSequence, *RefSeqFileName, *IndexFileName, *QueryFileName, *OutputPrefix, *vcfFileName, *alnFileName, *snpFileName, *indFileName, *svsFileName, *gpFileName, *GnuPlotPath;

bool LoadQueryFile()
{
	char ch;
	fstream file;
	string seq, str;
	QueryChr_t QueryChr = {"", ""};

	file.open(QueryFileName, ios_base::in);
	if (!file.is_open()) return false;
	else
	{
		while (!file.eof())
		{
			getline(file, str);
			if (str == "") continue;
			else if (str[0] == '>')
			{
				if (QueryChr.seq != "") QueryChrVec.push_back(QueryChr);
				QueryChr.name = str.substr(1); QueryChr.seq.clear();
				//fprintf(stderr, "\r\tGet sequence of %s...", QueryChr.name.c_str());
			}
			else QueryChr.seq.append(str);
		}
		if (QueryChr.seq != "") QueryChrVec.push_back(QueryChr);
		fprintf(stderr, "\tLoad the query sequences (%d chromosome%s)\n", (int)QueryChrVec.size(), (QueryChrVec.size() > 1 ? "s":""));
	}
	file.close();

	if ((iQueryChrNum = QueryChrVec.size()) == 0) return false;
	else return true;
}

bool CheckOutputPrefix()
{
	int i = 0, len;
	bool bRet = true;

	if (strcmp(OutputPrefix, "/dev/null") == 0) return true;
	
	len = strlen(OutputPrefix);
	for (i = 0; i < len;i++)
	{
		if (!isprint(OutputPrefix[i]) ||
			((int)OutputPrefix[i] >= 32 && (int)OutputPrefix[i] <= 44) ||
			((int)OutputPrefix[i] >= 58 && (int)OutputPrefix[i] <= 64) ||
			((int)OutputPrefix[i] >= 123 && (int)OutputPrefix[i] <= 127)) break;
		else continue;
	}
	if (i < len)
	{
		bRet = false;
		fprintf(stderr, "FatalError: Please specify a valid prefix name\n");
	}
	return bRet;
}

void InitializeOutputFiles()
{
	FILE *outFile;
	int len = strlen(OutputPrefix);

	alnFileName = new char[len + 5]; strcpy(alnFileName, OutputPrefix), strcpy(alnFileName + len, ".aln"); alnFileName[len + 4] = '\0'; outFile = fopen(alnFileName, "w"); fclose(outFile);
	vcfFileName = new char[len + 5]; strcpy(vcfFileName, OutputPrefix), strcpy(vcfFileName + len, ".vcf"); vcfFileName[len + 4] = '\0'; outFile = fopen(vcfFileName, "w"); fclose(outFile);
	if (bShowSubstitution)
	{
		snpFileName = new char[len + 5]; strcpy(snpFileName, OutputPrefix), strcpy(snpFileName + len, ".snp"); snpFileName[len + 4] = '\0'; outFile = fopen(snpFileName, "w"); fclose(outFile);
	}
	if (bShowIndel)
	{
		indFileName = new char[len + 5]; strcpy(indFileName, OutputPrefix), strcpy(indFileName + len, ".ind"); indFileName[len + 4] = '\0'; outFile = fopen(indFileName, "w"); fclose(outFile);
	}
	if (GnuPlotPath != NULL)
	{
		gpFileName = new char[len + 4];  strcpy(gpFileName, OutputPrefix); gpFileName[len + 3] = '\0'; strcpy(gpFileName + len, ".gp");
	}
	//svsFileName = new char[len + 5]; strcpy(svsFileName, OutputPrefix), strcpy(svsFileName + len, ".svs"); svsFileName[len + 4] = '\0';outFile = fopen(svsFileName, "w"); fclose(outFile);
}

void FindGnuPlotPath()
{
	fstream file;
	stringstream ss;
	string fullpath, cmd, str, tmp;

	cmd = "/usr/bin/whereis gnuplot > GnuPlotPath"; system(cmd.c_str());
	file.open("GnuPlotPath"); getline(file, str); file.close();
	ss.clear(); ss.str(str);

	ss >> tmp >> fullpath;
	if (fullpath == "") GnuPlotPath = NULL;
	else
	{
		GnuPlotPath = new char[(int)fullpath.length() + 1];
		strcpy(GnuPlotPath, (const char*)fullpath.c_str());
	}
}

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	iThreadNum = 16;
	bDebugMode = false;
	MinSeedLength = 20;
	bShowSubstitution = false;
	bShowIndel = false;
	RefSequence = RefSeqFileName = IndexFileName = QueryFileName = OutputPrefix = NULL;

	for (i = 1; i < argc; i++)
	{
		parameter = argv[i];

		if (parameter == "-i") IndexFileName = argv[++i];
		else if (parameter == "-r") RefSeqFileName = argv[++i];
		else if (parameter == "-q") QueryFileName = argv[++i];
		else if (parameter == "-t")
		{
			if ((iThreadNum = atoi(argv[++i])) > 40)
			{
				fprintf(stderr, "Warning! Thread number is limited to 40!\n");
				iThreadNum = 40;
			}
		}
		else if (parameter == "-sub") bShowSubstitution = true;
		else if (parameter == "-ind") bShowIndel = true;
		else if (parameter == "-o") OutputPrefix = argv[++i];
		else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
		else fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
	}

	if ((IndexFileName == NULL && RefSeqFileName == NULL)|| QueryFileName == NULL)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "GenAlign v%s\n", VersionStr);
		fprintf(stderr, "Usage: %s [-i IndexFile Prefix / -r Reference file] -q QueryFile[Fasta]\n\n", argv[0]);
		fprintf(stderr, "Options: -t INT        number of threads [16]\n");
		fprintf(stderr, "         -o            Set the prefix of the output files\n");
		fprintf(stderr, "\n");
		exit(0);
	}
	if (OutputPrefix == NULL) OutputPrefix = (char*)"output";
	else if (CheckOutputPrefix() == false) exit(0);

	StartProcessTime = time(NULL);
	fprintf(stderr, "Step1. Load the two genome sequences...\n");

	if (LoadQueryFile() == false) fprintf(stderr, "Please check the query file: %s\n", QueryFileName), exit(0);

	if (IndexFileName != NULL && CheckBWAIndexFiles(IndexFileName)) RefIdx = bwa_idx_load(IndexFileName);
	else if (RefSeqFileName != NULL)
	{
		string cmd(argv[0]), prefix(RefSeqFileName); prefix.resize(prefix.find_first_of('.'));
		cmd = cmd.substr(0, cmd.find_last_of('/') + 1) + "bwa_index " + RefSeqFileName + " " + prefix ;
		system((char*)cmd.c_str()); RefIdx = bwa_idx_load(prefix.c_str());
	}
	else fprintf(stderr, "Please specify a reference genome\n"), exit(0);

	if (RefIdx == 0) fprintf(stderr, "\n\nError! Please check your input!\n");
	else
	{
		Refbwt = RefIdx->bwt;
		RestoreReferenceInfo();

		FindGnuPlotPath(); InitializeOutputFiles();
		GenomeComparison();
		bwa_idx_destroy(RefIdx);
		if (RefSequence != NULL) delete[] RefSequence;
	}
	fprintf(stderr, "\nDone! It took %lld seconds.\n\n\n", (long long)(time(NULL) - StartProcessTime));

	return 0;
}
