#include "structure.h"
#include <sys/stat.h>

bwt_t *Refbwt;
string cmd_line;
bwaidx_t *RefIdx;
time_t StartProcessTime;
vector<AlnBlock_t> AlnBlockVec;
vector<QueryChr_t> QueryChrVec;
const char* VersionStr = "1.0.0";
bool bDebugMode, bDUPmode, bSensitive, bVCF, bShowPlot;
int QueryChrIdx, iThreadNum, iQueryChrNum, MaxIndelSize, MinSeedLength, MaxSeedLength, MinSeqIdy, MinClusterSize, MinAlnLength, OutputFormat = 0;
char *RefSequence, *RefSeqFileName, *IndexFileName, *QueryFileName, *OutputPrefix, *vcfFileName, *mafFileName, *alnFileName, *gpFileName, *GnuPlotPath;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "GenAlign v%s\n", VersionStr);
	fprintf(stderr, "Usage: %s [-i IndexFile Prefix / -r Reference file] -q QueryFile[Fasta]\n\n", program);
	fprintf(stderr, "Options: -t     INT     number of threads [%d]\n", iThreadNum);
	fprintf(stderr, "         -o     STR     Set the prefix of the output files [output]\n");
	fprintf(stderr, "         -sensitive     Sensitive mode [False]\n");
	fprintf(stderr, "         -fmt   INT     Set the output format 0:maf, 1:aln [%d]\n", OutputFormat);
	fprintf(stderr, "         -idy   INT     Set the minimal sequence identity (0-100) of a local alignment [%d]\n", MinSeqIdy);
	fprintf(stderr, "         -slen  INT     Set the minimal seed length [%d]\n", MinSeedLength);
	fprintf(stderr, "         -mlen  INT     Set the maximal seed length [%d]\n", MaxSeedLength);
	fprintf(stderr, "         -alen  INT     Set the minimal alignment length [%d]\n", MinAlnLength);
	fprintf(stderr, "         -ind   INT     Set the maximal indel size [%d]\n", MaxIndelSize);
	fprintf(stderr, "         -clr   INT     Set the minimal cluster size [%d]\n", MinClusterSize);
	fprintf(stderr, "         -dp            Output Dot-plots\n");
	fprintf(stderr, "         -dup           Duplication detection mode [false]\n");
	fprintf(stderr, "\n");
}

string TrimChromosomeName(string name)
{
	string str;
	int i, len = (int)name.length();

	for (i = 0; i < len; i++)
	{
		//if (isalnum(name[i]) == 0) break;
		if (name[i] == '|') name[i] = '-';
		else if (name[i] == ' ' || name[i] == '#' || name[i] == ':' || name[i] == '=' || name[i] == '\t') break;
	}
	return name.substr(0, i);
}

bool LoadQueryFile()
{
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
				QueryChr.name = TrimChromosomeName(str.substr(1)); QueryChr.seq.clear();
				//fprintf(stderr, "\r\tGet sequence of %s...", QueryChr.name.c_str());
			}
			else QueryChr.seq.append(str);
		}
		if (QueryChr.seq != "") QueryChrVec.push_back(QueryChr);
		fprintf(stderr, "\tLoad the query sequences (%d %s)\n", (int)QueryChrVec.size(), (QueryChrVec.size() > 1 ? "chromosomes":"chromosome"));
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

	mafFileName = alnFileName = vcfFileName = gpFileName = NULL;

	if (OutputFormat == 0)
	{
		mafFileName = new char[len + 5]; strcpy(mafFileName, OutputPrefix), strcpy(mafFileName + len, ".maf"); mafFileName[len + 4] = '\0'; outFile = fopen(mafFileName, "w"); fclose(outFile);
	}
	if (OutputFormat == 1)
	{
		alnFileName = new char[len + 5]; strcpy(alnFileName, OutputPrefix), strcpy(alnFileName + len, ".aln"); alnFileName[len + 4] = '\0'; outFile = fopen(alnFileName, "w"); fclose(outFile);
	}
	if (GnuPlotPath != NULL)
	{
		gpFileName = new char[len + 4];  strcpy(gpFileName, OutputPrefix); gpFileName[len + 3] = '\0'; strcpy(gpFileName + len, ".gp");
	}
	vcfFileName = new char[len + 5]; strcpy(vcfFileName, OutputPrefix), strcpy(vcfFileName + len, ".vcf"); vcfFileName[len + 4] = '\0'; outFile = fopen(vcfFileName, "w"); fclose(outFile);
}

void DestroyOutputFileNames()
{
	if (mafFileName != NULL) delete[] mafFileName;
	if (alnFileName != NULL) delete[] alnFileName;
	if (vcfFileName != NULL) delete[] vcfFileName;
	if (gpFileName != NULL) delete[] gpFileName;
}

void FindGnuPlotPath()
{
	int i;
	fstream file;
	stringstream ss;
	string fullpath, cmd, str, tmp;

	cmd = "/usr/bin/whereis gnuplot > GnuPlotPath"; i = system(cmd.c_str());
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
	bSensitive = false;
	bShowPlot = false;
	bDebugMode = false;
	bVCF = true;
	MinSeedLength = 20;
	MaxSeedLength = 100000000;
	MinClusterSize = 50;
	MinAlnLength = 200;
	MinSeqIdy = 30;
	MaxIndelSize = 35;
	RefSequence = RefSeqFileName = IndexFileName = QueryFileName = OutputPrefix = NULL;

	if (argc == 1 || strcmp(argv[1], "-h") == 0)
	{
		ShowProgramUsage(argv[0]);
		exit(0);
	}
	else if (strcmp(argv[1], "update") == 0)
	{
		i = system("git fetch; git merge origin/master master;make");
		exit(0);
	}
	else
	{
		cmd_line = argv[0];
		for (i = 1; i < argc; i++)
		{
			parameter = argv[i]; cmd_line += " " + parameter;

			if (parameter == "-i") IndexFileName = argv[++i];
			else if (parameter == "-r" &&  i + 1 < argc) RefSeqFileName = argv[++i];
			else if (parameter == "-q" && i + 1 < argc) QueryFileName = argv[++i];
			else if (parameter == "-dup") bDUPmode = true;
			else if (parameter == "-t" && i + 1 < argc)
			{
				if ((iThreadNum = atoi(argv[++i])) < 0)
				{
					fprintf(stderr, "Warning! Thread number should be greater than 0!\n");
					iThreadNum = 16;
				}
			}
			else if (parameter == "-slen" && i + 1 < argc)
			{
				MinSeedLength = atoi(argv[++i]);
				if (MinSeedLength < 10 || MinSeedLength > 20)
				{
					fprintf(stderr, "Warning! minimal seed length is between 10~20!\n");
					exit(0);
				}
			}
			else if (parameter == "-ind" && i + 1 < argc)
			{
				MaxIndelSize = atoi(argv[++i]);
				if (MaxIndelSize < 10 || MaxIndelSize > 100)
				{
					fprintf(stderr, "Warning! maximal indel size is between 10~100!\n");
					exit(0);
				}
			}
			else if (parameter == "-mlen" && i + 1 < argc) MaxSeedLength = atoi(argv[++i]);
			else if (parameter == "-sensitive") bSensitive = true;
			else if (parameter == "-idy" && i + 1 < argc) MinSeqIdy = atoi(argv[++i]);
			else if (parameter == "-alen" && i + 1 < argc) MinAlnLength = atoi(argv[++i]);
			else if (parameter == "-clr" && i + 1 < argc) MinClusterSize = atoi(argv[++i]);
			else if (parameter == "-dp") bShowPlot = true;
			else if (parameter == "-fmt" && i + 1 < argc) OutputFormat = atoi(argv[++i]);
			else if (parameter == "-o") OutputPrefix = argv[++i];
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			else fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
		}
	}
	if ((IndexFileName == NULL && RefSeqFileName == NULL)|| QueryFileName == NULL)
	{
		ShowProgramUsage(argv[0]);
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
		cmd = cmd.substr(0, cmd.find_last_of('/') + 1) + "bwt_index " + RefSeqFileName + " " + prefix ;
		i = system((char*)cmd.c_str()); RefIdx = bwa_idx_load(prefix.c_str());
	}
	else fprintf(stderr, "Please specify a reference genome\n"), exit(0);

	if (RefIdx == 0) fprintf(stderr, "\n\nError! Please check your input!\n");
	else
	{
		if (MaxSeedLength < MinSeedLength) MaxSeedLength = MinSeedLength;

		Refbwt = RefIdx->bwt;
		RestoreReferenceInfo();
		if (bSensitive) MinSeedLength = 10;
		if (bShowPlot) FindGnuPlotPath();
		InitializeOutputFiles();
		if (bDUPmode) dupDetection();
		else GenomeComparison();
		if (bVCF) fprintf(stderr, "\tOutput all identified sequence variants to %s...\n", vcfFileName), OutputSequenceVariants();
		bwa_idx_destroy(RefIdx);
		if (RefSequence != NULL) delete[] RefSequence;
		DestroyOutputFileNames();
	}
	fprintf(stderr, "\nDone! It took %lld seconds.\n\n\n", (long long)(time(NULL) - StartProcessTime));

	return 0;
}
