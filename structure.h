#include <iostream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <ctype.h>
#include <pthread.h>
#include <sys/time.h>

using namespace std;

typedef unsigned char ubyte_t;
typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
typedef unsigned short uint16_t;
typedef unsigned long long uint64_t;
typedef unsigned long long bwtint_t;

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	uint32_t cnt_table[256];
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

typedef struct
{
	bwtint_t x[3];
} bwtintv_t;

typedef struct
{
	int len;
	int freq;
	bwtint_t* LocArr;
} bwtSearchResult_t;

typedef struct
{
	char* name; // chromosome name
	int len; // chromosome length
	int64_t FowardLocation;
	int64_t ReverseLocation;
} Chromosome_t;

typedef struct
{
	bool bDir;
	int gPos;
	int ChromosomeIdx;
} Coordinate_t;

typedef struct
{
	int qPos;
	int64_t rPos;
} SeedFragPair_t;

typedef struct
{
	bool bSeed;
	int qPos; // query position
	int64_t rPos; // reference position
	int qLen; // query frag length
	int rLen; // ref frag length
	int64_t PosDiff;
	string aln1;
	string aln2;
} FragPair_t;

typedef struct
{
	int score;
	int aln_len;
	Coordinate_t coor;
	vector<FragPair_t> FragPairVec;
} AlnBlock_t;

typedef struct
{
	string name;  // chromosome name
	string seq; // chromosome sequence
} QueryChr_t;

// Global variables
extern bwt_t *Refbwt;
extern bwaidx_t *RefIdx;
extern const char* VersionStr;
extern time_t StartProcessTime;
extern map<int64_t, int> ChrLocMap;
extern vector<AlnBlock_t> AlnBlockVec;
extern vector<QueryChr_t> QueryChrVec;
extern unsigned char nst_nt4_table[256];
extern vector<Chromosome_t> ChromosomeVec;
extern bool bDebugMode, bLowSimilarity, bShowPlot;
extern int64_t GenomeSize, TwoGenomeSize, TotalAlignmentLength, LocalAlignmentNum;
extern int QueryChrIdx, iThreadNum, iQueryChrNum, iChromsomeNum, MinSeqIdy, MinSeedLength, MinAlnLength, MinClusterSize, OutputFormat;
extern char *RefSequence, *RefSeqFileName, *IndexFileName, *QueryFileName, *OutputPrefix, *vcfFileName, *mafFileName, *alnFileName, *gpFileName, *GnuPlotPath;

// bwt_index.cpp
extern void RestoreReferenceInfo();
extern void bwa_idx_destroy(bwaidx_t *idx);
extern bwaidx_t *bwa_idx_load(const char *hint);

// bwt_search.cpp
extern bwtSearchResult_t BWT_Search(string& seq, int start, int stop);

// GSAlign.cpp
extern void GenomeComparison();
extern void *GenerateFragAlignment(void *arg);
extern void RemoveOverlaps(vector<FragPair_t>& FragPairVec);
extern void ShowFragPairVec(vector<FragPair_t>& FragPairVec);
extern int CalAlnBlockScore(vector<FragPair_t>& FragPairVec);
extern void IdentifyNormalPairs(vector<FragPair_t>& FragPairVec);
extern bool CompByPosDiff(const FragPair_t& p1, const FragPair_t& p2);

// dupDetection.cpp
extern void dupDetection();

// KmerAnalysis.cpp
extern vector<uint32_t> CreateKmerVecFromReadSeq(int len, char* seq);
//extern vector<FragPair_t> PartitionNormalPair(int len1, char* frag1, int len2, char* frag2);

// GetData.cpp
extern void GetQueryGenomeSeq();
extern bool CheckBWAIndexFiles(string IndexPrefix);

// tools.cpp
extern void OutputMAF();
extern void OutputAlignment();
extern int CheckMemoryUsage();
extern void ShowSeedLocationInfo(int64_t MyPos);
extern Coordinate_t GenCoordinateInfo(int64_t rPos);
extern void SelfComplementarySeq(int len, char* seq);
extern void ShowFragPairVec(vector<FragPair_t>& FragPairVec);

// nw_alignment.cpp
extern void nw_alignment(int m, string& s1, int n, string& s2);
