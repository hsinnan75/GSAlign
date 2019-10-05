#include "structure.h"

//#define KmerSize 8
//#define KmerPower 0x3FFF
#define KmerSize 5
#define KmerPower 0xFF

const char* DNA = "ACGT";

uint32_t CreateKmerID(const char* seq, short pos)
{
	uint32_t i, id, end_pos = pos + KmerSize;

	for (id = 0, i = pos; i < end_pos; i++) id = (id << 2) + nst_nt4_table[(int)seq[i]];

	return id;
}

string DecodeKmerID(uint32_t id)
{
	string kmer;

	kmer.resize(KmerSize);
	for (int i = KmerSize - 1; i >= 0; i--)
	{
		kmer[i] = DNA[(uint8_t)(id & 3)];
		id = (id >> 2);
	}
	return kmer;
}

vector<uint32_t> CreateKmerVecFromReadSeq(int len, char* seq)
{
	vector<uint32_t> vec;
	uint32_t wid, count, head, tail;

	tail = 0; count = 0; head = 0;

	while (count < (uint32_t)KmerSize && tail < (uint32_t)len)
	{
		if (seq[tail++] != 'N') count++;
		else count = 0;
	}
	if (count == KmerSize) // found the first kmer
	{
		wid = CreateKmerID(seq, head);
		vec.push_back(wid);

		for (head += 1; tail < (uint32_t)len; head++, tail++)
		{
			if (seq[tail] != 'N')
			{
				wid = ((wid & KmerPower) << 2) + nst_nt4_table[(int)seq[tail]];
				vec.push_back(wid);
			}
			else
			{
				// find next kmer without 'N'
				count = 0; tail++;
				while (count < (uint32_t)KmerSize && tail < (uint32_t)len)
				{
					if (seq[tail++] != 'N') count++;
					else count = 0;
				}
				if (count == KmerSize)
				{
					wid = CreateKmerID(seq, head);
					vec.push_back(wid);
				}
				else break;
			}
		}
		sort(vec.begin(), vec.end());
	}
	return vec;
}

bool CalGapSimilarity(int qPos1, int qPos2, int64_t rPos1, int64_t rPos2)
{
	int8_t nt1, nt2;
	int64_t r, PD1, PD2;
	bool bSimilar = false;
	int q, idy, q_len, r_len;
	string query_frag, ref_frag;
	vector<unsigned int> KmerVec1, KmerVec2, vec;

	q_len = qPos2 - qPos1; r_len = rPos2 - rPos1;

	if ((PD1 = rPos1 - qPos1) == (PD2 = rPos2 - qPos2))
	{
		//printf("q_len=%d (%d-%d) r_len=%d (%lld-%lld)\n", q_len, qPos1, qPos2 - 1, r_len, rPos1, rPos2 - 1); fflush(stdout);
		for (idy = 0, q = qPos1, r = rPos1; q < qPos2; q++, r++)
		{
			nt1 = nst_nt4_table[(unsigned char)RefSequence[r]];
			nt2 = nst_nt4_table[(unsigned char)QueryChrVec[QueryChrIdx].seq[q]];
			if (nt1 == nt2 || nt1 == 4 || nt2 == 4) idy++;
		}
		//printf("linear scan similarity = %.4f\n", 1.0*idy / q_len);
		if (idy >= q_len*0.5) bSimilar = true;
	}
	if (!bSimilar && q_len <= MaxSeedGap && r_len <= MaxSeedGap)
	{
		query_frag = QueryChrVec[QueryChrIdx].seq.substr(qPos1, q_len);
		ref_frag.resize(r_len); strncpy((char*)ref_frag.c_str(), RefSequence + rPos1, r_len);
		KmerVec1 = CreateKmerVecFromReadSeq(q_len, (char*)query_frag.c_str());
		KmerVec2 = CreateKmerVecFromReadSeq(r_len, (char*)ref_frag.c_str());
		set_intersection(KmerVec1.begin(), KmerVec1.end(), KmerVec2.begin(), KmerVec2.end(), back_inserter(vec));
		if ((int)vec.size() > (q_len + r_len)*0.1) bSimilar = true;

		//if (!bSimilar)
		//{
		//	printf("CommonKmer# = %d (%.5f)\n", (int)vec.size(), 1.0*(int)vec.size() / (q_len + r_len));
		//	string str1, str2;
		//	str1 = ref_frag; str2 = query_frag; ksw2_alignment(r_len, str1, q_len, str2);
		//	printf("ksw2\n%s\n%s\n\n", str1.c_str(), str2.c_str());
		//	//str1 = ref_frag; str2 = query_frag; nw_alignment(str1, str2);
		//	//printf("nw\n%s\n%s\n\n", str1.c_str(), str2.c_str());
		//}
	}
	return bSimilar;
}
