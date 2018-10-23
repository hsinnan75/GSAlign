#include "structure.h"

#define KmerSize 8
#define KmerPower 0x3FFF

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

	tail = 0; count = 0;

	while (count < KmerSize && tail < len)
	{
		if (seq[tail++] != 'N') count++;
		else count = 0;
	}
	if (count == KmerSize) // found the first kmer
	{
		wid = CreateKmerID(seq, head);
		vec.push_back(wid);

		for (head += 1; tail < len; head++, tail++)
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
				while (count < KmerSize && tail < len)
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
