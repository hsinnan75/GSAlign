#include "structure.h"

#define MaxSeedFreq 100

#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

//unsigned char nst_nt4_table[256] = {
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
//};

#define __occ_aux4(bwt, b)											\
	((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
	 + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))
#define bwt_bwt(b, k) ((b)->bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])
#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

static inline int __occ_aux(bwtint_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
	bwtint_t n;
	uint32_t *p, *end;

	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (bwtint_t)(-1)) return 0;
	k -= (k >= bwt->primary); // because $ is not in bwt

	// retrieve Occ at k/OCC_INTERVAL
	n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
	p += sizeof(bwtint_t); // jump to the start of the first BWT cell

	// calculate Occ up to the last k/32
	end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
	for (; p < end; p += 2) n += __occ_aux((bwtint_t)p[0]<<32 | p[1], c);

	// calculate Occ
	n += __occ_aux(((bwtint_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
	if (c == 0) n -= ~k&31; // corrected for the masked bits

	return n;
}

void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
	bwtint_t x;
	uint32_t *p, tmp, *end;
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, 4 * sizeof(bwtint_t));
		return;
	}
	k -= (k >= bwt->primary); // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	memcpy(cnt, p, 4 * sizeof(bwtint_t));
	p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
	end = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4)); // this is the end point of the following loop
	for (x = 0; p < end; ++p) x += __occ_aux4(bwt, *p);
	tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
	x += __occ_aux4(bwt, tmp) - (~k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
	bwtint_t _k, _l;
	_k = k - (k >= bwt->primary);
	_l = l - (l >= bwt->primary);
	if (_l >> OCC_INTV_SHIFT != _k >> OCC_INTV_SHIFT || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		bwtint_t x, y;
		uint32_t *p, tmp, *endk, *endl;
		k -= (k >= bwt->primary); // because $ is not in bwt
		l -= (l >= bwt->primary);
		p = bwt_occ_intv(bwt, k);
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
		p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
		// prepare cntk[]
		endk = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4));
		endl = p + ((l>>4) - ((l&~OCC_INTV_MASK)>>4));
		for (x = 0; p < endk; ++p) x += __occ_aux4(bwt, *p);
		y = x;
		tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
		x += __occ_aux4(bwt, tmp) - (~k&15);
		// calculate cntl[] and finalize cntk[]
		for (; p < endl; ++p) y += __occ_aux4(bwt, *p);
		tmp = *p & ~((1U<<((~l&15)<<1)) - 1);
		y += __occ_aux4(bwt, tmp) - (~l&15);
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
}

static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t k) // compute inverse CSA
{
	bwtint_t x = k - (k > bwt->primary);
	x = bwt_B0(bwt, x);
	x = bwt->L2[x] + bwt_occ(bwt, k, x);
	return k == bwt->primary? 0 : x;
}

bwtint_t bwt_sa(bwtint_t k)
{
	bwtint_t sa = 0, mask = Refbwt->sa_intv - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(Refbwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + Refbwt->sa[k/Refbwt->sa_intv];
}

bwtSearchResult_t BWT_Search(string& seq, int start, int stop)
{
	uint8_t nt;
	int i, pos, p;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;

	p = (int)nst_nt4_table[(int)seq[start]];
	ik.x[0] = Refbwt->L2[p] + 1;
	ik.x[1] = Refbwt->L2[3 - p] + 1;
	ik.x[2] = Refbwt->L2[p + 1] - Refbwt->L2[p];

	bwtSearchResult.freq = 0; bwtSearchResult.len = 0; bwtSearchResult.LocArr = NULL;
	for (pos = start + 1; pos < stop; pos++)
	{
		if ((nt = nst_nt4_table[(int)seq[pos]]) > 3) break;// ambiguous base

		bwt_2occ4(Refbwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = Refbwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= Refbwt->primary && ik.x[1] + ik.x[2] - 1 >= Refbwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];

		i = 3 - nt;
		if (ok[i].x[2] == 0) break; // extension ends
		else ik = ok[i];
	}

	if ((bwtSearchResult.len = pos - start) < MinSeedLength) bwtSearchResult.freq = 0;
	else
	{
		if ((bwtSearchResult.freq = (int)ik.x[2]) <= MaxSeedFreq)
		{
			bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}
		else bwtSearchResult.freq = 0;
	}
	return bwtSearchResult;
}

FragPair_t Specific_BWT_Search(string& seq, int start, int stop, int64_t rPos1, int64_t rPos2)
{
	uint8_t nt;
	int64_t rPos;
	FragPair_t FraPair;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	int i, len, pos, p, freq;

	p = (int)nst_nt4_table[(int)seq[start]];
	ik.x[0] = Refbwt->L2[p] + 1;
	ik.x[1] = Refbwt->L2[3 - p] + 1;
	ik.x[2] = Refbwt->L2[p + 1] - Refbwt->L2[p];

	FraPair.qLen = FraPair.rLen = 0;
	for (pos = start + 1; pos < stop; pos++)
	{
		if ((nt = nst_nt4_table[(int)seq[pos]]) == 4) break;// ambiguous base

		bwt_2occ4(Refbwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = Refbwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= Refbwt->primary && ik.x[1] + ik.x[2] - 1 >= Refbwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];

		i = 3 - nt;
		if (ok[i].x[2] == 0) break; // extension ends
		else ik = ok[i];
	}
	if ((len = pos - start) < MinSeedLength) freq = 0;
	else
	{
		if ((freq = (int)ik.x[2]) <= MaxSeedFreq)
		{
			for (i = 0; i < freq; i++)
			{
				rPos = (int64_t)bwt_sa(ik.x[0] + i);
				if (rPos >= rPos1 && rPos < rPos2)
				{
					FraPair.qPos = start;
					FraPair.rPos = rPos;
					FraPair.PosDiff = rPos - start;
					FraPair.qLen = FraPair.rLen = len;
					ShowFragPair(FraPair);
					break;
				}
			}
		}
	}
	return FraPair;
}
