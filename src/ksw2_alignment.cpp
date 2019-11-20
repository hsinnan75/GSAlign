#include "structure.h"
#include <emmintrin.h>
#include <smmintrin.h>

//copyright: Heng Li

#define KSW_NEG_INF -0x40000000
//int8_t mat[25] = { 2, -4, -4, -4, 0, -4, 2, -4, -4, 0, -4, -4, 2, -4, 0, -4, -4, -4, 2, 0, 0, 0, 0, 0, 0 };
int8_t mat[25] = { 1, -1, -4, -4, 0, -4, 2, -4, -4, 0, -4, -4, 2, -4, 0, -4, -4, -4, 2, 0, 0, 0, 0, 0, 0 };

typedef struct {
	uint32_t max;
	int max_q, max_t;      // max extension coordinate
	int mqe, mqe_t;        // max score when reaching the end of query
	int mte, mte_q;        // max score when reaching the end of target
	int score;             // max score reaching both ends; may be KSW_NEG_INF
} ksw_extz_t;

void ksw_reset_extz(ksw_extz_t *ez)
{
	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
	ez->max = 0, ez->score = ez->mqe = ez->mte = KSW_NEG_INF;
}

string ksw_backtrack(const uint8_t *p, const int *off, const int *off_end, int n_col, int i0, int j0)
{
	int i = i0, j = j0, r, state = 0;
	uint32_t tmp;
	string cigar;

	while (i >= 0 && j >= 0) { // at the beginning of the loop, _state_ tells us which state to check
		int force_state = -1;

		r = i + j;
		if (i < off[r]) force_state = 2;
		if (off_end && i > off_end[r]) force_state = 1;
		tmp = force_state < 0 ? p[r * n_col + i - off[r]] : 0;

		if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
		else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
		if (state == 0) state = tmp & 7; // TODO: probably this line can be merged into the "else if" line right above; not 100% sure
		if (force_state >= 0) state = force_state;
		if (state == 0)
		{
			cigar.push_back('M');
			--i, --j; // match
		}
		else if (state == 1 || state == 3)
		{
			cigar.push_back('D');
			--i; // deletion
		}
		else
		{
			cigar.push_back('I');
			--j; // insertion
		}
	}
	if (i >= 0)
	{
		cigar.append(string().assign(i + 1, 'D'));
	}
	if (j >= 0)
	{
		cigar.append(string().assign(j + 1, 'I'));
	}
	return cigar;
}

string ksw_extz2_sse(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, int8_t q, int8_t e, int w, ksw_extz_t *ez)
{
	string cigar;

#define __dp_code_block1 \
	z = _mm_add_epi8(_mm_load_si128(&s[t]), qe2_); \
	xt1 = _mm_load_si128(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(xt1, 15);                   /* tmp <- x[r-1][t+15] */ \
	xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); /* xt1 <- x[r-1][t-1..t+14] */ \
	x1_ = tmp; \
	vt1 = _mm_load_si128(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(vt1, 15);                   /* tmp <- v[r-1][t+15] */ \
	vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); /* vt1 <- v[r-1][t-1..t+14] */ \
	v1_ = tmp; \
	a = _mm_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
	ut = _mm_load_si128(&u[t]);                      /* ut <- u[t..t+15] */ \
	b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */

#define __dp_code_block2 \
	z = _mm_max_epu8(z, b);                          /* z = max(z, b); this works because both are non-negative */ \
	z = _mm_min_epu8(z, max_sc_); \
	_mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm_store_si128(&v[t], _mm_sub_epi8(z, ut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	z = _mm_sub_epi8(z, q_); \
	a = _mm_sub_epi8(a, z); \
	b = _mm_sub_epi8(b, z);

	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr;
	//int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	__m128i q_, qe2_, zero_, flag1_, flag2_, flag8_, flag16_, sc_mch_, sc_mis_, m1_, max_sc_;
	__m128i *u, *v, *x, *y, *s, *p = 0;

	ksw_reset_extz(ez);
	if (m <= 0 || qlen <= 0 || tlen <= 0) return cigar;

	zero_ = _mm_set1_epi8(0);
	q_ = _mm_set1_epi8(q);
	qe2_ = _mm_set1_epi8((q + e) * 2);
	flag1_ = _mm_set1_epi8(1);
	flag2_ = _mm_set1_epi8(2);
	flag8_ = _mm_set1_epi8(0x08);
	flag16_ = _mm_set1_epi8(0x10);
	sc_mch_ = _mm_set1_epi8(mat[0]);
	sc_mis_ = _mm_set1_epi8(mat[1]);
	m1_ = _mm_set1_epi8(m - 1); // wildcard
	max_sc_ = _mm_set1_epi8(mat[0] + (q + e) * 2);

	if (w < 0) w = tlen > qlen ? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = ((w + 1 < tlen ? (w + 1 < qlen ? w + 1 : qlen) : tlen) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;
	//for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
	//	max_sc = max_sc > mat[t] ? max_sc : mat[t];
	//	min_sc = min_sc < mat[t] ? min_sc : mat[t];
	//}
	//if (-min_sc > 2 * (q + e)) return cigar; // otherwise, we won't see any mismatches

	mem = (uint8_t*)calloc(tlen_ * 6 + qlen_ + 1, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, s = y + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;

	H = (int32_t*)malloc(tlen_ * 16 * 4);
	for (t = 0; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;

	mem2 = (uint8_t*)malloc(((qlen + tlen - 1) * n_col_ + 1) * 16);
	p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
	off = (int*)malloc((qlen + tlen - 1) * sizeof(int) * 2);
	off_end = off + qlen + tlen - 1;

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, v1;
		uint8_t *qrr = qr + (qlen - 1 - r), *u8 = (uint8_t*)u, *v8 = (uint8_t*)v;
		__m128i x1_, v1_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r - wr + 1) >> 1) st = (r - wr + 1) >> 1; // take the ceil
		if (en >(r + wl) >> 1) en = (r + wl) >> 1; // take the floor

		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en)
				x1 = ((uint8_t*)x)[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			else x1 = v1 = 0; // not calculated; set to zeros
		}
		else x1 = 0, v1 = r ? q : 0;
		if (en >= r) ((uint8_t*)y)[r] = 0, u8[r] = r ? q : 0;
		// loop fission: set scores first
		for (t = st0; t <= en0; t += 16) {
			__m128i sq, st, tmp, mask;
			sq = _mm_loadu_si128((__m128i*)&sf[t]);
			st = _mm_loadu_si128((__m128i*)&qrr[t]);
			mask = _mm_or_si128(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(st, m1_));
			tmp = _mm_cmpeq_epi8(sq, st);
			tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
			tmp = _mm_andnot_si128(mask, tmp);
			_mm_storeu_si128((__m128i*)((uint8_t*)s + t), tmp);
		}
		// core loop
		x1_ = _mm_cvtsi32_si128(x1);
		v1_ = _mm_cvtsi32_si128(v1);
		st_ = st / 16, en_ = en / 16;

		__m128i *pr = p + r * n_col_ - st_;
		off[r] = st, off_end[r] = en;
		for (t = st_; t <= en_; ++t) {
			__m128i d, z, a, b, xt1, vt1, ut, tmp;
			__dp_code_block1;
			d = _mm_and_si128(_mm_cmpgt_epi8(a, z), flag1_); // d = a > z? 1 : 0
			z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
			tmp = _mm_cmpgt_epi8(b, z);
			d = _mm_blendv_epi8(d, flag2_, tmp);             // d = b > z? 2 : d
			__dp_code_block2;
			tmp = _mm_cmpgt_epi8(a, zero_);
			_mm_store_si128(&x[t], _mm_and_si128(tmp, a));
			d = _mm_or_si128(d, _mm_and_si128(tmp, flag8_));  // d = a > 0? 0x08 : 0
			tmp = _mm_cmpgt_epi8(b, zero_);
			_mm_store_si128(&y[t], _mm_and_si128(tmp, b));
			d = _mm_or_si128(d, _mm_and_si128(tmp, flag16_)); // d = b > 0? 0x10 : 0
			_mm_store_si128(&pr[t], d);
		}
		int32_t max_H, max_t;
		// compute H[], max_H and max_t
		if (r > 0) {
			int32_t HH[4], tt[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
			__m128i max_H_, max_t_, qe_;
			max_H = H[en0] = en0 > 0 ? H[en0 - 1] + u8[en0] - qe : H[en0] + v8[en0] - qe; // special casing the last element
			max_t = en0;
			max_H_ = _mm_set1_epi32(max_H);
			max_t_ = _mm_set1_epi32(max_t);
			qe_ = _mm_set1_epi32(q + e);
			for (t = st0; t < en1; t += 4) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
				__m128i H1, tmp, t_;
				H1 = _mm_loadu_si128((__m128i*)&H[t]);
				t_ = _mm_setr_epi32(v8[t], v8[t + 1], v8[t + 2], v8[t + 3]);
				H1 = _mm_add_epi32(H1, t_);
				H1 = _mm_sub_epi32(H1, qe_);
				_mm_storeu_si128((__m128i*)&H[t], H1);
				t_ = _mm_set1_epi32(t);
				tmp = _mm_cmpgt_epi32(H1, max_H_);
				max_H_ = _mm_blendv_epi8(max_H_, H1, tmp);
				max_t_ = _mm_blendv_epi8(max_t_, t_, tmp);
			}
			_mm_storeu_si128((__m128i*)HH, max_H_);
			_mm_storeu_si128((__m128i*)tt, max_t_);
			for (i = 0; i < 4; ++i)
				if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
			for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
				H[t] += (int32_t)v8[t] - qe;
				if (H[t] > max_H)
					max_H = H[t], max_t = t;
			}
		}
		else H[0] = v8[0] - qe - qe, max_H = H[0], max_t = 0; // special casing r==0
															  // update ez
		if (en0 == tlen - 1 && H[en0] > ez->mte)
			ez->mte = H[en0], ez->mte_q = r - en;
		if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
			ez->mqe = H[st0], ez->mqe_t = st0;
		if (r == qlen + tlen - 2 && en0 == tlen - 1)
			ez->score = H[tlen - 1];

		last_st = st, last_en = en;
		}
	cigar = ksw_backtrack((uint8_t*)p, off, off_end, n_col_ * 16, tlen - 1, qlen - 1);
	//printf("cigar=%s\n", cigar.c_str());
	free(mem); free(H); free(mem2); free(off);
	mem = mem2 = NULL; off = NULL; H = NULL;

	return cigar;
}

void ksw2_alignment(int m, string& s1, int n, string& s2)
{
	ksw_extz_t ez;
	int i, p, len;
	string cigar, aln1, aln2;

	uint8_t *str1 = (uint8_t*)malloc(m), *str2 = (uint8_t*)malloc(n);
	for (i = 0; i < m; ++i) str1[i] = nst_nt4_table[(uint8_t)s1[i]];
	for (i = 0; i < n; ++i) str2[i] = nst_nt4_table[(uint8_t)s2[i]];

	cigar = ksw_extz2_sse(m, str1, n, str2, 5, 2, 1, -1, &ez);
	free(str1); free(str2); str1 = str2 = NULL;

	len = (int)cigar.length();
	for (p = 0, i = len - 1; i >= 0; i--,p++)
	{
		switch (cigar[i])
		{
		case 'D': s1.insert(s1.begin() + p, 1, '-'); break;
		case 'I': s2.insert(s2.begin() + p, 1, '-'); break;
		}
	}
}
