#include "structure.h"

void nw_alignment(string& s1, string& s2)
{
	nw_aligner_t *nw = needleman_wunsch_new();
	alignment_t *result = alignment_create(256);

	needleman_wunsch_align(s1.c_str(), s2.c_str(), &scoring, nw, result);

	s1 = result->result_a; s2 = result->result_b;

	needleman_wunsch_free(nw); alignment_free(result);
}
