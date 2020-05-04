#include "structure.h"

vector<Variant_t> VarVec;
int iSNV, iInsertion, iDeletion;

bool CompByVariantPos(const Variant_t& p1, const Variant_t& p2)
{
	if (p1.chr_idx == p2.chr_idx) return p1.pos < p2.pos;
	else return p1.chr_idx < p2.chr_idx;
}

void VariantIdentification()
{
	int64_t rPos;
	Variant_t Variant;
	string frag1, frag2;
	int i, qPos, aln_len, ind_len;
	vector<FragPair_t>::iterator FragPairIter;

	for (vector<AlnBlock_t>::iterator ABiter = AlnBlockVec.begin(); ABiter != AlnBlockVec.end(); ABiter++)
	{
		if (ABiter->bDup) continue;

		Variant.chr_idx = ABiter->coor.ChromosomeIdx;
		Variant.query_idx = QueryChrIdx;

		for (FragPairIter = ABiter->FragPairVec.begin(); FragPairIter != ABiter->FragPairVec.end(); FragPairIter++)
		{
			if (!FragPairIter->bSeed)
			{
				if (FragPairIter->qLen == 0 && FragPairIter->rLen == 0) continue;
				else if (FragPairIter->qLen == 0) // delete
				{
					iDeletion++;
					frag1.resize(FragPairIter->rLen + 1); strncpy((char*)frag1.c_str(), RefSequence + FragPairIter->rPos - 1, FragPairIter->rLen + 1);
					Variant.type = 2;
					Variant.pos = GenCoordinateInfo(FragPairIter->rPos - 1).gPos;
					Variant.ref_frag = frag1;
					Variant.alt_frag.resize(1); Variant.alt_frag[0] = QueryChrVec[QueryChrIdx].seq[FragPairIter->qPos - 1];
					VarVec.push_back(Variant);
				}
				else if (FragPairIter->rLen == 0) // insert
				{
					iInsertion++;
					frag2.resize(FragPairIter->qLen + 1); strncpy((char*)frag2.c_str(), QueryChrVec[QueryChrIdx].seq.c_str() + FragPairIter->qPos - 1, FragPairIter->qLen + 1);
					Variant.type = 1;
					Variant.pos = GenCoordinateInfo(FragPairIter->rPos - 1).gPos;
					Variant.ref_frag.resize(1); Variant.ref_frag[0] = RefSequence[FragPairIter->rPos - 1];
					Variant.alt_frag = frag2;
					VarVec.push_back(Variant);
				}
				else if (FragPairIter->qLen == 1 && FragPairIter->rLen == 1) // substitution
				{
					if (nst_nt4_table[(int)FragPairIter->aln1[0]] != nst_nt4_table[(int)FragPairIter->aln2[0]] && nst_nt4_table[(int)FragPairIter->aln2[0]] != 4)
					{
						iSNV++;
						Variant.type = 0;
						Variant.pos = GenCoordinateInfo(FragPairIter->rPos).gPos;
						Variant.ref_frag = FragPairIter->aln1;
						Variant.alt_frag = FragPairIter->aln2;
						VarVec.push_back(Variant);
					}
				}
				else
				{
					//fprintf(stdout, "ref=%s\nqry=%s\n", FragPairIter->aln2.c_str(), FragPairIter->aln1.c_str());
					for (aln_len = (int)FragPairIter->aln1.length(), rPos = FragPairIter->rPos, qPos = FragPairIter->qPos, i = 0; i < aln_len; i++)
					{
						if (FragPairIter->aln1[i] == '-') // insert
						{
							iInsertion++;
							ind_len = 1; while (FragPairIter->aln1[i + ind_len] == '-') ind_len++;
							frag2 = QueryChrVec[QueryChrIdx].seq.substr(qPos - 1, ind_len + 1);
							Variant.type = 1;
							Variant.pos = GenCoordinateInfo(rPos - 1).gPos;
							Variant.ref_frag.resize(1); Variant.ref_frag[0] = frag2[0];
							Variant.alt_frag = frag2;
							VarVec.push_back(Variant);

							qPos += ind_len; i += ind_len - 1;
						}
						else if (FragPairIter->aln2[i] == '-') // delete
						{
							iDeletion++;
							ind_len = 1; while (FragPairIter->aln2[i + ind_len] == '-') ind_len++;
							frag1.resize(ind_len + 2); frag1[ind_len + 1] = '\0';
							strncpy((char*)frag1.c_str(), RefSequence + rPos - 1, ind_len + 1);

							Variant.type = 2;
							Variant.pos = GenCoordinateInfo(rPos - 1).gPos;
							Variant.ref_frag = frag1;
							Variant.alt_frag.resize(1); Variant.alt_frag[0] = frag1[0];
							VarVec.push_back(Variant);

							rPos += ind_len; i += ind_len - 1;
						}
						else if (nst_nt4_table[(int)FragPairIter->aln1[i]] != nst_nt4_table[(int)FragPairIter->aln2[i]])// substitute
						{
							if (nst_nt4_table[(int)FragPairIter->aln2[i]] != 4)
							{
								iSNV++;
								Variant.type = 0;
								Variant.pos = GenCoordinateInfo(rPos).gPos;
								Variant.ref_frag.resize(1); Variant.ref_frag[0] = FragPairIter->aln1[i];
								Variant.alt_frag.resize(1); Variant.alt_frag[0] = FragPairIter->aln2[i];
								VarVec.push_back(Variant);
							}
							rPos++; qPos++;
						}
						else
						{
							rPos++; qPos++;
						}
					}
				}
			}
		}
	}
}

void OutputSequenceVariants()
{
	FILE *outFile;
	const char *MutType[3] = { "SUBSTITUTE", "INSERT", "DELETE" };

	sort(VarVec.begin(), VarVec.end(), CompByVariantPos);

	iSNV = iInsertion = iDeletion = 0;

	outFile = fopen(vcfFileName, "w");
	fprintf(outFile, "##fileformat=VCFv4.1\n");
	if (IndexFileName != NULL) fprintf(outFile, "##reference=%s\n", IndexFileName);
	else fprintf(outFile, "##reference=%s\n", RefSeqFileName);
	fprintf(outFile, "##source=GSAlign %s\n", VersionStr);
	fprintf(outFile, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"The type of allele, either SUBSTITUTE, INSERT, or DELETE.\">\n");
	for (int i = 0; i < iChromsomeNum; i++) fprintf(outFile, "##contig=<ID=%s,length=%d>\n", ChromosomeVec[i].name, ChromosomeVec[i].len);
	fprintf(outFile, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n");
	for (vector<Variant_t>::iterator iter = VarVec.begin(); iter != VarVec.end(); iter++)
	{
		fprintf(outFile, "%s\t%d\t.\t%s\t%s\t100\t*\tTYPE=%s\n", ChromosomeVec[iter->chr_idx].name, iter->pos, (char*)iter->ref_frag.c_str(), (char*)iter->alt_frag.c_str(), MutType[iter->type]);
	}
	std::fclose(outFile);
}
