#include <stdio.h>
#include <string.h>

int bwa_idx_build(const char *fa, const char *prefix);
int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s Ref_File[ex. ref.fa] Prefix[ex. MyRef]\n", argv[0]);
	}
	else bwa_idx_build(argv[1], argv[2]);

	return 0;
}
