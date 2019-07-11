#!/bin/bash
#test alignment
echo
echo "Test -- Align two E.coli strains with 4 threads"
echo "Command=./GSAlign -t 4 -r test/ecoli.fa -q test/ecoli.mut -o test/output"
echo
./GSAlign -t 4 -r test/ecoli.fa -q test/ecoli.mut -o test/output
echo
echo "[End of test]"
