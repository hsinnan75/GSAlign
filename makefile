.KEEP_STAT:

all: align main

align:
		make -C seq-align
		
main:		
		make -C src && mv src/GSAlign .

index:
		make -C BWT_Index && mv BWT_Index/bwt_index .
		
