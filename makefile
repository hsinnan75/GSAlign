.KEEP_STAT:

all: main index

main:		
		make -C src && mv src/GSAlign .

index:
		make -C BWT_Index && mv BWT_Index/bwt_index .
		
