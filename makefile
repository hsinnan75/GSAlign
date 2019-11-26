.KEEP_STAT:

all: 		main bwt_index

bin:
		mkdir -p bin

main: bwt_index
		$(MAKE) -C src
		mv -f src/GSAlign bin

bwt_index: bin
		$(MAKE) -C src/BWT_Index
		mv -f src/BWT_Index/$@ bin
		
clean:
		rm -f bin/GSAlign bin/bwt_index
		$(MAKE) clean -C src
		$(MAKE) clean -C src/BWT_Index
