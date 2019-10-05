.KEEP_STAT:

all: 		main bwt_index

main:
		$(MAKE) -C src
		cp -f src/GSAlign .

bwt_index:
		$(MAKE) -C src/BWT_Index
		cp -f src/BWT_Index/$@ .
		
clean:
		rm -f GSAlign bwt_index
		$(MAKE) clean -C src
		$(MAKE) clean -C src/BWT_Index
