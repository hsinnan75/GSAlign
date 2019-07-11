.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64 -msse4.1 -fPIC
INCS		= -I ../seq-align/libs -I ../seq-align/src
LIBS		= -L ../seq-align/libs/string_buffer -L ../seq-align/src
LINK		= -lalign -lstrbuf -lpthread -lz -lm -lstdc++
SOURCE		= main.cpp GetData.cpp GSAlign.cpp ProcessCandidateAlignment.cpp KmerAnalysis.cpp nw_alignment.cpp tools.cpp SeqVariant.cpp DotPloting.cpp bwt_index.cpp bwt_search.cpp ksw2_alignment.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main

main:		$(OBJECT) ../seq-align/src/libalign.a
			$(Compiler) $(FLAGS) $(INCS) $(LIBS) $(OBJECT) -o GSAlign $(LINK)
%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) $(INCS) $(LIBS) -c $<
clean:
		rm -f *.o *~
