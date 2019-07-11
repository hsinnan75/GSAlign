.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64 -msse4.1 -fPIC
LINK		= -lpthread -lz -lm -lstdc++
SOURCE		= main.cpp GetData.cpp GSAlign.cpp ProcessCandidateAlignment.cpp KmerAnalysis.cpp tools.cpp SeqVariant.cpp DotPloting.cpp bwt_index.cpp bwt_search.cpp ksw2_alignment.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o GSAlign $(LINK)
%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $<
clean:
		rm -f *.o *~
