.KEEP_STAT:

all: main

CXX		= g++
CFLAGS		= -Wall -D NDEBUG -O3 -m64 -msse4.1 -fPIC -I${PREFIX}/include 
LDFLAGS		= -lpthread -lz -lm -lstdc++ -L${PREFIX}/lib 
SOURCE		= main.cpp GetData.cpp GSAlign.cpp ProcessCandidateAlignment.cpp KmerAnalysis.cpp tools.cpp SeqVariant.cpp DotPloting.cpp bwt_index.cpp bwt_search.cpp ksw2_alignment.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main

main:		$(OBJECT)
			$(CXX) $(CFLAGS) $(OBJECT) -o GSAlign $(LDFLAGS)
%.o:		%.cpp $(HEADER)
			$(CXX) $(CFLAGS) -c $<
clean:
		rm -f *.o *~
