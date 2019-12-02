.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -I${PREFIX}/include -Wall -D NDEBUG -O3 -m64 -msse4.2 -mpopcnt -fPIC
LIB		= -lpthread -lstdc++ -lz -lm
SOURCE		= main.cpp GetData.cpp GSAlign.cpp ProcessCandidateAlignment.cpp KmerAnalysis.cpp tools.cpp SeqVariant.cpp DotPloting.cpp bwt_index.cpp bwt_search.cpp ksw2_alignment.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o GSAlign BWT_Index/libbwa.a $(LIB) 
clean:
		rm -f *.o *~
