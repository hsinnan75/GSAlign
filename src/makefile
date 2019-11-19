.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -I${PREFIX}/include -L${PREFIX}/lib -Wall -D NDEBUG -O3 -m64 -msse4.2 -mpopcnt -fPIC
LIB		= -lm -lpthread -lstdc++
SOURCE		= main.cpp GetData.cpp GSAlign.cpp ProcessCandidateAlignment.cpp KmerAnalysis.cpp tools.cpp SeqVariant.cpp DotPloting.cpp bwt_index.cpp bwt_search.cpp ksw2_alignment.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o GSAlign $(LIB) 
clean:
		rm -f *.o *~
