.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64 -msse4.1
LIB		= -lm -lpthread
SOURCE		= main.cpp GetData.cpp GSAlign.cpp nw_alignment.cpp ksw2_alignment.cpp edlib_alignment.cpp tools.cpp bwt_index.cpp bwt_search.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o GSAlign $(LIB)
%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $<
clean:
		rm -f *.o *~

eva:		Evaluation.cpp
		$(Compiler) $(FLAGS) Evaluation.cpp -o eva

