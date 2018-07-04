.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64
LIB		= -lm -lpthread
SOURCE		= main.cpp GetData.cpp GSAlign.cpp dupDetection.cpp KmerAnalysis.cpp nw_alignment.cpp tools.cpp bwt_index.cpp bwt_search.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

all:		main index

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o GSAlign $(LIB)
%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $<

index:
		make -C BWT_Index && mv BWT_Index/bwt_index .
		
clean:
		rm -f *.o *~

eva:		Evaluation.cpp
		$(Compiler) $(FLAGS) Evaluation.cpp -o eva

