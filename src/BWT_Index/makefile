CC=		gcc
CFLAGS=		-I${PREFIX}/include -Wall -O3
LDFLAGS=	-L${PREFIX}/lib -lm -lz
AR=ar
LOBJS=		utils.o bwt.o bntseq.o QSufSort.o bwt_gen.o bwtindex.o
PROG=		bwt_index

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bwt_index:libbwa.a main.o
		$(CC) $(CFLAGS) main.o -o $@ -L. -lbwa $(LDFLAGS)

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f *.o $(PROG) *~ *.a
