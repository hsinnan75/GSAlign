LIBS_PATH=libs

#ifdef DEBUG
#	OPT = -O0 --debug -g -ggdb
#else
#	OPT = -O3
#endif

CFLAGS = -Wall -Wextra -std=c99 -O3
OBJFLAGS = -fPIC
LINKFLAGS = -lalign -lstrbuf -lpthread -lz

INCS=-I $(LIBS_PATH) -I src
LIBS=-L $(LIBS_PATH)/string_buffer -L src
LINK=-lalign -lstrbuf -lpthread -lz

# Compile and bundle all non-main files into library
SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:.c=.o)

all: src/libalign.a

# Build libraries only if they're downloaded
src/libalign.a: $(OBJS)
	[ -d libs/string_buffer ] && cd libs && $(MAKE)
	ar -csru src/libalign.a $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(OBJFLAGS) $(INCS) -c $< -o $@

#needleman_wunsch: tools/nw_cmdline.c src/libalign.a
#	$(CC) -o needleman_wunsch $(SRCS) $(TGTFLAGS) $(INCS) $(LIBS) tools/nw_cmdline.c $(LINKFLAGS)

clean:
	rm -rf src/*.o src/libalign.a
