osname := $(shell uname -s)
cpuarch := $(shell uname -m)

CC-CYGWIN_NT-6.1-WOW64-i686 := c:/mingw/bin/mingw32-gcc
CC-Linux-x86_64 := gcc

CC     := $(CC-$(osname)-$(cpuarch))
CFLAGS  = -shared -fPIC
LDFLAGS =
SOURCES = find_index.c smallest_index.c filter.c
OBJECTS = $(SOURCES:.c=.o)

all: $(SOURCES) $(OBJECTS)

clean:
	rm *.dylib

.c.o:
	$(CC) $(CFLAGS) $< -o $@
	mv $@ $@.dylib
