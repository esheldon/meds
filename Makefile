CC=gcc
LD=gcc
AR=ar

prefix := /usr/local

CFLAGS=-std=gnu99 -Wall -Werror -O2
ARFLAGS=rcs

SRCDIR=./src

LIB_SOURCES = $(SRCDIR)/meds.c

TEST_SOURCES = $(SRCDIR)/test.c
TEST_SPEED_SOURCES = $(SRCDIR)/test-speed.c
GET_CUTOUT_SOURCES= $(SRCDIR)/meds-get-cutout.c
MAKE_INPUT_SOURCES = $(SRCDIR)/make-meds-input.c

ALL_SOURCES = $(LIB_SOURCES) \
			  $(TEST_SOURCES) \
			  $(TEST_SPEED_SOURCES) \
			  $(MAKE_INPUT_SOURCES) \
			  $(GET_CUTOUT_SOURCES)

LIB_OBJECTS=$(patsubst %.c,%.o,$(LIB_SOURCES)) 
TEST_OBJECTS=$(patsubst %.c,%.o,$(TEST_SOURCES)) 
TEST_SPEED_OBJECTS=$(patsubst %.c,%.o,$(TEST_SPEED_SOURCES)) 
MAKE_INPUT_OBJECTS=$(patsubst %.c,%.o,$(MAKE_INPUT_SOURCES)) 
GET_CUTOUT_OBJECTS=$(patsubst %.c,%.o,$(GET_CUTOUT_SOURCES)) 

# these installed
LIB_BASE=libmeds.a
MI_BASE=make-meds-input
GET_CUTOUT_BASE=meds-get-cutout
LIB = $(SRCDIR)/$(LIB_BASE)
HEADER = $(SRCDIR)/meds.h
MAKE_INPUT_PROG = $(SRCDIR)/$(MI_BASE)
GET_CUTOUT_PROG=$(SRCDIR)/$(GET_CUTOUT_BASE)

# note order
TEST_LINKFLAGS=-L$(SRCDIR) -lmeds -lcfitsio -lm
MAKE_INPUT_LINKFLAGS=-lcfitsio -lm
GET_CUTOUT_LINKFLAGS=-L$(SRCDIR) -lmeds -lcfitsio -lm


# just for tests
TEST_PROG = $(SRCDIR)/test
TEST_SPEED_PROG = $(SRCDIR)/test-speed

PROGS=$(TEST_PROG) $(TEST_SPEED_PROG) $(MAKE_INPUT_PROG)
DEPFILE=$(SRCDIR)/.depend

default: all

depend: $(DEPFILE)

$(DEPFILE): $(ALL_SOURCES)
	$(CC) $(CFLAGS) -MM $^ > $(DEPFILE);

-include $(DEPFILE)


install: $(LIB) $(MAKE_INPUT_PROG)
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/bin
	cp $(LIB) $(prefix)/lib/
	cp $(HEADER) $(prefix)/include/

	cp $(MAKE_INPUT_PROG) $(prefix)/bin/$(MI_BASE)
	chmod a+x $(prefix)/bin/$(MI_BASE)

	cp $(GET_CUTOUT_PROG) $(prefix)/bin
	chmod a+x $(prefix)/bin/$(GET_CUTOUT_BASE)


all: $(TEST_PROG) $(TEST_SPEED_PROG) $(MAKE_INPUT_PROG) $(GET_CUTOUT_PROG)

lib: $(LIB)
	
$(LIB): $(LIB_OBJECTS)
	$(AR) $(ARFLAGS) $(LIB) $(LIB_OBJECTS)

$(TEST_PROG): $(LIB) $(TEST_OBJECTS)
	$(LD) -o $@  $(TEST_OBJECTS) $(TEST_LINKFLAGS)

$(TEST_SPEED_PROG): $(LIB) $(TEST_SPEED_OBJECTS)
	$(LD) -o $@  $(TEST_SPEED_OBJECTS) $(TEST_LINKFLAGS)

$(MAKE_INPUT_PROG): $(MAKE_INPUT_OBJECTS)
	$(LD) -o $@ $(MAKE_INPUT_OBJECTS) $(MAKE_INPUT_LINKFLAGS) 

$(GET_CUTOUT_PROG): $(GET_CUTOUT_OBJECTS)
	$(LD) -o $@ $(GET_CUTOUT_OBJECTS) $(GET_CUTOUT_LINKFLAGS) 

clean:
	rm -f $(SRCDIR)/*.o $(LIB) $(PROGS) $(DEPFILE)
