
# compilation and linkage flags.
CC=gcc
CFLAGS=-Wall -O2 -D_GNU_SOURCE
LDFLAGS=-lm

# installation flags.
DESTDIR=
PREFIX=$(DESTDIR)/usr
BINDIR=$(PREFIX)/bin
MANDIR=$(PREFIX)/man/man1

# output binaries.
BIN=pca-ellipses pca-ellipsoids pca-dendrogram pca-bootstrap pca-distances
BIN+= pca-overlap pca-stats pca-rand pca-maps
BINS=$(addprefix bin/,$(BIN))

# output manual pages.
MANS=$(addprefix man/,$(addsuffix .1,$(BINS)))

# intermediate 'non-binary' object code files.
LIBOBJ=pca-utils.o pca-utils-math.o pca-utils-rand.o pca-utils-stat.o
LIBOBJ+= pca-utils-list.o pca-utils-dist.o pca-utils-tree.o pca-utils-draw.o
LIBOBJ+= pca-utils-getopt.o

# full object file paths.
LIBOBJS=$(addprefix src/,$(LIBOBJ))
BINOBJS=$(addsuffix .o,$(BINS))

# all object code files.
OBJS=$(LIBOBJS) $(BINOBJS)

# specification of compilation suffixes.
.SUFFIXES: .c .o

# default target.
all: $(OBJS) $(BINS)

# source compilation target.
.c.o:
	@echo " CC $^"
	@$(CC) $(CFLAGS) -I./src -c $^ -o $@

# binary linkage target.
$(BINS): $(OBJS)
	@echo " LD $@"
	@$(CC) $(CFLAGS) $(LIBOBJS) $(addsuffix .o,$@) $(LDFLAGS) -o $@

# installation target.
install: all
	@for dir in $(BINDIR) $(MANDIR); do \
	   if [ ! -d $${dir} ]; then \
	      echo " MKDIR $${dir}"; \
	      install -d -m 755 -o root -g root $${dir}; \
	   fi; \
	 done
	@for bin in $(BINS); do \
	   echo " INSTALL $${bin}"; \
	   install -m 755 -o root -g root $${bin} $(BINDIR)/; \
	 done
	@for man in $(MANS); do \
	   echo " INSTALL $${man}"; \
	   install -m 644 -o root -g root $${man} $(MANDIR)/; \
	 done

# intermediate file removal target.
clean:
	@echo " CLEAN"
	@rm -f $(BINS) $(OBJS)

# recompilation target.
again: clean all

# fixme code search target.
fixme:
	@grep -RHni fixme *.[ch1] || \
	 echo "No FIXMEs found. Good work, stud."

# line counting target.
lines:
	@wc -l *.[ch1]

