AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
	-I$(MY_INSTALL)/include \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
  libAnaDST.la

pkginclude_HEADERS = \
	AnaEventGen.h \
	AnaPythiaEventGen.h \
	AnaSvtxVertex.h

libAnaDST_la_LDFLAGS = \
  -L$(libdir) \
	-I$(MY_INSTALL)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

libAnaDST_la_LIBADD = \
  -lfun4all \
  -lg4detectors \
  -lphg4hit \
  -lg4jets \
  -lg4eval \
  -lphool \
  -lgenfit2 \
  -lgenfit2exp \
  -lPHGenFit \
	-lphhepmc

if ! MAKEROOT6
ROOT5_DICTS = \
	AnaEventGenDict.cc \
	AnaPythiaEventGenDict.cc \
	AnaSvtxVertexDict.cc
endif

libAnaDST_la_SOURCES = \
	$(ROOT5_DICTS) \
  AnaEventGen.cc \
  AnaPythiaEventGen.cc \
	AnaSvtxVertex.cc

install-exec-hook:
	@[ $$OSTYPE = darwin ] && ln -sf $(DESTDIR)$(libdir)/libAnaDST.dylib $(DESTDIR)$(libdir)/libAnaDST.so || true

BUILT_SOURCES = \
  testexternals.C

testexternals_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libAnaDST.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating CINT dictionaries from class headers.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
