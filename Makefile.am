
AM_CPPFLAGS  = `root-config --cflags` -I inc
AM_LDFLAGS   = `root-config --libs`
## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = analysis

# --- build my programs ---
analysis_SOURCES = src/analysis.cxx src/PDF.cxx src/Likelihood.cxx src/PValueTest.cxx src/PseudoExperiment.cxx \
	src/AtlasStyle.cxx src/TestStatMonitor.cxx src/PDFMonitor.cxx
analysis_LDADD   = -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -lm -ldl -lMinuit2

