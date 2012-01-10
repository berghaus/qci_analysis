
AM_CPPFLAGS  = `root-config --cflags` -I ../inc
AM_LDFLAGS   = `root-config --libs`
## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = analysis # testPE

# --- build my programs ---
analysisdir=../bin
analysis_SOURCES = analysis.cxx PDF.cxx Likelihood.cxx PValueTest.cxx PseudoExperiment.cxx \
	AtlasStyle.cxx TestStatMonitor.cxx PDFMonitor.cxx ControlFrame.cxx ControlFrameDict.cxx \
	CertaintyLevel.cxx
analysis_CPPFLAGS = `root-config --cflags` -I ../inc -pthread -fopenmp
analysis_LDADD   = -lboost_program_options -lGui -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -lm -ldl -lMinuit2

#testPE_SOURCES = PseudoExperiment_test.cxx PseudoExperiment.cxx PDF.cxx
#testPE_LDADD   = -lboost_test_exec_monitor -lboost_unit_test_framework -lGui -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -lm -ldl -lMinuit2


