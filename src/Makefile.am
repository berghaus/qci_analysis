
AM_CPPFLAGS  = `root-config --cflags` -I ../inc
AM_LDFLAGS   = -v -L /usr/lib/x86_64-linux-gnu -L/usr/local/lib/root
## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = analysis postProcess plotErrors pValueTest # testPE

# parallel analysis program
analysisdir = ../bin
analysis_SOURCES = analysis.cxx Prediction.cxx Likelihood.cxx PValueTest.cxx PseudoExperiment.cxx \
	AtlasStyle.cxx TestStatMonitor.cxx PredictionMonitor.cxx CertaintyLevel.cxx Effect.cxx
analysis_CPPFLAGS = `root-config --cflags` -I ../inc -pthread -fopenmp -I ../../applgrid/include/ `lhapdf-config --incdir`
analysis_LDADD = -lboost_program_options -lGui -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -lm -ldl -lMinuit2 -L../../applgrid/lib -lAPPLgrid -lfAPPLgrid -lLHAPDF

# post-processing/merging program
postProcessdir = ../bin
postProcess_SOURCES = postProcess.cxx Prediction.cxx TestStatMonitor.cxx Likelihood.cxx PseudoExperiment.cxx \
	PValueTest.cxx CertaintyLevel.cxx PostProcessCL.cxx AtlasStyle.cxx Effect.cxx
postProcess_CPPFLAGS = `root-config --cflags` -I ../inc -I ../../applgrid/include/ `lhapdf-config --incdir`
postProcess_LDADD = -lboost_program_options -lGui -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -lm -ldl -lMinuit2 -L../../applgrid/lib -lAPPLgrid -lfAPPLgrid -lLHAPDF


# tool to plot statistical errors on MC prediction program
plotErrorsdir = ../bin
plotErrors_SOURCES = plotErrors.cxx Prediction.cxx PseudoExperiment.cxx AtlasStyle.cxx PredictionMonitor.cxx Effect.cxx
plotErrors_CPPFLAGS = `root-config --cflags` -I ../inc  -I ../../applgrid/include/ `lhapdf-config --incdir` 
plotErrors_LDADD  = -lboost_program_options -lGui -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -lm -ldl -lMinuit2 -L../../applgrid/lib -lAPPLgrid -lfAPPLgrid -lLHAPDF

# --- build my programs ---
# p-value test only program
pValueTestdir = ../bin
pValueTest_SOURCES = pValueTest.cxx Prediction.cxx Likelihood.cxx PValueTest.cxx PseudoExperiment.cxx \
	AtlasStyle.cxx TestStatMonitor.cxx PredictionMonitor.cxx
pValueTest_CPPFLAGS = `root-config --cflags` -I ../inc
pValueTest_LDADD = -lboost_program_options -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -lm -ldl -lMinuit2


#testPE_SOURCES = PseudoExperiment_test.cxx PseudoExperiment.cxx Prediction.cxx
#testPE_LDADD   = -lboost_test_exec_monitor -lboost_unit_test_framework -lGui -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -lm -ldl -lMinuit2


