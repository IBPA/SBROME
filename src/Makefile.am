bin_PROGRAMS 		= sbrome
AUTOMAKE_OPTIONS 	= foreign
sbrome_SOURCES 		= 	DesignProblem.cpp DesignProblem.h \
						DynamicGeneCircuit.cpp DynamicGeneCircuit.h \
						Framework.cpp Framework.h \
						GeneCircuitGraph.cpp GeneCircuitGraph.h \
						GRN.h main.cpp \
						Module.cpp Module.h \
						PartDatabase.cpp PartDatabase.h \
						xmlParser.cpp xmlParser.h tcl.h
sbrome_LDADD 		= GA/libGA.a METIS/libmetis/libMETIS.a METIS/GKlib/libGK.a SQLite/libsqlite.a VFLIB/libvf.a -ldl -lpthread
sbrome_CPPFLAGS		= -O0
SUBDIRS 			= GA METIS SQLite VFLIB
