include Make.inc

# Switch for Brutus
JNI_DIR = /System/Library/Frameworks/JavaVM.framework/Headers/
#JNI_DIR = $(JAVA_HOME)/include/

# Switch for different OS
# Mac OS X
EXPMV = macosx_expomv
EXPMVBEAST = libexpoTree.jnilib

# Linux (configured for Brutus, using lapack and goto)
#EXPMV = brutus_expomv
#EXPMVBEAST = libexpoTree.so

# Depends on the structure of the BEAST 2 project
PHYLODYNAMICS_BUILD = ../../beast2-read-only/build/Phylodynamics/
BEAST2_BUILD = ../../beast2-read-only/build/Beast2-read-only/
LIB_DIR = $(PHYLODYNAMICS_BUILD)/lib


# ExactBDSIS.class should be produced by IntelliJ/Eclipse/Netbeans etc.
expotreelib: JavaExpoTree.c $(PHYLODYNAMICS_BUILD)/beast/phylodynamics/ExactBDSIS.class $(LIB_DIR)
	javah -classpath $(PHYLODYNAMICS_BUILD):$(BEAST2_BUILD) beast.phylodynamics.ExactBDSIS
	$(CC) -shared -o $(LIB_DIR)/$(EXPMVBEAST) JavaExpoTree.c -I$(JNI_DIR) -L. -l$(EXPMV) $(LDFLAGS) -lgsl -fpic

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

clean: 
	rm -f $(LIB_DIR)/*.jnilib
	rm -f $(LIB_DIR)/*.so
	rm -f beast_phylodynamics_*.h
	rmdir $(LIB_DIR)
