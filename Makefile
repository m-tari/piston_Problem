# Makefile for Wave

# use gcc compiler
CC   = gcc

############################
###      DIRECTORIES     ###
############################

# SDIR: directory containing the source files (.c)
SDIR   = src

# ODIR: directory containing the object files (.o)
ODIR   = obj

# HDIR: directory containing the header files (.h)
HDIR   = head

BINDIR = bin
LIBDIR = lib

# libraries
MKL_INCLUDE  = $(MKLROOT)/include
MKL_LIB      = $(MKLROOT)/lib/intel64

L_GENERAL    = -lpthread -lm -ldl -Wall
L_MKL        = -L$(MKL_LIB) -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LIBS_PISTON  = $(L_MKL) $(L_GENERAL)

I_HEAD       = -I$(HDIR)

EXE          = $(BINDIR)/piston

# list of flags to pass to the compilation command:
CFLAGS_WAVE  = -m64 -I$(MKL_INCLUDE) -I$(HDIR)

TARGETS      =  $(EXE)

############################
###        FILES         ###
############################

# functions library
SRC     = $(wildcard $(SDIR)/*.c)
HEAD    = $(wildcard $(HDIR)/*.h)
OBJ     = $(subst $(SDIR),$(ODIR),$(patsubst %.c,%.o,$(SRC))) 


############################
###   RULES DEFINITION   ###
############################

# Target to make
# Call "make piston" or "make" to use
.PHONY: piston
wave: CFLAGS = $(CFLAGS_WAVE)
wave: LIBS   = $(LIBS_PISTON)
wave: $(TARGETS)

# compilation
$(ODIR)/%.o: $(SDIR)/%.c  $(HEAD)
	$(CC) $(CFLAGS_WAVE) $(I_HEAD) -c -o $@ $<

# object files are prerequisites
$(EXE): $(OBJ) $(COM_OBJ)
	@$(CC) $(CFLAGS) -o $@ $(OBJ) $(LIBS)
	@echo "Compiled piston executable"


clean:
	rm -f $(ODIR)/*.o $(BINDIR)/piston	