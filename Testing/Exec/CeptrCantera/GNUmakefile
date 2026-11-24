# Builds the evaluator of a CEPTR-generated mechanism for the comparison with
# Cantera (Support/ceptr/tests/test_cantera.py):
#
#   make MECH_DIR=<directory with mechanism.H and mechanism.cpp>
#
# Only the AMReX headers are used: nothing is linked but the mechanism, so
# neither AMReX nor SUNDIALS need to be built. Everything is built in
# MECH_DIR (evaluate.ex), so that a mechanism generated in a temporary
# directory leaves no files here.

ifndef MECH_DIR
  $(error Set MECH_DIR to the directory of the generated mechanism)
endif

# Set here, not taken from the environment (a PELE_PHYSICS_HOME or CXXFLAGS
# exported for another build), but they can be set on the command line
EVALUATOR_DIR    := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
PELE_PHYSICS_HOME = $(abspath $(EVALUATOR_DIR)/../../..)
AMREX_HOME        = $(PELE_PHYSICS_HOME)/Submodules/amrex
PYTHON           ?= python3

# No fast-math: the comparison is to round-off
CXXFLAGS = -O2
# mechanism.cpp only provides the species names, the reaction map and CKKFKR
# here, but its sparsity routines instantiate the analytical Jacobian, which
# takes minutes to optimize
MECH_CXXFLAGS = -O0

PP_SRC = $(PELE_PHYSICS_HOME)/Source
AMREX_DEFINES = -DAMREX_SPACEDIM=1 -DBL_SPACEDIM=1 -DAMREX_USE_DOUBLE \
                -DBL_USE_DOUBLE -DAMREX_GPU_MAX_THREADS=0
DEFINES  = $(AMREX_DEFINES) -DUSE_FUEGO_EOS -DUSE_SIMPLE_TRANSPORT
INCLUDES = -I$(MECH_DIR) -I$(PP_SRC) -I$(PP_SRC)/Eos -I$(PP_SRC)/Transport \
           -I$(PP_SRC)/Utility/BlackBoxFunction -I$(PP_SRC)/Utility/Utilities \
           -I$(AMREX_HOME)/Src/Base -I$(AMREX_HOME)/Src/Base/Parser
COMPILE  = $(CXX) -std=c++20 $(DEFINES) $(INCLUDES)

.PHONY: all clean
all: $(MECH_DIR)/evaluate.ex

# AMReX's GNU make build generates it the same way
$(MECH_DIR)/AMReX_Config.H:
	$(PYTHON) $(AMREX_HOME)/Tools/libamrex/mkconfig.py \
	  --defines="$(AMREX_DEFINES)" --comp=gnu --allow-different-compiler=TRUE > $@

$(MECH_DIR)/main.o: $(EVALUATOR_DIR)main.cpp $(MECH_DIR)/mechanism.H $(MECH_DIR)/AMReX_Config.H
	$(COMPILE) $(CXXFLAGS) -c $< -o $@

$(MECH_DIR)/mechanism.o: $(MECH_DIR)/mechanism.cpp $(MECH_DIR)/mechanism.H $(MECH_DIR)/AMReX_Config.H
	$(COMPILE) $(MECH_CXXFLAGS) -c $< -o $@

$(MECH_DIR)/evaluate.ex: $(MECH_DIR)/main.o $(MECH_DIR)/mechanism.o
	$(CXX) -o $@ $^

clean:
	$(RM) $(addprefix $(MECH_DIR)/,evaluate.ex main.o mechanism.o AMReX_Config.H)
