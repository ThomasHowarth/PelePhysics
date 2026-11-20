# AMReX
DIM          = 2
PRECISION    = DOUBLE
PROFILE      = FALSE
VERBOSE      = FALSE
DEBUG        = FALSE

# Compiler
COMP       = gnu
FCOMP      = gfortran
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE

# PelePhysics
TINY_PROFILE = FALSE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    ?= ../../..

PELE_USE_KLU = FALSE

# Needs a heterogeneous (gas + surface) mechanism
Eos_Model       = Fuego
Chemistry_Model = ptcombust
Transport_Model = Simple

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
