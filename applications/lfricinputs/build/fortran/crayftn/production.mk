##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
##############################################################################
# Various things specific to production lfric_atm when using the
# Cray Fortran compiler.
##############################################################################

# lfric2um_initialise_um_mod calculates levels to include in fieldsfile.
# Rigorous setting required to prevent differences from fast-debug.
%lfric2um_initialise_um_mod.o: private FFLAGS_RISKY_OPTIMISATION = -O3 -hipa3 -hflex_mp=rigorous
