##############################################################################
# (c) Crown copyright xMet Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
##############################################################################
# Various things specific to fast-debug lfricinputs when using the
# Cray Fortran compiler.
# To override flags specified in LFRic Core, use
# FFLAGS_SAFE_OPTIMISATION, to add an additional flag to what is specified
# in Core, use FFLAGS_EXTRA. If you need to add debugging flags to a module,
# modify FFLAGS_DEBUG for the target module.
# All flag modification should be `private` so they are not inherited by
# submodules.
##############################################################################

# ==========================================================================
# NON-DEBUG FLAGS
# ==========================================================================

# ==========================================================================
# DEBUG FLAGS
# ==========================================================================
