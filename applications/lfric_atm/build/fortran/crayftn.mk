##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the Cray Fortran compiler.
##############################################################################

$(info Project specials for Cray compiler)

export FFLAGS_UM_PHYSICS = -s real64

# The lfric_atm app defines an extra set of debug flags for
# fast-debug. For this compiler use the same as the full-debug
# settings
FFLAGS_FASTD_INIT         = $(FFLAGS_INIT) 
FFLAGS_FASTD_RUNTIME      = $(FFLAGS_RUNTIME)

#to try and ease compile time for CCE on EXZ
ifeq ($(shell expr ${CRAYFTN_VERSION} \>= 015000000), 1)
    %bl_imp_alg_mod_psy.o %bl_imp_alg_mod_psy.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %bl_imp_alg_mod_psy.o %bl_imp_alg_mod_psy.mod: private FFLAGS_DEBUG = -G0
    %ukca_emiss_mode_mod.o %ukca_emiss_mode_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %ukca_emiss_mode_mod.o %ukca_emiss_mode_mod.mod: private FFLAGS_DEBUG = -G0
    %bl_exp_alg_mod_psy.o %bl_exp_alg_mod_psy.mod: private FFLAGS_SAFE_OPTIMISATION = -O0 
    %bl_exp_alg_mod_psy.o %bl_exp_alg_mod_psy.mod: privateFFLAGS_DEBUG = -G0
    %init_aerosol_fields_alg_mod_psy.o init_aerosol_fields_alg_mod_psy.mod: private FFLAGS_SAFE_OPTIMISATION = -O0 
    %init_aerosol_fields_alg_mod_psy.o init_aerosol_fields_alg_mod_psy.mod: private FFLAGS_DEBUG = -G0
    %jules_extra_kernel_mod.o %jules_extra_kernel_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0 
    %jules_extra_kernel_mod.o %jules_extra_kernel_mod.mod: private FFLAGS_DEBUG = -G0
 endif
$(info $(FFLAGS_SAFE_OPTIMISATION))
$(info $(FFLAGS_DEBUG))
