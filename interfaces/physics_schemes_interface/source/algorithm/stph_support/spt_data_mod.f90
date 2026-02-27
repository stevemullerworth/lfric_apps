!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module spt_data_mod

  use constants_mod,                 only : r_def, i_def
  use key_value_collection_mod,      only : key_value_collection_type
  use key_value_mod,                 only : abstract_value_type

  use extrusion_config_mod,          only : planet_radius
  use stochastic_physics_config_mod, only : spt_decorrelation_time, stph_n_max

  implicit none

  private
  public :: get_spt_data_from_collection

  !> Collection of data for Stochastic Physics SPT option

  type, extends(abstract_value_type), public :: spt_data_type

    private

    ! Spectral coefficients and power law
    real(r_def), allocatable :: spt_spectral_coeffc(:)
    real(r_def), allocatable :: spt_spectral_coeffs(:)
    real(r_def), allocatable :: spt_power_law(:)
    integer(i_def)           :: spectral_dimensions

  contains
    private

    procedure, public :: initialise
    procedure, public :: get_spt_spectral_coeffc
    procedure, public :: get_spt_spectral_coeffs
    procedure, public :: get_spt_power_law
    procedure, public :: get_stph_spectral_dim


    procedure, public :: set_spt_spectral_coeffc
    procedure, public :: set_spt_spectral_coeffs
    procedure, public :: set_spt_power_law

  end type spt_data_type

contains

  subroutine initialise(self, stph_spectral_dim, dt)

    class(spt_data_type), intent(inout) :: self

    integer(i_def), intent(in) :: stph_spectral_dim
    real(r_def),    intent(in) :: dt

    ! Scalars for power law & 1AR
    real(kind=r_def) :: beta, gamma, Lbeta, sigma, spt_alpha
    ! Physics scalars
    real(kind=r_def) :: mlcrcp
    ! iterators in for loops
    integer(i_def) :: n,n_row, m


    self%spectral_dimensions = stph_spectral_dim

    allocate(self%spt_spectral_coeffc(stph_spectral_dim))
    allocate(self%spt_spectral_coeffs(stph_spectral_dim))
    allocate(self%spt_power_law(stph_spectral_dim))

    ! set them to zero (invokes don't work for non-fields types)
    self%spt_spectral_coeffc = 0.0_r_def
    self%spt_spectral_coeffs = 0.0_r_def
    self%spt_power_law= 0.0_r_def

    !!!!!! 1.a compute power law

    ! compute alpha for temporal decorrelation,
    spt_alpha=1-exp(-dt/spt_decorrelation_time)

    ! Compute spatial decorrelation for 500km
    ! Maybe move 500 to namelist (hardwired in UM)
    Lbeta = 500e3_r_def !Original value eq 500e3
    beta  = 0.5_r_def*(Lbeta/planet_radius)**2.0_r_def

    ! Compute Tau, add zonal one by hand
    ! (as done in UM)
    gamma=1.0_r_def/(2.0_r_def-spt_alpha)
    do n=1,stph_n_max
      gamma=gamma+(2*n + 1)/(1-0.5_r_def*spt_alpha)*exp(-2*beta*n*(n+1))
    end do

    ! Compute sigma
    sigma=1.0_r_def/sqrt(12.0_r_def)

    ! Define the power law over all stph_spectral_dim
    ! for simplicity (thus repeating same values)
    ! on the single line
    n_row = 0
    do n = 1, stph_n_max
      n_row = n_row + n
      do m = 0,n
        self%spt_power_law( n_row + m )=exp(-beta*n*(n+1))/(sigma*sqrt(gamma))
      end do
    end do

  end subroutine initialise

  function get_spt_spectral_coeffc(self) result(spt_spectral_coeffc)

    class(spt_data_type), intent(inout) :: self

    real(r_def) :: spt_spectral_coeffc(self%spectral_dimensions)

    spt_spectral_coeffc = self%spt_spectral_coeffc

  end function get_spt_spectral_coeffc

  function get_spt_spectral_coeffs(self) result(spt_spectral_coeffs)

    class(spt_data_type), intent(inout) :: self

    real(r_def) :: spt_spectral_coeffs(self%spectral_dimensions)

    spt_spectral_coeffs = self%spt_spectral_coeffs

  end function get_spt_spectral_coeffs

  function get_spt_power_law(self) result(spt_power_law)

    class(spt_data_type), intent(inout) :: self

    real(r_def) :: spt_power_law(self%spectral_dimensions)

    spt_power_law = self%spt_power_law

  end function get_spt_power_law

  function get_stph_spectral_dim(self) result(stph_spectral_dim)

    class(spt_data_type), intent(inout) :: self

    integer(i_def) :: stph_spectral_dim

    stph_spectral_dim = self%spectral_dimensions

  end function get_stph_spectral_dim

  subroutine set_spt_spectral_coeffc(self, spectral_coeffc)

    class(spt_data_type), intent(inout) :: self

    real(r_def), intent(in) :: spectral_coeffc(:)

    self%spt_spectral_coeffc = spectral_coeffc

  end subroutine set_spt_spectral_coeffc

  subroutine set_spt_spectral_coeffs(self, spectral_coeffs)

    class(spt_data_type), intent(inout) :: self

    real(r_def), intent(in) :: spectral_coeffs(:)

    self%spt_spectral_coeffs = spectral_coeffs

  end subroutine set_spt_spectral_coeffs

  subroutine set_spt_power_law(self, power_law)

    class(spt_data_type), intent(inout) :: self

    real(r_def), intent(in) :: power_law(:)

    self%spt_power_law = power_law

  end subroutine set_spt_power_law

  !-----------------------------------------------------------------------------
  !> @brief Helper function to extract a spt_data_type object from a
  !>        key-value collection
  !> @param[in] collection The key-value collection to extract from
  !> @param[in] name       The name of the timestep method object to extract
  !> @return    timestep_method   The requested timestep object
  function get_spt_data_from_collection(collection, name) &
                                               result(spt_data)

  implicit none

    type(key_value_collection_type), intent(in) :: collection
    character(*),                    intent(in) :: name

    class(spt_data_type), pointer :: spt_data

    class(abstract_value_type), pointer :: abstract_value

    call collection%get_value(trim(name), abstract_value)
    select type(abstract_value)
      type is (spt_data_type)
        spt_data => abstract_value
    end select

  end function get_spt_data_from_collection

end module spt_data_mod
