!
! Module providing API to library routines
!
! TODO(Alex) Things to Add:
!> Dielectric constant, and if use screening
!> tam-dancoff approximation
!> Polarisation of the light 
!> Spectral or 2p-eigenvalues? Only spectrum atm. Spectral broadening 
!> Solver, although only CG atm
!> Neighbour cut-off. Have the library internally determine number of neighbours
!> per atom from a distance matrix 


!TODO(Alex) Rename as bse_api
module api

  use precision, only : wp
  use mpi_bse,   only : mpi_env_type
  use atoms,     only : atom_type, construct_molecule

  implicit none
  private 

  
  !> A BSE calculation: Should contain all data types and methods to set them
  !> for a BSE calculation. Should also contain run method
  type bse_type
     private
     !> MPI environment variables
     type(mpi_env_type) :: mpi_env
     !> Eigenvalues and vectors
     type(eigen_type) :: eigen
     !> Molecular system
     type(atom_type), allocatable :: molecule(:)

   contains
     ! Each of these must be called prior to running a calculation
     procedure :: set_mpi_env => set_bse_mpi_environment     
     procedure :: set_eigenstates => set_bse_eigenstates     !TODO(Alex) Write me
     procedure :: set_molecule => set_bse_atoms_container
     ! Should only need to take self, then access all the data 
     procedure :: run => run_bse  !TODO(Alex) Write me
     !> Get the energy of the lowest exciton state, the optical gap
     procedure :: get_optical_gap => get_bse_optical_gap     !TODO(Alex) Write me
     !> Get the two-particle absorption spectrum
     procedure :: get_spectrum => get_two_particle_spectrum  !TODO(Alex) Write me
  end type bse


  public :: bse_type

contains

  !> This should be a wrapper
  subroutine set_bse_mpi_environment(self, input_communicator)
    class(mpi_env_type), intent(inout) :: self
    !> MPI communicator of calling program. Note, this will be cast to
    ! MPI_COMM_TYPE in mpi_env%init 
    integer, intent(in) :: input_communicator

    call self%mpi_env%init(input_communicator = input_communicator)

  end subroutine set_bse_mpi_environment

  
  !> API wrapper for function construct_molecule
  !> Initialises molecule array and assigns a position, atomic number and symbol
  !> to each atom in the array
  subroutine set_bse_atoms_container(self, positions, atomic_numbers)
    class(atom_type), intent(inout) :: self
    real(wp), intent(in) :: positons(:,:)
    integer, intent(in) :: atomic_numbers(:)
    integer :: n_atoms 

    n_atoms = size(atomic_numbers)
    if (.not. allocated(self%molecule)) allocate(self%molecule(n_atoms))
    self%molecule = construct_molecule(positions, atomic_numbers)
 
  end subroutine set_bse_atoms_container
  

 
end module api
  
