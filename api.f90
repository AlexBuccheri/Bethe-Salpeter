module api

  use precision, only : wp 
  implicit none
  private 
  
  
  !TODO(Alex) Move these types elsewhere

  !> Use F2003 declaration for this
  character(len=3), dimension(n_elements) :: atomic_number_to_symbol = (/ /)
  
  !> Type containing eigenstates
  type eigen_type
     private
     !> Eigenvectors
     real(wp), allocatable :: values(:)
     !> Eigenvectors, stored columnwise. Must be real 
     real(wp), allocatable :: vectors(:,:)
  end type eigen

  !> Type Atom
  type atom_type
     private
     !> Position in angstrom 
     real(wp), dimension(3) :: position
     !> Species label
     character(len=3) :: symbol
     !> Atomic number
     integer :: an 
  end type atom_type

  !> Type molecule
  type atoms_type
     private 
     !> Vector of atoms 
     type(atom_type), allocatable :: molecule(:)

   contains
     !>TODO(Alex0 Maybe this isn't a typer, but a function in a module
     !> Ah, it can be a method of BSE -> Should be a wrapper and the implementation
     ! should be elsewhere for use 
     !> Given a set of positions and atomic numbers, construct a
     !> vector of atom
     procedure :: init => initialise_atoms

  end type atoms

  
  type homo_lumo_type

  end type homo_lumo_type

  
  type two_particle_basis_type
     !>
     integer :: n_occupied
     integer :: n_virtual
  end type two_particle_basis_type

  
  type atomic_orbital_basis_type

  end type atomic_orbital_basis_type

  !> See DL_POLY for this and build on it
  type mpi_env_type
     use mpi_f08
     !> MPI communicator with F08 types  
     type(MPI_Comm) :: comm

  end type mpi_env_type

  !> Dielectric constant, and if use screening
  !> tam-dancoff approximation
  !> Polarisation of the light 
  !> Spectral or 2p-eigenvalues? Only spectrum atm. Spectral broadening 
  !> Solver, although only CG atm
  !> Neighbour cut-off. Have the library internally determine number of neighbours
  !> per atom from a distance matrix 
  
  !-------------------------------
  
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
     !> Each of these must be called prior to running a calculation
     procedure :: set_mpi_env => set_bse_mpi_environment
     procedure :: set_eigenstates => set_bse_eigenstates
     procedure :: set_molecule => set_bse_atoms_container
     ! Should only need to take self, then access all the data 
     procedure :: run => run_bse
     !> Get the energy of the lowest exciton state, the optical gap
     procedure :: get_optical_gap => get_bse_optical_gap
     !> Get the two-particle absorption spectrum
     procedure :: get_spectrum => get_two_particle_spectrum
  end type bse


  public :: bse_type

contains

  !> This should be a wrapper
  subroutine set_bse_mpi_environment(self, input_communicator)
    class(mpi_env_type), intent(inout) :: self

    ! mpi_env should have an init
    call self%mpi_env%init(input_communicator = input_communicator)

  end subroutine set_bse_mpi_environment


  subroutine 
  
  
  !> For type atoms
  subroutine initialise_atoms(self, positions, atomic_numbers)

    class(atoms_type), intent(inout) :: self(:)
    real(wp),          intent(in)    :: positions(:,:)
    integer,           intent(in)    :: atomic_numbers(:)

    integer :: ia, n_atoms
    
    call assert(size(positions,1) == 3), message="")
    call assert(size(atomic_numbers) == size(positions,2), message="")

    n_atoms = size(atomic_numbers)
    allocate(self(n_atoms))
    
    do ia = 1, n_atoms
       self(ia)%position(:) = positions(:,ia)
       self(ia)%an = atomic_numbers(ia)
       self(ia)%symbol = atomic_number_to_symbol(atomic_numbers(ia))
    enddo
       
  end subroutine initialise_atoms

  

end module api
  
