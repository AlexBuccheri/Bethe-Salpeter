module api_bse
  use bse, only : 
  implicit none

contains

  ! Top level API for calling a BSE calculation


  !eigenvalues vector
  !eigenvectors 2D array
  !Nelectrons, NHoles - basis size for BSE Hamiltonian, OR enery cutoff for states
  !Precision
  ! molecule: atom_species and coordinates 
  ! single-particle basis: basis per atom
  ! Tamm-Dancoff approx Logical
  ! Compute-spectrum Logical - default true
  ! Compute 2p eigenvectors: excitation energies + eigenvectors
  !  Exciton wave functions, excitation energies, and potentially spectrum
  ! Screening function
  ! Will need to pass in MPI communicator

  
  
  !> Calculator for BSE 
  type, public :: bse_type

     type(results_type) :: results 
     
  end type bse_type

  
  !> 
  type, private :: results_type

  end type results_type

  
  !>
  subroutine run_bse()

    !Put data into types
    !Initialise stuff

    call bethe_salpeter(bse_calc%results)

    !Spectrum

    !
    
  end subroutine run_bse

  
  
  
end module api_bse
