module mpi_bse
  use kinds 
  implicit none

  private

  !> MPI type

  public, parameter, logical :: IO = .false. 
  
  
contains

  !> Initialise communicator
  !> Not likely to be needed
  subroutine initialise_mpi
    
  end subroutine initialise_mpi

  
  !> Finalise communicator


  
  subroutine set_IO()
    IO = rank == 0
  end subroutine set_IO
  
  
end module mpi_bse
