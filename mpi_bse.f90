module mpi_bse
  
  use kinds
  ! No proprocessing of MPI inclusion as this library is designed
  ! to be used with MPI 
  use mpif08 
  implicit none

  private

  !> MPI environment type
  !  See DL_POLY for this and build on it
  type mpi_env_type
     type(MPI_Comm) :: comm
     
     !> Read/write from this process
     logical :: IO = .false. 

   contains
     procedure :: init => mpi_initialise
     procedure :: destruct => mpi_destruct
     
  end type mpi_env_type


  
  
contains

  subroutine mpi_initialise(self, input_communicator)
    class(mpi_env_type), intent(inout) :: self 
    integer, intent(in), optional :: input_communicator

    if (.not. present(input_communicator)) then
       ! Initialise MPI_COMM_WORLD
    endif

    ! Duplicate comm

    ! Get rank, np, etc 
    
  end subroutine mpi_initialise
    
 
  subroutine mpi_destruct(self)
    
  end subroutine mpi_destruct


  

  
end module mpi_bse
