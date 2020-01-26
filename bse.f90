!********************************************************************
! Module for Bethe-Salpeter Calculations
!
! Compute resonant 2-particle Hamiltonian using orbital integrals 
! Matrix-vector multiplication
! Conjugate gradient
!********************************************************************

module bse
  !Modules used
  use var,     only: sp,dp,qp,Ne,Nh,Np,eigen_cnt,evcr,ortho_basis,Natoms,orb_bas,eps,TD,&
                     varb,emin=>emin_2p,emax=>emax_2p,dE=>delE_2p,Eg_SP=>Eg
  use parallel
  use exbnd,   only: construct_H2p,construct_H2p_v2,orthogonal_construct_H2p,&
                     orthogonal_construct_h2p_k2, serial_orth_construct_H2p
                     !TEST_H2p, construct_H2p_onsite
  !use kernel  - All required routines are in here instead
  use sort,    only: reallocate

  !Test module for BSE, not under version control
  !use kernels, only: orbital_integral_sys,resonant_h2p,construct_H2p_zeroth,&
  !                   test_H2p,dbtest_H2p

  !Interfaces
  !Matrix-vector product in TD approximation
  interface mvp_TD
     module procedure mvp_real_TD, mvp_cmplx_TD
  end interface

  !Matrix-vector product for full BSE Hamiltonian
  interface mvp_full
     module procedure mvp_real_full, mvp_cmplx_full
  end interface

  !Random array generation with desired symmetry/diagonals 
  interface trial_arrays
     module procedure trial_arrays_complex,trial_arrays_real
  end interface


  !Modular Variables
  integer rstart,rstop


contains


  !************************************************
  !Bethe-Salpeter Main Routine
  !Computes two-particle excitations
  !************************************************
  subroutine bethe_salpeter(eigen)
    use var, only:radius 
    implicit none


    !****************
    !Declarations
    !****************

    !Global/Modular Variables 
    integer      iter
    real(dp)     ethr,anorm
    complex(dp)  omega
         
    !Local Variables 
    integer      cnt,i,j,N2p,N,Nlt,irow,icol,VBT,CBB,Npts
    real(dp)     erange,wr,t1_local,t2_local,t21,t21_lg,orth,&
                 CG_t1_local,CG_t2_local,CG_time,sval,Eg_opt,Ee,Eh
    logical      onset
    
    !Local Arrays
    integer,     allocatable               :: CGitr(:),tl(:,:)
    real(dp),    allocatable               :: CGnorm(:),bse_time_lcl(:),bse_time_glb(:)

    !Global Arrays
    real(dp),    allocatable,   intent(in) :: eigen(:,:)
    real(dp),    allocatable               :: H2p(:),K2(:),e2(:),worb(:),H2p_s(:,:)
    complex(dp), allocatable               :: X(:),b(:)
    
    !Test Variables/Arrays
    integer      flag
    logical      cj
    real(dp)                 :: xr,xi
    real(dp),    allocatable :: H(:,:),Hpkd(:)
    complex(dp), allocatable :: A(:,:),Xid(:)
    

    !****************
    !Main Code
    !****************
    call cpu_time(t1_local)


    !Library
    ! Call read_data()
    ! Call initialise_data()
    !If 


    
    !Determine Ne and Nh for Np pair-products in the basis
    if(Ne==1 .and. Nh==1) call get_single_particle_basis(eigen,Np, Ne,Nh,tl)

    !Discard any eigenvectors that aren't utilised by the 2-particle basis
    call reduce_eigenstates(Nh,Ne,eigen_cnt,eigen,evcr,tl)

    !2-particle resonant Hamiltonian dimensions 
    N2p=Ne*Nh
  
    !Dimensions of full BSE matrix A
    N=2*N2p

    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1

    !Single-particle gap
    Eg_SP=eigen(1,CBB)-eigen(1,VBT)  

    !Inform user 
    if(IO)then
       write(*,'(X,A,I3,A,I3,A,I3,A)')'The pair-product basis of ',Np,' states is comprised of: ',&
                                       Nh,' holes and ',Ne,' electrons'
       write(*,*) 'Two-particle Hamiltonian basis: ',N
    endif

    !Distribute lower triangle of resonant H2p
    call distribute_LT(N2p, rstart,rstop)

    !Distribute lower triangular elements
    if(numprocs==1)then
       Nlt=0.5*N2p*(N2p+1)
    else
       Nlt=0
       do i=rstart,rstop
          Nlt=Nlt+i
       enddo
    endif

    !Resonant 2-particle Hamiltonian
    allocate(H2p(Nlt))
    H2p(:)=0._dp
   !write(*,*) 'Size of kernel held in memory: ',Nlt, rank


    !--------------------------------------------------------------------------
    !Kernel construction
    !Compute lower tiangle(s) of BSE kernel in an orthogonal orbital basis
    !Name H2p is misleading as diagonal transition energies are not included
    !--------------------------------------------------------------------------

    !Tamm-Dancoff Approximation. Turn off coupling submatrices
    if(TD .eqv. .true.)then
       if(IO) write(*,*) 'Using Tamm-Dancoff Approximation'
       !H2p = K_(vc)(v'c')
       call orthogonal_construct_H2p(H2p,rstart,rstop)

    !Full BSE calculation requires coupling submatrices
    elseif(TD .eqv. .false.)then
       if(IO) write(*,*) 'Full BSE calculation'
       allocate(K2(Nlt))
       K2(:)=0._dp
       !H2p = K_(vc)(v'c') and K2 = K(vc)(c'v')   
       call orthogonal_construct_H2p_K2(H2p,K2,rstart,rstop)
    endif
    if(IO) write(*,*)'Finished kernel construction'


    !Output the diagonal elements of the BSE kernel
    !call output_diagonal_elements(Ne,N2p,Np,tl,H2p)
    !deallocate(tl)    



    !--------------------------------------------------------------
    !Unused routine variations and test calls
    !--------------------------------------------------------------
    !Nonorthogonal implementation: Orbital interactions restricted to NN
    !Compute resonant 2-particle Hamiltonian using orbital integrals 
    !Returns lower triangle in packed form 
    !call construct_H2p(H2p,rstart,rstop)
    !if(IO) write(*,*)'Finished H2p construction'

    !Orbital interactions restricted to on-site only
    !call construct_H2p_onsite(H2p,rstart,rstop)
    !if(IO) write(*,*)'Finished onsite-H2p construction'

    !if(IO)write(*,*)'Call to TEST_H2p bse.f90 line 128.'
    !call TEST_H2p(H2p,rstart,rstop)

    !Two test routines for orbital integral
    !call coulomb_integral_NN_a2
    !call coulomb_integral_NN_a1_ALTr()

    !Cleaner loop structure but heavier memory requirement
    !call construct_H2p_v2(H2p,rstart,rstop)
    !if(IO) write(*,*)'Finished H2p construction'

    !Compute worb and H2p using stripped-down routine
    !call orbital_integral_sys(worb)
    !call resonant_H2p(H2p,rstart,rstop,worb)

    !Construct H2p in zeroth approximation (untested)
    !call construct_H2p_zeroth(H2p,rstart,rstop)

    !Test db_factor
    !call orbital_integral_sys(worb)
    !call dbtest_H2p(H2p,worb,rstart,rstop)
    !if(IO) write(*,*) 'bse stops line 109'
    !call stop_mpi
    !--------------------------------------------------------------
  

    !Right-hand side and 2-particle coefficients
    allocate(b(N), X(N))

    !Initial guess at X. Could also be random vector
    X(:)=0._dp

    !Construct RHS of Ax=b, where b=-/+<ic|e.r|iv>
    if(.not. varb) call dipole_matrix_elements(eigen,evcr, b)
    if(varb  )    call var_dipole_matrix_elements(eigen,evcr, b)

    !Threshold for convergence in X
    ethr=1.e-12


    !*********************************************************
    !Compute excitonic spectrum, looping over frequency input
    !*********************************************************
    
    !Energy range 
    erange=emax-emin
    
    !Energy spacing (must be <0.5*Img(w), where Img(w)==smearing)
    if(0.5*eps<dE)then
       dE=0.5*eps
       if(IO)write(*,*)'Energy spacing should be half the size of the broadening'
       if(IO)write(*,*)'Reset energy spacing to:',dE,'eV'
    endif
 
    !Number of points required to span energy range with spacing dE
    Npts=ceiling(erange/dE)+1

    !Absorption array
    allocate(e2(Npts))
    e2(:)=0._dp

    !CG convergence details arrays
    allocate(CGnorm(Npts),CGitr(Npts))
    CGitr=0
    CGnorm=0._dp
    call cpu_time(CG_t1_local)


    !Initialise variables used in finding 2p onset
    sval=0._dp
    onset=.false.

    !------------------------------
    !Static Broadening 
    !------------------------------
    if(.not. varb)then
       if(IO) write(*,*)'Static broadening applied to spectra'
       !Create 2-particle coefficients for [emin:emax]
       do i=1,Npts
          
          !Real Energy
          wr=emin+(i-1)*dE
          
          !Static broadening
          omega=cmplx(wr,eps)
          
          !Use Conjugate gradient to solve Ax=b, for energy omega
          if(TD)       call bcgsolve_parallel(N2p,H2p,eigen,omega,ethr,X,b,iter,anorm)        !, A,flag)
          if(.not. TD) call bcgsolve_parallel(N2p,H2p,eigen,omega,ethr,X,b,iter,anorm,K2=K2)  !, A,flag)
          
          !Store number of iterations + norm of the residuals per energy point
          CGitr(i)=iter
          CGnorm(i)=anorm
          
          !Absorption at energy wr
          e2(i) = aimag(dot_product(X(1:N2p),-b(1:N2p)) + &
                        dot_product(X(N2p+1:N),b(N2p+1:N)) )

          !Establish onset (optical gap)
          if(e2(i)>sval) sval=e2(i)
          if(e2(i)<sval .and. (onset .eqv. .false.))then
             if(IO) write(*,*) 'Two-particle onset:',emin+(i-2)*dE,'eV'
             Eg_opt=emin+(i-2)*dE  !i.e. wr for prior point i
             onset=.true.
          endif
      
       enddo
    endif

    !------------------------------
    !Variable Broadening 
    !------------------------------
    if(varb)then
       if(IO) write(*,*)'Variable broadening applied to spectra'
       !Create 2-particle coefficients for [emin:emax]
       do i=1,Npts
          
          !Real Energy
          wr=emin+(i-1)*dE
          
          !Complex energy input
          if(onset .eqv. .false.) omega=cmplx(wr,eps)
          !Variable broadening 
          if(onset .eqv. .true.)  omega=cmplx(wr,veps(wr,Eg_opt))
          
          !Use Conjugate gradient to solve Ax=b, for energy omega
          if(TD)       call bcgsolve_parallel(N2p,H2p,eigen,omega,ethr,X,b,iter,anorm)        !, A,flag)
          if(.not. TD) call bcgsolve_parallel(N2p,H2p,eigen,omega,ethr,X,b,iter,anorm,K2=K2)  !, A,flag)
          
          !Store number of iterations + norm of the residuals per energy point
          CGitr(i)=iter
          CGnorm(i)=anorm
          
          !Absorption at energy wr
          e2(i) = aimag(dot_product(X(1:N2p),-b(1:N2p)) + &
                        dot_product(X(N2p+1:N),b(N2p+1:N)) )
          
          !Establish onset (optical gap)
          !Also allows variable broadening to applied, starting at the onset
          if(e2(i)>sval) sval=e2(i)
          if(e2(i)<sval .and. (onset .eqv. .false.))then
             if(IO) write(*,*) 'Two-particle onset:',emin+(i-2)*dE,'eV'
             Eg_opt=emin+(i-2)*dE  !i.e. wr for prior point i
             onset=.true.
          endif
          
       enddo
    endif


    !Exciton binding energy
    if(IO)then
       write(*,*) 'Exciton binding energy:',(Eg_SP-Eg_opt)*1000._dp,' meV'
       write(*,*) 'Determined from the 2p spectrum with a grid spacing:',dE*1000._dp,' meV'
    endif

    !Ouput CG convergence details 
    call cpu_time(CG_t2_local)
    if(IO)then
       write(*,*) 'Avg number of CG iterations',(1._dp/dble(Npts))*sum(CGitr)
       write(*,*) 'Max number of CG iterations', maxval(CGitr)
       write(*,*) 'Avg CG convergence ',(1._dp/dble(Npts))*sum(CGnorm)
       write(*,*) 'Least converged CG residual norm', maxval(CGnorm)
    endif
    
    !Output spectrum
    if(IO) call output_2p_spectrum(emin,emax,dE,Npts,aimag(omega),e2)

    !Output useful values
    Eh=(eigen(1,VBT)-eigen(1,1))*1000._dp
    Ee=(eigen(1,VBT+Ne)-eigen(1,CBB))*1000._dp

    if(IO)then
       open(unit=100,file='opt.dat')
       write(100,'(A)') '#Diameter (Ang), Binding energy (meV), <Eh,Ee> (meV),  Eh,  Eh'
       write(100,'(F6.2,2X,4(F9.4,2X))') radius*2.,(Eg_SP-Eg_opt)*1000._dp, 0.5*(Eh+Ee),Eh,Ee
       close(100)
    endif

    !*********************************************************
    !Timings
    !*********************************************************

    !Time full routine
    call cpu_time(t2_local)

    !Local and global time arrays
    allocate(bse_time_lcl(numprocs),bse_time_glb(numprocs))
    bse_time_lcl(:)=0._dp
    bse_time_glb(:)=0._dp

    !Fill element of local time array, corresponding to process
    bse_time_lcl(rank+1)=(t2_local-t1_local)
    !Fill global time array with all local values
    call mpi_allreduce(bse_time_lcl,bse_time_glb,numprocs,MPI_DOUBLE_PRECISION,&
                       MPI_SUM,MPI_COMM_WORLD,MPIierr)

    !Sum local timings
    call mpi_allreduce((CG_t2_local-CG_t1_local),CG_time,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce((t2_local-t1_local),t21,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
                        MPI_COMM_WORLD,MPIierr)

    !Output timings
    if(IO)then
       write(*,'(A,X,F20.10,A)') ' Avg CG time =',CG_time/(dble(numprocs)*60._dp),' mins'
       write(*,'(A,X,F20.10,A)') ' Avg BSE runtime =',t21/(dble(numprocs)*60._dp),' mins'
       write(*,'(A,X,F20.10,A,X,I3)') ' Longest BSE runtime =',maxval(bse_time_glb)/60._dp,&
                                      ' mins on process',maxloc(bse_time_glb)-1
    endif

    !Local Timing
    !write(*,'(A,X,I3,A,F20.10,A)') 'Full BSE subroutine time for process',rank,'=',&
                                     !(t2-t1)/60._dp,' mins'


  end subroutine bethe_salpeter


  !*********************************************************
  !Vary the Lorentzian broadening parameter as a function
  !of the transition energy
  !
  !Inputs:
  !Et    Transition energy Ec-Ev
  !Eg    Absorption on-set
  !*********************************************************
  real(dp) function veps(Et,Eg)
    use var,only:dp
    implicit none

    !**************
    !Declarations 
    !**************
    real(dp),    intent(in) :: Et,Eg
    real(dp)  :: emin,emax,eps_min,eps_max,eps_range,f


    !**************
    !Main Routine
    !**************
    !Energy range over which the variable broadening is applied
    emin=Eg
    emax=Eg+2._dp       !4._dp
    !Max and min values for spectral broadening, over the spectral range [emin:emax]
    eps_min=0.025
    eps_max=0.25    !0.35

    !If transition energy exceeds emax, assign maximum value for broadening
    if(Et>emax)then
       veps=eps_max
       return
    endif

    !Else, apply variable broadening, which increases linearly with energy
    
    !Range for application of broadening 
    eps_range=eps_max-eps_min
    
    !Fraction too small number so 0/(emax-emin) /=NaN
    f=((Et-emin)/(emax-emin))+1.e-6

    !Variable broadening in the range[emin:emax]
    veps=eps_min+(f*eps_range)
    return
  end function veps




  !***************************************************************
  !Serial routine to construct H2p in unpacked form, using the
  !lower triangle, L. 
  !Current not including transition energies on the diagonals
  !***************************************************************
  subroutine distribute_full_H2p(Ne,Nh,L,eigen,omega, H2p_full)
    use var, only:dp,eigen_cnt
    implicit none
    
    !***********************
    !Declarations
    !***********************
    !Global integers
    integer,  intent(in)  :: Ne,Nh
    complex(dp), intent(in) :: omega

    !Global Arrays
    real(dp), allocatable,intent(in) :: L(:,:),eigen(:,:)
    complex(dp), allocatable ::            H2p_full(:,:)

    !Local Variables
    integer row,col,N2p,h_row_cnt,e_row_cnt,h_row,e_row,VBT,CBB

    !***********************
    !Main Routine
    !***********************
    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1

    !HALF of the Basis size of 2p Hamiltonian 
    N2p=Ne*Nh

    !Off-diagonal terms
    do row=1,Ne*Nh
       do col=1,row-1
          
          !Top left Quad
          H2p_full(row,col) = L(row,col)  !Lower triangle
          H2p_full(col,row) = L(row,col)  !Upper triangle 
          !Top Right Quad
          H2p_full(row,col+N2p) =  L(row,col) !Lower triangle
          H2p_full(col,row+N2p) =  L(row,col) !Upper triangle
          !Bottom Left Quad
          H2p_full(row+N2p,col) = -L(row,col) !Lower triangle
          H2p_full(col+N2p,row) = -L(row,col) !Upper triangle
          !Bottom Right Quad
          H2p_full(row+N2p,col+N2p) = -L(row,col) !Lower triangle
          H2p_full(col+N2p,row+N2p) = -L(row,col) !Upper triangle

       enddo
    enddo

    !Diagonal terms 
    do row=1,Ne*Nh

       !Counters for electron and hole row indices
       !of resonant sub-matrix
       h_row_cnt=((row-1)/Ne)+1
       e_row_cnt=row-((h_row_cnt-1)*Ne)
       !Hole and electron indices for rows of (h,e)
       h_row=(VBT+1) - h_row_cnt
       e_row= e_row_cnt + (CBB-1)

       !Top left
       H2p_full(row,row) =  L(row,row) + (eigen(1,e_row)-eigen(1,h_row)-omega)
       !Top right
       H2p_full(row,row+N2p) = L(row,row) 
       !Bottom left
       H2p_full(row+N2p,row) = -L(row,row)
       !Bottom right
       H2p_full(row+N2p,row+N2p) = -L(row,row) + (-eigen(1,e_row)+eigen(1,h_row)-omega)
    enddo

 
  end subroutine distribute_full_H2p


  !**************************************************************
  !Downside the eigenvector array such that any eigenvectors not
  !utilised by the two-particle basis are discarded
  !**************************************************************
  subroutine reduce_eigenstates(Nh,Ne,eigen_cnt,eigen,evcr,tl)
    use var,  only: dp,Natoms,orb_bas
    implicit none

    !****************
    !Declarations
    !****************

    !Global Variables
    integer,  intent(in)  :: Nh,Ne

    !Global Arrays
    integer,  allocatable :: eigen_cnt(:,:),tl(:,:)
    real(dp), allocatable :: eigen(:,:),evcr(:,:)

    !Local Variables
    integer  :: N,VBT,i,j,VBi,CBf,ivc

    !Local Arrays
    real(dp), allocatable :: tmp(:,:)


    !****************
    !Main
    !****************

    !Total number of orbital states
    N=Natoms*orb_bas

    !VBT
    VBT=eigen_cnt(1,1)

    !(Update iv,ic) pair basis index array
    do ivc=1,size(tl,1)
       tl(ivc,1)=tl(ivc,1)-VBT+Nh
       tl(ivc,2)=tl(ivc,2)-VBT+Nh
    enddo


    !Eigenvectors
    allocate(tmp(N,Nh+Ne))
    tmp(:,:)=0._dp
    
    VBi=VBT-Nh+1
    CBf=VBT+Ne
    tmp(1:N,1:Nh+Ne)=evcr(1:N,VBi:CBf)


    !Keep subset of eigenvectors for 2-particle basis
!!$    do i=1,Nh+Ne
!!$       j=i-VBT+Nh     !Erroneous j=i+VBT-Nh
!!$       tmp(:,i)=evcr(:,j)
!!$    enddo

    !Pass from tmp to evcr array
    deallocate(evcr)
    allocate(evcr(N,Nh+Ne))
    evcr(:,:)=tmp(:,:)
    deallocate(tmp)

    !Repeat for eigenvalues
    allocate(tmp(1,Nh+Ne))
!!$    tmp(:,:)=0._dp
!!$    do i=1,Nh+Ne
!!$       j=i-VBT+Nh 
!!$       tmp(1,i)=eigen(1,i)
!!$    enddo


    tmp(1,1:Nh+Ne)=eigen(1,VBi:CBf)

    !Pass from tmp to evcr array
    deallocate(eigen)
    allocate(eigen(1,Nh+Ne))
    eigen(:,:)=tmp(:,:)
    deallocate(tmp)

    !Update index for VBT state (== number of valence/hole states)
    eigen_cnt(1,1)=Nh
    !Update number of conduction/electron states
    eigen_cnt(1,2)=Ne


  end subroutine reduce_eigenstates



  !***************************************************************
  !Matrix-vector Product on full BSE Hamiltonian
  !Real inputs, complex output due to omega
  !Computes distributed matrix-vector product for Ax=b
  !          _        _
  !where A= |   H   K  | and H = (Ec-Ev)I + K
  !         |_ -K  -H _|
  ! Inputs
  ! H      H= Lower triangle of resonant 2-particle Hamiltonian
  !        in packed form (distributed)
  ! K      K = lower triangle of K_(vc)(c'v') in packed form
  ! X      Guess at two-particle coefficients (global)
  ! eigen  Tight-binding eigenvalues
  ! omega  Complex frequency-dependence
  !
  ! Outputs
  ! b      Vector product (global)
  !***************************************************************

  subroutine mvp_real_full(H,K2,N2p,X,eigen,omega,cj,b)
    use var, only: eigen_cnt
    implicit none


    !****************
    !Declarations
    !****************

    !Local Variables
    integer      irow,icol,i1,cnt,i,VBT,CBB,H_row_cnt,e_row_cnt,e_row,h_row,&
                 ik,N  
    real(dp)     factor
    complex(dp)  w

    !Global Variables
    integer,     intent(in) :: N2p
    complex(dp), intent(in) :: omega
    logical,     intent(in) :: cj

    !Global Arrays
    real(dp),    allocatable, intent(in)  :: eigen(:,:),H(:),K2(:)
    complex(dp), allocatable, intent(in)  :: X(:)
    complex(dp), allocatable              :: b(:) 

    !Local Arrays
    complex(dp), allocatable :: brow(:)


    !****************
    !Main Routine
    !****************

    !Finite systems
    ik=1

    !Dimensions of full BSE matrix
    N=2*N2p

    !Initialise outputs
    !A(:,:)=0._dp
    b(:)=0._dp

    !Local output vector 
    allocate(brow(N))
    brow(:)=0._dp


    !*****************************
    !Matrix-Vector Multiplication
    !*****************************
    !Ax=b
    if(cj .eqv. .false.) factor=1._dp
    !A^dagger x=b
    if(cj .eqv. .true. ) factor=-1._dp
  
    cnt=0
    do irow=rstart,rstop
       do icol=1,irow

          cnt=cnt+1
          !Lower triangle * X 
          !Top left and top right Quadrants of A:   TL         TR
          brow(irow)=brow(irow)+ ( H(cnt)*X(icol) + factor*K2(cnt)*X(icol+N2p))
          !Bottom left and bottom right Quadrants of A:      BL         BR
          brow(irow+N2p)=brow(irow+N2p)- (factor*K2(cnt)*X(icol) + H(cnt)*X(icol+N2p))
       enddo
       
       !Use symmetry to multiply elements above diagonal with X
       !Note, -1 in brow(1:irow-1) prevents double-counting diagonal elements from
       !upper triangular matrix
       
       !1st element of ith row in packed LT array
       i1=cnt-irow+1
       
       !Top left and top right Quadrants of A:         TL         TR       
       brow(1:irow-1)=brow(1:irow-1)+  ( H(i1:cnt-1)*X(irow) + factor*K2(i1:cnt-1)*X(irow+N2p)) 
       !Bottom left and bottom right Quadrants of A:                     BL         BR
       brow(N2p+1:irow+N2p-1)=brow(N2p+1:irow+N2p-1) -(factor*K2(i1:cnt-1)*X(irow) + &
                                                                       H(i1:cnt-1)*X(irow+N2p))

    enddo


    !*****************************************************************
    !MV product for Diagonal terms due to quasi-particle contribution
    !*****************************************************************
    if(.not. cj) w=omega
    if(cj)       w=conjg(omega)

    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1

    !Loop over diagonals
    do irow=rstart,rstop
       !Counters for electron and hole row indices
       !of resonant sub-matrix
       h_row_cnt=((irow-1)/Ne)+1
       e_row_cnt=irow-((h_row_cnt-1)*Ne)
       !Hole and electron indices for rows of (h,e)
       h_row=(VBT+1) - h_row_cnt
       e_row= e_row_cnt + (CBB-1)

       !Add quasi-particle contribution to MV product (distributed)
       !Resonant (TL quadrant)
       brow(irow)=brow(irow)+( (eigen(ik,e_row)-eigen(ik,h_row)-w)* X(irow) )
       !Anti-resonant (BR quadrant)
       brow(irow+N2p)=brow(irow+N2p)+( (-eigen(ik,e_row)+eigen(ik,h_row)-w)* X(irow+N2p)  )
    enddo

 
    !Sum elements of local output vector to give full vector
    call mpi_allreduce(brow,b,N,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    deallocate(brow)


  end subroutine mvp_real_full



  !***************************************************************
  !Matrix-vector Product on full BSE Hamiltonian
  !Complex H and K due to complex orbital basis
  !Computes distributed matrix-vector product for Ax=b
  !          _        _
  !where A= |   H   K  | and H = (Ec-Ev)I + K
  !         |_ -K* -H*_|
  ! Inputs
  ! H      H= Lower triangle of resonant 2-particle Hamiltonian
  !        in packed form (distributed)
  ! K      K = lower triangle of K_(vc)(c'v') in packed form
  ! X      Guess at two-particle coefficients (global)
  ! eigen  Tight-binding eigenvalues
  ! omega  Complex frequency-dependence
  !
  ! Outputs
  ! b      Vector product (global)
  !***************************************************************
  subroutine mvp_cmplx_full(H,K2,N2p,X,eigen,omega,cj,b)
    use var, only: eigen_cnt
    implicit none


    !****************
    !Declarations
    !****************

    !Local Variables
    integer      irow,icol,i1,cnt,i,VBT,CBB,H_row_cnt,e_row_cnt,e_row,h_row,&
                 ik,N  !,H_row_cnt2,e_row_cnt2,e_row2,h_row2
    real(dp)     factor
    complex(dp)  w

    !Global Variables
    integer,     intent(in) :: N2p
    complex(dp), intent(in) :: omega
    logical,     intent(in) :: cj

    !Global Arrays
    real(dp),    allocatable, intent(in)  :: eigen(:,:)
    complex(dp), allocatable, intent(in)  :: H(:),K2(:),X(:)
    complex(dp), allocatable              :: b(:) 

    !Local Arrays
    complex(dp), allocatable :: brow(:)


    !****************
    !Main Routine
    !****************

    !Finite systems
    ik=1

    !Dimensions of full BSE matrix
    N=2*N2p

    !Initialise outputs
    !A(:,:)=0._dp
    b(:)=0._dp

    !Local output vector 
    allocate(brow(N))
    brow(:)=0._dp

    if(IO)then
       write(*,*)'Written in same manner as prior complex product routine. Untested'
    endif

    !*****************************
    !Matrix-Vector Multiplication
    !*****************************
    
    !Ax=b 
    if(.not. cj)then
       cnt=0
       do irow=rstart,rstop
          do icol=1,irow
             cnt=cnt+1
             !Lower triangle * X 
             !Top left and top right Quadrants of A:   TL         TR
             brow(irow)=brow(irow)+ ( H(cnt)*X(icol) + K2(cnt)*X(icol+N2p))
             !Bottom left and bottom right Quadrants of A:      BL         BR
             brow(irow+N2p)=brow(irow+N2p)- (conjg(K2(cnt))*X(icol) + conjg(H(cnt))*X(icol+N2p))
          enddo
          
          !Use symmetry to multiply elements above diagonal with X
          !Note, -1 in brow(1:irow-1) prevents double-counting diagonal elements from
          !upper triangular matrix
          
          !1st element of ith row in packed LT array
          i1=cnt-irow+1
          
          !Top left and top right Quadrants of A:         TL         TR       
          brow(1:irow-1)=brow(1:irow-1)+  ( conjg(H(i1:cnt-1))*X(irow) + conjg(K2(i1:cnt-1))*X(irow+N2p)) 
          !Bottom left and bottom right Quadrants of A:                     BL         BR
          brow(N2p+1:irow+N2p-1)=brow(N2p+1:irow+N2p-1)-(K2(i1:cnt-1)*X(irow) + H(i1:cnt-1)*X(irow+N2p))
          
       enddo
    endif


    !A^dagger x=b
    if(cj)then
       cnt=0
       do irow=rstart,rstop
          do icol=1,irow
             cnt=cnt+1
             !Lower triangle * X 
             !Top left and top right Quadrants of A:   TL         TR
             brow(irow)=brow(irow)+ ( H(cnt)*X(icol) - conjg(K2(cnt))*X(icol+N2p))
             !Bottom left and bottom right Quadrants of A:      BL         BR
             brow(irow+N2p)=brow(irow+N2p)+ (K2(cnt)*X(icol) - conjg(H(cnt))*X(icol+N2p))
          enddo
          
          !Use symmetry to multiply elements above diagonal with X
          !Note, -1 in brow(1:irow-1) prevents double-counting diagonal elements from
          !upper triangular matrix
          
          !1st element of ith row in packed LT array
          i1=cnt-irow+1
          
          !Top left and top right Quadrants of A:         TL         TR       
          brow(1:irow-1)=brow(1:irow-1)+  ( conjg(H(i1:cnt-1))*X(irow) - K2(i1:cnt-1)*X(irow+N2p)) 
          !Bottom left and bottom right Quadrants of A:                     BL         BR
          brow(N2p+1:irow+N2p-1)=brow(N2p+1:irow+N2p-1)+(conjg(K2(i1:cnt-1))*X(irow) - H(i1:cnt-1)*X(irow+N2p))
          
       enddo
    endif


    !*****************************************************************
    !MV product for Diagonal terms due to quasi-particle contribution
    !*****************************************************************
    if(.not. cj) w=omega
    if(cj)       w=conjg(omega)

    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1

    !Loop over diagonals
    do irow=rstart,rstop
       !Counters for electron and hole row indices
       !of resonant sub-matrix
       h_row_cnt=((irow-1)/Ne)+1
       e_row_cnt=irow-((h_row_cnt-1)*Ne)
       !Hole and electron indices for rows of (h,e)
       h_row=(VBT+1) - h_row_cnt
       e_row= e_row_cnt + (CBB-1)

       !Add quasi-particle contribution to MV product (distributed)
       !Resonant (TL quadrant)
       brow(irow)=brow(irow)+( (eigen(ik,e_row)-eigen(ik,h_row)-w)* X(irow) )
       !Anti-resonant (BR quadrant)
       brow(irow+N2p)=brow(irow+N2p)+( (-eigen(ik,e_row)+eigen(ik,h_row)-w)* X(irow+N2p)  )
    enddo
 

    !Sum elements of local output vector to give full vector
    call mpi_allreduce(brow,b,N,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    deallocate(brow)

  end subroutine mvp_cmplx_full



  !***************************************************************
  !Computes distributed matrix-vector product for Ax=b
  !H is real but the product is complex due to omega
  !H is symmetric in real basis, therefore H_2p = (H_2p)^T
  !          _       _
  !where A= |  H      | and H = (Ec-Ev)I + K
  !         |_    -H _|
  ! Inputs
  ! H      H= Lower triangle of resonant 2-particle Hamiltonian
  !        in packed form (distributed)
  ! X      Guess at two-particle coefficients (global)
  ! eigen  Tight-binding eigenvalues
  ! omega  Complex frequency-dependence
  !
  ! Outputs
  ! b      Vector product (global)
  !***************************************************************
  subroutine mvp_real_TD(H,N2p,X,eigen,omega,cj,b)
    use var, only: eigen_cnt
    implicit none


    !****************
    !Declarations
    !****************

    !Local Variables
    integer      irow,icol,i1,cnt,i,VBT,CBB,H_row_cnt,e_row_cnt,e_row,h_row,&
                 ik,N  !,H_row_cnt2,e_row_cnt2,e_row2,h_row2
    real(dp)     factor
    complex(dp)  w

    !Global Variables
    integer,     intent(in) :: N2p
    complex(dp), intent(in) :: omega
    logical,     intent(in) :: cj

    !Global Arrays
    real(dp),    allocatable, intent(in)  :: eigen(:,:),H(:)
    complex(dp), allocatable, intent(in)  :: X(:)
    complex(dp), allocatable              :: b(:) 

    !Local Arrays
    complex(dp), allocatable :: brow(:)


    !****************
    !Main Routine
    !****************

    !Finite systems
    ik=1

    !Dimensions of full BSE matrix
    N=2*N2p

    !Initialise outputs
    !A(:,:)=0._dp
    b(:)=0._dp

    !Local output vector 
    allocate(brow(N))
    brow(:)=0._dp


    !*****************************
    !Matrix-Vector Multiplication
    !*****************************
    !Only used for TD approximation
    factor=0._dp
    !Would make more sense to remove terms multiplied by factor
    !Implement that once these changes are found to be stable 
    
    cnt=0
    do irow=rstart,rstop
       do icol=1,irow
          cnt=cnt+1
          !Lower triangle * X 
          !Top left and top right Quadrants of A:   TL         TR
          brow(irow)=brow(irow)+ H(cnt)*( X(icol)+ factor*X(icol+N2p))
          !Bottom left and bottom right Quadrants of A:      BL         BR
          brow(irow+N2p)=brow(irow+N2p)- H(cnt)*( factor*X(icol) +X(icol+N2p))
       enddo
       
       !Use Hermitian symmetry to multiply elements above diagonal with X
       !Note, -1 in brow(1:irow-1) prevents double-counting diagonal elements from
       !upper triangular matrix
       
       !1st element of ith row in packed LT array
       i1=cnt-irow+1
       
       !Top left and top right Quadrants of A:         TL         TR       
       brow(1:irow-1)=brow(1:irow-1)+  H(i1:cnt-1)*( X(irow)+factor*X(irow+N2p)) 
       !Bottom left and bottom right Quadrants of A:                     BL         BR
       brow(N2p+1:irow+N2p-1)=brow(N2p+1:irow+N2p-1)-H(i1:cnt-1)*( factor*X(irow)+X(irow+N2p))
       
    enddo


    !*****************************************************************
    !MV product for Diagonal terms due to quasi-particle contribution
    !*****************************************************************
    if(.not. cj) w=omega
    if(cj)       w=conjg(omega)

    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1

    !Loop over diagonals
    do irow=rstart,rstop
       !Counters for electron and hole row indices
       !of resonant sub-matrix
       h_row_cnt=((irow-1)/Ne)+1
       e_row_cnt=irow-((h_row_cnt-1)*Ne)
       !Hole and electron indices for rows of (h,e)
       h_row=(VBT+1) - h_row_cnt
       e_row= e_row_cnt + (CBB-1)

       !Add quasi-particle contribution to MV product (distributed)
       !Resonant (TL quadrant)
       brow(irow)=brow(irow)+( (eigen(ik,e_row)-eigen(ik,h_row)-w)* X(irow) )
       !Anti-resonant (BR quadrant)
       brow(irow+N2p)=brow(irow+N2p)+( (-eigen(ik,e_row)+eigen(ik,h_row)-w)* X(irow+N2p)  )
    enddo

 
    !Sum elements of local output vector to give full vector
    call mpi_allreduce(brow,b,N,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    deallocate(brow)

  end subroutine mvp_real_TD




  !***************************************************************
  !Computes distributed matrix-vector product for Ax=b
  !Complex H due to complex orbital basis 
  !H is Hermitian, therefore H_2p = H_2p^dagger
  !
  !          _       _
  !where A= |  H      | and H = (Ec-Ev)I + K
  !         |_    -H*_|
  ! Inputs
  ! H      H= Lower triangle of resonant 2-particle Hamiltonian
  !        in packed form (distributed)
  ! X      Guess at two-particle coefficients (global)
  ! eigen  Tight-binding eigenvalues
  ! omega  Complex frequency-dependence
  !
  ! Outputs
  ! b      Vector product (global)
  !***************************************************************
  subroutine mvp_cmplx_TD(H,N2p,X,eigen,omega,cj,b)
    use var, only: eigen_cnt
    implicit none


    !****************
    !Declarations
    !****************

    !Local Variables
    integer      irow,icol,i1,cnt,i,VBT,CBB,H_row_cnt,e_row_cnt,e_row,h_row,&
                 ik,N   !,H_row_cnt2,e_row_cnt2,e_row2,h_row2
    complex(dp)  w

    !Global Variables
    integer,     intent(in) :: N2p
    complex(dp), intent(in) :: omega
    logical,     intent(in) :: cj

    !Global Arrays
    real(dp),    allocatable, intent(in)  :: eigen(:,:)
    complex(dp), allocatable, intent(in)  :: H(:),X(:)
    complex(dp), allocatable              :: b(:)  !A(:,:)

    !Local Arrays
    complex(dp), allocatable :: brow(:)


    !****************
    !Main Routine
    !****************

    !Finite systems
    ik=1

    !Dimensions of full BSE matrix
    N=2*N2p

    !Initialise outputs
    !A(:,:)=0._dp
    b(:)=0._dp

    !Local output vector 
    allocate(brow(N))
    brow(:)=0._dp


    !*************************************************
    !Matrix-Vector Multiplication in TD Approximation
    !*************************************************
    !H has hermitian symmetry, therefore in the TD approximation 
    !H_2p = (H_2p)^dagger

    !Ax=b == A^daggerx=b
    cnt=0
    do irow=rstart,rstop
       do icol=1,irow
          cnt=cnt+1
          !Lower triangle * X 
          !Top left and top right Quadrants of A:  TL         
          brow(irow)=brow(irow)+ H(cnt)*X(icol)
          !Bottom left and bottom right Quadrants of A:        BR
          brow(irow+N2p)=brow(irow+N2p)-conjg(H(cnt))*X(icol+N2p)
       enddo
       
       !Use Hermitian symmetry to multiply elements above diagonal with X
       !Note, -1 in brow(1:irow-1) prevents double-counting diagonal elements from
       !upper triangular matrix
       
       !1st element of ith row in packed LT array
       i1=cnt-irow+1
       
       !Top left and top right Quadrants of A:         TL               
       brow(1:irow-1)=brow(1:irow-1)+ conjg(H(i1:cnt-1))*X(irow) 
       !Bottom left and bottom right Quadrants of A:     BR
       brow(N2p+1:irow+N2p-1)=brow(N2p+1:irow+N2p-1)-H(i1:cnt-1)*X(irow+N2p)
       
    enddo


    !*****************************************************************
    !MV product for Diagonal terms due to quasi-particle contribution
    !*****************************************************************
    if(.not. cj) w=omega
    if(cj)       w=conjg(omega)

    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1

    !Loop over diagonals
    do irow=rstart,rstop
       !Counters for electron and hole row indices
       !of resonant sub-matrix
       h_row_cnt=((irow-1)/Ne)+1
       e_row_cnt=irow-((h_row_cnt-1)*Ne)
       !Hole and electron indices for rowS OF  (h,e)
       h_row=(VBT+1) - h_row_cnt
       e_row= e_row_cnt + (CBB-1)

       !Add quasi-particle contribution to MV product (distributed)
       !Resonant (TL quadrant)
       brow(irow)=brow(irow)+( (eigen(ik,e_row)-eigen(ik,h_row)-w)* X(irow) )
       !Anti-resonant (BR quadrant)
       brow(irow+N2p)=brow(irow+N2p)+( (-eigen(ik,e_row)+eigen(ik,h_row)-w)* X(irow+N2p) )
    enddo

    !Sum elements of local output vector to give full vector
    call mpi_allreduce(brow,b,N,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    deallocate(brow)


  end subroutine mvp_cmplx_TD




  !********************************************************************
  !Symmetric matrix-vector product test routine
  !********************************************************************
  !Look at symmetric product: AX=b, then A^daggerX =b'
  !Add eigenvalue to diagonals: AX=b, then A^daggerX =b'
  !Try a few omega values on diagonals:  AX-wI=b, then A^daggerX-wI =b'
  !********************************************************************

  subroutine test_mvp_real(eigen)
    use var,     only: sp,dp,qp,Ne,Nh,eigen_cnt
    use parallel
    implicit none
    

    !****************
    !Declarations
    !****************
    
    !Local Variables 
    integer      N2p,N,Nlt,rstart,rstop,i,j
    complex(dp)  omega
    logical      cj,TD


    !Global/Modular Variables 
    integer      iter,flag
    real(dp)     ethr,anorm

    !Global Arrays
    real(dp),    allocatable,   intent(in) :: eigen(:,:)
    real(dp),    allocatable    :: H(:,:),H2p(:)
    complex(dp), allocatable    :: X(:),A(:,:),b(:),Xid(:)



    !****************
    !Main Routine
    !****************
    if(IO) write(*,*) 'Symmetric matrix mp testing'

    !2-particle resonant Hamiltonian dimensions 
    N2p=Ne*Nh

    !Dimensions of full BSE matrix A
    N=2*N2p

    !Dimension of packed, lower triangle
    Nlt=N2p*(0.5*(N2p-1)+1)

    !Distribute lower triangle of resonant H2p
    !Called from main routine instead

    !Transpose off
    !!cj=.false.
    !Tamm-Dancoff Approximation off
    TD=.false.
    !Frequency off
    omega=cmplx(4.5,0.001)  !0._dp
    
    !2-particle resonant matrix (real), packed form (real)
    !multiplying vector (complex) and 2-particle matrix (comple(
    allocate(H(N2p,N2p), H2p(Nlt), X(N),Xid(N), A(N,N))
    
    !Set up random, trial arrays with correct symmetry
    !noting that eigen and omega are currently not added to A
    call trial_arrays(H,H2p,N2p,A,Xid,omega,eigen)
    !!if(cj) A(:,:)=transpose(conjg(A(:,:)))
    
    
    !******************************************
    !Conjugate Gradient Test
    !******************************************
    ethr=1.e-12

    !RHS
    allocate(b(N))
    b=matmul(A,Xid)

    !Xid= vector used to produce RHS
    !X = CG attempt to reproduce/find Xid
    X(:)=0._dp

    !Testing CG with matmul
    omega=0._dp  !as A already has omega built into it
    flag=1
    call bcgsolve_parallel(N2p,H2p,eigen,omega,ethr,X,b,iter,anorm)  !, A,flag)

    !Write out convergence details
    if(IO)then
       write(*,*) 'CG using matmul:'
       write(*,*) 'Number of Iterations: ',iter
       write(*,*) 'Convergence (norm of residuals): ',anorm
    endif
    
    !Output X to file for CG with matmul
    do i=1,N
       if(IO)write(101,*) X(i)
    enddo
    
    !Testing CG with mvp
    X(:)=0._dp
    omega=cmplx(4.5,0.001)  
    flag=2
    write(*,*) rank, allocated(X)
    call bcgsolve_parallel(N2p,H2p,eigen,omega,ethr,X,b,iter,anorm)  !, A,flag)

    !Write out convergence details
    if(IO)then
       write(*,*) 'CG using parallel routine:'
       write(*,*) 'Number of Iterations: ',iter
       write(*,*) 'Convergence (norm of residuals): ',anorm
    endif

    !Output X to file (for CG with mvp)
    do i=1,N
       if(IO)write(201,*) X(i)
    enddo

    !Output Xid to file 
    do i=1,N
       if(IO)write(301,*) Xid(i)
    enddo




    !******************************************
    call mpi_barrier(MPI_COMM_WORLD,MPIierr)
    if(IO) write(*,*) 'CG test stops line 621, bse.f90'
    call stop_mpi


    !RHS
    allocate(b(N))
    b=0._dp
    
    !Compute and output matmul vector product
    if(IO)then
       b=matmul(A,X)
       write(*,*)'Matrix vector product with matmul'
       do i=1,N
          write(100,*) b(i)
       enddo
    endif
    
    !Matrix-vector product with my routine
    !Disable eigen contribution
    b=0._dp
    !Call needs changing to reflect new choice of routines 
    !call mvp_realH(H2p,N2p,X,eigen,omega,cj,TD,b)

    !Output b from mvp_real, on root process
    if(IO)then
       write(*,*)'Matrix vector product with mvp'
       do i=1,N
          write(200,*) b(i)
       enddo
    endif


  end subroutine test_mvp_real
 

  !************************************************************************
  !For a user input of total number of pair products in the BSE basis
  !the routine determines the number of occupied (Nh) and unoccupuied (Ne)
  !states systematically, such that pair products are built from the lowest
  !set of transitions
  !************************************************************************
  subroutine get_single_particle_basis(eigen,Np, Ne,Nh,tl)
    use var,      only: dp,Nev
    use sorting,  only: ORDVEC
    implicit none

    !*******************
    !Declarations
    !*******************

    !Global Variables
    integer, intent(in) :: Np
    integer, intent(out):: Ne,Nh

    !Global Arrays
    real(dp), intent(in), allocatable :: eigen(:,:)

    !Local Variables
    integer             :: ik,Nv,Nc,VBi,CBf,iv,ic,cnt,i,imin,imax,Nvp,Ncp
    real(dp)            :: tol,E_h,E_e

    !Local Arrays
    integer,         allocatable :: tl(:,:),tmp(:,:),indx(:)
    real(dp),        allocatable :: trnEn(:)


    !*******************
    !Main Routine
    !*******************
    !First part of routine Copies code from transition_details2 in spectra.f90

    !Finite system
    ik=1

    !Number of valence states.
    Nv=eigen_cnt(ik,1)
    !Number of conduction states 
    Nc=Nev-Nv
    
    !Set the number of valence and conduction states to sum over
    !to Np => One returns a total of Np**2 transitions but guarantees
    !that the lowest Np transtions are included

    !This is too much, particularly for the conduction states
    !Best way to do it would be to determine using the DOS
    !In the meantime, try 0.8Np for valence and 0.4Np for conduction
    !Or just kill the flags and see if any given run crashes 

    VBi=Nv-(int(0.8*Np)-1)
    CBf=Nv+1+(int(0.4*Np)-1)
    
    !Check enough states have been returned from eigensolver
    if(VBi<1)then
       if(IO)then
          write(*,*) 'Not enough valence states returned to compute pair-product basis.' 
          write(*,*) 'Code has stopped'
       endif
       call stop_mpi()
    endif
    if(Nc<int(0.4*Np))then
       if(IO)then
          write(*,*) 'Not enough conduction states returned to compute pair-product basis.' 
          write(*,*) 'Code has stopped'
       endif
       call stop_mpi()
    endif

    !Allocate arrays. If initialising, trnEn needs to be set
    !to a large transition energy such as 20 eV, else it might
    !get included in the lowest states when sorted 
    !allocate(trnEn(Np**2), tl(Np**2,2))

    !Number of valence and conduction state products considered
    Nvp=(Nv-VBi)+1
    Ncp=CBf-Nv
    allocate(trnEn(Nvp*Nc), tl(Nvp*Ncp,2))

    !**************************************************
    !Loop over the subset of lowest transitions
    !**************************************************
    cnt=0
    !Loop over selected valence band states
    do iv=VBi,Nv     
       !Loop over selected conduction band states
       do ic=Nv+1,CBf

          !Composite counter
          cnt=cnt+1

          !Transition energies
          trnEn(cnt)=(eigen(ik,ic)-eigen(ik,iv))

          !Transition Indices
          tl(cnt,1)=iv
          tl(cnt,2)=ic
          
       enddo
    enddo

     !*********************************************
    !Establish and store the lowest Np transitions
    !*********************************************
    !Index array
    allocate(indx(Nvp*Ncp))
    !Tolerance for values to be considered the same.
    tol=1.e-15
    !1=> Dimension of trnEn
    call ORDVEC(tol,1,Nvp*Ncp,trnEn,indx)
    !Returns ordered trnEn and an index to map the sorted values
    !to their unsorted array positions

    !Reorder tl
    allocate(tmp(Np,2))
    do i=1,Np
       tmp(i,:)=tl(indx(i),:)
       !if(IO) write(*,*) tmp(i,:)trnEn(indx(i))
    enddo
    deallocate(tl,indx)
    allocate(tl(Np,2))
    tl(:,:)=tmp(:,:)
    deallocate(tmp)


    !Largest state index (highest conduction state) in transitions
    imax=maxval(tl(:,2))
    !Number of conduction (electron) states
    Ne=imax-Nv

    !Lowest state index (largest hole index) in the transitions
    imin=minval(tl(:,1))
    !Number of valence (hole) states
    Nh=Nv-imin+1

    !Information on state energies 
    if(IO)then
       !Energy of max hole w.r.t. h1
       E_h=eigen(1,Nv)-eigen(1,imin)
       !Energy of max electron w.r.t. e1
       E_e=eigen(1,imax)-eigen(1,Nv+1)

       if(E_h<1.)then
          write(*,'(X,A,F7.2,A)')'The largest hole state in the BSE basis is ',    &
                                  E_h*1000.,' meV above h1'
       else
          write(*,'(X,A,F8.3,A)')'The largest hole state in the BSE basis is ',    &
                                  E_h,' eV above h1'
       endif
       if(E_e<1.)then
          write(*,'(X,A,F7.2,A)')'The largest electron state in the BSE basis is', &
                                  E_e*1000.,' meV above e1'
       else
          write(*,'(X,A,F8.3,A)')'The largest electron state in the BSE basis is', &
                                  E_e,' eV above e1'
       endif

       write(*,'(X,A,F6.3,A)')'Largest transition energy of a pair-product in the BSE basis is:'&
                              ,eigen(1,imax)-eigen(1,imin) ,' eV'
    endif


    !Return Nh and Ne
  end subroutine get_single_particle_basis



  !****************************************************************
  !Output the diagonal elements of the BSE kernel to a single file
  !That is K_(vc)(vc) i.e. (vc)=(v'c')
  !Note, K=H2p in rest of code
  !****************************************************************
  subroutine output_diagonal_elements(Ne,N2p,Np,tl,K)
    implicit none

    !******************
    !Declarations
    !******************

    !Global Variables
    integer,  intent(in) :: N2p,Ne,Np
    
    !Global Arrays
    integer,  allocatable, intent(in) :: tl(:,:)
    real(dp), allocatable, intent(in) :: K(:)
    
    !Local Variables
    integer   :: Nrow_l,kcnt,dcnt,row,col,ivc,iv,ic,h_row_cnt,e_row_cnt,VBT,i
    
    !Local Arrays
    integer(sp),  allocatable ::  row_size_local(:),row_size(:),displ(:)                  
    real(dp),     allocatable ::  K_diag_l(:),K_diag(:)   


    !******************
    !Main Routine
    !******************
    if(IO) write(*,'(X,A,I2,A)')'Outputting ',Np,' diagonal elements of the BSE kernel'

  
    !Number of rows per process
    Nrow_l=rstop-rstart+1

    !--------------------------------------------------
    !Store diagonal elements of packed LT kernel
    !--------------------------------------------------
    allocate(K_diag_l(Nrow_l))
    K_diag_l(:)=0._dp

    !Loop over diagonal elements of distributed LT kernel and store in new array
    kcnt=0
    dcnt=0
    do row=rstart,rstop
       col=row !Don't need to iterate col index when looping over diagonal

       !Composite index defines diagonal elements of packed LT array
       kcnt=kcnt+row
       !Fill local,digonal kernel using packed LT of kernel
       dcnt=dcnt+1
       K_diag_l(dcnt)=K(kcnt)

    enddo
    
   
    !--------------------------------------------------
    !Gather local arrays storing diagonal elements
    !--------------------------------------------------

    !Fill array which stores Nrow_l for all processes
    allocate(row_size_local(numprocs),row_size(numprocs))
    row_size_local(:)=0
    row_size_local(rank+1)=Nrow_l
    row_size(:)=0
    call mpi_allreduce(row_size_local,row_size,numprocs,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    deallocate(row_size_local)

    !Fill displacement index for use with allgatherv 
    allocate(displ(numprocs))
    !displ(:)=0
    !do i=2,numprocs
    !   displ(i)=displ(i-1)+row_size(i-1)
    !enddo

    !Use this definition of displ as code distributes rows starting at numprocs, numprocs-1,...,0
    !i.e. the opposite way round to how one would typically distribute 
    displ(:)=0
    do i=numprocs-1,1,-1
       displ(i)=displ(i+1)+row_size(i+1)
    enddo
    
    !Allocate global K_diag array
    allocate(K_diag(N2p))
    K_diag(:)=0._dp

    !Gather local arrays into single global array
    call mpi_allgatherv(K_diag_l, Nrow_l,          MPI_DOUBLE_PRECISION,  & !Local send
                        K_diag  , row_size, displ, MPI_DOUBLE_PRECISION,  & !Global receive
                        MPI_COMM_WORLD,MPIierr)

    !Free local array
    deallocate(K_diag_l)


    !--------------------------------------------------------------
    !Loop over the diagonal kernel elements in an order consistent
    !with the lowest-energy transitions, as these are thought to 
    !be the most physically relevant to the lowest exciton state
    !--------------------------------------------------------------
    VBT=eigen_cnt(1,1)
    !if(IO) write(*,*) 'VBT',VBT
    !call mpi_barrier(MPI_COMM_WORLD,MPIierr)
    
    !Open file and write header
    if(IO)then
       open(unit=001,file='kdiag.dat')
       write(001,*)'# iv,  ic,  diagonal,  K(eV)'
    endif

    do ivc=1,Np
       !Valence and conduction indices, consistent with the global
       !state index 1:Nev
       iv=tl(ivc,1)
       ic=tl(ivc,2)
       !if(IO) write(*,*) 'k routine:',iv,ic
       !Electron and hole indices
       h_row_cnt=VBT+1-iv
       e_row_cnt=ic-VBT
       !Composite index for 2-particle kernel
       row=e_row_cnt+(h_row_cnt-1)*Ne
       if(IO) write(001,'(3(I2,X),F12.8)') h_row_cnt,e_row_cnt,row,K_diag(row)
    enddo
    if(IO) close(001)

    !Free memory
    deallocate(K_diag)

 
  end subroutine output_diagonal_elements







  !*******************************************************
  !Construct trial H,X and A for matrix-vector testing
  !using correct symmetry properties
  !*******************************************************

  subroutine trial_arrays_complex(H,Hpkd,N2p,A,X,omega,eigen)
    use var,only: eigen_cnt
    implicit none

    !****************
    !Declarations
    !****************

    !Global Variables
    complex(dp), intent(in)  :: omega

    !Global Arrays
    real(dp),    allocatable, intent(in) :: eigen(:,:)
    complex(dp), allocatable :: H(:,:),A(:,:),X(:),Hpkd(:)

    !Local Variables
    integer      N,N2p,i,Nlt,Nlt_local,h_row_cnt,e_row_cnt,e_row,h_row,&
         VBT,CBB,cnt,j

    !Local Arrays
    real(dp),    allocatable :: C1(:,:),C2(:,:),C3(:),C4(:)


    !****************
    !Main Routine
    !****************


    !Dimensions of lower triangle of H
    Nlt=N2p*(0.5*(N2p-1)+1)

    !Dimensions of full matrix A
    N=2*N2p


    !Do on one process then broadcast to others
    if(IO)then

       !Allocate temporary arrays to hold random numbers
       allocate(C1(N2p,N2p),C2(N2p,N2p),C3(N),C4(N))

       !Fill H with random numbers
       call random_seed
       call random_number(C1)
       call random_number(C2)

       !Assign real + complex parts
       H=cmplx(C1,C2)

       !Use random complex matrix to construct a random Hermitian matrix
       H=0.5*(H+transpose(conjg(H)))

       !Set diagonals to (1,0). Img parts must be zero for Hermitian symmetry
       do i=1,N2p
          H(i,i)=cmplx(1.,0.)
       enddo

       !Multiplying vector X 
       call random_number(C3)
       call random_number(C4)
       X=cmplx(C3,C4)

    endif


    !Distribute H and X to other CPUs
    call mpi_bcast(H,N2p*N2p,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIierr)
    call mpi_bcast(X,N,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIierr)

    !Wait for all CPUs to get arrays before continuing
    call mpi_barrier(MPI_COMM_WORLD,MPIierr)

    !*******************************************
    !Define A 
    !*******************************************

    !Assign upper left quad with H
    A(1:N2p,1:N2p)=H(:,:)

    !Assign upper right quad with H
    A(1:N2p,N2p+1:N)=H(:,:)

    !Assign lower left quad with -(H)*
    A(N2p+1:N,1:N2p)=-1._dp*conjg(H(:,:))

    !Assign lower right quad with -(H)*
    A(N2p+1:N,N2p+1:N)=-1._dp*conjg(H(:,:))


!!$    !Add omega diagonal contribution to A
!!$    do i=1,N
!!$       A(i,i)=A(i,i)-omega
!!$    enddo

    !************************************************
    !Add quasi-particle contribution to the diagonals
    !************************************************
    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1

    do i=1,N2p
       !Counters for electron and hole row indices
       !of full matrix
       h_row_cnt=((i-1)/Ne)+1
       e_row_cnt=i-((h_row_cnt-1)*Ne)
       !Hole and electron indices for rows of full matrix
       h_row=(VBT+1) - h_row_cnt
       e_row= e_row_cnt + (CBB-1)

       !Add quasi-particle contributions to:
       !Resonant
       A(i,i)=A(i,i)+(eigen(1,e_row)-eigen(1,h_row))-omega
       !Anti-resonant
       A(i+N2p,i+N2p)=A(i+N2p,i+N2p)-(eigen(1,e_row)-eigen(1,h_row))-omega
    enddo



    !*******************************************
    !Define Hpkd (diagonals included later)
    !*******************************************

    !Construct packed lower triangle of H (allocated upper bound)
    !Trailing zeroes shouldn't be an issue as not looped over
    !in matrix-vector product
    cnt=0
    do i=rstart,rstop
       do j=1,i
          cnt=cnt+1
          Hpkd(cnt)=H(i,j)
       enddo
    enddo


  end subroutine trial_arrays_complex


  !*******************************************************
  !Construct trial H,X and A for matrix-vector testing
  !using correct symmetry properties (real)
  !*******************************************************

  subroutine trial_arrays_real(H,Hpkd,N2p,A,X,omega,eigen)
    use var,only: eigen_cnt
    implicit none

    !****************
    !Declarations
    !****************

    !Global Variables
    complex(dp), intent(in)  :: omega

    !Global Arrays
    real(dp),    allocatable, intent(in) :: eigen(:,:)
    real(dp),    allocatable  :: H(:,:),Hpkd(:)
    complex(dp), allocatable  :: A(:,:),X(:)

    !Local Variables
    integer      N,N2p,i,Nlt,Nlt_local,h_row_cnt,e_row_cnt,e_row,h_row,&
                 VBT,CBB,cnt,j

    !Local Arrays
    real(dp),    allocatable :: C1(:),C2(:)


    !****************
    !Main Routine
    !****************


    !Dimensions of lower triangle of H
    Nlt=N2p*(0.5*(N2p-1)+1)

    !Dimensions of full matrix A
    N=2*N2p


    !Do on one process then broadcast to others
    if(IO)then

       !Allocate temporary arrays to hold random numbers
       allocate(C1(N),C2(N))

       !Fill H with random numbers
       call random_seed
       call random_number(H)

       !Use random real matrix to construct symmetric matrix
       H=0.5*(H+transpose(H))

       !Set diagonals to 1
       do i=1,N2p
          H(i,i)=1._dp
       enddo

       !Multiplying random vector X (complex)
       call random_number(C1)
       call random_number(C2)
       X=cmplx(C1,C2)

    endif


    !Distribute H and X to other CPUs
    call mpi_bcast(H,N2p*N2p,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIierr)
    call mpi_bcast(X,N,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIierr)

    !Wait for all CPUs to get arrays before continuing
    call mpi_barrier(MPI_COMM_WORLD,MPIierr)


    !*******************************************
    !Define complex A using real H
    !*******************************************
    !Initialise A
    A(:,:)=cmplx(0._dp,0._dp)

    !Assign upper left quad with H
    A(1:N2p,1:N2p)=H(:,:)

    !Assign upper right quad with H
    A(1:N2p,N2p+1:N)=H(:,:)

    !Assign lower left quad with -H
    A(N2p+1:N,1:N2p)=-1._dp*H(:,:)

    !Assign lower right quad with -H
    A(N2p+1:N,N2p+1:N)=-1._dp*H(:,:)


!!$    !Add omega diagonal contribution to A
!!$    do i=1,N
!!$       A(i,i)=A(i,i)-omega
!!$    enddo

    !************************************************
    !Add quasi-particle contribution to the diagonals
    !************************************************
    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1

    do i=1,N2p
       !Counters for electron and hole row indices
       !of full matrix
       h_row_cnt=((i-1)/Ne)+1
       e_row_cnt=i-((h_row_cnt-1)*Ne)
       !Hole and electron indices for rows of full matrix
       h_row=(VBT+1) - h_row_cnt
       e_row= e_row_cnt + (CBB-1)

       !Add quasi-particle contributions to:
       !Resonant
       A(i,i)=A(i,i)+(eigen(1,e_row)-eigen(1,h_row))-omega
       !Anti-resonant
       A(i+N2p,i+N2p)=A(i+N2p,i+N2p)-(eigen(1,e_row)-eigen(1,h_row))-omega
    enddo


    !*******************************************
    !Define Hpkd (diagonals included later)
    !*******************************************

    !Construct packed lower triangle of H (allocated upper bound)
    !Trailing zeroes shouldn't be an issue as not looped over
    !in matrix-vector product
    cnt=0
    do i=rstart,rstop
       do j=1,i
          cnt=cnt+1
          Hpkd(cnt)=H(i,j)
       enddo
    enddo


  end subroutine trial_arrays_real




  !*********************************************************************************
  !
  !   Iterative solution of the linear system:
  !
  !                 ( A - wI ) * x = b         
  !
  !   where A is the complex 2-particle BSE Hamiltonian, w is a complex energy
  !   required for solution to converge, and x and b are complex vectors
  !
  !   Uses biorthogonal conjugate gradient method as described in
  !   D. A. H. Jacobs, A generalization of the conjugate-gradient method to
  !   solve complex systems, IMA Journal of Numerical Analysis 6, 447 (1986).
  !   Most subsequent relevant CG methods are based on this paper.
  !
  !   Original author: Felix Giustino, Berkeley Aug 2008
  !   Adapted by Alex Buccheri,        Oxford   Aug 2014
  !
  !*********************************************************************************
  !Inputs
  ! ethr      Threshold for convergence in norm of the residual
  ! N2p       Dimensions of resonant 2-particle H (mvp uses it to produce N)
  ! H2p       Lower triangle of resonant 2-particle H in packed form, for use in 
  !           distributed MV product (minus eigenvalue contribution to diagonals)
  ! K2        Optional coupling matrix for full BSE calculations 
  ! eval      Quasi-particle Eigenvalues (for use in mvp)
  ! Omega     or w, is a small, complex energy which is required for matrices with
  !           this specific symmetry to converge with biconjugate gradient
  ! b         Solution/right-hand side of Ax=b 
  ! x         Inital guess at x vector. Currently set to 0 in this routine
  !
  !Outputs
  ! x         Iterative solution for x. Represents 2-particle coefficients
  !*********************************************************************************

  subroutine bcgsolve_parallel(N2p,H2p,eval,omega,ethr,x,b,iter,anorm, K2)  !, Am,flag)
    use var, only: cgmaxitr
    implicit none

    !************************
    !Declarations
    !************************

    !Modular variables
    integer,     intent(in) :: N2p
    real(dp),    intent(in) :: ethr
    complex(dp), intent(in) :: omega

    !Local Variables
    integer              :: iter, ndmx, ndim, N, i
    real(dp)             :: anorm     
    complex(dp)          :: beta, alpha,a, c,ZDOTC
    logical              :: conv

    complex(dp), parameter :: cone = (1.d0, 0.d0)
    complex(dp), parameter :: czero = (0.d0, 0.d0)

    !Global/modular arrays
    real(dp),    allocatable, intent(in)    :: eval(:,:),H2p(:)
    real(dp),    allocatable, intent(in),   optional ::  K2(:)
    complex(dp), allocatable, intent(in)    :: b(:) !,H2p(:)
    complex(dp), allocatable, intent(inout) :: x(:)

    !Local Work Arrays
    complex(dp), allocatable :: r  (:), q  (:), p  (:), pold  (:)
    complex(dp), allocatable :: rt (:), qt (:), pt (:), ptold (:)
    complex(dp), allocatable :: rp (:), rtp (:)

    !Test quantities
    !integer         flag,j
    !complex(dp),    allocatable :: Am(:,:)



    !****************
    !Main Routine
    !****************
    !if(IO) write(*,*)'Conjugate gradient computing two-particle coefficients'

    !Dimensions of full, 2-particle matrix A (of Ax=b)
    N=2*N2p

    !Allocate conjugate gradient work arrays
    ndmx = n
    ndim = n
    allocate ( r (ndmx), q (ndmx), p (ndmx), pold (ndmx) )    
    allocate ( rt(ndmx), qt(ndmx), pt(ndmx), ptold(ndmx) )    
    allocate ( rp(ndmx), rtp(ndmx) )  
    r(:)    =0._dp
    rt(:)   =0._dp
    rp(:)   =0._dp
    rtp(:)  =0._dp
    p(:)    =0._dp
    pt(:)   =0._dp
    pold(:) =0._dp
    ptold(:)=0._dp
    q(:)    =0._dp
    qt(:)   =0._dp


    !**************************************************
    !Do conjugate gradient until convergence is reached 
    !or cgmaxitr is reached
    !**************************************************
    iter = 0
    anorm=0._dp
    conv=.false.

    do while ( iter.lt.cgmaxitr)
       !Iterate counter
       iter = iter + 1

       !Initialise variables for first iteration of CG
       !No reason for this to be in the while loop
       if (iter .eq. 1) then

          !Residual + transpose
          ! r  = b - A*x
          ! rt = bt - xtA^H 
          ! rt == conjg ( r )
          !On 1st iteration, if x=0 => r0=b
          !p0=r0 when no preconditioner is present


          !Compute Am.x and store in r
          !if(flag==1) r=matmul(Am,x)
          !Compute r=Am.x using H2p
          !if(flag==2)
          if(TD)       call mvp_TD(  H2p,   N2p,x,eval,omega,.false.,r)
          if(.not. TD) call mvp_full(H2p,K2,N2p,x,eval,omega,.false.,r)

          !p <- r
          p=r

          !Update r=Am.x   like so:   r <- r-b  where r was Ax
          !r  <--  alpha*b + r where alpha= (-1,0)
          call ZAXPY (ndim, -cone, b, 1, r, 1)
          !r  <--  -cone*r
          call ZSCAL (ndim, -cone, r, 1)

          !Transpose of Residual
          rt = conjg ( r )

          ! p  = inv(M) * r
          ! pt = conjg ( p )

          !p  <--  r
          call ZCOPY (ndim, r , 1, p , 1)
          pt = conjg ( p )

       endif


       !Check convergence, where ZDOTC=dot_product(r,r)
       !If |r| < ethr, exit 
       anorm = sqrt ( abs ( ZDOTC (ndim, r, 1, r, 1)  ) )
       if (anorm.lt.ethr) conv = .true.
       if ( conv ) exit
       !if(IO) write(*,*) iter,anorm



       !**************************************************************
       !Do conjugate gradient steps
       !Not sure why q and qt aren't within the if statement below
       !**************************************************************

       !**************************************************************
       ! Calculate q = A * p and qt = A^\dagger * pt
       !**************************************************************

       !  q  = A * p 
       !  qt = A^\dagger * pt

       ! calculate q=A.p
       !if(flag==1) q=matmul(Am,p)
       !if(flag==2) 
       !call mvp(H2p,N2p,p,eval,omega,.false.,TD,q)
       if(TD)       call mvp_TD(  H2p,   N2p,p,eval,omega,.false.,q)
       if(.not. TD) call mvp_full(H2p,K2,N2p,p,eval,omega,.false.,q)


       ! calculate qt=A^\dagger * pt
       !if(flag==1) qt=matmul(transpose(conjg(Am)),pt)
       !if(flag==2) 
       !call mvp(H2p,N2p,pt,eval,omega,.true.,TD,qt)
       if(TD)       call mvp_TD(  H2p,   N2p,pt,eval,omega,.true.,qt)
       if(.not. TD) call mvp_full(H2p,K2,N2p,pt,eval,omega,.true.,qt)



       if ( .not.conv) then

          !   rp  = r preconditioned [ inv(M) * r ]   
          !   rtp = r tidla preconditioned [ inv(M) * rt ]          
          !   Predconditioner not being used =>
          !   rp == r and rtp == rt 

          !No preconditioner used, so let rp <-- r. 
          call ZCOPY (ndim, r, 1, rp, 1)


          ! alpha = <rt|rp>/<pt|q>
          ! [ the denominator is stored for subsequent use in beta ]          
          a = dot_product(rt,rp)        !ZDOTC (ndim, rt, 1, rp, 1) 
          c = dot_product(pt,q)          !ZDOTC (ndim, pt, 1, q, 1) 
          alpha = a / c 


          !  x  = x  + alpha        * p
          !  r  = r  - alpha        * q 
          !  rt = rt - conjg(alpha) * qt
          call ZAXPY (ndim,  alpha,        p , 1, x , 1)
          call ZAXPY (ndim, -alpha,        q , 1, r , 1)
          call ZAXPY (ndim, -conjg(alpha), qt, 1, rt, 1)


          !  rp  = inv(M) * r
          !  rtp = inv(M) * rt          
          call ZCOPY (ndim, r , 1, rp, 1)
          call ZCOPY (ndim, rt, 1, rtp, 1)


          !  beta = - <qt|rp>/<pt|q>
          a = dot_product(qt,rp)    !ZDOTC (ndim, qt, 1, rp, 1)
          beta = - a / c 


          ! Update p and pt
          ! pold  = p
          ! ptold = pt
          call ZCOPY (ndim, p  , 1, pold  , 1)
          call ZCOPY (ndim, pt , 1, ptold , 1)

          !  p  = rp  +       beta  * pold
          !  pt = rtp + conjg(beta) * ptold
          call ZCOPY (ndim, rp , 1, p , 1)
          call ZCOPY (ndim, rtp, 1, pt, 1)
          call ZAXPY (ndim,       beta,  pold , 1, p , 1)
          call ZAXPY (ndim, conjg(beta), ptold, 1, pt, 1)          
       endif

    enddo
    !End CG self-consistent loop

    !Write out convergence details
    !if(IO)then
    !   write(*,*) 'CG Details:'
    !   write(*,*) 'Number of Iterations: ',iter
    !   write(*,*) 'Convergence (norm of residuals): ',anorm
    !endif

    !Write X to screen
    !if(IO)then 
    !   write(*,*)'Best-fit X'
    !   do i=1,N
    !      write(*,*) x(i)
    !   enddo
    !endif

    !Free memory
    deallocate ( r , q , p , pold, rt, qt, pt, ptold)


    return
  end subroutine bcgsolve_parallel



  !**********************************************************************
  !Essentially the same as absorption_spectrum2 but only loops
  !over Ne electrons and Nh holes, consistent with the 2-particle basis 
  !and returns the dipole transition elements in a vector
  !**********************************************************************
  !Inputs:
  ! eigen     Single-particle Eigenvalues
  ! evcr      Real Eigenvectors
  !
  !Output:
  ! b         Transition dipole matrix elements for use in BSE
  !           b(iv ic) = <ic|e.r|iv>, stored in vector  
  !**********************************************************************

  subroutine dipole_matrix_elements(eigen,evcr, b)
    use var,      only: dp,Natoms,orb_bas,eigen_cnt,NTT,NumNN,NN_list,Ne,Nh,delE_2p,&
                        eps,range
    use density,  only: lorentzian 
    use spectra,  only: PackedOrbital_TransMatrix
    implicit none

    !*************
    !Declarations
    !*************

    !Global Variables
    integer      midpnt

    !Global Arrays
    integer,     allocatable              :: Torb_indx(:,:)
    real(dp),    allocatable, intent(in)  :: eigen(:,:),evcr(:,:)
    real(dp),    allocatable              :: lor(:)
    complex(dp), allocatable              :: T_orb(:),b(:)

    !Local Variables
    integer      ik,N,VBT,CBB,VBf,CBf,ic,iv,Nv,icmp,istart,istop,i,j,N2p,&
                 cnt,ih,Npts,gi,gf,gmid,xi
    real(dp)     Emax,Emin,Erange,delE
    complex(dp)  T,T_local

    !Local Arrays 
    real(dp),    allocatable :: spectrum(:)

    !*********************************************************
    !Main Routine
    !*********************************************************

    !Inform user that routine has been called
    if(IO) write(*,*) 'Computing dipole transition matrix elements for RHS'

    !Number of orbital states in the system
    N=Natoms*orb_bas

    !Finite systems only
    ik=1

    !Absorption spectrum range determined by Ne and Nh
    !No need for determining the energy range 
    
    !Initialise b
    b(:)=0._dp


    !***********************************************************************
    !Valence (hole) and conduction (electron) state indices
    !***********************************************************************

    !In 2-particle basis, the hole index starts at the VBT then counts 
    !down (i.e. to higher VB states)

    !VBT and CBB indices
    VBT=eigen_cnt(ik,1)
    CBB=VBT+1

    !Final valence state (hole)
    VBf=VBT-Nh+1

    !Final conduction state (electron)
    CBf=CBB+Ne-1

   
    !***********************************************************************
    !Produce single-particle specturm in this routine for direct comparison 
    !with 2-particle spectrum
    !***********************************************************************
    !Use same spacing as 2p-spectrum
    delE=delE_2p

    !Produce Lorentzian profile
    call lorentzian(delE,eps,range,midpnt,lor)

    !Max energy transition
    Emax=eigen(ik,CBf)-eigen(ik,VBf)

    !Min energy SP transition = Bandgap - 0.5eV
    Emin=eigen(ik,CBB)-eigen(ik,VBT)-0.5
    
    !Define energy range
    Erange=Emax-Emin

    !With energy range covered by the 2-particle basis states, and
    !energy spacing set in file, determine the number of idscrete points
    !at which to store the spectrum 
    Npts=ceiling(Erange/delE) +1

    !Allocate spectrum array size
    allocate(spectrum(Npts)) 
    spectrum(:)=0._dp


    !**********************************************************
    !Transition Matrix in orbital basis
    !**********************************************************

    !Orbital transition matrix is diagonal in this approximation
    !hence dimension is exact
    if(NTT==1) allocate(T_orb(N), Torb_indx(N,2))

    !Include all on-site transitions. This is an upper bound
    !-Atomic selection rules apply
    if(NTT==2) allocate(T_orb(N*orb_bas*orb_bas),Torb_indx(N*orb_bas*orb_bas,2))

    !Include transitions between adjacent sites. Upper bound
    !as not all atoms have NumNN nearest neighbours
    if(NTT==3)then
       allocate(T_orb(N*NumNN*orb_bas*orb_bas),Torb_indx(N*NumNN*orb_bas*orb_bas,2))
    endif

    !Generate packed orbital transition dipole matrix   
    call PackedOrbital_TransMatrix(T_orb,Torb_indx)


    !********************************************************************
    !Transition Matrix in single-particle basis
    !********************************************************************
 

    !NOTE: valence and conduction loop ordering 
    !must be consistent with BSE routine 


    cnt=0
    !do iv=VBT,VBf,iv-1  Loop won't work for some reason
    !Loop over selected valence band states
    do ih=1,Nh
       
       !Valence index
       iv=VBT-(ih-1)

       !Loop over selected conduction band states
       do ic=CBB,CBf

          !Output electron and hole indices
          !write(*,*) iv,ic

          !Composite (iv,ic) index
          cnt=cnt+1


          !***********************************************************
          !Lorentzian Energy Profile
          !***********************************************************
          
          !Convert the transition energy, (E_c-E_v), to its corresponding index
          !w.r.t. the total number of points used to sample the energy range (Gnpts)
          gmid=anint( (eigen(ik,ic)-eigen(ik,iv)-Emin)/delE )
          
          !Indices for Lorentzian
          gi = gmid-(midpnt-1)
          gf = gmid+(midpnt-1)
          
          !***********************************************
          !Prevent out-of-bounds due to Lorentzian tails
          !***********************************************
          !Define 1st value of Lorentzian profile
          xi=1
          !If any Lorentzian points falls outside of the global
          !energy range, cap the global value at the maximum value of the energy array
          if(gf > Npts) gf = Npts
          !If the Lorentzian totally exceeds the max energy value, ignore it
          if(gi > Npts) cycle
          !If the minimum value of gi is less than 1, ignore it
          !If any value of gi is less than 1, cut it
          if(gi<1) then
             xi=abs(gi)+1
             gi=1
          endif


          !********************************************
          !Sum the transition matrix for given (iv,ic)
          !********************************************

          !Distribute loop indices
          call distribute_loop(size(T_orb),istart,istop)
          T_local=0._dp
          T=0._dp

          !Distributed loop
          do icmp=istart,istop  
             !Retrieve global row and column indices for the orbital transition matrix
             !Row index
             j=Torb_indx(icmp,1)
             !Column index
             i=Torb_indx(icmp,2)

             !Sum the states which comprise the Transition matrix element:
             !<ic|e.r|iv> = sum(ia,iorb,ja,jorb)[c*_ic;(ja,jorb) c_iv;(ia,iorb) <ja,jorb|e.r|ia,iorb>]
             !T_local=T_local +( (eigen(ik,ic)-eigen(ik,iv))*evcr(j,ic)*evcr(i,iv)*T_orb(icmp) )
             T_local=T_local +( evcr(j,ic)*evcr(i,iv)*T_orb(icmp) )

          enddo

          !T <- Sum[ T_local ] 
          call mpi_allreduce(T_local,T,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIierr)


          !Store transition matrix element in single-particle basis: b(1:N2p) = -<ic|e.r|iv>
          !for use in RHS of 2-particle equation
          b(cnt)=-T

          !Single-particle Spectrum
          !spectrum(gi:gf) = spectrum(gi:gf) + &
          !                  ((eigen(ik,ic)-eigen(ik,iv)) * (T*conjg(T)) * lor(xi:gf-gi+1)) 
          !spectrum(gi:gf) = spectrum(gi:gf) + (T*conjg(T)*lor(xi:gf-gi+1))
          spectrum(gi:gf) = spectrum(gi:gf) + (real(T)*real(T)*lor(xi:gf-gi+1))

        

       enddo
    enddo

    !Output un-normalised SP spectrum
    if(IO)then
       open(unit=002,file='absp_CdSe.dat')
       do i=1,Npts
          write(002,*) Emin+(i-1)*delE,spectrum(i)
       enddo
       close(002)
       write(*,*) 'absp_CdSe.dat written to file'
    endif

    !Resonant Hamiltonian basis
    N2p=Ne*Nh

    !Fill lower half of the vector and multiply by -1 
    !b(N2p+1:2*N2p) = <ic|e.r|iv>
    b(N2p+1:2*N2p)=-1._dp*b(1:N2p)


  end subroutine dipole_matrix_elements



  !Equivalent to routine above, except Lorentizan broadening varies as a
  !function of the transition energy 
  subroutine var_dipole_matrix_elements(eigen,evcr, b)
    use var,      only: dp,Natoms,orb_bas,eigen_cnt,NTT,NumNN,NN_list,Ne,Nh,delE_2p,&
                        eps,range
    use density,  only: lorentzian 
    use spectra,  only: PackedOrbital_TransMatrix
    implicit none

    !*************
    !Declarations
    !*************
    !Global Variables
    integer      midpnt

    !Global Arrays
    integer,     allocatable              :: Torb_indx(:,:)
    real(dp),    allocatable, intent(in)  :: eigen(:,:),evcr(:,:)
    real(dp),    allocatable              :: lor(:)
    complex(dp), allocatable              :: T_orb(:),b(:)

    !Local Variables
    integer      ik,N,VBT,CBB,VBf,CBf,ic,iv,Nv,icmp,istart,istop,i,j,N2p,&
                 cnt,ih,Npts,gi,gf,gmid,xi
    real(dp)     Emax,Emin,Erange,Eg,delE
    complex(dp)  T,T_local

    !Local Arrays 
    real(dp),    allocatable :: spectrum(:)

    !*********************************************************
    !Main Routine
    !*********************************************************

    !Inform user that routine has been called
    if(IO) write(*,*) 'Computing dipole transition matrix elements for RHS'
    !Number of orbital states in the system
    N=Natoms*orb_bas
    !Finite systems only
    ik=1
    !Initialise b
    b(:)=0._dp
    !VBT and CBB indices
    VBT=eigen_cnt(ik,1)
    CBB=VBT+1
    !Final valence state (hole)
    VBf=VBT-Nh+1
    !Final conduction state (electron)
    CBf=CBB+Ne-1
   
    !***********************************************************************
    !Produce single-particle specturm in this routine for direct comparison 
    !with 2-particle spectrum
    !***********************************************************************

    !Absorption on-set(==band-gap at single-particle level)
    Eg=eigen(1,CBB)-eigen(1,VBT)
    !Max energy transition
    Emax=eigen(ik,CBf)-eigen(ik,VBf)
    !Min energy SP transition = Bandgap - 0.5eV
    Emin=eigen(ik,CBB)-eigen(ik,VBT)-0.5
    !Define energy range
    Erange=Emax-Emin
    !Same spacing as used by 2p-spectrum  
    delE=delE_2p
    !Number of points in energy range
    Npts=ceiling(Erange/delE) +1
    !Allocate spectrum array size
    allocate(spectrum(Npts)) 
    spectrum(:)=0._dp

    !**********************************************************
    !Transition Matrix in orbital basis
    !**********************************************************
    if(NTT==1) allocate(T_orb(N), Torb_indx(N,2))
    if(NTT==2) allocate(T_orb(N*orb_bas*orb_bas),Torb_indx(N*orb_bas*orb_bas,2))
    if(NTT==3)then
       allocate(T_orb(N*NumNN*orb_bas*orb_bas),Torb_indx(N*NumNN*orb_bas*orb_bas,2))
    endif
    !Generate packed orbital transition dipole matrix   
    call PackedOrbital_TransMatrix(T_orb,Torb_indx)

    !********************************************************************
    !Transition Matrix in single-particle basis
    !********************************************************************
    !NOTE: valence and conduction loop ordering 
    !must be consistent with BSE routine 
    cnt=0
    do ih=1,Nh
       iv=VBT-(ih-1)
       do ic=CBB,CBf

          !Composite (iv,ic) index
          cnt=cnt+1
          !Produce variably-broadened Lorentzian profile with limits
          call lorentzian(delE,veps(eigen(ik,ic)-eigen(ik,iv),Eg),range,midpnt,lor)
          gmid=anint( (eigen(ik,ic)-eigen(ik,iv)-Emin)/delE )
          gi = gmid-(midpnt-1)
          gf = gmid+(midpnt-1)
          xi=1
          if(gf > Npts) gf = Npts
          if(gi > Npts) cycle
          if(gi<1) then
             xi=abs(gi)+1
             gi=1
          endif

          !Sum the transition matrix for given (iv,ic)
          call distribute_loop(size(T_orb),istart,istop)
          T_local=0._dp
          T=0._dp
          do icmp=istart,istop  
             j=Torb_indx(icmp,1)
             i=Torb_indx(icmp,2)
             T_local=T_local +( evcr(j,ic)*evcr(i,iv)*T_orb(icmp) )
          enddo
          call mpi_allreduce(T_local,T,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIierr)

          !For use in 2p-spectrum
          b(cnt)=-T
          !Single-particle Spectrum
          spectrum(gi:gf) = spectrum(gi:gf) + (real(T)*real(T)*lor(xi:gf-gi+1))

       enddo
    enddo

    !Output un-normalised SP spectrum
    if(IO)then
       open(unit=002,file='absp_CdSe.dat')
       do i=1,Npts
          write(002,*) Emin+(i-1)*delE,spectrum(i)
       enddo
       close(002)
       write(*,*) 'absp_CdSe.dat written to file'
    endif

    !Resonant Hamiltonian basis
    N2p=Ne*Nh
    b(N2p+1:2*N2p)=-1._dp*b(1:N2p)

  end subroutine var_dipole_matrix_elements





  !*******************************************
  !Dynamic File/header for 2p-spectrum
  !*******************************************
  !Energy Range        [emin:max]
  !Energy Spacing      dE
  !# sampling points   Npts
  !Broadening          wi
  !Absorption(w)       e2
  !*******************************************
  subroutine output_2p_spectrum(emin,emax,dE,Npts,wi,e2)
    use var, only: epol,Natoms,nx,ny,nz,system,radius,delEa,delEb,delEc,delEd,&
                   shift_hyb,compounds,ic,Ne,Nh,Lx,Ly,Lz,NumML
    implicit none

    
    !*************
    !Declarations
    !*************

    !Modular Variables
    integer,  intent(in) :: Npts
    real(dp), intent(in) :: emin,emax,dE,wi

    !Input Array
    real(dp),        allocatable,  intent(in) :: e2(:)

    !Local Variables
    integer      ::  i,ios
    character    ::  nx_lab*3,ny_lab*3,nz_lab*3,sys*12,Natoms_lab*5,filename*70,&
                     dat_typ*7

    !*************
    !Main Routine
    !*************

    !Create strings from numerical variables
    write(nx_lab,'(I3)') nx 
    write(ny_lab,'(I3)') ny 
    write(nz_lab,'(I3)') nz
    write(Natoms_lab,'(I5)') Natoms

    !Data 
    dat_typ='2p_absp'
    sys='finite'


    !**********************************************
    !Determine file name based on system parameters
    !**********************************************

    !Arbitrary Cell
    if(system == 'arb')then
       !Filename
       filename= dat_typ//'_'//trim(adjustl(sys))//'_'//trim(adjustl(nx_lab))&
                 //','//trim(adjustl(ny_lab))//','//trim(adjustl(nz_lab))//'.dat'
    endif
    
    !Sphere or nanoplatelet
    if(system=='sph' .or. system=='npl')then
       !Filename
       filename= dat_typ//'_'//trim(adjustl(sys))//'_'//trim(adjustl(system))//'_'//&
                 trim(adjustl(Natoms_lab))// '_'//trim(adjustl(compounds(ic)))//'.dat'
    endif


    !**********************
    !Write to file
    !**********************
    
    !Open File
    open(unit=001, iostat=ios, file=trim(adjustl(filename)) )
    if(ios /= 0) write(*,*)  trim(adjustl(filename)),' can''t be created/openned'

    
    !***********
    !Header
    !***********

    !System specific
    if(system=='arb')then
       write(001,'(A,X,A4,X,A3,X,I4,X,A)') '#', compounds(ic),system,Natoms,'atoms'
       if(shift_hyb .eqv. .true.)write(001,'(A,X,4(F12.6,X))')  '# Passivation',delEa,delEb,delEc,delEd
    elseif(system=='sph')then
       write(001,'(A,X,A4,X,A3,X,I4,X,A,X,A,F5.2)') &
            '#', compounds(ic),system,Natoms,'atoms','diameter(nm):',(0.2*radius)
       if(shift_hyb .eqv. .true.)write(001,'(A,X,4(F13.6,X))')  '# Passivation',delEa,delEb,delEc,delEd
    elseif(system=='npl')then
       write(001,'(A,X,A4,X,A3,X,I4,X,A,X,A,I2)') &
            '#', compounds(ic),system,Natoms,'atoms','Thickness(ML):',int(NumML)
       write(001,'(A,X,3(F6.2,X))') 'Nanoplatelet (x,y,z) dimensions in nm:',Lx*0.1,Ly*0.1,Lz*0.1
       if(shift_hyb .eqv. .true.)write(001,'(A,X,4(F13.6,X))')  '# Passivation',delEa,delEb,delEc,delEd
    endif

    !Applicable to all spectra
    write(001,'(2(A,F6.3,X),A,3(F4.1))') '# Energy spacing',dE,'broadening',wi,'Polarisation vector:',epol(:)
    write(001,'(A,I3,X,A,I3)') '# Two-particle Basis. Ne ',Ne,'Nh ',Nh
    
    
    !**********************
    !Write spectrum to file
    !**********************
    do i=1,Npts
       write(001,*) emin+(i-1)*dE,e2(i)
    enddo
    close(001)
    if(IO) write(*,*) 'Written spectrum to: ', filename

  end subroutine output_2p_spectrum





  !End of Module
end module bse
