!**************************************************
!Module for constructing the BSE Hamiltonian 
!**************************************************
module kernel
  !--------------------------
  !Modules used
  !--------------------------
  use var,     only: dp
  use parallel  

  !--------------------------
  !Modular variables & arrays
  !--------------------------
  implicit none
  private
  integer   ::   N,VBT
  real(dp)  ::   ieps,unt
  integer,  allocatable  :: atom_species(:)
  real(dp), allocatable  :: U_c(:,:),U_x(:,:),c_vi(:),c_ci(:),c_vj(:),c_cj(:)

  !---------------------
  !Interfaces
  !---------------------
  !interface construct_orthogonal_H2p
  !   module procedure  construct_orthogonal_H2p_packed,construct_orthogonal_H2p_2D
  !end interface construct_orthogonal_H2p

  !------------------------------------------
  !Subroutines 
  !------------------------------------------
  public :: construct_orthogonal_H2p_packed

contains 
  

  !******************************************************************
  !Wrapper for Hamiltonian construction calls 
  !******************************************************************
  subroutine construct_orthogonal_H2p_packed(H2p,K2,TD,rstart,rstop)
    implicit none
    !-----------------
    !Declarations
    !-----------------
    integer,               intent(in)    :: rstart,rstop
    logical,               intent(in)    :: TD
    real(dp), allocatable, intent(inout) :: H2p(:),K2(:)
    integer                              :: i,Nlt

    !-----------------
    !Main Routine
    !-----------------
    !Set up eigenvector arrays, dielectric function and read in orbital integrals
    call initialise_2p_workspace()

    !Number of elements in packed lower triangle
    !Elements per row == row index
    Nlt=0
    do i=rstart,rstop
       Nlt=Nlt+i
    enddo

    !Resonant two-particle Hamiltonian
    allocate(H2p(Nlt))
    H2p(:)=0._dp

    !In both instances, the diagonal transition energies are excluded
    !Tamm-Dancoff Approximation. Turn off coupling submatrices
    !H2p = K_(vc)(v'c')
    if(TD)then
       call construct_orthogonal_H2p(H2p,rstart,rstop)
    endif

    !Full BSE calculation requires coupling submatrices
    !H2p = K_(vc)(v'c') and K2 = K(vc)(c'v')
    if(.not. TD)then
       allocate(K2(Nlt))
       K2(:)=0._dp
       call construct_orthogonal_H2p_K2(H2p,K2,rstart,rstop)
    endif

  end subroutine construct_orthogonal_H2p_packed



  !******************************************************************
  !Set up modular variables/arrays used in two-particle calculations
  !******************************************************************
  subroutine initialise_2p_workspace()
    use var, only: Natoms,orb_bas,eigen_cnt,ic,Eg,Ne,Nh,atom_type,system
    implicit none
    integer :: ia
  
    !Dimensions of the system (==total number of orbital states)
    !in the absence of spin
    N=Natoms*orb_bas

    !VBT
    VBT=eigen_cnt(1,1)

    !Read in tabulated, two-centre orbital integrals (in eV)
    call orbital_integral_energies(U_c,U_x)
    
    !Dielectric screening
    if(IO .and. ic/=19)then
       write(*,*)'No dielectric constant stored for this material.'
       write(*,*)'Code has stopped'
       call stop_mpi()
    endif

    !Bulk, high-frequency dielectric constant for CdSe (Bulk static 9.7)
    !used for NPL in absence of dielectric function 
    if(ic==19 .and. system=='npl')then  
       ieps=(1._dp/6.2)  
       if(IO) write(*,*) 'Screening with bulk eps_inf' 

    !High-frequency dielectric as a function of NC diameter (enters via band-gap)
    elseif(ic==19 .and. system=='sph')then
       ieps=(1._dp/eps_inf_dot(Eg))
    endif

    !Units. e**2/(4 pi eps_0 a_0) = 1 Ha
    !Convert atomic positions to Bohr and Hartree to eV
    !unt=(0.5291772109217*27.2113850560)
    unt=14.399644849251

    !1D arrays to hold TB coefficients
    allocate(c_vi(Nh),c_vj(Nh),c_ci(Ne),c_cj(Ne))
    c_vi(:)=0._dp
    c_vj(:)=0._dp
    c_ci(:)=0._dp
    c_cj(:)=0._dp
    
    !In U arrays, ia=1 -> anion and ia=2 -> cation which is the opposite
    !of the globally defined atom_species. Either need to swap the block 
    !ordering in U arrays or redefine atom_species.
    !Simpler to redefine atom_species as a local array
    allocate(atom_species(Natoms))
    atom_species(:)=0
    do ia=1,Natoms
       !Cation = 2
       if(atom_type(ia)==1) atom_species(ia)=2
       !Anion = 1
       if(atom_type(ia)==-1)atom_species(ia)=1
    enddo
    
  end subroutine initialise_2p_workspace



  !**************************************************
  !Size-dependent, high-frequency dielectric constant 
  !See PHYSICAL REVIEW B, VOLUME 63, 195318
  !Electron-hole correlations in semiconductor 
  !quantum dots with tight-binding wave functions
  !**************************************************
  real(dp) function eps_inf_dot(Eg)
    use var, only:ic
    implicit none
    real(dp),   intent(in) :: Eg
    real(dp) :: Eg_b,eps_inf_b,delta
    !User Flag
    if(ic/=19 .and. IO) write(*,*) 'Constants in eps_inf function only defined for CdSe.'
    !CdSe Values
    !TB bulk single-particle gap 
    Eg_b=1.90
    !Bulk, high-frequency dielectric constant
    eps_inf_b=6.2
    !Energy of second peak in bulk abs spectrum w.r.t. Eg_b
    delta=5.72-Eg_b
    eps_inf_dot= 1._dp + (eps_inf_b-1._dp)*((Eg_b+delta)**2./(Eg+delta)**2.)
    if(IO) write(*,*) 'Size-dependent eps_inf = ',eps_inf_dot
    return
  end function eps_inf_dot


  
  !**********************************************************************
  !Construct the 2-particle kernel in an orthogonal, atomic orbital basis
  !**********************************************************************
  subroutine construct_orthogonal_H2p(H2p,rstart,rstop)
    use var, only: Natoms,NumNN,orb_bas,eigen_cnt,r,Ne,Nh,al,evcr,NN_list,ic,atom_type,bond_length,Eg
    implicit none

    !******************
    !Declarations
    !******************
    !Global Variables
    integer,      intent(in)          :: rstart,rstop

    !Global Arrays
    real(dp),     allocatable,        intent(inout) :: H2p(:)

    !Local Variables
    integer  ::    N,VBT,ia,ja,ia_NN,ia_l,ja_l,iorb,jorb,i,j,i_l,j_l,&
                   Hcnt,row,col,h_row_cnt,h_col_cnt,e_row_cnt,e_col_cnt, cnt
    real(dp) ::    t1_local,t2_local,Uc,Ux,d_2N,denom

 
    !******************
    !Main Routine
    !******************
    call cpu_time(t1_local)
    if(IO)then
       write(*,*) 'Using Tamm-Dancoff Approximation'
       write(*,*) 'Constructing H2p in an orthogonal basis'
    endif
 
  
    !****************************************************************
    !Integrals between on-site and NN terms explicitly computed 
    !Monopole term of multipole approximation is used for Coulomb-like 
    !integrals such that they reduce to 1/|r_ia - r_ja|
    !Exchange-like terms ignored beyond NN
    !****************************************************************

    !Sum on-site and NN terms 
    do ia=1,Natoms
       !Local atom index 
       ia_l=atom_species(ia)

       !Sum over NN of atom ia
       do ia_NN=1,NumNN
          ja=NN_list(ia,ia_NN)
          if(ja==0) cycle
          ja_l=atom_species(ja)
     

          do iorb=1,orb_bas
             !Orbital index for a system comprised of 1 anion and 1 cation
             i_l=iorb+(ia_l-1)*orb_bas
             !Global orbital index
             i=iorb+(ia-1)*orb_bas
             !Assign TB coefficients to 1D arrays 
             c_vi(:) = evcr(i, VBT:VBT-Nh+1:-1) 
             c_ci(:) = evcr(i, VBT+1:VBT+Ne)    

             do jorb=1,orb_bas
                j_l=jorb+(ja_l-1)*orb_bas
                j=jorb+(ja-1)*orb_bas
                c_vj(:) = evcr(j, VBT:VBT-Nh+1:-1) 
                c_cj(:) = evcr(j, VBT+1:VBT+Ne)    
                
                !Integrals between orbitals and bare Coulomb potential
                Uc=U_c(i_l,j_l)
                !Only on-site Coulomb integral explicitly calculated
                Ux=U_x(i_l,j_l)

 
                !Counter for packed matrix H2p
                Hcnt=0
                !Row loop distributed => H2p is distributed
                do row=rstart,rstop
                   
                   !Counters for electron and hole row indices
                   h_row_cnt=((row-1)/Ne)+1
                   e_row_cnt=row-((h_row_cnt-1)*Ne)
                   
                   !Loop over lower triangle only
                   do col=1,row
                    
                      !Iterate H2p packed index counter
                      Hcnt=Hcnt+1
                      !Counters for electron and hole col indices
                      h_col_cnt=((col-1)/Ne)+1
                      e_col_cnt=col-((h_col_cnt-1)*Ne)
                      
                      !Packed Two-particle Hamiltonian (technically just
                      !the kernel elements)
                      H2p(Hcnt) = H2p(Hcnt) + &
		    
                                  c_vi(h_row_cnt)*c_ci(e_row_cnt)*c_vj(h_col_cnt)* &
                                  c_cj(e_col_cnt)*(2._dp*Uc-ieps*Ux)+              &

                                  c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vi(h_col_cnt)* &
                                  c_cj(e_col_cnt)*(2._dp*Ux-ieps*Uc)+              &

                                  c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vj(h_col_cnt)* &
                                  c_ci(e_col_cnt)*(2._dp-ieps)*Ux                                         
           
                   enddo
                enddo
                !Shut loops over (vc)(v'c')

             enddo
          enddo
       enddo
    enddo


    !Free Memory
    deallocate(U_c,U_x)

    !2nd NN distance - tolerance. Shouldn't need it but just to be safe
    d_2N=al/sqrt(2._dp)-0.05*bond_length(ic)

    !Continue summing terms in H2p
    !Beyond NN, only retain terms multiplied by Uc=(ii|v|jj)
    do ia=1,Natoms
       do ja=1,Natoms
          !|r-r'|
          denom=sqrt(dot_product(r(ia,:)-r(ja,:),r(ia,:)-r(ja,:)))
          !If atoms ia and ja are NN or same atom, cycle
          if(denom<d_2N) cycle 
          !Inter-atomic integral (ii|v|jj) =~ 1/|r_ia-r_ja|
          Uc=unt/denom
    
          do iorb=1,orb_bas
             i=iorb+(ia-1)*orb_bas
             c_vi(:) = evcr(i, VBT:VBT-Nh+1:-1) 
             c_ci(:) = evcr(i, VBT+1:VBT+Ne)    

             do jorb=1,orb_bas
                j=jorb+(ja-1)*orb_bas
                c_vj(:) = evcr(j, VBT:VBT-Nh+1:-1) 
                c_cj(:) = evcr(j, VBT+1:VBT+Ne)    


                Hcnt=0
                do row=rstart,rstop
                   h_row_cnt=((row-1)/Ne)+1
                   e_row_cnt=row-((h_row_cnt-1)*Ne)
                   
                   do col=1,row                    
                      Hcnt=Hcnt+1
                      h_col_cnt=((col-1)/Ne)+1
                      e_col_cnt=col-((h_col_cnt-1)*Ne)
                      
                      !Packed Two-particle Hamiltonian (technically just
                      !the kernel elements) 
                      H2p(Hcnt) = H2p(Hcnt) +  &
                                  ( 2._dp*c_vi(h_row_cnt)*c_ci(e_row_cnt)*c_vj(h_col_cnt)*c_cj(e_col_cnt)  &
                                    -ieps*c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vi(h_col_cnt)*c_cj(e_col_cnt) )*Uc

                   enddo
                enddo
                
             enddo
          enddo
       enddo
    enddo

   
    !Output routine runtime
    call cpu_time(t2_local)
    if(IO)write(*,'(A,F20.10,A)') 'Time taken to construct H2p = ',(t2_local-t1_local)/60._dp,' mins'


  end subroutine construct_orthogonal_H2p




  !**********************************************************************
  !Construct the 2-particle kernels in an orthognal, atomic orbital basis
  !H2p=K_(vc)(v'c') and K2=K_(vc)(c'v'), required for the definition of
  !the full BSE Hamiltonian
  !**********************************************************************
  subroutine construct_orthogonal_H2p_K2(H2p,K2,rstart,rstop)
    use var, only: Natoms,orb_bas,NumNN,eigen_cnt,r,Ne,Nh,al,evcr,NN_list,ic,atom_type,bond_length,Eg
    implicit none

    !******************
    !Declarations
    !******************

    !Global Variables
    integer,      intent(in)          :: rstart,rstop

    !Global Arrays
    real(dp),     allocatable,        intent(inout) :: H2p(:),K2(:)
  
    !Local Variables
    integer  ::    ia,ja,ia_NN,ia_l,ja_l,iorb,jorb,i,j,i_l,j_l,&
                   Hcnt,row,col,h_row_cnt,h_col_cnt,e_row_cnt,e_col_cnt,cnt
    real(dp) ::    t1_local,t2_local,Uc,Ux,d_2N,denom,c1,c2,c3

 
    !******************
    !Main Routine
    !******************
    call cpu_time(t1_local)
    if(IO)then
       write(*,*) 'Full BSE calculation'
       write(*,*) 'Constructing H2p and K2 in an orthogonal basis'
    endif


    !****************************************************************
    !Integrals between on-site and NN terms explicitly computed 
    !Monopole term of multipole approximation is used for Coulomb-like 
    !integrals such that they reduce to 1/|r_ia - r_ja|
    !Exchange-like terms ignored beyond NN
    !****************************************************************

    !Sum on-site and NN terms 
    do ia=1,Natoms
       !Local atom index 
       ia_l=atom_species(ia)

       !Sum over NN of atom ia
       do ia_NN=1,NumNN
          ja=NN_list(ia,ia_NN)
          if(ja==0) cycle
          ja_l=atom_species(ja)
     
          do iorb=1,orb_bas
             !Orbital index for a system comprised of 1 anion and 1 cation
             i_l=iorb+(ia_l-1)*orb_bas
             !Global orbital index
             i=iorb+(ia-1)*orb_bas
             !Assign TB coefficients to 1D arrays 
             c_vi(:) = evcr(i, VBT:VBT-Nh+1:-1) 
             c_ci(:) = evcr(i, VBT+1:VBT+Ne)    

             do jorb=1,orb_bas
                j_l=jorb+(ja_l-1)*orb_bas
                j=jorb+(ja-1)*orb_bas
                c_vj(:) = evcr(j, VBT:VBT-Nh+1:-1) 
                c_cj(:) = evcr(j, VBT+1:VBT+Ne)    
                
                !Integrals between orbitals and bare Coulomb potential
                Uc=U_c(i_l,j_l)
                !Only on-site Coulomb integral explicitly calculated
                Ux=U_x(i_l,j_l)
 
                !Counter for packed matrix H2p
                Hcnt=0
                !Row loop distributed => H2p is distributed
                do row=rstart,rstop
                   
                   !Counters for electron and hole row indices
                   h_row_cnt=((row-1)/Ne)+1
                   e_row_cnt=row-((h_row_cnt-1)*Ne)
                   
                   !Loop over lower triangle only
                   do col=1,row
                    
                      !Iterate H2p packed index counter
                      Hcnt=Hcnt+1
                      !Counters for electron and hole col indices
                      h_col_cnt=((col-1)/Ne)+1
                      e_col_cnt=col-((h_col_cnt-1)*Ne)
                      
                      !Same set of coefficients in both kernels, just multiplied
                      !by different orbital integrals, hence assign here
                      c1= c_vi(h_row_cnt)*c_ci(e_row_cnt)*c_vj(h_col_cnt)*c_cj(e_col_cnt)
                      c2= c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vi(h_col_cnt)*c_cj(e_col_cnt)
                      c3= c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vj(h_col_cnt)*c_ci(e_col_cnt)

                      !Packed Two-particle Hamiltonian 
                      !(technically just the kernel elements) K_(vc)(v'c')
                      H2p(Hcnt) = H2p(Hcnt) + &
                                  c1*(2._dp*Uc-ieps*Ux)+   &
                                  c2*(2._dp*Ux-ieps*Uc)+   &
                                  c3*(2._dp   -ieps   )*Ux                                      
 
                      !Second Coupling matrix K_(vc)(c'v')
                      K2(Hcnt)  = K2(Hcnt)  + &
                                  c1*(2._dp*Uc-ieps*Ux)   +  &
                                  c2*(2._dp   -ieps   )*Ux+  &
                                  c3*(2._dp*Ux-ieps*Uc) 

                      !Terms * Coulomb integral only
                      !H2p(Hcnt) = H2p(Hcnt) + c1*(2._dp*Uc) + c2*(-ieps*Uc) 
                      !K2(Hcnt)  = K2(Hcnt)  + c1*(2._dp*Uc) + c3*(-ieps*Uc) 
                      !Terms * Exchange integrals  only
                      !H2p(Hcnt) = H2p(Hcnt) + c1*(-ieps*Ux)+ c2*(2._dp*Ux)+ c3*(2._dp-ieps)*Ux           
                      !K2(Hcnt)  = K2(Hcnt)  + c1*(-ieps*Ux)+ c2*(2._dp-ieps)*Ux+c3*(2._dp*Ux) 


                   enddo
                enddo
                !Shut loops over (vc)(v'c')

             enddo
          enddo
       enddo
    enddo


    !Free Memory
    deallocate(U_c,U_x)

    !2nd NN distance - tolerance. Shouldn't need it but just to be safe
    d_2N=al/sqrt(2._dp)-0.05*bond_length(ic)

    !Continue summing terms in H2p
    !Beyond NN, only retain terms multiplied by Uc=(ii|v|jj)
    do ia=1,Natoms
       do ja=1,Natoms
          !|r-r'|
          denom=sqrt(dot_product(r(ia,:)-r(ja,:),r(ia,:)-r(ja,:)))
          !If atoms ia and ja are NN or same atom, cycle
          if(denom<d_2N) cycle 
          !Inter-atomic integral (ii|v|jj) =~ 1/|r_ia-r_ja|
          Uc=unt/denom
    
          do iorb=1,orb_bas
             i=iorb+(ia-1)*orb_bas
             c_vi(:) = evcr(i, VBT:VBT-Nh+1:-1) 
             c_ci(:) = evcr(i, VBT+1:VBT+Ne)    

             do jorb=1,orb_bas
                j=jorb+(ja-1)*orb_bas
                c_vj(:) = evcr(j, VBT:VBT-Nh+1:-1) 
                c_cj(:) = evcr(j, VBT+1:VBT+Ne)    

                Hcnt=0
                do row=rstart,rstop
                   h_row_cnt=((row-1)/Ne)+1
                   e_row_cnt=row-((h_row_cnt-1)*Ne)
                   
                   do col=1,row                    
                      Hcnt=Hcnt+1
                      h_col_cnt=((col-1)/Ne)+1
                      e_col_cnt=col-((h_col_cnt-1)*Ne)
                      
                      !Same set of coefficients in both kernels, just multiplied
                      !by different orbital integrals, hence assign here
                      c1= c_vi(h_row_cnt)*c_ci(e_row_cnt)*c_vj(h_col_cnt)*c_cj(e_col_cnt)
                      c2= c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vi(h_col_cnt)*c_cj(e_col_cnt)
                      c3= c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vj(h_col_cnt)*c_ci(e_col_cnt)

                      !Only retain terms in the kernel multiplied by Uc 

                      !Packed Two-particle Hamiltonian 
                      !(technically just the kernel elements) K_(vc)(v'c')
                      H2p(Hcnt) = H2p(Hcnt) +  (2._dp*c1-ieps*c2)*Uc
      
                      !Second Coupling matrix K_(vc)(c'v')
                      K2(Hcnt)  = K2(Hcnt)  +  (2._dp*c1-ieps*c3)*Uc

                   enddo
                enddo
                
             enddo
          enddo
       enddo
    enddo

   
    !Output routine runtime
    call cpu_time(t2_local)
    if(IO)write(*,'(A,F20.10,A)') 'Time taken to construct H2p = ',(t2_local-t1_local)/60._dp,' mins'


  end subroutine construct_orthogonal_H2p_K2


  !****************************************************************************
  !Read in and store on-site and NN Coulomb and exchange energies for use
  !in Bethe-Salpeter calculations. All integrals computed in the two-centre
  !approximation using a real basis of single-zeta functions from Siesta
  !****************************************************************************
  subroutine orbital_integral_energies(U_c,U_x)
    use var,     only: ns,orb_bas
    implicit none

    !********************
    !Declarations
    !********************

    !Modular Arrays
    real(dp), allocatable :: U_c(:,:),U_x(:,:)

    !Local Variables
    integer   ::  M,Nt,i,j,k,ierr
    real(dp)  ::  val
    character ::  fname1*50,fname2*50


    !********************
    !Main Routine
    !********************

    !One doesn't need to know grid details but can generate them if required:
    !Know the grid integers are (5,5,5), so use that to generate the number 
    !of points and grid spacing 
    !call gen_NN_rspacing(npts_r,rspacing,grid,iptr,atm_pntr,xyz_ptr,atom_species_NN)
    !deallocate(grid,iptr,atm_pntr,xyz_ptr,atom_species_NN)

    !NN array dimensions 
    M=ns*orb_bas

    !Number of elements in upper triangle of matrix
    Nt=int(0.5*M*(M+1))

    !Allocate arrays. Both are symmetric but easier (and faster?)
    !implementation to fill the whole thing
    allocate(U_c(M,M),U_x(M,M))
    U_c(:,:)=0._dp
    U_x(:,:)=0._dp


    !-----------------------
    !Read from file
    !-----------------------
    if(IO)then
       !Tabulated Coulomb and exchange integrals (energy in eV)
       fname1='parameter_files/orbital_integrals/coul_int.dat'
       fname2='parameter_files/orbital_integrals/exc_int.dat'

       !Open file
       open(unit=100,file=trim(adjustl(fname1)),iostat=ierr)
       if(ierr/=0 .and. IO) write(*,*) trim(adjustl(fname1)),' can''t be found'

       !Skip two lines of header
       read(100,*)
       read(100,*) 
     
       !Read in data 
       do k=1,Nt
          read(100,*) i,j,val
          U_c(i,j)=val
          U_c(j,i)=val
       enddo
       close(100)

       !Exchange integral
       open(unit=101,file=trim(adjustl(fname2)),iostat=ierr)
       if(ierr/=0 .and. IO) write(*,*) trim(adjustl(fname2)),' can''t be found'
       
       read(101,*) 
       read(101,*)  
       
       do k=1,Nt
          read(101,*) i,j,val
          U_x(i,j)=val
          U_x(j,i)=val
       enddo
       close(101)

    endif

    !Send arrays to other CPUs from master
    call mpi_barrier(MPI_COMM_WORLD,MPIierr)
    call mpi_bcast(U_c,M*M,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIierr)
    call mpi_bcast(U_x,M*M,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIierr)


    !Return U_c and U_x
  end subroutine orbital_integral_energies




  !******************************************************
  !Construct resonant Hamiltonian and store in 2D array
  !******************************************************
  subroutine construct_orthogonal_H2p_2D(H2p,rstart,rstop)
    use var, only: Natoms,NumNN,orb_bas,r,Ne,Nh,al,evcr,NN_list,ic,bond_length
    implicit none

    !******************
    !Declarations
    !******************
    !Global Variables
    integer,      intent(in)          :: rstart,rstop

    !Global Arrays
    real(dp),     allocatable,        intent(inout) :: H2p(:,:)

    !Local Variables
    integer  ::    N2p,N,VBT,ia,ja,ia_NN,ia_l,ja_l,iorb,jorb,i,j,i_l,j_l,&
                   Hcnt,row,col,h_row_cnt,h_col_cnt,e_row_cnt,e_col_cnt, cnt
    real(dp) ::    t1_local,t2_local,Uc,Ux,d_2N,denom

 
    !******************
    !Main Routine
    !******************
    call cpu_time(t1_local)
    if(IO)then
       write(*,*) 'Using Tamm-Dancoff Approximation'
       write(*,*) 'Constructing H2p in an orthogonal basis'
    endif
 

    !Resonant two-particle Hamiltonian dimensions
    N2p=Ne*Nh
  
    !****************************************************************
    !Integrals between on-site and NN terms explicitly computed 
    !Monopole term of multipole approximation is used for Coulomb-like 
    !integrals such that they reduce to 1/|r_ia - r_ja|
    !Exchange-like terms ignored beyond NN
    !****************************************************************

    !Sum on-site and NN terms 
    do ia=1,Natoms
       !Local atom index 
       ia_l=atom_species(ia)

       !Sum over NN of atom ia
       do ia_NN=1,NumNN
          ja=NN_list(ia,ia_NN)
          if(ja==0) cycle
          ja_l=atom_species(ja)
     

          do iorb=1,orb_bas
             !Orbital index for a system comprised of 1 anion and 1 cation
             i_l=iorb+(ia_l-1)*orb_bas
             !Global orbital index
             i=iorb+(ia-1)*orb_bas
             !Assign TB coefficients to 1D arrays 
             c_vi(:) = evcr(i, VBT:VBT-Nh+1:-1) 
             c_ci(:) = evcr(i, VBT+1:VBT+Ne)    

             do jorb=1,orb_bas
                j_l=jorb+(ja_l-1)*orb_bas
                j=jorb+(ja-1)*orb_bas
                c_vj(:) = evcr(j, VBT:VBT-Nh+1:-1) 
                c_cj(:) = evcr(j, VBT+1:VBT+Ne)    
                
                !Integrals between orbitals and bare Coulomb potential
                Uc=U_c(i_l,j_l)
                !Only on-site Coulomb integral explicitly calculated
                Ux=U_x(i_l,j_l)

                !Row loop not distributed
                do row=1,N2p
                   !Counters for electron and hole row indices
                   h_row_cnt=((row-1)/Ne)+1
                   e_row_cnt=row-((h_row_cnt-1)*Ne)
                   
                   !Loop over lower triangle only
                   do col=1,row
                      !Counters for electron and hole col indices
                      h_col_cnt=((col-1)/Ne)+1
                      e_col_cnt=col-((h_col_cnt-1)*Ne)
                      
                      !Packed Two-particle Hamiltonian (technically just
                      !the kernel elements)
                      H2p(row,col) = H2p(row,col) + &
		    
                                  c_vi(h_row_cnt)*c_ci(e_row_cnt)*c_vj(h_col_cnt)* &
                                  c_cj(e_col_cnt)*(2._dp*Uc-ieps*Ux)+              &

                                  c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vi(h_col_cnt)* &
                                  c_cj(e_col_cnt)*(2._dp*Ux-ieps*Uc)+              &

                                  c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vj(h_col_cnt)* &
                                  c_ci(e_col_cnt)*(2._dp-ieps)*Ux                                         
           
                   enddo
                enddo
                !Shut loops over (vc)(v'c')

             enddo
          enddo
       enddo
    enddo


    !Free Memory
    deallocate(U_c,U_x)

    !2nd NN distance - tolerance. Shouldn't need it but just to be safe
    d_2N=al/sqrt(2._dp)-0.05*bond_length(ic)

    !Continue summing terms in H2p
    !Beyond NN, only retain terms multiplied by Uc=(ii|v|jj)
    do ia=1,Natoms
       do ja=1,Natoms
          !|r-r'|
          denom=sqrt(dot_product(r(ia,:)-r(ja,:),r(ia,:)-r(ja,:)))
          !If atoms ia and ja are NN or same atom, cycle
          if(denom<d_2N) cycle 
          !Inter-atomic integral (ii|v|jj) =~ 1/|r_ia-r_ja|
          Uc=unt/denom
    
          do iorb=1,orb_bas
             i=iorb+(ia-1)*orb_bas
             c_vi(:) = evcr(i, VBT:VBT-Nh+1:-1) 
             c_ci(:) = evcr(i, VBT+1:VBT+Ne)    

             do jorb=1,orb_bas
                j=jorb+(ja-1)*orb_bas
                c_vj(:) = evcr(j, VBT:VBT-Nh+1:-1) 
                c_cj(:) = evcr(j, VBT+1:VBT+Ne)    

                do row=1,N2p
                   h_row_cnt=((row-1)/Ne)+1
                   e_row_cnt=row-((h_row_cnt-1)*Ne)
                   
                   do col=1,row                    
                      h_col_cnt=((col-1)/Ne)+1
                      e_col_cnt=col-((h_col_cnt-1)*Ne)
                      
                      !Two-particle Hamiltonian (technically just
                      !the kernel elements) 
                      H2p(row,col) = H2p(row,col) +  &
                                  ( 2._dp*c_vi(h_row_cnt)*c_ci(e_row_cnt)*c_vj(h_col_cnt)*c_cj(e_col_cnt)  &
                                    -ieps*c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vi(h_col_cnt)*c_cj(e_col_cnt) )*Uc

                   enddo
                enddo
                
             enddo
          enddo
       enddo
    enddo

   
    !Output routine runtime
    call cpu_time(t2_local)
    if(IO)write(*,'(A,F20.10,A)') 'Time taken to construct H2p = ',(t2_local-t1_local)/60._dp,' mins'


  end subroutine construct_orthogonal_H2p_2D



  !*******************************************
  !Construct exciton wave function 
  !*******************************************
  subroutine exciton_wavefunction_TDA()
    implicit none
  end subroutine exciton_wavefunction_TDA






end module kernel
  
