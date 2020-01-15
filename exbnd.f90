module exbnd

  !************************************************************
  !Module: exbnd 'exciton binding energy'
  !Routines to compute the Coulomb and exchange integrals which
  !will give the exciton binding energy.
  !************************************************************  
  use var,   only: dp
  use parallel
  use sort,  only: reallocate
  use setup, only: snpts,sspacing,sgrid, satm_pntr, siptr, sxyz_ptr, sys_pntr, gsindex,&
                   setup_system_grid, setup_arb_system_grid,setup_BSEtest_system_grid
  implicit none


  !*******************************
  !Modular Types and Declarations
  !*******************************

  !NN grid arrays
  integer,     allocatable ::  NNnpts_r(:),NNiptr(:,:),NNatm_pntr(:),NNxyz_ptr(:,:,:),&
                               atom_species_NN(:)
  real(dp),    allocatable ::  NNrspacing(:,:),NNgrid(:,:)




contains
  !************************************************************
  !Subroutines and functions
  !************************************************************

  

  !****************************************************************************
  !Precompute Coulomb integrals between atomic orbitals in the NN approximation
  !****************************************************************************

  subroutine coulomb_integral_NN_a1()
    use var,   only: dp,qp,Natoms,orb_bas,evc,atom_species,d,Nexsp,eigen_cnt,NN_list,NumNN,&
                     r,time_CoulInt,al,NNdp_list,atom_type
    use grids, only: gen_NN_grid, gen_NN_rspacing
    use smat,  only: NN_disp
    use sort,  only: reallocate
    implicit none


    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer       N,Nm,iro,ir,jro,jr,  ia,ja,ka,la,  n1,n2,n3,n4,n12,n34 , iorb,jorb,korb,lorb,&
                  n1_l,n2_l,n3_l,n4_l,Natm,ia_l,ja_l,ka_l,la_l,n1234,NOI      
    real(dp)      dV,integral

    !Modular arrays
    integer,      allocatable         :: orb_grd_pntr(:,:,:),nnz(:,:)
    real(dp),     allocatable         :: orb_overlap(:,:,:)
    real(dp),     allocatable         :: worb(:)

    !Local Arrays
    real(dp),     allocatable         :: r1(:),r2(:),r_local(:,:)



    !***************************
    !Main Routine
    !***************************

    !Number of atoms in NN system
    Natm=NumNN
    !Atomic positions in NN system
    allocate(r_local(d,Natm))
    !r_local in units of Angstroms 
    r_local = 0.25*al*reshape( (/ 0.,  0.,  0., &   
                                  1.,  1.,  1., &
                                  1., -1., -1., &
                                 -1.,  1., -1., &
                                 -1., -1.,  1./), (/d,Natm/))




    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)
    
    !All arguments here are defined in the module
    !call gen_NN_grid(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)
    !This NN grid uses the same spacing as the system grid
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)

    !Volume element associated with this grid spacing
    dV=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Construct orb_overlap matrix on NN grid
    call construct_orb_overlap(orb_overlap,orb_grd_pntr,nnz)

    !Dimensions of the system
    N=Natoms*orb_bas

    !Maximum number of orbital interactions for two centres, in NN regime
    NOI=N*(NumNN*orb_bas)

    !Orbital integral, allocated with upper bounds
!!$    allocate( worb(ceiling(0.5*(NOI**2))) )
!!$    worb(:)=0._dp

    !Allocate spatial work arrays, r(r1) and r'(r2)
    allocate(r1(d),r2(d))


    !***************************
    !Loop over n12 and n34>n12
    !***************************
    !Open File
    open(unit=001,file='w1.dat')

    !Initialise 8-fold composite index
    n1234=0

    !First composite index over orbitals 1&2 at defined at spatial coordinate r (ir)
    do n12=1,NOI


       !Return indices warpped in composite index n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)
       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle


       !Global n1 and n2 values
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)


       !Second Composite index over orbitals 3&4 at defined at spatial coordinate r' (jr)
       do n34=n12,NOI

          !Return indices warpped in composite index n34
          call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
          !If ka doesn't have a neighbour la, cycle loop
          if(la==0)cycle


          !Global n3 and n4 values
          n3=korb+((ka-1)*orb_bas)
          n4=lorb+((la-1)*orb_bas)


          !Iterate 8-fold coposite index
          n1234=n1234+1

!!$          !Cycle if either ia or ka are cations. 
!!$          if(atom_type(ia)==1 .or. atom_type(ka)==1) cycle
          
          
          !*******************************************************
          !Replace integral at r=r' with analytic approx
          !*******************************************************
          !Integrate over points which overlap only
          !iro=Composite spatial index for overlapping points
          r1(:)=0._dp
          r2(:)=0._dp
          integral=0._dp
    
          do iro=1,nnz(n1_l,n2_l)
             !Composite spatial index for NN grid
             ir=orb_grd_pntr(n1_l,n2_l,iro)
             !Spatial variable 1 (r), centred at atom ia 
             r1(:)= NNgrid(:,ir)-r_local(:,ia_l)+r(ia,:)
             
             do jro=1,nnz(n3_l,n4_l)
                !Composite spatial index for NN grid
                jr=orb_grd_pntr(n3_l,n4_l,jro)
                !Spatial variable 2 (r'), centred at atom ka
                r2(:)= NNgrid(:,jr)-r_local(:,ka_l)+r(ka,:)
                
                !Cycle if r=r'
                if(sum(r1-r2)==0._dp)cycle
                !if(sqrt(dot_product(r1-r2,r1-r2))==0._dp)cycle

                !Integrand 
                integral=integral + (orb_overlap(n1_l,n2_l,iro)*orb_overlap(n3_l,n4_l,jro)/&
                                    sqrt(dot_product(r1-r2,r1-r2)) )

             enddo
          enddo
          !*******************************************************


          !Output integral 
           write(001,'(8(I4,X),F20.10)') ia,iorb,ja,jorb, ka,korb,la,lorb, integral*dV*dV

          !Store orbital Coulomb integral
!!$          worb(n1234) = integral*dV*dV

       enddo
    enddo

    close(001)

    !Reallocate worb
!!$    call reallocate(worb,n1234)

    !Return worb(n1234), noting that only the upper triangle is filled. 
  end subroutine coulomb_integral_NN_a1





  !****************************************************************************
  !Unwrap the individual atomic and orbital indices from the composite orbital 
  !index, nij
  !Reads more clearly in a subroutine than in the calling routine
  !****************************************************************************

  subroutine indices(n12, n1_local,n2_local,ia_local,ja_local, ia,ja,iorb,jorb)
    use var, only:NN_list,NumNN,orb_bas,NNdp_list,atom_type
    implicit none
    
    integer, intent(in)  :: n12
    integer, intent(out) :: ia,ja,iorb,jorb,n1_local,n2_local,ia_local,ja_local
    integer              :: n2NN,iNN,n1,n2

    !Orbital Index 1
    n1=((n12-1)/(NumNN*orb_bas))+1
    !NN of Orbital index 1 (slightly erroneously called n2NN)
    n2NN=n12-((n1-1)*NumNN*orb_bas)

    !Atomic index 1
    ia=((n1-1)/orb_bas)+1
    !Orbital index 1 
    iorb=n1-((ia-1)*orb_bas)  
       
    !NN of atomic index 1 
    iNN=((n2NN-1)/orb_bas)+1
       
    !Atomic index 2
    ja=NN_list(ia,iNN)
    !If ia doesn't have neighbour at iNN, cycle
    !Moved to calling routine
    !Orbital index 2 (must be retrieved this way due to loop structure)
    jorb=mod(n2NN-1,orb_bas)+1
    !Use jorb and ja to retrieve n2
    !n2=jorb+((ja-1)*orb_bas)


    !Using the retrieved indices and the atom_species, 
    !relate (ia,iorb) and (ja,jorb) to orbital indices consistent with
    !the NN loop structure used to construct orb_overlap
    !See ham.f90 construction of S for exactly the same concept

!!$    !Anion
!!$    if(atom_type(ia)==-1)then
!!$       ia_local=1
!!$       ja_local=NNdp_list(ia,iNN)
!!$    
!!$    !Cation (atom_type(ia)==1)
!!$    else
!!$       !Reverse atomic indices for cation-centred
!!$       ia_local=NNdp_list(ia,iNN)
!!$       ja_local=1
!!$    endif


    !Anion-centred
    if(atom_type(ia)==-1)then
       ia_local=1
       ja_local=NNdp_list(ia,iNN)
    !Cation-centred. Interatomic 
    elseif(atom_type(ia)==1 .and. ia/=ja)then
       ia_local=NNdp_list(ia,iNN)
       ja_local=1
    !Cation-centred. On-site. This is the extra statement w.r.t. above
    elseif(atom_type(ia)==1 .and. ia==ja)then
       ia_local=2
       ja_local=2
    endif

  
    !Local orbital index 1
    n1_local=iorb+((ia_local-1)*orb_bas)
    !Local orbital index 2
    n2_local=jorb+((ja_local-1)*orb_bas)

  end subroutine indices


  !---------------------------------------------------------
  !Unwrap global orbital index n1=(ia,iorb) where ia
  !is an atomic index and iorb is a local orbital index
  !Also map global orbital n1 to NN system which contains
  !just two atoms, therefore 10 orbitals
  !---------------------------------------------------------
  subroutine index(n1, n1_local,ia_local, ia,iorb)
    use var, only:NN_list,NumNN,orb_bas,NNdp_list,atom_type
    implicit none
    
    integer, intent(in)  :: n1
    integer, intent(out) :: ia,iorb,n1_local,ia_local
    
    !Atomic index 
    ia=((n1-1)/orb_bas)+1
    !Orbital index 1 
    iorb=n1-((ia-1)*orb_bas)  
    
    !Anion on-site
    if(atom_type(ia)==-1)then
       ia_local=1
       !Cation-centred. On-site. 
    elseif(atom_type(ia)==1)then
       ia_local=2
    endif
    
    !Local orbital index 1
    n1_local=iorb+((ia_local-1)*orb_bas)
    
  end subroutine index




  !****************************************************************************
  !Testing routine for construct_orb_overlap
  !Call from main
  !****************************************************************************

  subroutine test_orb_overlap()
    use var,     only: NumNN,orb_bas,Natoms,NNdp_list,atom_type,NN_list
    use grids,   only: gen_NN_rspacing
    use smat,    only: NN_disp,NN_smatrix
    implicit none


    !*******************************
    !Declarations
    !*******************************

    !Local Variables
    integer        n12,N,n1_l,n2_l,ia,ja,iorb,jorb,n1,n2,iro, i,l,jNN,m,row_cnt,col_cnt,j,&
                   ja_local,ia_local,i_l,j_l,ia_l,ja_l,ir
    real(16)       integral,dV

    !Local arrays
    real(dp),     allocatable  :: S(:,:),S2(:,:)
    complex(dp),  allocatable  :: S_NN(:,:)
    
    !Modular arrays
    integer,      allocatable  :: orb_grd_pntr(:,:,:),nnz(:,:)
    real(dp),     allocatable  :: orb_overlap(:,:,:)


    !*******************************
    !Main Routine
    !*******************************


    !Set up NN grid. Natoms_grid=NumNN
    !call gen_NN_grid(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)




    !***************************************************************************
    !Generate Overlap matrix using standard approach but utilising this NN grid
    !****************************************************************************

    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)
    
    !Construct NN S matrix (allocated in smat.f90)
    !Modified the routine to use the same NN grid as 'orb_overlap'
    call NN_smatrix(S_NN)

    !Allocate S
    N=Natoms*orb_bas
    allocate(S(N,N))
    S(N,N)=0._dp


    do i=1,Natoms
       do l=1,orb_bas
          do jNN=1,NumNN

             !Skip if NN doesn't exist
             j=NN_list(i,jNN)
             if(j==0)cycle
             
             do m=1,orb_bas

                !Shouldn't need local-global mapping for super cells
                !nor element checker for finite system

                !Overlap counters for S:
                row_cnt= l + ((i-1)*(orb_bas))
                col_cnt= m + ((j-1)*(orb_bas))

                !Map atomic index j to local NN grid index
                !Required for the indices used in S_NN
                ja_local=NNdp_list(i,jNN)

                
                !i will always relate to the central atom in S_NN, 
                ia_local=1
                      
                !Orbital indices for S_NN
                if(atom_type(i)==-1)then
                   i_l= l + ((ia_local-1)*(orb_bas))
                   j_l= m + ((ja_local-1)*(orb_bas))
                   !If cation, NNdp_list will give local atomic indices relative to cation centre
                   !Interchange them to make them correct in context to S_NN.  
                elseif(atom_type(i)==1 .and. i/=j)then
                   i_l= l + ((ja_local-1)*(orb_bas))
                   j_l= m + ((ia_local-1)*(orb_bas))
                !On-site cation orbitals. Required for debugging, when explicitly
                !comparing computed integrals instead of forcing terms to 1/0 
                !via orthonormality
                elseif(atom_type(i)==1 .and. i==j)then
                   i_l= l + orb_bas
                   j_l= m + orb_bas
                endif

                
                !Assign overlap matrix element
                S(row_cnt,col_cnt)=real(S_NN(i_l,j_l))
                !write(100,'(4(I3,X),F12.6,X)') i,j,l,m,S(row_cnt,col_cnt)
                
             enddo
          enddo
       enddo
    enddo

    !Output all zero elements
!!$    do ia=1,Natoms
!!$       do ja=1,Natoms
!!$          do iorb=1,orb_bas
!!$             do jorb=1,orb_bas
!!$
!!$                i=iorb+((ia-1)*orb_bas)
!!$                j=jorb+((ja-1)*orb_bas)
!!$
!!$                if( S(i,j)==0._dp)then
!!$                   write(100,*) ia,ja, S(i,j)
!!$                endif
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
  



    !*****************************************************************
    !Approach 2. 
    !*****************************************************************

    !NN grid that uses the spacing of the system grid 
    write(*,*) 'Generating NN grid with spacing of system grid'
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)

    !Volume element associated with this grid spacing
    dV=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Construct orb_overlap matrix on NN grid
    write(*,*)'Constructing orbital product array'
    call construct_orb_overlap(orb_overlap,orb_grd_pntr,nnz)

    !Overlap Matrix constructed from this routine
    allocate(S2(N,N))
    S2(N,N)=0._dp

    !Open file
    open(unit=002, file='s2.dat')
    open(unit=001, file='s.dat')    

    !****************************************************************************
    !Compute integral equivalent to the overlap integral for (n1,n2)
    !****************************************************************************

    !4-fold Composite index. n12 = (n1,n2) = (ia,iorb,ja,jorb)
    do n12=1,N*NumNN*orb_bas

       !Return indices warpped in composite index n12
       !Local n1 and n2 indices. Global (ia,iorb, ja,jorb) indices
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)

       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle

       !Global n1 and n2 values
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)


       !*******************************
       !Integral
       !*******************************
       integral=0._dp
       !Loop over Non-zero values
       do iro=1,nnz(n1_l,n2_l)
          ir=orb_grd_pntr(n1_l,n2_l,iro)
          !Sum integral
          integral=integral+orb_overlap(n1_l,n2_l,ir)
!!$          integral=integral+orb_overlap(n1_l,n2_l,iro)
       enddo
       !Mulitply by volume element
       S2(n1,n2)=integral*dV
       !*******************************
       
       !Write both overlap matrices to different files
!!$       if(integral/=0._dp)then
          write(001,*) n1,n2, S(n1,n2)
          write(002,*) n1,n2, S2(n1,n2)
!!$       endif
    enddo

    !Close files
    close(001)
    close(002)

    write(*,*) 'written out non-zero integrals to machine precision'

  end subroutine test_orb_overlap




  !*****************************************************************
  !Check integral{ n12(r) n34(r') dr dr'  against the product of 
  !overlap matrices S(n1,n2)*S(n3,n4)
  !*****************************************************************

  subroutine test2_orb_overlap()
    use var,     only: NumNN,orb_bas,Natoms,NNdp_list,atom_type,NN_list
    use grids,   only: gen_NN_rspacing
    use smat,    only: NN_disp,NN_smatrix
    implicit none


    !*******************************
    !Declarations
    !*******************************

    !Local Variables
    integer        n12,N,n1_l,n2_l,ia,ja,iorb,jorb,n1,n2,iro, i,l,jNN,m,row_cnt,col_cnt,j,&
                   ja_local,ia_local,i_l,j_l,n34,n3_l,n4_l,ka,la,n3,n4,korb,lorb,jro,ka_l,la_l,&
                   ia_l,ja_l,jr_l,ir_l
    real(16)       integral,dV

    !Local arrays
    real(dp),     allocatable  :: S(:,:),S2(:,:)
    complex(dp),  allocatable  :: S_NN(:,:)
    
    !Modular arrays
    integer,      allocatable  :: orb_grd_pntr(:,:,:),nnz(:,:)
    real(dp),     allocatable  :: orb_overlap(:,:,:)



    !*******************************
    !Main Routine
    !*******************************



    !***************************************************************************
    !Generate overlap matrix S(n1,n2)
    !****************************************************************************

    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)
    
    !Construct NN S matrix (allocated in smat.f90)
    !Modified the routine to use the same NN grid as 'orb_overlap'
    call NN_smatrix(S_NN)

    !Allocate S
    N=Natoms*orb_bas
    allocate(S(N,N))
    S(:,:)=0._dp


    do i=1,Natoms
       do l=1,orb_bas
          do jNN=1,NumNN

             !Skip if NN doesn't exist
             j=NN_list(i,jNN)
             if(j==0)cycle
             
             do m=1,orb_bas

                !Shouldn't need local-global mapping for super cells
                !nor element checker for finite system

                !Overlap counters for S:
                row_cnt= l + ((i-1)*(orb_bas))
                col_cnt= m + ((j-1)*(orb_bas))

       
                !Map atomic index j to local NN grid index
                !Required for the indices used in S_NN
                ja_local=NNdp_list(i,jNN)

                
                !i will always relate to the central atom in S_NN, 
                ia_local=1
                      
                !Orbital indices for S_NN
                if(atom_type(i)==-1)then
                   i_l= l + ((ia_local-1)*(orb_bas))
                   j_l= m + ((ja_local-1)*(orb_bas))
                   !If cation, NNdp_list will give local atomic indices relative to cation centre
                   !Interchange them to make them correct in context to S_NN.  
                elseif(atom_type(i)==1 .and. i/=j)then
                   i_l= l + ((ja_local-1)*(orb_bas))
                   j_l= m + ((ia_local-1)*(orb_bas))
                !On-site cation orbitals. Required for debugging, when explicitly
                !comparing computed integrals instead of forcing terms to 1/0 
                !via orthonormality
                elseif(atom_type(i)==1 .and. i==j)then
                   i_l= l + orb_bas
                   j_l= m + orb_bas
                endif
                
                !Assign overlap matrix element
                S(row_cnt,col_cnt)=real(S_NN(i_l,j_l))
                !write(100,'(4(I3,X),F12.6,X)') i,j,l,m,S(row_cnt,col_cnt)
                
             enddo
          enddo
       enddo
    enddo




    !*****************************************************************
    !Check integral{ n12(r) n34(r') dr dr'
    !*****************************************************************

    !NN grid that uses the spacing of the system grid 
    write(*,*) 'Generating NN grid with spacing of system grid'
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)

    !Volume element associated with this grid spacing
    dV=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Construct orb_overlap matrix on NN grid
    write(*,*)'Constructing orbital product array'
    call construct_orb_overlap(orb_overlap,orb_grd_pntr,nnz)

    !Overlap Matrix constructed from this routine
    allocate(S2(N*NumNN*orb_bas,N*NumNN*orb_bas))
    S2(:,:)=0._dp

    !Open files
    open(unit=001, file='s.dat')    
    open(unit=002, file='s2.dat')



    !****************************************************************************
    !Compute integral{ n12(r) n34(r') } dr dr'
    !equivalent to S(n1,n2)*S(n3,n4)
    !Only looping over upper triangular 
    !****************************************************************************

    !4-fold Composite index. n12 = (n1,n2) = (ia,iorb,ja,jorb)
    do n12=1,N*NumNN*orb_bas

       !Return indices warpped in composite index n12
       !Local n1 and n2 indices. Global (ia,iorb, ja,jorb) indices
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)

       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle

       !Global n1 and n2 values
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)


       !2nd Composite Loop
       do n34=n12,N*NumNN*orb_bas

          !Indices
          call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
          if(la==0)cycle
          n3=korb+((ka-1)*orb_bas)
          n4=lorb+((la-1)*orb_bas)


          !*******************************
          !Integral over r,r'
          !*******************************
          integral=0._dp

          do iro=1,nnz(n1_l,n2_l)
             do jro=1,nnz(n3_l,n4_l)
                jr_l=orb_grd_pntr(n3_l,n4_l,jro)
                ir_l=orb_grd_pntr(n1_l,n2_l,iro)
                integral=integral+( orb_overlap(n1_l,n2_l,ir_l)*orb_overlap(n3_l,n4_l,jr_l))
             enddo
          enddo

          !Mulitply by volume elements
          S2(n12,n34)=integral*dV*dV
          !*******************************
          

          !Write both overlap matrices to different files
!!$          if(integral/=0._dp)then
             write(001,*) n12,n34, S(n1,n2)*S(n3,n4)
             write(002,*) n12,n34, S2(n12,n34)
!!$          endif

       enddo
    enddo
    !Shut loop over n12


    !Close files
    close(001)
    close(002)
    
    write(*,*)'Output non-zero values only, to improve comparison'

  end subroutine test2_orb_overlap





  !****************************************************************************
  !Construct orb_overlap matrix on NN grid
  !****************************************************************************
  !Return array containing all permutations of overlaps between two NN orbitals
  !Return an index array 'orb_grd_pntr' which indexs points of overlap only
  !****************************************************************************
  
  subroutine construct_orb_overlap(orb_overlap,orb_grd_pntr,nnz)
    
    !Global variables and arrays
    use var,     only: NumNN,orb_bas
    !Orbital functions/profiles
    use orbital, only: radial_function,orbital_profile 
    !Place orbitals on grid
    use smat,    only: place_orbitals 

    implicit none


    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer                 Natm,Ntotal,Nsub    
    integer                 ia,iorb,ja,jorb,n1,n2,ir,iro,ix,iy,iz,i,j
    integer                 x0,y0,z0,xf,yf,zf       !Unused

    !Global Arrays    
    integer,                allocatable :: midpnt(:)
    type(radial_function),  allocatable :: &
                            sorb(:,:,:),pxorb(:,:,:),pyorb(:,:,:),pzorb(:,:,:),sstorb(:,:,:)

    !Modular arrays
    integer,      allocatable, intent(out) :: orb_grd_pntr(:,:,:),nnz(:,:)
    real(dp),     allocatable, intent(out) :: orb_overlap(:,:,:)
    real(dp),     allocatable              :: orbdummy(:,:,:)
    
    !Local Arrays
    complex(dp),            allocatable :: orbital1(:,:,:),orbital2(:,:,:)
    


    !****************
    !Main Routine
    !****************
    !NN grid set up in calling routine and available to the module

    !Generate 3D orbital profiles
    call orbital_profile(NNrspacing,midpnt,sorb,pxorb,pyorb,pzorb,sstorb)

    !Number of atoms in the grid = 1 central atom + 4 NN
    Natm=NumNN

    !Total number of orbital states in NN subsystem
    Nsub=Natm*orb_bas

    !Total number of points in NN grid
    Ntotal=product(NNnpts_r)

    !Orbitals 1 and 2, filled per iteration. Dimensions of the  NN grid
    allocate(orbital1(NNnpts_r(1),NNnpts_r(2),NNnpts_r(3)),orbital2(NNnpts_r(1),NNnpts_r(2),NNnpts_r(3)))

    !orb_overlap has dimensions of the NN grid and contains the overlap between orbitals n1 and n2
    !3rd dimension allocated with upper bound 
    allocate(orb_overlap(Nsub,Nsub,Ntotal))

    !Spatial index, mapping the composite NN grid index to the composite overlap index  
    !Allocated with an upper bound
    allocate(orb_grd_pntr(Nsub,Nsub,Ntotal), nnz(Nsub,Nsub))

    orb_overlap(:,:,:)=0.
    orb_grd_pntr(:,:,:)=0
    nnz(:,:)=0

    !*****************************************************************
    !Loop over all atoms in a NN grid to create an array containing all 
    !orbital overlaps/products
    !Order is consistent with loops in orbital integral (important!)
    !*****************************************************************

    do ia=1,Natm
       do iorb=1,orb_bas

          !Composite state index 1
          n1=iorb+((ia-1)*orb_bas)

          do ja=1,Natm
             do jorb=1,orb_bas

                !Composite state index 2
                n2=jorb+((ja-1)*orb_bas)

                !Place appropriate orbitals on grid:orbital1 and orbital2 returned.
                call place_orbitals(orbital1,orbital2,NNatm_pntr,NNiptr,&
                     ia,ja,iorb,jorb,midpnt,sorb,pxorb,pyorb,pzorb,sstorb,&
                     x0,y0,z0,xf,yf,zf,atom_species_NN) 

         
                !Construct an array containing the overlap between two orbitals 
                !Loop over all points in the NNgrid. Can then easily define a composite
                !counter ir and map that to the overlap spatial index iro

                iro=0
                ir=0
  
                do ix=1,NNnpts_r(1)
                   do iy=1,NNnpts_r(2)
                      do iz=1,NNnpts_r(3)

                         !Composite index for NN grid
                         ir=ir+1

                         !Store for all space, not just where the product is non-zero
                         !orb_overlap(n1,n2,ir)=real(orbital1(ix,iy,iz))+real(orbital2(ix,iy,iz))

                         !If the product of orbital1 and orbital2 is non-zero
                         !store the it
                         if( real(orbital1(ix,iy,iz)*orbital2(ix,iy,iz))/=0._dp)then
                            !Composite index for spatial overlap
                            iro=iro+1
                            !Array contain overlap between two orbitals
                            !Summation for visualisation ONLY:
                            !orb_overlap(n1,n2,iro)=real(orbital1(ix,iy,iz))+real(orbital2(ix,iy,iz))
                            !Product for ALL computations
                            orb_overlap(n1,n2,iro)=real(orbital1(ix,iy,iz))*real(orbital2(ix,iy,iz))
                            !Map points of overlap from three spatial indices to one spatial index
                            orb_grd_pntr(n1,n2,iro)=ir
                         endif

                      enddo
                   enddo
                enddo

                !Store number of non-zero values per (n1,n2)
                nnz(n1,n2)=iro                            

             enddo
          enddo

       enddo
    enddo


    !Allocate correct dimension using the overlap
    !as stored in nnz(n1,n2)
    call reallocate(orb_grd_pntr,Nsub,Nsub,maxval(nnz))

    !Storing all spatial points for orb_overlap, hence don't reallocate
    call reallocate(orb_overlap,Nsub,Nsub,maxval(nnz))

    !Return now to avoid writing to cube file
    return


    !**********************************
    !Cube output
    !**********************************

    !If I want to visualise the full orbitals, not just the region of
    !overlap, need to replace this code segment (above)

    !if( real(orbital1(ix,iy,iz)*orbital2(ix,iy,iz))/=0._dp)then
       !Composite index for spatial overlap
    !   iro=iro+1
       !Array contain overlap between two orbitals
    !   orb_overlap(n1,n2,iro)=real(orbital1(ix,iy,iz))*real(orbital2(ix,iy,iz))
       !Map points of overlap from three spatial indices to one spatial index
    !   orb_grd_pntr(n1,n2,iro)=ir
    !endif
    
    !with this: commented-out if statement and sum the orbitals, not take the product
   
    ! !!$if( real(orbital1(ix,iy,iz)*orbital2(ix,iy,iz))/=0._dp)then
       !Composite index for spatial overlap
    !   iro=iro+1
       !Array contain overlap between two orbitals
    !   orb_overlap(n1,n2,iro)=real(orbital1(ix,iy,iz))+real(orbital2(ix,iy,iz))
       !Map points of overlap from three spatial indices to one spatial index
    !   orb_grd_pntr(n1,n2,iro)=ir
    !!!$endif
    
  
    allocate(orbdummy(NNnpts_r(1),NNnpts_r(2),NNnpts_r(3)))

    !Loop over orbital states n1 and n2
    do ia=1,1
       do iorb=1,orb_bas
          
          !Composite state index 1
          n1=iorb+((ia-1)*orb_bas)

          do ja=1,Natm
             do jorb=1,orb_bas

                !Composite state index 2
                n2=jorb+((ja-1)*orb_bas)

                !Loop over x,y,z dimensions to fill dummy array for (n1,n2)
                ir=0
                orbdummy(:,:,:)=0.
                do ix=1,NNnpts_r(1)
                   do iy=1,NNnpts_r(2)
                      do iz=1,NNnpts_r(3)
                         !Composite index
                         ir=ir+1
                         orbdummy(ix,iy,iz)=orb_overlap(n1,n2,ir)
                      enddo
                   enddo
                enddo

                !Output the overlap to a cube file
                call cube_output(orbdummy,n1,n2)

             enddo
          enddo
       enddo
    enddo

    !Return orb_overlap and orb_grd_pntr
  end subroutine construct_orb_overlap




  !*******************************************************************
  !Output NN overlaps to cube file
  !Requires dummy array as orb_overlap only has one spatial dimension
  !not 3.
  !*******************************************************************
  subroutine cube_output(orbdummy,n1,n2)
    use var,only:dp,d,NumNN ,al  
    implicit none


    !**********************
    !Declarations
    !**********************

    !Both are modular arrays 
    !integer,     allocatable, intent(in) :: NNnpts_r(:)
    !real(dp),    allocatable, intent(in) :: NNrspacing(:,:)

    !Modular variables
    integer,     intent(in) :: n1,n2

    !Modular arrays
    real(dp),    allocatable, intent(in) :: orbdummy(:,:,:)

    !Local arrays
    real(dp),    allocatable             :: r_local(:,:),r0(:)

    !Local variables
    integer             :: Natoms,atnum,iatm,rw_len,Nrowz,rmdr,rowz,jz,jmax,iz1,iz2
    integer             :: i,ix,iy,iz
    real(dp)            :: Bohr
    character           :: n1_lab*5, n2_lab*5



    !**************************************
    !Main Routine
    !**************************************

    !Convert Angstroms to Bohr
    Bohr = 1./0.5291772109217

    !Label the output with orbital and site index
    write(n1_lab,'(I2)')   n1
    write(n2_lab,'(I2)')   n2

    !Open file to write orbital iorb, at atom ia to. 
    open(unit=002,file='orb_'//trim(adjustl(n1_lab))//'_'//trim(adjustl(n2_lab))//'.cube')



    !**************************************
    !Hard-coded variables for NN grid
    !**************************************
    !Number of atoms
    Natoms=NumNN

    !Starting point for NN grid: x0,y0,z0. See grids.f90/gen_NN_grid
    !Converted to Bohr
    allocate(r0(d))
    r0 = Bohr*al*(/-1.,-1.,-1./)

    !Atomic positions 
    allocate(r_local(Natoms,d))
    r_local = 0.25*al*transpose( reshape( (/ 0.,  0.,  0., &   
                                             1.,  1.,  1., &
                                             1., -1., -1., &
                                            -1.,  1., -1., &
                                            -1., -1.,  1./), (/d,Natoms/)))
 
    !Convert to Bohr
    r_local(:,:) = Bohr*r_local(:,:)


    !**************************************
    !Beginning of output to cube file
    !**************************************

    !Lines 1 & 2 are optional comments
    write(002,'(A,X,I5,X,A,X,I5)') 'Orbital cube file for overlap:',n1,' ',n2
    !Second line is blank
    write(002,*)

    !Line 3. Natoms with minus sign infront: Seems to be important if outputted data is in Bohr
    !Origin of the voxel space == grid. Cube file regenerates it from its origin,
    !starting vector, number of points and the grid (==lattice) vectors.
    write(002,'(I5,3(F12.6))') -1*Natoms, r0(:)


    !Lines 4,5,6: Number of data points along the lattice vectors r1, r2 and r3. 
    !as defined by rspacing(i,:).
    do i=1,d
       write(002,'(I5,3(F12.6))')  NNnpts_r(i), NNrspacing(i,:)*Bohr
    enddo

    !Next block of lines give the atomic number, dble(atomic number) and the atomic positions 
    !Atomic numbers are dummy values as they're not important
    do iatm=1,Natoms
       if(iatm==1) atnum=1
       if(iatm/=1) atnum=2
       write(002,'(I5,4(F12.6))') atnum, dble(atnum), r_local(iatm,:)
    enddo

    !Number of (molecular) orbitals evaluated per file and the index for the orbital
    write(002,'(3I5)') 1,n1,n2


    !Remaining lines = Volumetric data: The value of |wav|^2 at each discrete point on the grid 
    !Row length fixed by cube file format
    rw_len=6

    !Number of rows of 6 required for all NNnpts_r(3) values
    !Nrowz= NNnpts_r(3)/(rw_len-1)
    Nrowz=int(NNnpts_r(3)/rw_len)+1

    !If NNnpts_r(3) isn't divisible by 6, how many values will be put on the last row:
    rmdr=mod(NNnpts_r(3),rw_len)
    !Equivalent expression for rmdr:
    !rmdr=NN npts_r(3)-(int(NNnpts_r(3)/rw_len)*rw_len)

    do ix=1,NNnpts_r(1)
       do iy=1,NNnpts_r(2)
          do rowz=1,Nrowz    

             !Column counter for each row
             if(rowz<Nrowz)  jmax=rw_len
             if(rowz==Nrowz) jmax=rmdr

             !Form of the iz counter
             !iz=jz+((rowz-1)*rw_len)
             iz1=1+((rowz-1)*rw_len)
             iz2=jmax+((rowz-1)*rw_len)

             !Write the pdf_x1,y1,z1, pdf_x1,y1,z2,.. one row at a time, 
             !6 wav pdfs per row
             write(002,'(6E13.5)') orbdummy(ix,iy,iz1:iz2)

          enddo
       enddo
    enddo

    close(002)


    !Inform user of output
    write(*,'(A,I2,A,I2,A)') 'Ouputted overlap ',n1,' ',n2,' to file.'

  end subroutine cube_output



  !****************************************************************************
  !Construct orb_overlap matrix on NN grid
  !On-site approximation means that only |phi(r-r)ia)|^2 is computed 
  !=> 5 terms for the anion and five terms for the cation. 
  !****************************************************************************
  !Return array containing all permutations of overlaps between two NN orbitals
  !Return an index array 'orb_grd_pntr' which indexs points of overlap only
  !****************************************************************************

  subroutine construct_onsite_orb_products(orb_overlap,orb_grd_pntr,nnz)
    
    !Global variables and arrays
    use var,     only: NumNN,orb_bas
    !Orbital functions/profiles
    use orbital, only: radial_function,orbital_profile 
    !Place orbitals on grid
    use smat,    only: place_orbitals 

    implicit none


    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer                 Natm,Ntotal,Norb    
    integer                 ia,iorb,ja,jorb,n1,n2,ir,iro,ix,iy,iz,i,j
    integer                 x0,y0,z0,xf,yf,zf       !Unused but required

    !Global Arrays    
    integer,                allocatable :: midpnt(:)
    type(radial_function),  allocatable :: &
                            sorb(:,:,:),pxorb(:,:,:),pyorb(:,:,:),pzorb(:,:,:),sstorb(:,:,:)

    !Modular arrays
    integer,      allocatable, intent(out) :: orb_grd_pntr(:,:),nnz(:)
    real(dp),     allocatable, intent(out) :: orb_overlap(:,:)
    real(dp),     allocatable              :: orbdummy(:,:)
    
    !Local Arrays
    complex(dp),            allocatable :: orbital1(:,:,:),orbital2(:,:,:)
    


    !****************
    !Main Routine
    !****************
    !NN grid set up in calling routine and available to the module

    !Generate 3D orbital profiles
    call orbital_profile(NNrspacing,midpnt,sorb,pxorb,pyorb,pzorb,sstorb)

    !Number of atoms in the grid = Number of atomic species
    Natm=2

    !Total number of points in NN grid
    Ntotal=product(NNnpts_r)

    !Total number of orbital overlaps 
    Norb=orb_bas*Natm

    !Orbitals 1 and 2, filled per iteration. Dimensions of the  NN grid
    allocate(orbital1(NNnpts_r(1),NNnpts_r(2),NNnpts_r(3)),orbital2(NNnpts_r(1),NNnpts_r(2),NNnpts_r(3)))

    !orb_overlap has dimensions of the NN grid and contains on-site anion or cation function:
    !|phi(r-r(ia))|^2
    allocate(orb_overlap(Norb,Ntotal))

    !Spatial index, mapping the composite NN grid index to the composite overlap index  
    !Allocated with an upper bound
    allocate(orb_grd_pntr(Norb,Ntotal), nnz(Norb))

    orb_overlap(:,:)=0._dp
    orb_grd_pntr(:,:)=0
    nnz(:)=0


    !*****************************************************************
    !Loop over all atoms in a NN grid to create an array containing all 
    !orbital overlaps/products
    !Order is consistent with loops in orbital integral (important!)
    !*****************************************************************

    !Loop over anion and cation
    do ia=1,Natm
       do iorb=1,orb_bas

          !Composite state index 1
          n1=iorb+((ia-1)*orb_bas)

          !Place appropriate orbitals on grid: orbital1 
          !Second argument is equivalent and therefore a dummy
          call place_orbitals(orbital1,orbital2,NNatm_pntr,NNiptr,ia,ia,iorb,iorb,&
                              midpnt,sorb,pxorb,pyorb,pzorb,sstorb,&
                              x0,y0,z0,xf,yf,zf,atom_species_NN) 

         
          !Construct an array containing the overlap between two orbitals 
          !Loop over all points in the NNgrid. 
          iro=0
          ir=0
          
          do ix=1,NNnpts_r(1)
             do iy=1,NNnpts_r(2)
                do iz=1,NNnpts_r(3)
                   
                   !Composite index for NN grid
                   ir=ir+1

                   !If the product of |orbital1|^2 is non-zero store the it
                   if( real(orbital1(ix,iy,iz)*conjg(orbital1(ix,iy,iz)))/=0._dp)then
                      !Composite index for spatial overlap
                      iro=iro+1
                      !Product for ALL unique orbital overlaps
                      orb_overlap(n1,iro)=real(orbital1(ix,iy,iz)*conjg(orbital1(ix,iy,iz)))
                      !Map points of overlap from three spatial indices to one spatial index
                      orb_grd_pntr(n1,iro)=ir
                   endif

                enddo
             enddo
          enddo
          
          !Store number of non-zero values per orbital
          nnz(n1)=iro                            
          
       enddo
    enddo
    

    !Allocate correct dimension using the overlap
    !as stored in nnz(n1,n2)
    call reallocate(orb_grd_pntr,Norb,maxval(nnz))

    call reallocate(orb_overlap,Norb,maxval(nnz))


    !Return orb_overlap and orb_grd_pntr
  end subroutine construct_onsite_orb_products










  !*****************************************************************************
  !Test Routine 1. Approach 2
  !Integrate the orbital overlaps/products to reproduce the 
  !Overlap matrix elements. Do this on a system grid
  !
  !Details:
  !  Set up the system grid
  !  Generate orbital product array (keeping spatial indices)
  !  Compute integral{ n12(r) } dr where n12(r)=nprod(nij,ix,iy,iz)
  !  instead of orb_overlap(ni,nj,ir)
  !*****************************************************************************

  subroutine test_orb_construct_app2()
    use var,   only: dp,NumNN,orb_bas,Natoms,atom_species
    use grids, only: gen_NN_rspacing
    implicit none

    !******************
    !Declarations
    !******************
    
    !Local Variables
    integer   Ntotal,ir,jr,n12,N,ix,iy,iz, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb,&
              n1,n2,iro
    real(dp)  integral,dr

    !Modular Arrays
    integer,  allocatable :: grd_pntr(:,:),nnz(:)
    real(dp), allocatable :: nprod(:,:,:,:),S(:,:),S2(:,:)



    !******************
    !Main Routine
    !******************

    !Produce overlap matrix for comparison to output
    write(*,*)'Generating S matrix'
    call generate_Smatrix(S)
 
    !Construct NN grid (with spacing = system grid)
    !Declared within this module. 
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)
    
    !Construct system grid with vacuum region included
    !Declared within setup module:
    !snpts,sspacing,sgrid, satm_pntr, siptr, sxyz_ptr, sys_pntr,setup_system_grid
    !call setup_system_grid()
    call setup_arb_system_grid()

    !NN grid and system grid volume element
    dr=sspacing(1,1)*sspacing(2,2)*sspacing(3,3)

    !Construct orbital products for all n12, stored in an array with
    !dimensions of the system grid. 
    write(*,*) 'Generating orbital product on system grid'
    call orbitalproduct_sys(nprod)

    !Array indexing spatial points where nprod is non-zero
    call nz_indexing(nprod, grd_pntr,nnz)

    !Construct S2(n1,n2) = integral{ nprod(nij,r) } dr
    !Total number of orbitals in system
    N=Natoms*orb_bas

    !Total number of system grid points
    Ntotal=product(snpts)

    !Overlap Matrix
    allocate(S2(N,N))
    S2(:,:)=0._dp

    !Open files
    open(unit=001,file='s.dat')
    open(unit=002,file='s2.dat')


    !*********************************************
    !Loop over all orbital products
    !*********************************************
    write(*,*) 'Integrating orbital product'
    do n12=1,N*NumNN*orb_bas
          
       !Unwrap n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)

       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle

       !Global n1 and n2 values
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)

       !*********************************************
       !Integral
       integral=0._dp
       !Integrate over spatial variable r
!!$       do ir=1,Ntotal  
       do iro=1,nnz(n12)
          
          !Index mapping
          ir=grd_pntr(n12,iro)
       
          !Map composite spatial index to (x,y,z) indices
          ix=siptr(ir,1)
          iy=siptr(ir,2)
          iz=siptr(ir,3)  

          !Sum integral to produce overlap matrix for orbitals n12=(n1,n2)
          integral=integral+nprod(n12,ix,iy,iz)
       enddo
       !*********************************************

       !Overlap Matrix
       S2(n1,n2)=integral*dr

       !Output overlap matrix  '(I3,X,I3,X,F12.6)'
       write(001,*) n1,n2, S(n1,n2)
       write(002,*) n1,n2, S2(n1,n2)

    enddo

    !Close files
    close(001)
    close(002)
    write(*,*)'Test routine finished'
  end subroutine test_orb_construct_app2





  !*****************************************************************************
  !Test Routine 2. Approach 2
  !System grid approach
  !Compute integral { n12(r) n34(r') } dr dr' = S(n1,n2)*S(n3,n4)
  !More sensible loop structure than prior attempt (removed from code)
  !*****************************************************************************

  subroutine test2_orb_construct_app2()
    use var,   only: dp,NumNN,orb_bas,Natoms,atom_species
    use grids, only: gen_NN_rspacing
    implicit none

    !******************
    !Declarations
    !******************
    
    !Local Variables
    integer   Ntotal,ir,jr,n12,n34,N,ix,iy,iz, n1_l,n2_l,n3_l,n4_l,ia_l,ja_l,ka_l,la_l,ia,ja,ka,&
              la,iorb,jorb,korb,lorb,n1,n2,n3,n4,jx,jy,jz,iro,jro
    real(dp)  integral,dr

    !Modular Arrays
    integer,  allocatable :: grd_pntr(:,:), nnz(:)
    real(dp), allocatable :: nprod(:,:,:,:),S(:,:),S2(:,:)



    !******************
    !Main Routine
    !******************

    !Produce overlap matrix for comparison to output
    write(*,*)'Generating S matrix'
    call generate_Smatrix(S)
 
    !Construct system grid with vacuum region included
    call setup_arb_system_grid()

    !NN grid and system grid volume element
    dr=sspacing(1,1)*sspacing(2,2)*sspacing(3,3)

    !Construct orbital products for all n12, stored in an array with
    !dimensions of the system grid. 
    write(*,*) 'Generating orbital product on system grid'
    call orbitalproduct_sys(nprod)

    !Total number of orbitals in system
    N=Natoms*orb_bas

    !Total number of system grid points
    Ntotal=product(snpts)

    !Array indexing spatial points where nprod is non-zero
    call nz_indexing(nprod, grd_pntr,nnz)


    !*********************************************
    !Construct S2(n1,n2)
    !*********************************************
    write(*,*) 'Step 1. Integrating S2(n1,n2)'

    !Overlap Matrix 1
    allocate(S2(N,N))
    S2(:,:)=0._dp


    do n12=1,N*NumNN*orb_bas
          
       !Unwrap n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)

       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle

       !Global n1 and n2 values
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)

       !*********************************************
       !Integral
       integral=0._dp
       !Integrate over spatial variable r
!!$       do ir=1,Ntotal
       do iro=1,nnz(n12)
       
          !Map non-zero space onto full space
          ir=grd_pntr(n12,iro)

          !Map composite spatial index to (x,y,z) indices
          ix=siptr(ir,1)
          iy=siptr(ir,2)
          iz=siptr(ir,3)  

          !Sum integral to produce overlap matrix for orbitals n12=(n1,n2)
          integral=integral+nprod(n12,ix,iy,iz)
       enddo
       !*********************************************

       !Overlap Matrix
       S2(n1,n2)=integral*dr

    enddo


    !Open files
    open(unit=001,file='s.dat')
    open(unit=002,file='s2.dat')


    !*****************************************************************
    !Construct integral{ n34(r') } dr' * S2(n1,n2) = S(n1,n2)S(n3,n4) 
    !A bit of a trivial test but an intermediate step nonetheless
    !*****************************************************************
    
    !Composite Index 1
    do n12=1,N*NumNN*orb_bas

       !Indices w.r.t. n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)
       if(ja==0)cycle
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)


       !Composite Index 2
       do n34=n12,N*NumNN*orb_bas

          !Indices w.r.t. n34
          call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
          if(la==0)cycle
          n3=korb+((ka-1)*orb_bas)
          n4=lorb+((la-1)*orb_bas)


          !*******************************
          !Integrate over r'
          !*******************************
          integral=0._dp
!!$          do jr=1,Ntotal
          do jro=1,nnz(n34)

             jr=grd_pntr(n34,jro)

             !Map composite spatial index to (x,y,z) indices
             jx=siptr(jr,1)
             jy=siptr(jr,2)
             jz=siptr(jr,3)  
             
             !Sum integral to produce overlap matrix for orbitals n34=(n3,n4)
             integral=integral + nprod(n34,jx,jy,jz)
          enddo


          !Output values
          !Product of overlap matrices
          write(001,*)  n12,n34, S(n1,n2)*S(n3,n4)
          !Computed integrals: S3(n3,n4)*S2(n1,n2) where S3(n3,n4)=integral*dr
          write(002,*)  n12,n34, (integral*dr)*S2(n1,n2)

       enddo
    enddo


    !Close files
    close(001)
    close(002)

    write(*,*)'App2.Test2. Routine finished'
  end subroutine test2_orb_construct_app2

 






  !************************************************
  !Subroutine to generate the full overlap matrix 
  !(with no contribution from exponentials)
  !************************************************

  subroutine generate_Smatrix(S)

    use var,     only: NumNN,orb_bas,Natoms,NNdp_list,atom_type,NN_list
    use smat,    only: NN_disp,NN_smatrix
    implicit none

    !*******************************
    !Declarations
    !*******************************

    !Local Variables
    integer        n12,N,n1_l,n2_l,ia,ja,iorb,jorb,n1,n2,iro, i,l,jNN,m,row_cnt,col_cnt,j,&
                   ja_local,ia_local,i_l,j_l,n34,n3_l,n4_l,ka,la,n3,n4,korb,lorb,jro,ka_l,la_l,&
                   ia_l,ja_l
    real(16)       integral,dV

    !Local arrays
    complex(dp),  allocatable  :: S_NN(:,:)
    
    !Modular arrays
    real(dp),     allocatable,  intent(out)  :: S(:,:)


    !***************************************************************************
    !Generate overlap matrix S(n1,n2)
    !****************************************************************************

    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)
    
    !Construct NN S matrix (allocated in smat.f90)
    !Modified the routine to use the same NN grid as 'orb_overlap'
    call NN_smatrix(S_NN)

    !Allocate S
    N=Natoms*orb_bas
    allocate(S(N,N))
    S(N,N)=0._dp


    do i=1,Natoms
       do l=1,orb_bas
          do jNN=1,NumNN

             !Skip if NN doesn't exist
             j=NN_list(i,jNN)
             if(j==0)cycle
             
             do m=1,orb_bas

                !Shouldn't need local-global mapping for super cells
                !nor element checker for finite system

                !Overlap counters for S:
                row_cnt= l + ((i-1)*(orb_bas))
                col_cnt= m + ((j-1)*(orb_bas))

       
                !Map atomic index j to local NN grid index
                !Required for the indices used in S_NN
                ja_local=NNdp_list(i,jNN)

                
                !i will always relate to the central atom in S_NN, 
                ia_local=1
                      
                !Orbital indices for S_NN
                if(atom_type(i)==-1)then
                   i_l= l + ((ia_local-1)*(orb_bas))
                   j_l= m + ((ja_local-1)*(orb_bas))
                   !If cation, NNdp_list will give local atomic indices relative to cation centre
                   !Interchange them to make them correct in context to S_NN.  
                elseif(atom_type(i)==1 .and. i/=j)then
                   i_l= l + ((ja_local-1)*(orb_bas))
                   j_l= m + ((ia_local-1)*(orb_bas))
                !On-site cation orbitals. Required for debugging, when explicitly
                !comparing computed integrals instead of forcing terms to 1/0 
                !via orthonormality
                elseif(atom_type(i)==1 .and. i==j)then
                   i_l= l + orb_bas
                   j_l= m + orb_bas
                endif

                !Assign overlap matrix element
                S(row_cnt,col_cnt)=real(S_NN(i_l,j_l))
                !write(100,'(4(I3,X),F12.6,X)') i,j,l,m,S(row_cnt,col_cnt)
                
             enddo
          enddo
       enddo
    enddo

    !Return Overlap matrix S
  end subroutine generate_Smatrix




  
  !***************************************************************************
  !Define an array with dimensions of the system which stores all
  !orbital PRODUCTS per nij
  !
  !Details: 
  !  Assumes calling routine has set up NN grid. Available in module.
  !  Assumes calling routine has set up sys grid. Available in module.
  !  Sets up orbital profiles and places them on NN grid arrays
  !  Maps the product of orbitals n1 and n2 onto an array with 
  !  dimensions of the system grid
  !
  !Output
  !nprod(1:N*NumNN*orb_bas, ix,iy,iz)
  !***************************************************************************
 
  subroutine orbitalproduct_sys(nprod)
    use var,     only: dp,d,orb_bas,NumNN,Natoms,ns,atom_species
    use orbital, only: radial_function, orbital_profile
    use smat,    only: place_orbitals
    use evec,    only: output_wav_cube
    implicit none

    !**************
    !Declarations
    !**************

    !Local Variables
    integer      n12,N,ix1,ix2,iy1,iy2,iz1,iz2,ixm_NN,iym_NN,izm_NN,ixm_sys,&
                 iym_sys,izm_sys,x0,y0,z0,xf,yf,zf,n1_l,n2_l,ia_l,ja_l,ia,ja,&
                 iorb,jorb,gx,gy,gz,n1,n2

    !Local Arrays
    complex(dp), allocatable :: orbital1(:,:,:),orbital2(:,:,:)

    !Modular Arrays
    real(dp),    allocatable,intent(out) :: nprod(:,:,:,:) 
    real(dp),     allocatable            :: orbdummy(:,:,:)

    !Orbital arrays
    integer,                allocatable :: midpnt(:)
    type(radial_function),  allocatable :: &
                            sorb(:,:,:),pxorb(:,:,:),pyorb(:,:,:),pzorb(:,:,:),sstorb(:,:,:)


    !**************
    !Main Routine
    !**************

    !Total number of orbitals in system
    N=Natoms*orb_bas

    !Generate 3D orbital profiles.
    !sspacing is available from module and is same for system and NN grids
    call orbital_profile(sspacing,midpnt,sorb,pxorb,pyorb,pzorb,sstorb)

!!$    write(*,*) 'Testing with square functions. See line 1520 in exbnd/orbitalproduct_sys'
!!$    sorb(:,:,:)%an=1._dp
!!$    pxorb(:,:,:)%an=1._dp
!!$    pyorb(:,:,:)%an=1._dp
!!$    pzorb(:,:,:)%an=1._dp
!!$    sstorb(:,:,:)%an=1._dp
!!$    sorb(:,:,:)%cat=1._dp
!!$    pxorb(:,:,:)%cat=1._dp
!!$    pyorb(:,:,:)%cat=1._dp
!!$    pzorb(:,:,:)%cat=1._dp
!!$    sstorb(:,:,:)%cat=1._dp


    !Orbital Product with dimensions of system 
    !Stupidly big array that won't scale beyond tests
    allocate(nprod(N*NumNN*orb_bas,snpts(1),snpts(2),snpts(3)))
    nprod(:,:,:,:)=0._dp

    !Orbitals 1 and 2, filled per iteration. Dimensions of the system grid
    allocate(orbital1(snpts(1),snpts(2),snpts(3)),orbital2(snpts(1),snpts(2),snpts(3)))

    !Dummy array for cube outputting
    !allocate(orbdummy(snpts(1),snpts(2),snpts(3)))
    !orbdummy(:,:,:)=0._dp


    !For each orbital pair product, store the product on a grid
    !with dimensions of the system
    do n12=1,N*NumNN*orb_bas


       !Return indices warpped in composite index n12
       !Local n1 and n2 indices. Global (ia,iorb, ja,jorb) indices
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)

       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle

       !Global n1 and n2 values
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)

       !Place orbitals on orbital1 and orbital2, defined with dimensions of the full system
       !This avoids the need to map a second time, below.
       !As such, out of bounds also shouldn't need to be considered
       call place_orbitals(orbital1,orbital2,satm_pntr,siptr,&
                           ia,ja,iorb,jorb,midpnt,sorb,pxorb,pyorb,pzorb,sstorb,&
                           x0,y0,z0,xf,yf,zf,atom_species) 


       !Fill orbital product at correct points in the system grid
       nprod(n12,:,:,:)= real(orbital1(:,:,:))*real(orbital2(:,:,:))


       !Dummy array for visualisation
       !orbdummy(:,:,:)= real(orbital1(:,:,:))*real(orbital2(:,:,:))
       !Cube file outputting
       !orbdummy(:,:,:)=real(orbital1(:,:,:))
       !On-site cube files
       !if(n1==n2)then
       !   call cube_system_output(orbdummy,n1,n2) 
       !endif

    enddo



    !Return nprod
  end subroutine orbitalproduct_sys







  !*****************************************************************************************
  !Check on-site anion and cation integrals. Should give an indication of which routine is
  !mixing the cations up
  !*****************************************************************************************
  subroutine check_onsite
   use var,only:dp,d
    use orbital, only: radial_function, orbital_profile
    implicit none


    !******************************
    !Declarations
    !******************************

    !Local Variables
    integer                  ib,ia,ix1,ix2,iy1,iy2,iz1,iz2,ix,iy,iz,ir,Ntotal
    real(dp)                 integral,dr

    !Local Arrays
    real(dp),               allocatable :: orbital1(:,:,:)

    !Orbital arrays
    integer,                allocatable :: midpnt(:)
    type(radial_function),  allocatable :: &
                            sorb(:,:,:),pxorb(:,:,:),pyorb(:,:,:),pzorb(:,:,:),sstorb(:,:,:)


    !******************************
    !Main Routine
    !******************************

    !Construct system grid with vacuum region included
    call setup_arb_system_grid()


    !Generate 3D orbital profiles.
    !sspacing is available from module and is same for system and NN grids
    call orbital_profile(sspacing,midpnt,sorb,pxorb,pyorb,pzorb,sstorb)

    !Allocate grid function    
    allocate(orbital1(snpts(1),snpts(2),snpts(3)))



    !******************************
    !Place orbital at atom ia
    !******************************
    !Atom 
    ia=5

    !Point to atom 'iatm' in grid array
    ib=satm_pntr(ia)

    !Orbital 1, centred at ia, defined on sampling grid
    ix1= siptr(ib,1)-(midpnt(1)-1)
    ix2= siptr(ib,1)+(midpnt(1)-1)
    iy1= siptr(ib,2)-(midpnt(2)-1)
    iy2= siptr(ib,2)+(midpnt(2)-1)
    iz1= siptr(ib,3)-(midpnt(3)-1)
    iz2= siptr(ib,3)+(midpnt(3)-1)


    !*****************
    !Do integrals
    !*****************

    !Volume element
    dr=sspacing(1,1)*sspacing(2,2)*sspacing(3,3)

    !Total number of grid points
    Ntotal=product(snpts)

    !Anion s
    orbital1(:,:,:)=0._dp
    orbital1(ix1:ix2, iy1:iy2, iz1:iz2)=sorb(:,:,:)%an
    !Integrate
    integral=0._dp
    do ir=1,Ntotal
       !Map composite spatial index to (x,y,z) indices
       ix=siptr(ir,1)
       iy=siptr(ir,2)
       iz=siptr(ir,3)  
       !Sum integral to produce overlap matrix for orbitals n12=(n1,n2)
       integral=integral+(orbital1(ix,iy,iz)*orbital1(ix,iy,iz))
    enddo
    integral=integral*dr
    write(*,*) 'Integral for Anion s',integral


    !Anion px
    orbital1(:,:,:)=0._dp
    orbital1(ix1:ix2, iy1:iy2, iz1:iz2)=pxorb(:,:,:)%an
    !Integrate
    integral=0._dp
    do ir=1,Ntotal
       !Map composite spatial index to (x,y,z) indices
       ix=siptr(ir,1)
       iy=siptr(ir,2)
       iz=siptr(ir,3)  
       !Sum integral to produce overlap matrix for orbitals n12=(n1,n2)
       integral=integral+(orbital1(ix,iy,iz)*orbital1(ix,iy,iz))
    enddo
    integral=integral*dr
    write(*,*) 'Integral for Anion px',integral


    !Cation s
    orbital1(:,:,:)=0._dp
    orbital1(ix1:ix2, iy1:iy2, iz1:iz2)=sorb(:,:,:)%cat
    !Integrate
    integral=0._dp
    do ir=1,Ntotal
       !Map composite spatial index to (x,y,z) indices
       ix=siptr(ir,1)
       iy=siptr(ir,2)
       iz=siptr(ir,3)  
       !Sum integral to produce overlap matrix for orbitals n12=(n1,n2)
       integral=integral+(orbital1(ix,iy,iz)*orbital1(ix,iy,iz))
    enddo
    integral=integral*dr
    write(*,*) 'Integral for cation s',integral


    !Cation px
    orbital1(:,:,:)=0._dp
    orbital1(ix1:ix2, iy1:iy2, iz1:iz2)=pxorb(:,:,:)%cat
    !Integrate
    integral=0._dp
    do ir=1,Ntotal
       !Map composite spatial index to (x,y,z) indices
       ix=siptr(ir,1)
       iy=siptr(ir,2)
       iz=siptr(ir,3)  
       !Sum integral to produce overlap matrix for orbitals n12=(n1,n2)
       integral=integral+(orbital1(ix,iy,iz)*orbital1(ix,iy,iz))
    enddo
    integral=integral*dr
    write(*,*) 'Integral for Cation px',integral


  end subroutine check_onsite








  !*************************************************************************************
  !Routine to output function on system grid defined with the cubic unit cell 
  !+ vacuum layer, only. Origin defined with this specifically in mind
  !In future, modify this routine so all the required arrays can be passed as arguments
  !*************************************************************************************

  subroutine cube_system_output(orbdummy,n1,n2)  
    use var,only:dp,d,NumNN,al,r,Natoms,atom_type  
    implicit none


    !**********************
    !Declarations
    !**********************

    !Both are modular arrays 
    !integer,     allocatable, intent(in) :: snpts(:)
    !real(dp),    allocatable, intent(in) :: sspacing(:,:)

    !Modular variables
    integer,     intent(in) :: n1,n2

    !Modular arrays
    real(dp),    allocatable, intent(in) :: orbdummy(:,:,:)

    !Local arrays
    real(dp),    allocatable             :: origin(:)

    !Local variables
    integer             :: atnum,iatm,rw_len,Nrowz,rmdr,rowz,jz,jmax,iz1,iz2,nvx
    integer             :: i,ix,iy,iz
    real(dp)            :: Bohr,vac
    character           :: n1_lab*5, n2_lab*5




    !**************************************
    !Main Routine
    !**************************************

    !Convert Angstroms to Bohr
    Bohr = 1./0.5291772109217

    !Label the output with orbital and site index
    write(n1_lab,'(I2)')   n1
    write(n2_lab,'(I2)')   n2

    !Open file to write orbital iorb, at atom ia to. 
    open(unit=002,file='orb_'//trim(adjustl(n1_lab))//'_'//trim(adjustl(n2_lab))//'.cube')


    !********************************
    !Atomic positions from var.f90.
    !Grid from grids.f90
    !********************************
    !Define origin explicitly 
    allocate(origin(d))

    !Vacuum layer for arbitrary supercells = NN cut-off
    vac=al/sqrt(8._dp)
    
    !Number of points required to extend a vacuum to >= NN cutoff
    !for a supercell spacing 'rspacing'
    !Relies on same spacing for all dimensions
    nvx=ceiling(vac/sspacing(1,1))

    !Origin of grid with vacuum layer = (0,0,0) - vacuum layer
    origin(:)= (-nvx*sspacing(1,:)) + (-nvx*sspacing(2,:)) +(-nvx*sspacing(3,:))

    !Convert to Bohr
    origin(:)=origin(:)*Bohr



    !**************************************
    !Beginning of output to cube file
    !**************************************

    !Lines 1 & 2 are optional comments
    write(002,'(A,X,I5,X,A,X,I5)') 'Orbital cube file for overlap:',n1,' ',n2
    !Second line is blank
    write(002,*)

    !Line 3. Natoms with minus sign infront: Seems to be important if outputted data is in Bohr
    !Origin of the voxel space == grid. Cube file regenerates it from its origin,
    !starting vector, number of points and the grid (==lattice) vectors.
    write(002,'(I5,3(F12.6))') -1*Natoms, origin(:)


    !Lines 4,5,6: Number of data points along the lattice vectors r1, r2 and r3. 
    !as defined by rspacing(i,:). -1*npts_r(i) indicates Angstroms didn't work.
    do i=1,d
       write(002,'(I5,3(F12.6))')  snpts(i), sspacing(i,:)*Bohr
    enddo

    !Next block of lines give the atomic number, dble(atomic number) and the atomic positions 
    !Atomic numbers are dummy values as they're not important
    do iatm=1,Natoms
       if(atom_type(iatm)==1  ) atnum=1
       if(atom_type(iatm)==-1 ) atnum=2
       write(002,'(I5,4(F12.6))') atnum, dble(atnum), r(iatm,:)*Bohr
    enddo

    !Number of (molecular) orbitals evaluated per file and the index for the orbital
    write(002,'(3I5)') 1,n1,n2


    !Remaining lines = Volumetric data: The value of |wav|^2 at each discrete point on the grid 
    !Row length fixed by cube file format
    rw_len=6

    !Number of rows of 6 required for all npts_r(3) values
    !Nrowz= npts_r(3)/(rw_len-1)
    Nrowz=int(snpts(3)/rw_len)+1

    !If npts_r(3) isn't divisible by 6, how many values will be put on the last row:
    rmdr=mod(snpts(3),rw_len)
    !Equivalent expression for rmdr:
    !rmdr= npts_r(3)-(int(npts_r(3)/rw_len)*rw_len)

    do ix=1,snpts(1)
       do iy=1,snpts(2)
          do rowz=1,Nrowz    

             !Column counter for each row
             if(rowz<Nrowz)  jmax=rw_len
             if(rowz==Nrowz) jmax=rmdr

             !Form of the iz counter
             !iz=jz+((rowz-1)*rw_len)
             iz1=1+((rowz-1)*rw_len)
             iz2=jmax+((rowz-1)*rw_len)

             !Write the pdf_x1,y1,z1, pdf_x1,y1,z2,.. one row at a time, 
             !6 wav pdfs per row
             write(002,'(6E13.5)') orbdummy(ix,iy,iz1:iz2)

          enddo
       enddo
    enddo

    close(002)
 
    !Inform user of output
    write(*,'(A,I2,A,I2,A)') 'Ouputted overlap ',n1,' ',n2,' to file.'

  end subroutine cube_system_output





  !*****************************************************************************
  !Approach 2. System grid approach
  !Full routine to compute Coulomb integral.
  !
  !Compute integral { n12(r) n34(r')/|r-r'| } dr dr' in two steps:
  !1) V_H(n12,r') = integral{ n12(r)/|r-r'| } dr
  !2) w_orb(n1234)= integral{ n34(r') V_H(r') } dr'
  !
  !Once debugging has finished, modify to return w_orb
  !
  !*****************************************************************************

  subroutine coulomb_integral_NN_a2()
    use var,   only: dp,d,NumNN,orb_bas,Natoms,atom_species,NNdp_list,all_pnts,pi
    use grids, only: gen_NN_rspacing,return_grid_xyz,return_grid_xyz_arb
    use smat,  only: NN_disp
    use sort,  only: reallocate
    implicit none


    !******************
    !Declarations
    !******************
    
    !Local Variables
    integer   Ntotal,ir,jr,n12,n34,N,ix,iy,iz, n1_l,n2_l,n3_l,n4_l,ia_l,ja_l,ka_l,la_l,ia,ja,ka,&
              la,iorb,jorb,korb,lorb,n1,n2,n3,n4,jx,jy,jz,NOI,n1234,iro,jro, istart,istop
    real(dp)  integral,dr,x,y,z,orbij,r_sph
    character fname*12,rank_lab*3


    !Modular Arrays
    integer , allocatable :: grd_pntr(:,:), nnz(:)
    real(dp), allocatable :: nprod(:,:,:,:),V_H(:,:),r_i(:),r_j(:)
    real(dp), allocatable :: worb(:)

    !Local Arrays
    integer,  allocatable :: tmp(:)

    !******************
    !Main Routine
    !******************

    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    !Not utilised here but required by indices subroutine
    call NN_disp(NNdp_list)
 
    !For testing: Construct system grid with vacuum region included
    !call setup_arb_system_grid()
    !write(*,*)'Should be used with return_grid_xyz_arb'

    !Construct system grid for non-arbitrary shapes
    call setup_BSEtest_system_grid()
    write(*,*)'Should be used with return_grid_xyz'
    !Only compares to the NN grid approach as all aux
    !grid arrays are stored with all supercell points
    deallocate(sys_pntr,gsindex)
    

    !NN grid and system grid volume element. dr == dr'
    dr=sspacing(1,1)*sspacing(2,2)*sspacing(3,3)

    !Construct orbital products for all n12, stored in an array with
    !dimensions of the system grid. 
    if(IO) write(*,*) 'Generating orbital product on system grid'
    call orbitalproduct_sys(nprod)

    !Total number of orbitals in system
    N=Natoms*orb_bas

    !Total number of system grid points
    Ntotal=product(snpts)

    !Maximum number of orbital interactions for two centres
    NOI=N*NumNN*orb_bas

    !Generate array which indexes non-zero spatial points in nprod 
    call nz_indexing(nprod, grd_pntr,nnz)

    !Spatial Variable Arrays
    allocate(r_i(d),r_j(d))

    !Radius of cubic volume element, given as the length from 
    !the centre to the corner (see Brillouin approximation)
    r_sph=sqrt(3._dp)*0.5*sspacing(1,1)



    if(IO) write(*,*) 'App2'

    !*********************************************
    !Construct Hartree Potential V_H(n12,r')
    !*********************************************
   if(IO) write(*,*) 'Constructing Hartree Potential'

    !Allocate memory (massive relative amount)
    allocate(V_H(NOI,Ntotal))
    V_H(:,:)=0._dp

    !Shouldn't matter if I loop over jr or n12 first
    !Gone with n12 as => less evaluation of indices
       
   !Loop over composite index 1
   do n12=1,NOI
       
       !Unwrap n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)
       
       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle
       
       !Global n1 and n2 values
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)

       !Loop over all spatial variable r'
       do jr=1,Ntotal

          !Unwrap composite index for system grid: jr->(jx,jy,jz)
          call comp_to_ijk(jr,snpts(1),snpts(2),snpts(3),jx,jy,jz)

          !Return x,y,z variables associated with spatial indices
          call return_grid_xyz(jx,jy,jz,sspacing,snpts, x,y,z)
          r_j(:)=(/x,y,z/)


          !*********************************************
          !Integral
          integral=0._dp
          !Initialise orbital product
          orbij=0._dp

          !Integrate over spatial variable r
          do iro=1,nnz(n12)

             !Map non-zero product space to full space
             ir=grd_pntr(n12,iro)
                
             !Map composite spatial index to (x,y,z) indices
             ix=siptr(ir,1)
             iy=siptr(ir,2)
             iz=siptr(ir,3)  

             !Return x,y,z variables associated with spatial indices
             call return_grid_xyz(ix,iy,iz,sspacing,snpts, x,y,z)
             r_i(:)=(/x,y,z/)
             
             !Avoid integrating over singularities (r=r')
             !Can do below as ir and jr on exactly the same grid
             if(ir==jr)then
                !Store 'product of orbitals' value
                orbij=nprod(n12,ix,iy,iz)
                !write(200,'(2(I6,X),7(F12.6,X))') ir,jr,r_i(:),r_j(:),orbij
                cycle
             endif

             !Sum integral to produce Hartree potential: integral{ orbs(n12,r)/|r-r'| } dr
             integral=integral+ (nprod(n12,ix,iy,iz)/&
                                 sqrt(dot_product(r_i(:)-r_j(:),r_i(:)-r_j(:))) )

          enddo
          !*********************************************

          !Define Hartree Potential
          ! + Brillouin Approximation for when r=r'
          V_H(n12,jr)=(integral*dr)+(orbij*2._dp*pi*r_sph**2)
          !write(100,*) n12,jr,integral*dr

       enddo
       !Shut loop over r'
    enddo
    !Shut loop over n12

 

!!$    write(*,*) 'Writing out the Hartree Potential'
!!$    jr=0
!!$    do ix=1,snpts(1)
!!$       !write(100,*)
!!$       do iy=1,snpts(2)
!!$          do iz=1,snpts(3)
!!$
!!$             !Composite spatial counter
!!$             jr=jr+1
!!$            
!!$             !Not an integral: Just summing contributions of all n12 at point jr
!!$             integral=0._dp
!!$             do n12=1,NOI
!!$                integral=integral+V_H(n12,jr)
!!$             enddo
!!$             
!!$             !Output summed potential at point jr
!!$             write(100,*) sgrid(:,jr),integral
!!$             !write(100,'(4(f20.10,X))') sgrid(:,jr),integral
!!$
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    write(*,*) 'exbnd stops line 2126'
!!$    stop



    !*****************************************************************
    !Compute 2nd integral: orbs(n34,r')V_H(n12,r') dr'
    !*****************************************************************
    write(*,*)'Computing Coulomb orbital integral'
    
    !Allocate worb with an upper bound (Some atoms won't have 4 neighbours)
!!$    allocate( worb(ceiling(0.5*(NOI**2))) )
!!$    worb(:)=0._dp


    !Open file
    open(unit=001,file='w2.dat')

    !Initialise 8-fold composite index
    n1234=0


    !Composite Index 1
    do n12=1,NOI

       !Unwrap composite index n12=(ia,iorb,ja,jorb)
       !plus local indices for use in NN_smatrix (not needed in this routine)
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)

       !If ja, neighbour of atom ia, is zero cycle. 
       if(ja==0)cycle

       !n12=(n1,n2)
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)
       

       !Composite Index 2
       do n34=n12,NOI
          
          !Indices w.r.t. n34
          call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
          if(la==0)cycle
          n3=korb+((ka-1)*orb_bas)
          n4=lorb+((la-1)*orb_bas)
          
          !8-fold composite index
          n1234=n1234+1

          !*******************************
          !Integrate over r'
          !*******************************
          integral=0._dp

          do jro=1,nnz(n34)
             
             !Map non-zero product space to full space
             jr=grd_pntr(n34,jro)

             !Map composite spatial index r'=(x',y',z')
             jx=siptr(jr,1)
             jy=siptr(jr,2)
             jz=siptr(jr,3)  
             
             !Sum integral{ orbs(n34,r')V_H(n12,r') } dr'
             integral=integral + (nprod(n34,jx,jy,jz)*V_H(n12,jr))
          enddo

          !Store orbital Coulomb integral
!!$          worb(n1234)=integral*dr
          
          !Output to file for debugging
          write(001,*) n12,n34, integral*dr


       enddo
    enddo

    close(001)
    
    !Reallcate worb with appropriate size
!!$    call reallocate(worb,n1234)

    write(*,*)'Full Routine. Approach 2 finished'
  end subroutine coulomb_integral_NN_a2

 

  !*****************************************************************************
  !Approach 2 is very slow in the 1st loop when computing the full expression
  !sum_{n12,r'} integral{ n12(r)/|r-r'| } dr 
  !Integral can only be non-zero when numerator is non-zero, hence tabulate
  !all spatial ir indices for which n12(r) is non-zero
  !*****************************************************************************

  subroutine nz_indexing(nprod, grd_pntr,nnz)
    use var,   only: dp,Natoms,orb_bas,NumNN
    use sort,  only: reallocate
    implicit none

    !******************
    !Declarations
    !******************

    !Local Variables
    integer   N,NOI,Ntotal,nij,cnt,ix,iy,iz,ir

    !Modular Arrays
    integer,  allocatable, intent(out) :: grd_pntr(:,:),nnz(:)
    real(dp), allocatable, intent(in)  :: nprod(:,:,:,:)

    

    !******************
    !Main Routine
    !******************

    !Total number of orbitals in system
    N=Natoms*orb_bas

    !Maximum number of orbital interactions for two centres
    NOI=N*NumNN*orb_bas

    !Total number of system grid points
    Ntotal=product(snpts)

    !Allocate arrays. 2nd index of grd_pntr is an upper bound
    allocate(grd_pntr(NOI,int(0.5*Ntotal)),nnz(NOI))


    !Loop over composite index 1
    do nij=1,NOI
       cnt=0
       do ir=1,Ntotal
          ix=siptr(ir,1)
          iy=siptr(ir,2)
          iz=siptr(ir,3)     
          if(nprod(nij,ix,iy,iz)/=0._dp)then
             cnt=cnt+1
             grd_pntr(nij,cnt)=ir
          endif
       enddo
       nnz(nij)=cnt
    enddo

    !Reallocate 2nd index of grd_pntr
    call reallocate(grd_pntr,NOI,maxval(nnz))


  end subroutine nz_indexing
  


  !*****************************************************************************
  !|Map spatial index of NN grid to spatial index of system grid, using common
  !atomic centre indices
  !*****************************************************************************

  subroutine MapSpatialIndices_NN_sys(NNatm_pntr,satm_pntr,NNiptr,siptr,ia,ia_l,ix1,iy1,iz1)
    implicit none

    !Modular variables
    integer, intent(in)  :: ia,ia_l
    integer, intent(out) :: ix1,iy1,iz1
  
    !Global Arrays
    integer, intent(in), allocatable ::NNatm_pntr(:),satm_pntr(:),NNiptr(:,:),siptr(:,:)

    !Starting indices for system grid = 
    !Atomic centre on system grid - (atom centre on NN grid - 1)
    ix1=siptr(satm_pntr(ia),1) - (NNiptr(NNatm_pntr(ia_l),1)-1)
    iy1=siptr(satm_pntr(ia),2) - (NNiptr(NNatm_pntr(ia_l),2)-1)
    iz1=siptr(satm_pntr(ia),3) - (NNiptr(NNatm_pntr(ia_l),3)-1)
    
  end subroutine MapSpatialIndices_NN_sys



 !*****************************************************************************
  !|Map spatial index of NN grid to spatial index of system grid, using common
  !atomic centre indices and WITHOUT the need for siptr or NNiptr
  !*****************************************************************************

  subroutine MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,ia,ia_l,ix1,iy1,iz1)
    implicit none

    !Modular variables
    integer, intent(in)  :: ia,ia_l
    integer, intent(out) :: ix1,iy1,iz1

    !Local Variables
    integer  ix_s,iy_s,iz_s,ix_NN,iy_NN,iz_NN
  
    !Global Arrays
    integer, intent(in), allocatable ::NNatm_pntr(:),satm_pntr(:),NNnpts_r(:),snpts(:)

    !Convert composite spatial index locating atom ia into three spatial indices 
    !On system grid
    call comp_to_ijk(satm_pntr(ia),snpts(1),snpts(2),snpts(3), ix_s,iy_s,iz_s)
    !On NN grid
    call comp_to_ijk(NNatm_pntr(ia_l),NNnpts_r(1),NNnpts_r(2),NNnpts_r(3), ix_NN,iy_NN,iz_NN)

    !Starting indices for system grid = 
    !Atomic centre on system grid - (atom centre on NN grid - 1)
    ix1 = ix_s-(ix_NN-1)
    iy1 = iy_s-(iy_NN-1)
    iz1 = iz_s-(iz_NN-1)
    
  end subroutine MapSpatialIndices_NN_sys2




  !*****************************************************************************
  !Orbital Coulomb Integral using NN grid (1st approach) but constructs |r-r'|
  !in a different manner 
  !
  !Can use coulomb_integral_NN_a2() and stop after the Hartree construction
  !to provide a direct comparison
  !*****************************************************************************

  subroutine coulomb_integral_NN_a1_ALTr()
    use var,   only: sp,dp,qp,Natoms,orb_bas,evc,atom_species,d,Nexsp,eigen_cnt,NN_list,NumNN,&
                     r,time_CoulInt,al,NNdp_list,atom_type,pi,all_pnts
    use grids, only: gen_NN_grid, gen_NN_rspacing, return_grid_xyz,return_grid_xyz_arb
    use smat,  only: NN_disp
    use sort,  only: reallocate
    implicit none

    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer       N,Nm,iro,ir,jro,jr,  ia,ja,ka,la,  n1,n2,n3,n4,n12,n34 , iorb,jorb,korb,lorb,&
                  n1_l,n2_l,n3_l,n4_l,Natm,ia_l,ja_l,ka_l,la_l,n1234,NOI,ix1,iy1,iz1,ir_l,&
                  ix_l,iy_l,iz_l,Ntotal,ix,iy,iz,ixy_l, jx1,jy1,jz1,jr_l,jxy_l,jx_l,jy_l,jz_l,&
                  jx,jy,jz,istart,istop,Nreq
    real(dp)      dr,integral,x,y,z,denom,r_sph,orbij
    character     fname*12,rank_lab*3

    !Modular arrays
    integer,      allocatable         :: orb_grd_pntr(:,:,:),nnz(:,:)
    real(dp),     allocatable         :: orb_overlap(:,:,:)
    real(dp),     allocatable         :: worb(:)

    !Local Arrays
    integer,      allocatable         :: tmp(:)
    real(dp),     allocatable         :: V_H(:)
    real(dp),     allocatable         :: r_i(:),r_j(:)



    !**********************************
    !Main Routine
    !**********************************

    if(IO) write(*,*)'New Approach'

    !Dimensions of the system
    N=Natoms*orb_bas

    !Maximum number of orbital interactions for two centres, in NN regime
    NOI=N*(NumNN*orb_bas)
   
    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)

    !All arguments here are defined in the module
    !This NN grid uses the same spacing as the system grid 
    !- important for mapping between the system and NN grids
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)
    !Arrays not utilised so deallocate immediately
    deallocate(NNgrid,NNxyz_ptr)

    !Construct orb_overlap matrix on NN grid
    call construct_orb_overlap(orb_overlap,orb_grd_pntr,nnz)
    !Index array for NN grid not utilised after orbital products are constructed
    deallocate(NNiptr)


    !********************************************************
    !Set up system grid: 
    !Grid spacing, number of points and grid points 
    !that sample atomic centres
    !********************************************************

    !For testing: Construct system grid with vacuum region included
    !call setup_arb_system_grid()
    !write(*,*) 'Note, use return_grid_xyz_arb with setup_arb_system_grid'
    !deallocate(sgrid,siptr,sys_pntr)

    !Construct systme grid for non-arbitrary shapes
    !return_grid_xyz should be used with these calls
    !Two equivalent approaches:

    !Approach 1
    !Generate grid arrays with dimensions of the supercell
    !No auxilliary grid arrays are allocated, except atm_pntr
    all_pnts=.true.
    call setup_system_grid()

    !Approach 2.
    !Generate grid arrays with dimensions of the system+vacuum <
    !supercell
!!$    all_pnts=.false.
!!$    call setup_system_grid()
!!$    !Deallocate unused auxilliary grid arrays
!!$    deallocate(sgrid,siptr,sys_pntr)
!!$    !Loops below are only set up to work when atm_pntr is
!!$    !defined with points from the full supercell, not the smaller
!!$    !system+vacuum, so map sys+vac points in atm_pntr to supercell
!!$    !points
!!$    allocate(tmp(Natoms))
!!$    do ia=1,Natoms
!!$       tmp(ia)=gsindex(satm_pntr(ia))
!!$    enddo
!!$    satm_pntr(:)=0
!!$    satm_pntr(:)=tmp(:)
!!$    deallocate(tmp)

    !********************************************************

    !Total number of points in system grid
    !Note, can't define Ntotal this way if using other system grid routine
    !as number of points stored < product(snpts)
    Ntotal=product(snpts)

    !Volume element associated with this grid spacing
    dr=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Radius of cubic volume element, given as the length from 
    !the centre to the corner (see Brillouin approximation)
    r_sph=sqrt(3._dp)*0.5*NNrspacing(1,1)

 
    

    !*************************
    !Array Allocation
    !*************************

    !Position vectors for Coulomb Potential
    allocate(r_i(d),r_j(d))

    !Hartree Potential (Should reduce to single precision) 
    allocate(V_H(Ntotal))
    V_H(:)=0._dp

    !Orbital integrals stored in a packed, 1D array
    !Reduce to single precision?
    !allocate(worb(NOI**2/8._dp))

    !Open file
    open(unit=001,file='w1.dat')




    !****************************************************************************
    ! 1)   Compute Hartree Potential: V_H(r')=  integral{ orb(n12,r)/|r-r'| } dr
    ! 2)   Compute integral{ orb(n34,r') V_H(n12,r')} dr'      
    !****************************************************************************

    !Initalise full composite index
    n1234=0


    !First composite index over orbitals 1&2 defined with spatial coordinate r (ir)
    do n12=1,NOI

       !Return indices in composite index n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)
       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle
       

       !Map local NN grid spatial indices to system grid indices using atom ia 
       !on system grid and corresponding atom ia_l on NN grid
       call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,ia,ia_l,&
                                      ix1,iy1,iz1)

       
       !Loop over all points in system grid 
       do jr=1,Ntotal

          !Unwrap composite index for system grid: jr->(jx,jy,jz)
          call comp_to_ijk(jr,snpts(1),snpts(2),snpts(3),jx,jy,jz)

          !Return x,y,z variables associated with spatial indices
          call return_grid_xyz(jx,jy,jz,sspacing,snpts, x,y,z)
          r_j(:)=(/x,y,z/)


          !*********************************************
          !Integral over r
          !*********************************************
          !Initialise integral 
          integral=0._dp
          !Initialise orbital product
          orbij=0._dp

          !Loop over non-zero orbital products defined on NN grid
          do iro=1,nnz(n1_l,n2_l)

             !Map non-zero product space to NN grid space
             ir_l=orb_grd_pntr(n1_l,n2_l,iro)

             !Unwrap composite index for NN grid space: ir_l->(ix_l,iy_l,iz_l)
             call comp_to_ijk(ir_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                              ix_l,iy_l,iz_l)

             !Iterate system grid indices using (ix_l,iy_l,iz_l)
             ix=ix_l+(ix1-1)
             iy=iy_l+(iy1-1)
             iz=iz_l+(iz1-1)

             !Return x,y,z variables associated with spatial indices
             call return_grid_xyz(ix,iy,iz,sspacing,snpts, x,y,z)
             r_i(:)=(/x,y,z/)

             !Composite spatial counter for system grid
             ir=iz+( (iy+((ix-1)*snpts(2))-1)*snpts(3))

             !If |r-r'|=0
             if(ir==jr)then
                !Store 'product of orbitals' value
                orbij=orb_overlap(n1_l,n2_l,iro)
                !write(100,'(2(I6,X),7(F12.6,X))') ir,jr,r_i(:),r_j(:),orbij
                cycle
             endif

             !Sum Hartree Potential
             integral=integral+ (orb_overlap(n1_l,n2_l,iro)/&
                                 sqrt(dot_product(r_i(:)-r_j(:),r_i(:)-r_j(:))))
                  
          enddo
          !Shut integral over ir
          !*********************************************

          !Store Hartree potential at spatial point r' (jr)
          ! + Brillouin Approximation for when r=r'
          V_H(jr)=(integral*dr)+(orbij*2._dp*pi*r_sph**2)

       end do
       !Shut loop over jr



       !Loop over 2nd composite index n34
       do n34=n12,NOI

          !Return indices in composite index n34 
          call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
          !If ia doesn't have a neighbour la, cycle loop
          if(la==0)cycle

          !Map local NN grid spatial indices to system grid indices using atom ka 
          !on system grid and corresponding atom ka_l on NN grid
          call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,&
                                         ka,ka_l,jx1,jy1,jz1)


          !*****************************************
          !Integral over r'
          !*****************************************
          integral=0._dp

          !Compute integral over spatial variable r'
          do jro=1,nnz(n3_l,n4_l)

             !Map non-zero product space to NN grid space
             jr_l=orb_grd_pntr(n3_l,n4_l,jro)

             !Unwrap composite index for NN grid space
             call comp_to_ijk(jr_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                  jx_l,jy_l,jz_l)

             !Iterate system grid indices using (jx_l,jy_l,jz_l)
             jx=jx_l+(jx1-1)
             jy=jy_l+(jy1-1)
             jz=jz_l+(jz1-1)

             !Composite spatial counter for system grid
             jr=jz+( (jy+((jx-1)*snpts(2))-1)*snpts(3))

             !Integral
             integral=integral+ (orb_overlap(n3_l,n4_l,jro)*V_H(jr))

          enddo
          !*****************************************

          !Contiguous Composite index for all orbital interactions
          n1234=n1234+1

          !Store orbital integral (mathematically dr' but dr=dr')
          !worb(n1234)=integral*dr
          write(001,*) n12,n34,integral*dr

       enddo
    enddo
    !Shut loops over n12 and n34

    !Shut file
    close(001)

    !Reallocate worb with correct size
    !May be too large to reallocate
    !Better to write to memory, deallocate space, allocate
    !correct space and read back in?
    !call reallocate(worb,n1234)



    !Note on Brillouin Approximation for when r=r'
    !Cubic element approximated as spherical, with same radius.
    !Analytic integral for spherical element given as 4pi*(1/|r|)r**2dr, 
    !multiplied by orbital, which is assumed to be  slowly varying over 
    !this range, hence comes outside the integral
    !For Hartree potential, for a given point jr, there can only be one ir
    !where jr=ir, hence approximation can be added when storing V_H


    !End of the routine. Return worb
  end subroutine coulomb_integral_NN_a1_ALTr




  !*******************************************************************
  !Extract the indices from a generic, single composite index
  !Specific to when the composite index is comprised of three indices
  !Inputs:
  !ixyz       Composite Index
  !nx,ny,nz   Limit of of ix,iy,iz
  !
  !Outputs:
  !ix,iy,iz   Each of the three indices which comprise the composite
  !*******************************************************************
  subroutine comp_to_ijk(ixyz,nx,ny,nz, ix,iy,iz)
    implicit none

    !Inputs/Outputs
    integer, intent(in)  :: ixyz, nx,ny,nz
    integer, intent(out) :: ix,iy,iz
    !Local Variable
    integer              :: ixy

    !Conversion 
    ixy=((ixyz-1)/nz)+1
    iz=ixyz-((ixy-1)*nz)
    ix=((ixy-1)/ny)+1
    iy=ixy-((ix-1)*ny)

  end subroutine comp_to_ijk





 !*****************************************************************
  !Compute integral{ n12(r) n34(r') dr dr' 
  !*****************************************************************

  subroutine testr_app1()
    use var,     only: NumNN,orb_bas,Natoms,NNdp_list,atom_type,NN_list,d,al,r
    use grids,   only: gen_NN_rspacing
    use smat,    only: NN_disp
    implicit none


    !*******************************
    !Declarations
    !*******************************

    !Local Variables
    integer        n12,N,n1_l,n2_l,ia,ja,iorb,jorb,n1,n2,iro, i,l,jNN,m,row_cnt,col_cnt,j,&
                   ja_local,ia_local,i_l,j_l,n34,n3_l,n4_l,ka,la,n3,n4,korb,lorb,jro,ka_l,la_l,&
                   ia_l,ja_l,Natm,ir,NOI
    real(16)       integral,dV

    !Local arrays
    real(dp),     allocatable  :: r1(:),r_local(:,:)
    
    !Modular arrays
    integer,      allocatable  :: orb_grd_pntr(:,:,:),nnz(:,:)
    real(dp),     allocatable  :: orb_overlap(:,:,:)



    !*******************************
    !Main Routine
    !*******************************

    !Number of atoms in NN system
    Natm=NumNN
    !Atomic positions in NN system
    allocate(r_local(d,Natm))
    !r_local in units of Angstroms 
    r_local = 0.25*al*reshape( (/ 0.,  0.,  0., &   
                                  1.,  1.,  1., &
                                  1., -1., -1., &
                                 -1.,  1., -1., &
                                 -1., -1.,  1./), (/d,Natm/))


    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)
    
    !NN grid that uses the spacing of the system grid 
    write(*,*) 'Generating NN grid with spacing of system grid'
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)

    !Volume element associated with this grid spacing
    dV=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Construct orb_overlap matrix on NN grid
    write(*,*)'Constructing orbital product array'
    call construct_orb_overlap(orb_overlap,orb_grd_pntr,nnz)

    !Dimensions of the system
    N=Natoms*orb_bas

    !Maximum number of orbital interactions for two centres, in NN regime
    NOI=N*(NumNN*orb_bas)

    !Position vector r
    allocate(r1(d))
    r1(:)=0._dp


    !****************************************************************************
    !Compute integral{ n12(r)/|r| } dr 
    !****************************************************************************
    !Open file
    open(unit=001, file='w1.dat')    


    !4-fold Composite index. n12 = (n1,n2) = (ia,iorb,ja,jorb)
    do n12=1,NOI

       !Return indices warpped in composite index n12
       !Local n1 and n2 indices. Global (ia,iorb, ja,jorb) indices
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)

       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle

       !Global n1 and n2 values
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)

       !*******************************
       !Integral over r,r'
       !*******************************
       integral=0._dp
       
       do iro=1,nnz(n1_l,n2_l)
          !Composite spatial index for NN grid, indicating points where the product is non-zero
          ir=orb_grd_pntr(n1_l,n2_l,iro)
          !Vector r
          r1(:)= NNgrid(:,ir)-r_local(:,ia_l)+r(ia,:)

          !Integral
          if(r1(1)==0. .and. r1(2)==0. .and. r1(3)==0.) cycle
          integral=integral+( orb_overlap(n1_l,n2_l,iro)/sqrt(dot_product(r1(:),r1(:))) )
       enddo

       !Mulitply by volume element and write out
       write(001,'(4(I4,X),F20.10)') ia,iorb,ja,jorb,integral*dV
          
    enddo
    !Shut loop over n12sub


    !Close file
    close(001)
    

  end subroutine testr_app1



 !*****************************************************************
  !Compute integral{ n12(r) n34(r') dr dr' Approach 2
  !*****************************************************************
  subroutine testr_app2
    use var,   only: dp,NumNN,orb_bas,Natoms,atom_species,NNdp_list
    use grids, only: gen_NN_rspacing
    use smat,  only: NN_disp
    implicit none

    !******************
    !Declarations
    !******************
    
    !Local Variables
    integer   Ntotal,ir,jr,n12,N,ix,iy,iz, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb,&
              n1,n2,iro
    real(dp)  integral,dr

    !Modular Arrays
    integer,  allocatable :: grd_pntr(:,:),nnz(:)
    real(dp), allocatable :: nprod(:,:,:,:)



    !******************
    !Main Routine
    !******************


    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)
     
    !Construct system grid with vacuum region included
    call setup_arb_system_grid()

    !NN grid and system grid volume element
    dr=sspacing(1,1)*sspacing(2,2)*sspacing(3,3)

    !Construct orbital products for all n12, stored in an array with
    !dimensions of the system grid. 
    write(*,*) 'Generating orbital product on system grid'
    call orbitalproduct_sys(nprod)

    !Array indexing spatial points where nprod is non-zero
    call nz_indexing(nprod, grd_pntr,nnz)

    !Construct S2(n1,n2) = integral{ nprod(nij,r) } dr
    !Total number of orbitals in system
    N=Natoms*orb_bas

    !Total number of system grid points
    Ntotal=product(snpts)



    !*********************************************
    !Loop over all orbital products
    !*********************************************
    !Open file
    open(unit=001,file='w2.dat')

    do n12=1,N*NumNN*orb_bas
          
       !Unwrap n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)

       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle

       !Global n1 and n2 values
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)

       !*********************************************
       !Integral
       integral=0._dp
       !Integrate over spatial variable r
!!$       do ir=1,Ntotal
       do iro=1,nnz(n12)
          
          !Index mapping
          ir=grd_pntr(n12,iro)
       
          !Map composite spatial index to (x,y,z) indices
          ix=siptr(ir,1)
          iy=siptr(ir,2)
          iz=siptr(ir,3)  

          !Sum integral to produce overlap matrix for orbitals n12=(n1,n2)
          if(sgrid(1,ir)==0. .and. sgrid(2,ir)==0. .and.sgrid(3,ir)==0.)cycle
          integral = integral+(nprod(n12,ix,iy,iz)/&
                     sqrt(dot_product(sgrid(:,ir),sgrid(:,ir))))
       enddo
       !*********************************************

       !Output to file
       write(001,'(4(I4,X),F20.10)') ia,iorb,ja,jorb,integral*dr

    enddo

    !Close files
    close(001)

  end subroutine testr_app2

 
  !********************************************************************
  !Input real eigenvector array. Output electron and hole eigenvector
  !arrays, ordered consistently with the loop order of the two-particle
  !Hamiltonian 
  !********************************************************************
  subroutine eigenvector_order(evcr, evc_e,evc_h)
    use var, only:dp,Ne,Nh,Natoms,orb_bas,eigen_cnt
    implicit none

    !******************
    !Declarations
    !******************
    !Global Arrays
    real(dp), allocatable :: evcr(:,:),evc_e(:,:),evc_h(:,:)
 
    !Local variables
    integer   :: N,i,ie,VBT,CBB


    !******************
    !Main Routine
    !******************

    !Number of orbital states in system
    N=Natoms*orb_bas
    
    !Allocate arrays
    allocate(evc_e(Ne,N),evc_h(Nh,N))
    evc_e(:,:)=0._dp
    evc_h(:,:)=0._dp

    !Valence band top and conduction band bottom
    VBT=eigen_cnt(1,1)
    CBB=VBT+1

    !Fill electron array
    !i= electron index of H2p 
    do i=1,Ne
       !Map to single-particle index of evcr
       ie=i+CBB-1
       !Fill eigenvector array
       evc_e(i,1:N)=evcr(1:N,ie)
    enddo

    !Fill hole array
    do i=1,Nh
       ie=VBT-i+1
       evc_h(i,1:N)=evcr(1:N,ie)
    enddo

    !Return eigenvector arrays and free memeory
    !deallocate(evcr)
    if(IO)write(*,*)'Note, not deallocating evcr for testing, as it''s required by spectral routine.'
 
  end subroutine eigenvector_order




  !Store revelant electron and hole TB coefficients in 1D arrays
  subroutine organise_eigenvectors(evcr, c_e,c_h)
    use var, only: Ne,Nh,eigen_cnt,Natoms,orb_bas

    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer  ::  N,VBT,cnt,ia,iorb,iv,ic,i,ie

    !Modular Arrays
    real(dp),allocatable ::  c_h(:),c_e(:),evcr(:,:)


    !*****************
    !Main Routine
    !*****************
 
    !Total number of orbitals
    N=Natoms*orb_bas
    allocate(c_h(Nh*N),c_e(Ne*N))
    c_h(:)=0._dp
    c_e(:)=0._dp

    !VBT
    VBT=eigen_cnt(1,1)

    !Holes, ordered from VBT to VBT-Nh+1
    cnt=0
    do ia=1,Natoms
       do iorb=1,orb_bas
          do iv=1,Nh
             !Global orbital index
             i=iorb+((ia-1)*orb_bas)
             !Global single-particle index
             ie=VBT+1-iv
             !Index for 1D coefficient array
             cnt=cnt+1
             c_h(cnt)=evcr(i,ie)
          enddo
       enddo
    enddo

    !Electrons
    cnt=0
    do ia=1,Natoms
       do iorb=1,orb_bas
          do ic=1,Ne
             !Global orbital index
             i=iorb+((ia-1)*orb_bas)
             !Global single-particle index
             ie=VBT+ic
             !Index for 1D coefficient array
             cnt=cnt+1
             c_e(cnt)=evcr(i,ie)
          enddo
       enddo
    enddo

    !Will need to deallocate evcr and change the sp absorption definition
    !used by BSE routine

  end subroutine organise_eigenvectors





  !*****************************************************************************
  ! Construct H2p.
  !
  ! Computes Hartree potential and orbital Coulomb integral in parallel, using 
  ! orbitals defined on a NN grid.
  ! Computes Orbital Coulomb Integral
  ! Uses orbital integral to construct the lower triangle of 2-particle 
  ! Hamiltonian, in packed form, distributed across numproc processes
  !*****************************************************************************

  subroutine construct_H2p(H2p,rstart,rstop)
    use var,   only: sp,dp,qp,Natoms,orb_bas,evcr,atom_species,d,Nexsp,eigen_cnt,NN_list,NumNN,&
                     r,time_CoulInt,al,NNdp_list,atom_type,pi,all_pnts,eigen_cnt,Ne,Nh,ic,static_eps
    use grids, only: gen_NN_grid, gen_NN_rspacing, return_grid_xyz,return_grid_xyz_arb
    use smat,  only: NN_disp
    use sort,  only: reallocate
    implicit none

    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer       N,Nm,iro,ir,jro,jr,  ia,ja,ka,la,  n1,n2,n3,n4,n12,n34 , iorb,jorb,korb,lorb,&
                  n1_l,n2_l,n3_l,n4_l,Natm,ia_l,ja_l,ka_l,la_l,n1234,NOI,ix1,iy1,iz1,ir_l,&
                  ix_l,iy_l,iz_l,Ntotal,ix,iy,iz,ixy_l, jx1,jy1,jz1,jr_l,jxy_l,jx_l,jy_l,jz_l,&
                  jx,jy,jz,istart,istop,Nreq,Ntotal_lcl,jrcnt,i,n34_start,n34_stop,cnt,VBT,CBB
    integer       row,col,Nh_mx,Ne_mx,h_row_cnt,e_row_cnt,h_col_cnt,e_col_cnt,h_row,e_row,&
                  h_col,e_col,Hcnt,orb_cnt,j
    real(dp)      dr,integral,x,y,z,denom,r_sph,orbij,inv_eps, c1,c3,c6,V_1234 
    real(dp)      t1_local,t2_local,t21,t21_lg, VH_time,VH_time_lg,VH_time_local,VH_t1,VH_t2, &
                  worb_time, worb_time_lg, worb_time_local,worb_t1,worb_t2, &
                  H2p_time, H2p_time_lg, H2p_time_local,H2p_t1, H2p_t2
    real(qp)      unt
    character     fname*12,rank_lab*3
     

    !Modular Variables
    integer,      intent(in)          :: rstart,rstop
    integer       jr_start,jr_stop

    !Modular arrays
    integer,      allocatable         :: orb_grd_pntr(:,:,:),nnz(:,:)
    real(dp),     allocatable         :: orb_overlap(:,:,:)
    real(dp),     allocatable         :: worb(:),worb_local(:)
    real(dp),     allocatable,        intent(inout) :: H2p(:)
    !complex(dp), allocatable,        intent(inout) :: H2p(:)

    !Local Arrays
    integer,      allocatable         :: tmp(:)
    integer(sp),  allocatable         :: displ(:),rec_size(:),displ2(:),rec2(:)
    real(dp),     allocatable         :: V_H(:),V_H_local(:)
    real(dp),     allocatable         :: r_i(:),r_j(:),c_n1_h(:),c_n2_e(:),c_n2_h(:), &
                                         c_n3_h(:),c_n3_e(:),c_n4_e(:),evc_e(:,:),evc_h(:,:)
    real(dp),     allocatable         :: H2p_tg(:),VH_tg(:),H2p_tl(:),VH_tl(:)




 
    !**********************************
    !Main Routine
    !**********************************
    if(IO) write(*,*)'Construct_H2p'
    call cpu_time(t1_local)

    !Perturbative Test
    !if(IO) write(*,*) 'Only outputting K=-W from construct_H_2p'

    !Dimensions of the system
    N=Natoms*orb_bas

    !Maximum number of orbital interactions for two centres, in NN regime
    NOI=N*(NumNN*orb_bas)
    if(IO) write(*,*)'Max number of orbital interactions',NOI
   
    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1
    
    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)

    !All arguments here are defined in the module
    !This NN grid uses the same spacing as the system grid 
    !Important for mapping between the system and NN grids
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)

    !Arrays not utilised so deallocate immediately
    deallocate(NNgrid,NNxyz_ptr)

    !Construct orb_overlap matrix on NN grid
    call construct_orb_overlap(orb_overlap,orb_grd_pntr,nnz)
    !Index array for NN grid not utilised after orbital products are constructed
    deallocate(NNiptr)

    !Generate grid arrays with dimensions of the supercell
    !No auxilliary grid arrays are allocated, except atm_pntr
    all_pnts=.true.
    !Note, use return_grid_xyz with setup_system_grid
    call setup_system_grid()

    !For testing: Construct system grid with vacuum region included
    !call setup_arb_system_grid()
    !if(IO) write(*,*) 'Note, use return_grid_xyz_arb with setup_arb_system_grid'
    !deallocate(sgrid,siptr,sys_pntr)

    !Total number of points in system grid
    Ntotal=product(snpts)
    if(IO) write(*,*) 'Total number of points in Hartree potential',Ntotal

    !Volume element associated with this grid spacing
    dr=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Radius of cubic volume element, given as the length from 
    !the centre to the corner (see Brillouin approximation)
    r_sph=sqrt(3._dp)*0.5*NNrspacing(1,1)

    !Position vectors for Coulomb Potential
    allocate(r_i(d),r_j(d))

    !Index arrays for the two distributed loops
    allocate(rec_size(numprocs), displ(numprocs))
    allocate(rec2(numprocs), displ2(numprocs))

    !Units. e**2/(4 pi eps_0 a_0) = 1 Ha
    !Convert atomic positions to Bohr
    !Convert Hartree to eV  
    unt=(0.5291772109217*27.2113850560)

    !1/static dielectric constant 
    if(IO)write(*,*)'Static dielectric constant read from file:',static_eps
    inv_eps=(1._dp/static_eps)


    !******************************************************************************
    ! 1)   Compute Hartree Potential: V_H(r')=  integral{ orb(n12,r)/|r-r'| } dr
    !      per orbital product n12
    ! 2)   Compute integral{ orb(n34,r') V_H(n12,r')} dr'  per orbital product n12
    ! 3)   Compute elements of packed, distributed lower triangle of resonant, 
    !      2-particle matrix H2p
    !******************************************************************************

    !Distribute spatial index jr over supercell grid 
    !rec_size = size(V_H_local) for each process
    !disp = Index of V_H_local w.r.t. V_H
    call distribute_loop2(1,Ntotal, jr_start,jr_stop,rec_size,displ)

    !Number of grid points per CPU
    Ntotal_lcl=jr_stop-jr_start+1

    !Allocate global and local Hartree potentials
    allocate(V_H_local(Ntotal_lcl),V_H(Ntotal))
    V_H(:)=0._dp
    V_H_local(:)=0._dp

    !Distribute over worb n34=[1,NOI]. 
    !This has been moved out of the loop structure as the limits
    !are now fixed (couldn't get db_factor to work)
    call distribute_loop2(1,NOI, n34_start,n34_stop, rec2,displ2)
    allocate(worb(NOI),worb_local(n34_stop-n34_start+1))

    !1D arrays for holding tight-binding coefficients w.r.t 
    !orbital indices n1,n2,n3 and n4
    !allocate(c_n1_h(Nh),c_n2_e(Ne),c_n2_h(Nh), c_n3_h(Nh),c_n3_e(Ne),c_n4_e(Ne))
    
    !Reorder eigenvector arrays so that their inner index
    !corresponds to the inner-most loops in construction of H2p
    !!!!call eigenvector_order(evcr, evc_e,evc_h)



    !**********************************************************************
    ! Section 1
    !
    !**********************************************************************
    !Timings
    !Timing V_H construction
    VH_time_local=0._dp
    !Timing worb construction
    worb_time_local=0._dp
    !Timing H2p construction
    H2p_time_local=0._dp
   


    !First composite index over orbitals 1&2 defined with spatial coordinate r (ir)
    do n12=1,NOI

       !Time worb construction 
       call cpu_time(worb_t1)


       !Return indices in composite index n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)
       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle

       !Global orbital indices n1 and n2 (used in eigenvectors below)
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)
       
       !Map local NN grid spatial indices to system grid indices using atom ia 
       !on system grid and corresponding atom ia_l on NN grid
       call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,ia,ia_l,&
                                      ix1,iy1,iz1)



       !Loop over all points in system grid. Distributed
       call cpu_time(VH_t1)
       jrcnt=0
       V_H_local(:)=0._dp
     
       do jr=jr_start,jr_stop

          !Index for V_H_local
          jrcnt=jrcnt+1

          !Unwrap composite index for system grid: jr->(jx,jy,jz)
          call comp_to_ijk(jr,snpts(1),snpts(2),snpts(3),jx,jy,jz)

          !Return x,y,z variables associated with spatial indices
          call return_grid_xyz(jx,jy,jz,sspacing,snpts, x,y,z)
          r_j(:)=(/x,y,z/)


          !*********************************************
          !Integral over r
          !*********************************************
          !Initialise integral 
          integral=0._dp
          !Initialise orbital product
          orbij=0._dp

          !Loop over non-zero orbital products defined on NN grid
          do iro=1,nnz(n1_l,n2_l)

             !Map non-zero product space to NN grid space
             ir_l=orb_grd_pntr(n1_l,n2_l,iro)

             !Unwrap composite index for NN grid space: ir_l->(ix_l,iy_l,iz_l)
             call comp_to_ijk(ir_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                              ix_l,iy_l,iz_l)

             !Iterate system grid indices using (ix_l,iy_l,iz_l)
             ix=ix_l+(ix1-1)
             iy=iy_l+(iy1-1)
             iz=iz_l+(iz1-1)

             !Return x,y,z variables associated with spatial indices
             call return_grid_xyz(ix,iy,iz,sspacing,snpts, x,y,z)
             r_i(:)=(/x,y,z/)

             !Composite spatial counter for system grid
             ir=iz+( (iy+((ix-1)*snpts(2))-1)*snpts(3))

             !If |r-r'|=0
             if(ir==jr)then
                !Store 'product of orbitals' value
                orbij=orb_overlap(n1_l,n2_l,iro)
                !write(100,'(2(I6,X),7(F12.6,X))') ir,jr,r_i(:),r_j(:),orbij
                cycle
             endif

             !Sum Hartree Potential
             integral=integral+ (orb_overlap(n1_l,n2_l,iro)/&
                                 sqrt(dot_product(r_i(:)-r_j(:),r_i(:)-r_j(:))))
                                 
                  
          enddo
          !Shut integral over ir
          !*********************************************

          !Store local Hartree potential at spatial point r' (jr), orbital pair index n12
          ! + Brillouin Approximation for when r=r'
          V_H_local(jrcnt)=(integral*dr)+(orbij*2._dp*pi*r_sph**2)

       end do
       !Shut loop over jr

       
       !Gather the local Hartree potential per process into V_H, available to 
       !all process
       call mpi_allgatherv(V_H_local,Ntotal_lcl,MPI_DOUBLE_PRECISION,     & !Local send
                           V_H, rec_size, displ,MPI_DOUBLE_PRECISION,     & !Global receive
                           MPI_COMM_WORLD,MPIierr)

       !Deallocate local Hartree potentials to temporarily save memory
       !deallocate(V_H_local)

       !Time loop over Hartree potential, for all n12 
       call cpu_time(VH_t2)
       VH_time_local=VH_time_local+(VH_t2-VH_t1)



       !**********************************************************************
       ! Section 2
       !
       !**********************************************************************
       !Distribute loop over 2nd composite index
!!$       call distribute_loop2(n12,NOI, n34_start,n34_stop, rec2,displ2)
!!$       !Allocate local and global orbital integral arrays (per n12)
!!$       !Both must be defined here, as they're size is a function fo n12
!!$       allocate(worb(NOI-n12+1),worb_local(n34_stop-n34_start+1))
!!$       worb_local(:)=0._dp
!!$       worb(:)=0._dp   


       !Initialise orbital integrals
       worb_local(:)=0._dp
       worb(:)=0._dp

       !Initialise index for worb_local 
       cnt=0


       !Loop over 2nd composite index n34. Distributed 
       do n34=n34_start,n34_stop

          !Iterate w_local index
          !Must be defined before cycle, such that elements of worb_local
          !are consistent with the parallelisation (wastes some memory as a result)
          cnt=cnt+1

          !Return indices in composite index n34 
          call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
          !If ia doesn't have a neighbour la, cycle loop
          if(la==0)cycle

          !Map local NN grid spatial indices to system grid indices using atom ka 
          !on system grid and corresponding atom ka_l on NN grid
          call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,&
                                         ka,ka_l,jx1,jy1,jz1)


          !*****************************************
          !Integral over r'
          !*****************************************
          integral=0._dp

          !Compute integral over spatial variable r'
          do jro=1,nnz(n3_l,n4_l)

             !Map non-zero product space to NN grid space
             jr_l=orb_grd_pntr(n3_l,n4_l,jro)

             !Unwrap composite index for NN grid space
             call comp_to_ijk(jr_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                              jx_l,jy_l,jz_l)

             !Iterate system grid indices using (jx_l,jy_l,jz_l)
             jx=jx_l+(jx1-1)
             jy=jy_l+(jy1-1)
             jz=jz_l+(jz1-1)

             !Composite spatial counter for system grid
             jr=jz+( (jy+((jx-1)*snpts(2))-1)*snpts(3))

             !Integral
             integral=integral+ (orb_overlap(n3_l,n4_l,jro)*V_H(jr))

          enddo
          !*****************************************
      
          !Store orbital integral (mathematically dr' but dr=dr')
          worb_local(cnt)=integral*dr
          
       enddo
       !Shut loop over n34

       !Gather the local orbital integral  per process into worb
       call mpi_allgatherv(worb_local, n34_stop-n34_start+1 ,MPI_DOUBLE_PRECISION,     & !Local send
                           worb,       rec2, displ2         ,MPI_DOUBLE_PRECISION,     & !Global receive
                           MPI_COMM_WORLD,MPIierr)

       !Free local worb
       !deallocate(worb_local)

       !Allocate local Hartree potential, ready for the next iteration of n12
       !allocate(V_H_local(Ntotal_lcl))
       
       !Convert worb to eV
       worb(:)=worb(:)*unt

       !Time construction of worb, for all n12 
       call cpu_time(worb_t2)
       worb_time_local=worb_time_local+(worb_t2-worb_t1)
       




       !**********************************************************************
       ! Section 3
       ! Construct lower triangle of distributed resonant H2p (packed)
       !**********************************************************************
       call cpu_time(H2p_t1)

  
       !Loop over 1/8th of orbital products, per index n12
       !Retrieve all elements with 8 sets of coefficients
!!$       orb_cnt=0
!!$       do n34=n12,NOI
!!$          !Index for worb
!!$          !orb_cnt=orb_cnt+1


       do n34=1,NOI      
          !Return indices in composite index n34 
          call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
          !If ia doesn't have a neighbour la, cycle loop
          if(la==0)cycle
          !Define n3 and n4 orbital indices over all orbital states in the system
          n3=korb+((ka-1)*orb_bas)
          n4=lorb+((la-1)*orb_bas)

          !Assign (n1,n2|V|n3,n4) to variable, removing array access in loops below
          V_1234=worb(n34)


          !*********************************************************************************
          !Counter for packed matrix H2p
          Hcnt=0
          
          !Row loop distributed => H2p is distributed
          do row=rstart,rstop
                  
             !Counters for electron and hole row indices
             h_row_cnt=((row-1)/Ne)+1
             e_row_cnt=row-((h_row_cnt-1)*Ne)
             
             !Hole and electron indices for row (h,e)
             h_row=(VBT+1) - h_row_cnt
             e_row= e_row_cnt + (CBB-1)
             
             !Pass 2D arrays to variables (variable indices indicate order 
             !that the coefficients appear in the kernel)
             c1=evcr(n1,h_row)
             c3=evcr(n2,e_row)
             c6=evcr(n3,e_row)

           
             !Loop over lower triangle only
             do col=1,row
                
                !Iterate H2p packed index counter
                Hcnt=Hcnt+1
                
                !Counters for electron and hole col indices
                h_col_cnt=((col-1)/Ne)+1
                e_col_cnt=col-((h_col_cnt-1)*Ne)
                
                !Hole and electron indices for col (h'e')
                h_col=(VBT+1) - h_col_cnt
                e_col= e_col_cnt + (CBB-1)
                
                !Kernel in orbital basis. 2v-W = 
                ! c_n1_h c_n4_e' *(2 c_n2_e c_n3_h' - eps c_n2_h' c_n3_e) (12|v|34)
                !Coefficients common to both direct and exchange have been factorised        
                H2p(Hcnt) = H2p(Hcnt) +  c1*evcr(n4,e_col)*V_1234* &
                                      (2._dp*c3*evcr(n3,h_col)-inv_eps*evcr(n2,h_col)*c6)
    
        

             enddo
             !*********************************************************************************
          enddo
          !End loops over rows and columns of H2p
          !**********************************************************************
       enddo
       !Ends loop over moved segment (loop over n34)
       !deallocate(worb)

       !Time construction of H2p, for all n12 
       call cpu_time(H2p_t2)
       H2p_time_local=H2p_time_local+(H2p_t2-H2p_t1)

    enddo
    !Shut loop over n12 

    !Perturbative Test
    !if(IO) write(*,*) 'Only outputting K=-W from construct_H_2p'

    !Deallocate Hartree potentials
    deallocate(V_H_local,V_H)

    !Deallocate orbital integrals
    if(allocated(worb)) deallocate(worb)
    if(allocated(worb_local)) deallocate(worb_local)




    !***************************************************************
    !Routine Timing
    !Important timings are VH and H2p constructions
    !***************************************************************
    !Finish timing subroutine
    call cpu_time(t2_local)

    !Create arrays with only the local run time (_tl) and arrays
    !to contain run times per process (_tg)
    allocate(H2p_tl(numprocs),VH_tl(numprocs),H2p_tg(numprocs),VH_tg(numprocs))
    H2p_tg(:)=0._dp
    VH_tg(:)=0._dp
    H2p_tl(:)=0._dp
    VH_tl(:)=0._dp

    !Fill local run time into element corresponding to process rank
    H2p_tl(rank+1)=H2p_time_local
    VH_tl(rank+1)=VH_time_local

    !Sum local array times such that each process has a global array that contains run times
    !for all processes
    call mpi_allreduce(H2p_tl,H2p_tg,numprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce(VH_tl,VH_tg,numprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)

    !Sum time taken for each section of routine (time in seconds by default)
    !Redundant. Can be done from summing returned arrays above
    call mpi_allreduce(VH_time_local,VH_time,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce(worb_time_local,worb_time,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce(H2p_time_local,H2p_time,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce((t2_local-t1_local),t21,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)


    !Output average times and return in mins
    if(IO)then
       !VH
       write(*,'(A,X,F20.10,A)') 'Avg time taken to construct V_H =',VH_time/(dble(numprocs)*60._dp),' mins'
       write(*,'(A,X,F20.10,A,X,I3)') 'Longest time taken to construct V_H =',maxval(VH_tg)/60._dp,&
                                      ' mins on process',maxloc(VH_tg)-1
       write(*,'(A,X,F20.10,A)')   'The standard deviation for VH construction time =',stnd_dv(VH_tg)/60._dp,'mins'

       !H2p
       write(*,'(A,X,F20.10,A)') 'Avg time taken to construct H2p =',H2p_time/(dble(numprocs)*60._dp),' mins'
       write(*,'(A,X,F20.10,A,X,I3)') 'Longest time taken to construct H2p =',maxval(H2p_tg)/60._dp,&
                                      ' mins on process',maxloc(H2p_tg)-1
       write(*,'(A,X,F20.10,A)')   'The standard deviation for H2p construction time =',stnd_dv(H2p_tg)/60._dp,'mins'
    endif


    !More Global timings
    !write(*,'(A,X,F20.10,A)') 'Avg time taken to construct worb =',worb_time/(dble(numprocs)*60._dp),' mins'
    !write(*,'(A,X,F20.10,A)') 'Avg time taken to run routine "construct_H2p" =',t21/(dble(numprocs)*60._dp),' mins'

    !Local Timings
    !write(*,'(A,I3,X,A,X,F20.10,A)') 'Time taken to construct V_H on process',rank,'=',VH_time_local/60._dp,' mins'
    !write(*,'(A,I3,X,A,X,F20.10,A)') 'Time taken to construct worb on process',rank,'=',worb_time_local/60._dp,' mins'
    !write(*,'(A,I3,X,A,X,F20.10,A)') 'Time taken to construct H2p on process',rank,'=',H2p_time_local/60._dp,' mins'
    !write(*,'(A,I3,X,A,X,F20.10,A)') 'Time taken to run construct_H2p on process',rank,'=',(t2_local-t1_local)/60._dp,' mins'


    !End of the routine. 
  end subroutine construct_H2p



  !********************************************************************
  !Function to compute standard deviation of data set x of size(x)
  !********************************************************************
  real(dp) function stnd_dv(x)
    implicit none

    !**************
    !Declarations
    !**************
    !Input data array
    real(dp),  allocatable, intent(in) :: x(:)

    !Local variables
    integer   i,N
    real(dp)  mn,rs

    !**************
    !Main Routine
    !**************
    !Size of data set
    N=size(x)

    !Mean
    mn=(1._dp/dble(N))*sum(x)

    !Residuals squared
    rs=0._dp
    do i=1,N
       rs=rs+(x(i)-mn)**2.
    enddo
    
    !Standard Deviation
    stnd_dv=sqrt( (1._dp/dble(N))*rs )

    return
  end function stnd_dv




  !*****************************************************************************
  ! Construct H2p. Version 2
  !
  ! Computes Hartree potential and orbital Coulomb integral in parallel, using 
  ! orbitals defined on a NN grid.
  ! Computes Orbital Coulomb Integral
  ! Uses orbital integral to construct the lower triangle of 2-particle 
  ! Hamiltonian, in packed form, distributed across numproc processes
  !*****************************************************************************

  subroutine construct_H2p_v2(H2p,rstart,rstop)
    use var,   only: sp,dp,qp,Natoms,orb_bas,evcr,atom_species,d,Nexsp,eigen_cnt,NN_list,NumNN,&
                     r,time_CoulInt,al,NNdp_list,atom_type,pi,all_pnts,Ne,Nh
    use grids, only: gen_NN_grid, gen_NN_rspacing, return_grid_xyz,return_grid_xyz_arb
    use smat,  only: NN_disp
    use sort,  only: reallocate
    implicit none

    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer       N,Nm,iro,ir,jro,jr,  ia,ja,ka,la,  n1,n2,n3,n4,n12,n34 , iorb,jorb,korb,lorb,&
                  n1_l,n2_l,n3_l,n4_l,Natm,ia_l,ja_l,ka_l,la_l,n1234,NOI,ix1,iy1,iz1,ir_l,&
                  ix_l,iy_l,iz_l,Ntotal,ix,iy,iz,ixy_l, jx1,jy1,jz1,jr_l,jxy_l,jx_l,jy_l,jz_l,&
                  jx,jy,jz,istart,istop,Nreq,Ntotal_lcl,jrcnt,i,n34_start,n34_stop,VBT,CBB
    integer       row,col,Nh_mx,Ne_mx,h_row_cnt,e_row_cnt,h_col_cnt,e_col_cnt,h_row,e_row,&
                  h_col,e_col,Hcnt,orb_cnt,Norb
    real(dp)      dr,integral,x,y,z,denom,r_sph,orbij,t1,t2,VH_time_local,VH_t1,VH_t2,worb_t1,worb_t2,&
                  H2p_t1,H2p_t2
    real(qp)      unt
    character     fname*12,rank_lab*3



    !Modular Variables
    integer,      intent(in)          :: rstart,rstop
    integer       jr_start,jr_stop

    !Modular arrays
    integer,      allocatable         :: orb_grd_pntr(:,:,:),nnz(:,:)
    real(dp),     allocatable         :: orb_overlap(:,:,:)
    real(dp),     allocatable         :: worb(:),worb_local(:)
    real(dp),     allocatable,        intent(inout) :: H2p(:)

    !Local Arrays
    integer,      allocatable         :: tmp(:)
    integer(sp),  allocatable         :: displ(:),rec_size(:),displ2(:),rec2(:)
    real(dp),     allocatable         :: V_H(:),V_H_local(:)
    real(dp),     allocatable         :: r_i(:),r_j(:)


 
    !**********************************
    !Main Routine
    !**********************************
    if(IO) write(*,*)'Construct_H2p_v2'
    call cpu_time(t1)

    !Dimensions of the system
    N=Natoms*orb_bas

    !Maximum number of orbital interactions for two centres, in NN regime
    NOI=N*(NumNN*orb_bas)
    if(IO) write(*,*)'Max number of orbital interactions',NOI
   
    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1
    
    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)

    !All arguments here are defined in the module
    !This NN grid uses the same spacing as the system grid 
    !Important for mapping between the system and NN grids
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)

    !Arrays not utilised so deallocate immediately
    deallocate(NNgrid,NNxyz_ptr)

    !Construct orb_overlap matrix on NN grid
    call construct_orb_overlap(orb_overlap,orb_grd_pntr,nnz)
    !Index array for NN grid not utilised after orbital products are constructed
    deallocate(NNiptr)

    !Generate grid arrays with dimensions of the supercell
    !No auxilliary grid arrays are allocated, except atm_pntr
    all_pnts=.true.
    call setup_system_grid()

    !For testing: Construct system grid with vacuum region included
    !call setup_arb_system_grid()
    !if(IO) write(*,*) 'Note, use return_grid_xyz_arb with setup_arb_system_grid'
    !deallocate(sgrid,siptr,sys_pntr)

    !Total number of points in system grid
    Ntotal=product(snpts)
    if(IO) write(*,*) 'Total number of points in Hartree potential',Ntotal

    !Volume element associated with this grid spacing
    dr=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Radius of cubic volume element, given as the length from 
    !the centre to the corner (see Brillouin approximation)
    r_sph=sqrt(3._dp)*0.5*NNrspacing(1,1)

    !Position vectors for Coulomb Potential
    allocate(r_i(d),r_j(d))

    !Index arrays for the two distributed loops
    allocate(rec_size(numprocs), displ(numprocs))
    !allocate(rec2(numprocs), displ2(numprocs))

    !Units. e**2/(4 pi eps_0 a_0) = 1 Ha
    !Convert atomic positions to Bohr
    !Convert Hartree to eV  
    unt=(0.5291772109217*27.2113850560)




    !******************************************************************************
    ! 1)   Compute Hartree Potential: V_H(r')=  integral{ orb(n12,r)/|r-r'| } dr
    !      per orbital product n12
    ! 2)   Compute integral{ orb(n34,r') V_H(n12,r')} dr'  per orbital product n12
    !******************************************************************************

    !Distribute spatial index jr over supercell grid 
    !rec_size = size(V_H_local) for each process
    !disp = Index of V_H_local w.r.t. V_H
    call distribute_loop2(1,Ntotal, jr_start,jr_stop,rec_size,displ)

    !Number of grid points per CPU
    Ntotal_lcl=jr_stop-jr_start+1

    !Allocate global and local Hartree potentials
    allocate(V_H_local(Ntotal_lcl),V_H(Ntotal))
    V_H(:)=0._dp
    V_H_local(:)=0._dp

    !Orbital Integral allocated with upper bound
    !Norb=int(0.5*dble(NOI)**2)
    !Set to max upper bound, as was getting seg fault for the above
    Norb=NOI*NOI
    allocate(worb(Norb))
    worb(:)=0._dp
    !Open file to write worb to
    !if(IO) open(unit=100,file='worb.dat')





    !**********************************************************************
    ! Loop Structure to construct orbital integral worb=(n1n2|n3n4)
    !
    !**********************************************************************
    if(IO) write(*,*) 'Started orbital construction'
    call cpu_time(worb_t1)

    !Initialise composite orbital index
    n1234=0
    !Timing V_H construction
    VH_time_local=0._dp


    !First composite index over orbitals 1&2 defined with spatial coordinate r (ir)
    do n12=1,NOI

       !Return indices in composite index n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)
       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle
       
       !Map local NN grid spatial indices to system grid indices using atom ia 
       !on system grid and corresponding atom ia_l on NN grid
       call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,ia,ia_l,&
                                      ix1,iy1,iz1)

       
       !Loop over all points in system grid. Distributed
       call cpu_time(VH_t1)
       jrcnt=0
       V_H_local(:)=0._dp

       do jr=jr_start,jr_stop

          !Index for V_H_local
          jrcnt=jrcnt+1

          !Unwrap composite index for system grid: jr->(jx,jy,jz)
          call comp_to_ijk(jr,snpts(1),snpts(2),snpts(3),jx,jy,jz)

          !Return x,y,z variables associated with spatial indices
          call return_grid_xyz(jx,jy,jz,sspacing,snpts, x,y,z)
          r_j(:)=(/x,y,z/)


          !*********************************************
          !Integral over r
          !*********************************************
          !Initialise integral 
          integral=0._dp
          !Initialise orbital product
          orbij=0._dp

          !Loop over non-zero orbital products defined on NN grid
          do iro=1,nnz(n1_l,n2_l)

             !Map non-zero product space to NN grid space
             ir_l=orb_grd_pntr(n1_l,n2_l,iro)

             !Unwrap composite index for NN grid space: ir_l->(ix_l,iy_l,iz_l)
             call comp_to_ijk(ir_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                              ix_l,iy_l,iz_l)

             !Iterate system grid indices using (ix_l,iy_l,iz_l)
             ix=ix_l+(ix1-1)
             iy=iy_l+(iy1-1)
             iz=iz_l+(iz1-1)

             !Return x,y,z variables associated with spatial indices
             call return_grid_xyz(ix,iy,iz,sspacing,snpts, x,y,z)
             r_i(:)=(/x,y,z/)

             !Composite spatial counter for system grid
             ir=iz+( (iy+((ix-1)*snpts(2))-1)*snpts(3))

             !If |r-r'|=0
             if(ir==jr)then
                !Store 'product of orbitals' value
                orbij=orb_overlap(n1_l,n2_l,iro)
                cycle
             endif

             !Sum Hartree Potential
             integral=integral+ (orb_overlap(n1_l,n2_l,iro)/&
                                 sqrt(dot_product(r_i(:)-r_j(:),r_i(:)-r_j(:))))
                  
          enddo
          !Shut integral over ir
          !*********************************************

          !Store local Hartree potential at spatial point r' (jr), orbital pair index n12
          ! + Brillouin Approximation for when r=r'
          V_H_local(jrcnt)=(integral*dr)+(orbij*2._dp*pi*r_sph**2)
          
       end do
       !Shut loop over jr

       !Gather the local Hartree potential per process into V_H, available to 
       !all process
       call mpi_allgatherv(V_H_local,Ntotal_lcl,MPI_DOUBLE_PRECISION,     & !Local send
                           V_H, rec_size, displ,MPI_DOUBLE_PRECISION,     & !Global receive
                           MPI_COMM_WORLD,MPIierr)

       !Deallocate local Hartree potentials to temporarily save memory
       !Stopped doing these to see if it affects speed
       !deallocate(V_H_local)

       !Time loop over Hartree potential, for all n12 
       call cpu_time(VH_t2)
       VH_time_local=VH_time_local+(VH_t2-VH_t1)




       !**********************************************************************
       ! Section 2. USe Hartree potential to construct 
       ! (n1n2|n3n4)= int{ orbs(n34,r') VH(n12,r')} dr'
       !
       !**********************************************************************   
       
       !Loop over 2nd composite index n34. Serial
       do n34=1,NOI

          !Return indices in composite index n34 
          call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
          !If ia doesn't have a neighbour la, cycle loop
          if(la==0)cycle

          !Map local NN grid spatial indices to system grid indices using atom ka 
          !on system grid and corresponding atom ka_l on NN grid
          call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,&
                                         ka,ka_l,jx1,jy1,jz1)

          !Iterate worb index=(n12,n34)
          n1234=n1234+1

          !*****************************************
          !Integral over r'
          !*****************************************
          integral=0._dp

          !Compute integral over spatial variable r'
          do jro=1,nnz(n3_l,n4_l)

             !Map non-zero product space to NN grid space
             jr_l=orb_grd_pntr(n3_l,n4_l,jro)

             !Unwrap composite index for NN grid space
             call comp_to_ijk(jr_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                  jx_l,jy_l,jz_l)

             !Iterate system grid indices using (jx_l,jy_l,jz_l)
             jx=jx_l+(jx1-1)
             jy=jy_l+(jy1-1)
             jz=jz_l+(jz1-1)

             !Composite spatial counter for system grid
             jr=jz+( (jy+((jx-1)*snpts(2))-1)*snpts(3))

             !Integral
             integral=integral+ (orb_overlap(n3_l,n4_l,jro)*V_H(jr))

          enddo
          !*****************************************
      
          !Store orbital integral 
          worb(n1234)=integral*dr

          !if(IO) write(100,*) worb(n1234)
          !if(IO) write(100,'(2(I3,X),I7,X,F20.10)') n12,n34,n1234,worb(n1234)
       enddo
       !Shut loop over n34

       !Allocate local Hartree potential, ready for the next iteration of n12
       !allocate(V_H_local(Ntotal_lcl))
       
    enddo
    !Shut loop over n12 
    !if(IO) close(100)


    !****************************
    !Tidy up memory
    !****************************
    !COULD REALLOCATE WORB WITH CORRECT SIZE HERE

    !Deallocate Hartree potentials
    deallocate(V_H_local,V_H)

    !Convert worb into eV
    worb(:)=unt*worb(:)

    !Inform user of averga time to 1) Compute V_H.  2)Compute worb
    call cpu_time(worb_t2)
    write(*,'(A,I3,X,A,X,F20.10,A)') 'Time taken to construct V_H on process',rank,'=',VH_time_local/60._dp,' mins'
    write(*,'(A,I3,X,A,X,F20.10,A)') 'Time taken to construct worb on process',rank,'=',(worb_t2-worb_t1)/60._dp,' mins'

 



    !***************************************************************
    !Section 3
    !Construct H_2p in same routine
    !***************************************************************
    call cpu_time(H2p_t1)
    if(IO) write(*,*)'Constructing H2p using v2'

    !Initialise packed Hamiltonian index
    Hcnt=0

    !Initalise H2p
    H2p(:)=0._dp

    !Row loop distributed => H2p is distributed
    do row=rstart,rstop
               
       !Counters for electron and hole row indices
       h_row_cnt=((row-1)/Ne)+1
       e_row_cnt=row-((h_row_cnt-1)*Ne)

       !Hole and electron indices for row (h,e)
       h_row=(VBT+1) - h_row_cnt
       e_row= e_row_cnt + (CBB-1)

       
       
       !Loop over lower triangle only
       do col=1,row
          
          !Iterate H2p packed index counter
          Hcnt=Hcnt+1
          
          !Counters for electron and hole col indices
          h_col_cnt=((col-1)/Ne)+1
          e_col_cnt=col-((h_col_cnt-1)*Ne)
          
          !Hole and electron indices for col (h'e')
          h_col=(VBT+1) - h_col_cnt
          e_col= e_col_cnt + (CBB-1)
          
       
          !Loop over orbital composite indices
          n1234=0
          do n12=1,NOI
             !Return indices in composite index n12
             call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)
             !If ia doesn't have a neighbour ja, cycle loop
             if(ja==0)cycle
             !Orbital indices n1 and n2 
             n1=iorb+((ia-1)*orb_bas)
             n2=jorb+((ja-1)*orb_bas)
             
             do n34=1,NOI
                !Return indices in composite index n34 
                call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
                !If ia doesn't have a neighbour la, cycle loop
                if(la==0)cycle
                !Define n3 and n4 orbital indices over all orbital states in the system
                n3=korb+((ka-1)*orb_bas)
                n4=lorb+((la-1)*orb_bas)

                !Composite index for orbital integral - (n12,n34)
                n1234=n1234+1
                
                !Sum element H_(he),(h'e') for matrix A
                !Real eigenvectors
!!$                H2p(Hcnt) = H2p(Hcnt) + &
!!$                   ( (evcr(n1,e_row))*evcr(n2,e_col)*evcr(n3,h_row)*(evcr(n4,h_col)) + &
!!$                     (evcr(n2,e_row))*evcr(n1,e_col)*evcr(n3,h_row)*(evcr(n4,h_col)) + &
!!$                     (evcr(n1,e_row))*evcr(n2,e_col)*evcr(n4,h_row)*(evcr(n3,h_col)) + &
!!$                     (evcr(n2,e_row))*evcr(n1,e_col)*evcr(n4,h_row)*(evcr(n3,h_col)) + &
!!$                     (evcr(n3,e_row))*evcr(n4,e_col)*evcr(n1,h_row)*(evcr(n2,h_col)) + &
!!$                     (evcr(n4,e_row))*evcr(n3,e_col)*evcr(n1,h_row)*(evcr(n2,h_col)) + &
!!$                     (evcr(n3,e_row))*evcr(n4,e_col)*evcr(n2,h_row)*(evcr(n1,h_col)) + &
!!$                     (evcr(n4,e_row))*evcr(n3,e_col)*evcr(n2,h_row)*(evcr(n1,h_col)) ) &
!!$                     *dbfactor(n1,n2,n3,n4)*worb(n1234)

                !Expression for setting n34=1,NOI in two loops and removing dp_factor
                !Kernel given by K=2v-W
                H2p(Hcnt) = H2p(Hcnt) + &
                           ( -1._dp*(evcr(n1,e_row))*evcr(n2,e_col)*evcr(n3,h_row)*evcr(n4,h_col) +&
                              2._dp*(evcr(n1,e_row))*evcr(n2,h_row)*evcr(n3,e_col)*evcr(n4,h_col))*worb(n1234)

             enddo
          enddo
          !******************************
          !Return element 'Hcnt' of H2p
          !******************************
   
       enddo
    enddo
    !Shut loops over rows and columns
    
    
    call cpu_time(H2p_t2)
    call cpu_time(t2)
    write(*,'(A,I3,X,A,X,F20.10,A)') 'Time taken to construct H2p on process',rank,'=',(H2p_t2-H2p_t1)/60._dp,' mins'
    write(*,'(A,I3,X,A,X,F20.10,A)') 'Time taken to run construct_H2p_v2 on process',rank,'=',(t2-t1)/60._dp,' mins'

       
  end subroutine construct_H2p_v2




  !*******************************************************************************************
  ! Construct H2p_onsite.
  !
  ! Computes Hartree potential and orbital Coulomb integral in parallel, using 
  ! orbitals defined on a NN grid.
  ! Computes Orbital Coulomb Integral
  ! Uses orbital integral to construct the lower triangle of 2-particle 
  ! Hamiltonian, in packed form, distributed across numproc processes
  ! 
  ! Approximation:
  ! Constructs two-particle resonant Hamiltonian but only uses on-site terms for the 4-point
  ! orbital integrals, reducing them to 2-point quantities.
  ! Most extreme TB approximation but should result in ~400 speed-up in Hamiltonian construction
  !*******************************************************************************************

  subroutine construct_H2p_onsite(H2p,rstart,rstop)
    use var,   only: sp,dp,qp,Natoms,orb_bas,evcr,atom_species,d,Nexsp,eigen_cnt,NN_list,NumNN,&
                     r,time_CoulInt,al,NNdp_list,atom_type,pi,all_pnts,eigen_cnt,Ne,Nh,ic,static_eps
    use grids, only: gen_NN_grid, gen_NN_rspacing, return_grid_xyz,return_grid_xyz_arb
    use smat,  only: NN_disp
    use sort,  only: reallocate
    implicit none

    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer       N,Nm,iro,ir,jro,jr,  ia,ja,ka,la,  n1,n2,n3,n4,n12,n34 , iorb,jorb,korb,lorb,&
                  n1_l,n2_l,ia_l,ja_l, ix1,iy1,iz1, ir_l,Natm,&
                  ix_l,iy_l,iz_l,Ntotal,ix,iy,iz,ixy_l, jx1,jy1,jz1,jr_l,jxy_l,jx_l,jy_l,jz_l,&
                  jx,jy,jz,istart,istop,Nreq,Ntotal_lcl,jrcnt,i,cnt,VBT,CBB
    integer       row,col,Nh_mx,Ne_mx,h_row_cnt,e_row_cnt,h_col_cnt,e_col_cnt,h_row,e_row,&
                  h_col,e_col,Hcnt,orb_cnt,j
    real(dp)      dr,integral,x,y,z,denom,r_sph,orbij,inv_eps 
    real(dp)      t1_local,t2_local,t21,t21_lg, VH_time,VH_time_lg,VH_time_local,VH_t1,VH_t2, &
                  worb_time, worb_time_lg, worb_time_local,worb_t1,worb_t2, &
                  H2p_time, H2p_time_lg, H2p_time_local,H2p_t1, H2p_t2
    real(qp)      unt
    character     fname*12,rank_lab*3
     

    !Modular Variables
    integer,      intent(in)          :: rstart,rstop
    integer       jr_start,jr_stop

    !Modular arrays
    integer,      allocatable         :: orb_grd_pntr(:,:),nnz(:)
    real(dp),     allocatable         :: orb_overlap(:,:)
    real(dp),     allocatable         :: worb(:)
    real(dp),     allocatable,        intent(inout) :: H2p(:)

    !Local Arrays
    integer,      allocatable         :: tmp(:)
    integer(sp),  allocatable         :: displ(:),rec_size(:),displ2(:),rec2(:)
    real(dp),     allocatable         :: V_H(:),V_H_local(:)
    real(dp),     allocatable         :: r_i(:),r_j(:),c_n1_e(:),c_n2_e(:), &
                                         c_n1_h(:),c_n2_h(:)
    real(dp),     allocatable         :: H2p_tg(:),VH_tg(:),H2p_tl(:),VH_tl(:)


 
    !**********************************
    !Main Routine
    !**********************************
    if(IO) write(*,*)'Construct_H2p'
    call cpu_time(t1_local)

    !Dimensions of the system
    N=Natoms*orb_bas

    !Maximum number of orbital interactions, when one only retains on-site 
    !terms
    if(IO) write(*,*)'Max number of orbital interactions',N
   
    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1
    
    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)

    !All arguments here are defined in the module
    !This NN grid uses the same spacing as the system grid 
    !Important for mapping between the system and NN grids
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)
    deallocate(NNgrid,NNxyz_ptr)

    !Construct orb_overlap matrix on NN grid
    call construct_onsite_orb_products(orb_overlap,orb_grd_pntr,nnz)
    deallocate(NNiptr)

    !Generate grid arrays with dimensions of the supercell
    !No auxilliary grid arrays are allocated, except atm_pntr
    all_pnts=.true.
    !Note, use return_grid_xyz with setup_system_grid
    call setup_system_grid()

    !Total number of points in system grid
    Ntotal=product(snpts)
    if(IO) write(*,*) 'Total number of points in Hartree potential',Ntotal

    !Volume element associated with this grid spacing
    dr=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Radius of cubic volume element, given as the length from 
    !the centre to the corner (see Brillouin approximation)
    r_sph=sqrt(3._dp)*0.5*NNrspacing(1,1)

    !Position vectors for Coulomb Potential
    allocate(r_i(d),r_j(d))

    !Index arrays for the two distributed loops
    allocate(rec_size(numprocs), displ(numprocs))
    allocate(rec2(numprocs), displ2(numprocs))

    !Units. e**2/(4 pi eps_0 a_0) = 1 Ha
    !Convert atomic positions to Bohr
    !Convert Hartree to eV  
    unt=(0.5291772109217*27.2113850560)

    !1/static dielectric constant 
    if(IO)write(*,*)'Static dielectric constant read from file:',static_eps
    inv_eps=(1._dp/static_eps)



    !******************************************************************************
    ! 1)   Compute Hartree Potential: V_H(r')=  integral{ orb(n12,r)/|r-r'| } dr
    !      per orbital product n12
    ! 2)   Compute integral{ orb(n34,r') V_H(n12,r')} dr'  per orbital product n12
    ! 3)   Compute elements of packed, distributed lower triangle of resonant, 
    !      2-particle matrix H2p
    !******************************************************************************

    !Distribute spatial index jr over supercell grid 
    !rec_size = size(V_H_local) for each process
    !disp = Index of V_H_local w.r.t. V_H
    !Note, this only works for 1D arrays
    call distribute_loop2(1,Ntotal, jr_start,jr_stop,rec_size,displ)

    !Number of grid points per CPU
    Ntotal_lcl=jr_stop-jr_start+1

    !Allocate global and local Hartree potentials
    allocate(V_H_local(Ntotal_lcl),V_H(Ntotal))
    V_H(:)=0._dp
    V_H_local(:)=0._dp

    !Allocate worb, for n2=[1,N].
    allocate(worb(N))  

    !1D temporary arrays for holding tight-binding coefficients w.r.t 
    !orbital indices n1,n2
    allocate(c_n1_e(Ne),c_n2_e(Ne),c_n1_h(Nh),c_n2_h(Nh))




    !**********************************************************************
    ! Section 1
    !
    !**********************************************************************
    !Timings
    !Timing V_H construction
    VH_time_local=0._dp
    !Timing worb construction
    worb_time_local=0._dp
    !Timing H2p construction
    H2p_time_local=0._dp
   


    !First index over orbitals n1==n2 defined with spatial coordinate r (ir)
    do n1=1,N
       !Time worb construction 
       call cpu_time(worb_t1)     

       !Return indices              
       call index(n1, n1_l,ia_l, ia,iorb)     
       
       !Map local NN grid spatial indices to system grid indices using atom ia 
       !on system grid and corresponding atom ia_l on NN grid
       call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,ia,ia_l,&
                                      ix1,iy1,iz1)


       !Loop over all points in system grid. Distributed
       call cpu_time(VH_t1)
       jrcnt=0
       V_H_local(:)=0._dp
     
       do jr=jr_start,jr_stop

          !Index for V_H_local
          jrcnt=jrcnt+1

          !Unwrap composite index for system grid: jr->(jx,jy,jz)
          call comp_to_ijk(jr,snpts(1),snpts(2),snpts(3),jx,jy,jz)

          !Return x,y,z variables associated with spatial indices
          call return_grid_xyz(jx,jy,jz,sspacing,snpts, x,y,z)
          r_j(:)=(/x,y,z/)


          !*********************************************
          !Integral over r
          !*********************************************
          !Initialise integral 
          integral=0._dp
          !Initialise orbital product
          orbij=0._dp

          !Loop over non-zero orbital products defined on NN grid
          do iro=1,nnz(n1_l)

             !Map non-zero product space to NN grid space
             ir_l=orb_grd_pntr(n1_l,iro)

             !Unwrap composite index for NN grid space: ir_l->(ix_l,iy_l,iz_l)
             call comp_to_ijk(ir_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                              ix_l,iy_l,iz_l)

             !Iterate system grid indices using (ix_l,iy_l,iz_l)
             ix=ix_l+(ix1-1)
             iy=iy_l+(iy1-1)
             iz=iz_l+(iz1-1)

             !Return x,y,z variables associated with spatial indices
             call return_grid_xyz(ix,iy,iz,sspacing,snpts, x,y,z)
             r_i(:)=(/x,y,z/)

             !Composite spatial counter for system grid
             ir=iz+( (iy+((ix-1)*snpts(2))-1)*snpts(3))

             !If |r-r'|=0
             if(ir==jr)then
                !Store 'product of orbitals' value
                orbij=orb_overlap(n1_l,iro)
                cycle
             endif

             !Sum Hartree Potential
             integral=integral+ (orb_overlap(n1_l,iro)/&
                                 sqrt(dot_product(r_i(:)-r_j(:),r_i(:)-r_j(:))))
                                      
          enddo
          !Shut integral over ir
          !*********************************************

          !Store local Hartree potential at spatial point r' (jr), orbital index n1
          ! + Brillouin Approximation for when r=r'
          V_H_local(jrcnt)=(integral*dr)+(orbij*2._dp*pi*r_sph**2)

       end do
       !Shut loop over jr

       
       !Gather the local Hartree potential per process into V_H, available to 
       !all process
       call mpi_allgatherv(V_H_local,Ntotal_lcl,MPI_DOUBLE_PRECISION,     & !Local send
                           V_H, rec_size, displ,MPI_DOUBLE_PRECISION,     & !Global receive
                           MPI_COMM_WORLD,MPIierr)

       !Time loop over Hartree potential, for all n12 
       call cpu_time(VH_t2)
       VH_time_local=VH_time_local+(VH_t2-VH_t1)



       !**********************************************************************
       ! Section 2
       !
       !**********************************************************************
       !Initialise orbital integrals
       worb(:)=0._dp

       !Loop over 2nd orbital index
       do n2=1,N

          !Return indices:
          call index(n2, n2_l,ja_l, ja,jorb)
         

          !Map local NN grid spatial indices to system grid indices using atom ka 
          !on system grid and corresponding atom ka_l on NN grid
          call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,&
                                         ja,ja_l,jx1,jy1,jz1)


          !*****************************************
          !Integral over r'
          !*****************************************
          integral=0._dp

          !Compute integral over spatial variable r'
          do jro=1,nnz(n2_l)

             !Map non-zero product space to NN grid space
             jr_l=orb_grd_pntr(n2_l,jro)

             !Unwrap composite index for NN grid space
             call comp_to_ijk(jr_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                              jx_l,jy_l,jz_l)

             !Iterate system grid indices using (jx_l,jy_l,jz_l)
             jx=jx_l+(jx1-1)
             jy=jy_l+(jy1-1)
             jz=jz_l+(jz1-1)

             !Composite spatial counter for system grid
             jr=jz+( (jy+((jx-1)*snpts(2))-1)*snpts(3))

             !Integral
             integral=integral+ (orb_overlap(n2_l,jro)*V_H(jr))

          enddo
          !*****************************************
      
          !Store orbital integral (mathematically dr' but dr=dr')
          worb(n2)=integral*dr
          
       enddo
       !Shut loop over n2

       !Convert worb units to eV
       worb(:)=worb(:)*unt

       !Time construction of worb, for all n1 
       call cpu_time(worb_t2)
       worb_time_local=worb_time_local+(worb_t2-worb_t1)
       


       !**********************************************************************
       ! Section 3
       ! Construct lower triangle of distributed resonant H2p (packed)
       !**********************************************************************
       call cpu_time(H2p_t1)

       !-----------------------------------------------
       !Access ALL TB coefficients once per (v,v',c,c')
       !-----------------------------------------------

       !Set up 1D arrays to store tight-binding coefficients with Nh hole states
       !+ Ne electron states included, for orbitals n1 and n2
       c_n1_e(:)=evcr(n1, VBT+1:VBT+Ne)     !erow
       c_n1_h(:)=evcr(n1, VBT:VBT-Nh+1:-1)  !hrow

   
       !Loop over index n2, per iteration of index n1
       do n2=1,N
  
          !Set up 1D arrays to store tight-binding coefficients with Nh hole states
          !+ Ne electron states included, for orbital n2
          c_n2_e(:)=evcr(n2, VBT+1:VBT+Ne)    !ecol
          c_n2_h(:)=evcr(n2, VBT:VBT-Nh+1:-1) !hcol

          !----------------------------------------------------------------------------

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
                
               
                !BSE kernel with factorised coefficients
                H2p(Hcnt) = H2p(Hcnt) + &
                                        c_n1_e(e_row_cnt)*c_n2_h(h_col_cnt)*worb(n2)*&
                                       ( 2._dp*c_n1_h(h_row_cnt)*c_n2_e(e_col_cnt)&
                                         -inv_eps*c_n1_e(e_col_cnt)*c_n2_h(h_row_cnt))
                


             enddo
          enddo
          !End loops over rows and columns of H2p
       enddo
       !Ends loop over n2

       !Time construction of H2p, for all n12 
       call cpu_time(H2p_t2)
       H2p_time_local=H2p_time_local+(H2p_t2-H2p_t1)

    enddo
    !Shut loop over n12 


    !Deallocate Hartree potentials
    deallocate(V_H_local,V_H)
    !Deallocate orbital integrals
    if(allocated(worb)) deallocate(worb)




    !***************************************************************
    !Routine Timing
    !Important timings are VH and H2p constructions
    !***************************************************************
    !Finish timing subroutine
    call cpu_time(t2_local)

    !Create arrays with only the local run time (_tl) and arrays
    !to contain run times per process (_tg)
    allocate(H2p_tl(numprocs),VH_tl(numprocs),H2p_tg(numprocs),VH_tg(numprocs))
    H2p_tg(:)=0._dp
    VH_tg(:)=0._dp
    H2p_tl(:)=0._dp
    VH_tl(:)=0._dp

    !Fill local run time into element corresponding to process rank
    H2p_tl(rank+1)=H2p_time_local
    VH_tl(rank+1)=VH_time_local

    !Sum local array times such that each process has a global array that contains run times
    !for all processes
    call mpi_allreduce(H2p_tl,H2p_tg,numprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce(VH_tl,VH_tg,numprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)

    !Sum time taken for each section of routine (time in seconds by default)
    !Redundant. Can be done from summing returned arrays above
    call mpi_allreduce(VH_time_local,VH_time,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce(worb_time_local,worb_time,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce(H2p_time_local,H2p_time,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce((t2_local-t1_local),t21,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)


    !Output average times and return in mins
    if(IO)then
       !VH
       write(*,'(A,X,F20.10,A)') 'Avg time taken to construct V_H =',VH_time/(dble(numprocs)*60._dp),' mins'
       write(*,'(A,X,F20.10,A,X,I3)') 'Longest time taken to construct V_H =',maxval(VH_tg)/60._dp,&
                                      ' mins on process',maxloc(VH_tg)-1
       write(*,'(A,X,F20.10,A)')   'The standard deviation for VH construction time =',stnd_dv(VH_tg)/60._dp,'mins'

       !H2p
       write(*,'(A,X,F20.10,A)') 'Avg time taken to construct H2p =',H2p_time/(dble(numprocs)*60._dp),' mins'
       write(*,'(A,X,F20.10,A,X,I3)') 'Longest time taken to construct H2p =',maxval(H2p_tg)/60._dp,&
                                      ' mins on process',maxloc(H2p_tg)-1
       write(*,'(A,X,F20.10,A)')   'The standard deviation for H2p construction time =',stnd_dv(H2p_tg)/60._dp,'mins'
    endif


    !End of the routine. 
  end subroutine construct_H2p_onsite






































  !*****************************************************************************
  ! TEST H2p.
  ! 1.2nm BSE run crashed with a value of jr that was too large for V_H(jr) in 
  ! section 2. Very concering => could be an issue with mapping. 
  ! Could alternatively be a grid truncation issue and the mapping is actually ok.
  ! Hasn't occurred in any other instance
  !*****************************************************************************


  subroutine TEST_H2p(H2p,rstart,rstop)
    use var,   only: sp,dp,qp,Natoms,orb_bas,evcr,atom_species,d,Nexsp,eigen_cnt,NN_list,NumNN,&
                     r,time_CoulInt,al,NNdp_list,atom_type,pi,all_pnts,eigen_cnt,Ne,Nh,ic,static_eps
    use grids, only: gen_NN_grid, gen_NN_rspacing, return_grid_xyz,return_grid_xyz_arb
    use smat,  only: NN_disp
    use sort,  only: reallocate
    implicit none

    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer       N,Nm,iro,ir,jro,jr,  ia,ja,ka,la,  n1,n2,n3,n4,n12,n34 , iorb,jorb,korb,lorb,&
                  n1_l,n2_l,n3_l,n4_l,Natm,ia_l,ja_l,ka_l,la_l,n1234,NOI,ix1,iy1,iz1,ir_l,&
                  ix_l,iy_l,iz_l,Ntotal,ix,iy,iz,ixy_l, jx1,jy1,jz1,jr_l,jxy_l,jx_l,jy_l,jz_l,&
                  jx,jy,jz,istart,istop,Nreq,Ntotal_lcl,jrcnt,i,n34_start,n34_stop,cnt,VBT,CBB
    integer       row,col,Nh_mx,Ne_mx,h_row_cnt,e_row_cnt,h_col_cnt,e_col_cnt,h_row,e_row,&
                  h_col,e_col,Hcnt,orb_cnt,j
    real(dp)      dr,integral,x,y,z,denom,r_sph,orbij,inv_eps 
    real(dp)      t1_local,t2_local,t21,t21_lg, VH_time,VH_time_lg,VH_time_local,VH_t1,VH_t2, &
                  worb_time, worb_time_lg, worb_time_local,worb_t1,worb_t2, &
                  H2p_time, H2p_time_lg, H2p_time_local,H2p_t1, H2p_t2
    real(qp)      unt
    character     fname*12,rank_lab*3
     

    !Modular Variables
    integer,      intent(in)          :: rstart,rstop
    integer       jr_start,jr_stop

    !Modular arrays
    integer,      allocatable         :: orb_grd_pntr(:,:,:),nnz(:,:)
    real(dp),     allocatable         :: orb_overlap(:,:,:)
    real(dp),     allocatable         :: worb(:),worb_local(:)
    real(dp),     allocatable,        intent(inout) :: H2p(:)

    !Local Arrays
    integer,      allocatable         :: tmp(:)
    integer(sp),  allocatable         :: displ(:),rec_size(:),displ2(:),rec2(:)
    real(dp),     allocatable         :: V_H(:),V_H_local(:)
    real(dp),     allocatable         :: r_i(:),r_j(:),c_n1_e(:),c_n2_e(:),c_n3_e(:), &
                                         c_n2_h(:),c_n3_h(:),c_n4_h(:)
    real(dp),     allocatable         :: H2p_tg(:),VH_tg(:),H2p_tl(:),VH_tl(:)


 
    !**********************************
    !Main Routine
    !**********************************
    if(IO) write(*,*)'Test_H2p'
    call cpu_time(t1_local)

    !Dimensions of the system
    N=Natoms*orb_bas

    !Maximum number of orbital interactions for two centres, in NN regime
    NOI=N*(NumNN*orb_bas)
    if(IO) write(*,*)'Max number of orbital interactions',NOI
   
    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1
    
    !Generate NNdp_list which indicates which of the 4 
    !neighbours ja is to ia: Establish the displacement vector 
    !for each NN of all atoms in system.
    call NN_disp(NNdp_list)

    !All arguments here are defined in the module
    !This NN grid uses the same spacing as the system grid 
    !Important for mapping between the system and NN grids
    call gen_NN_rspacing(NNnpts_r,NNrspacing,NNgrid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)

    !Arrays not utilised so deallocate immediately
    deallocate(NNgrid,NNxyz_ptr)

    !Construct orb_overlap matrix on NN grid
    call construct_orb_overlap(orb_overlap,orb_grd_pntr,nnz)
    !call construct_FULL_orb_overlap(orb_overlap)
    !Index array for NN grid not utilised after orbital products are constructed
    deallocate(NNiptr)

    !Generate grid arrays with dimensions of the supercell
    !No auxilliary grid arrays are allocated, except atm_pntr
    all_pnts=.true.
    !Note, use return_grid_xyz with setup_system_grid
    call setup_system_grid()

    !For testing: Construct system grid with vacuum region included
    !call setup_arb_system_grid()
    !if(IO) write(*,*) 'Note, use return_grid_xyz_arb with setup_arb_system_grid'
    !deallocate(sgrid,siptr,sys_pntr)

    !Total number of points in system grid
    Ntotal=product(snpts)
    if(IO) write(*,*) 'Total number of points in Hartree potential',Ntotal

    !Volume element associated with this grid spacing
    dr=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Radius of cubic volume element, given as the length from 
    !the centre to the corner (see Brillouin approximation)
    r_sph=sqrt(3._dp)*0.5*NNrspacing(1,1)

    !Position vectors for Coulomb Potential
    allocate(r_i(d),r_j(d))

    !Index arrays for the two distributed loops
    allocate(rec_size(numprocs), displ(numprocs))
    allocate(rec2(numprocs), displ2(numprocs))

    !Units. e**2/(4 pi eps_0 a_0) = 1 Ha
    !Convert atomic positions to Bohr
    !Convert Hartree to eV  
    unt=(0.5291772109217*27.2113850560)

 


    !******************************************************************************
    ! 1)   Compute Hartree Potential: V_H(r')=  integral{ orb(n12,r)/|r-r'| } dr
    !      per orbital product n12
    ! 2)   Compute integral{ orb(n34,r') V_H(n12,r')} dr'  per orbital product n12
    ! 3)   Compute elements of packed, distributed lower triangle of resonant, 
    !      2-particle matrix H2p
    !******************************************************************************

    !Distribute spatial index jr over supercell grid 
    !rec_size = size(V_H_local) for each process
    !disp = Index of V_H_local w.r.t. V_H
    call distribute_loop2(1,Ntotal, jr_start,jr_stop,rec_size,displ)

    !Number of grid points per CPU
    Ntotal_lcl=jr_stop-jr_start+1

    !Allocate global and local Hartree potentials
    allocate(V_H_local(Ntotal_lcl),V_H(Ntotal))
    V_H(:)=0._dp
    V_H_local(:)=0._dp

    !Distribute over worb n34=[1,NOI]. 
    !This has been moved out of the loop structure as the limits
    !are now fixed (couldn't get db_factor to work)
    call distribute_loop2(1,NOI, n34_start,n34_stop, rec2,displ2)
    allocate(worb(NOI),worb_local(n34_stop-n34_start+1))




    !**********************************************************************
    ! Section 1
    !
    !**********************************************************************
    !Timings
    !Timing V_H construction
    VH_time_local=0._dp
    !Timing worb construction
    worb_time_local=0._dp
    !Timing H2p construction
    H2p_time_local=0._dp
   


    !First composite index over orbitals 1&2 defined with spatial coordinate r (ir)
    do n12=1,NOI

       !Time worb construction 
       call cpu_time(worb_t1)


       !Return indices in composite index n12
       call indices(n12, n1_l,n2_l,ia_l,ja_l,ia,ja,iorb,jorb)
       !If ia doesn't have a neighbour ja, cycle loop
       if(ja==0)cycle

       !Global orbital indices n1 and n2 (used in eigenvectors below)
       n1=iorb+((ia-1)*orb_bas)
       n2=jorb+((ja-1)*orb_bas)
       
       !Map local NN grid spatial indices to system grid indices using atom ia 
       !on system grid and corresponding atom ia_l on NN grid
       call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,ia,ia_l,&
                                      ix1,iy1,iz1)



       !Loop over all points in system grid. Distributed
       call cpu_time(VH_t1)
       jrcnt=0
       V_H_local(:)=0._dp
     
       do jr=jr_start,jr_stop

          !Index for V_H_local
          jrcnt=jrcnt+1

          !Unwrap composite index for system grid: jr->(jx,jy,jz)
          call comp_to_ijk(jr,snpts(1),snpts(2),snpts(3),jx,jy,jz)

          !Return x,y,z variables associated with spatial indices
          call return_grid_xyz(jx,jy,jz,sspacing,snpts, x,y,z)
          r_j(:)=(/x,y,z/)


          !*********************************************
          !Integral over r
          !*********************************************
          !Initialise integral 
          integral=0._dp
          !Initialise orbital product
          orbij=0._dp

          !Loop over non-zero orbital products defined on NN grid
          do iro=1,nnz(n1_l,n2_l)
          !do ir_l=1,product(NNnpts_r)

             !Map non-zero product space to NN grid space
             ir_l=orb_grd_pntr(n1_l,n2_l,iro)

             !Unwrap composite index for NN grid space: ir_l->(ix_l,iy_l,iz_l)
             call comp_to_ijk(ir_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                              ix_l,iy_l,iz_l)

             !Iterate system grid indices using (ix_l,iy_l,iz_l)
             ix=ix_l+(ix1-1)
             iy=iy_l+(iy1-1)
             iz=iz_l+(iz1-1)

             !Return x,y,z variables associated with spatial indices
             call return_grid_xyz(ix,iy,iz,sspacing,snpts, x,y,z)
             r_i(:)=(/x,y,z/)

             !Composite spatial counter for system grid
             ir=iz+( (iy+((ix-1)*snpts(2))-1)*snpts(3))

             !If |r-r'|=0
             if(ir==jr)then
                !Store 'product of orbitals' value
                orbij=orb_overlap(n1_l,n2_l,iro)
                !orbij=orb_overlap(n1_l,n2_l,ir_l)
                !write(100,'(2(I6,X),7(F12.6,X))') ir,jr,r_i(:),r_j(:),orbij
                cycle
             endif

             !Sum Hartree Potential
             integral=integral+ (orb_overlap(n1_l,n2_l,iro)/&
                                 sqrt(dot_product(r_i(:)-r_j(:),r_i(:)-r_j(:))))
                                 
                  
          enddo
          !Shut integral over ir
          !*********************************************

          !Store local Hartree potential at spatial point r' (jr), orbital pair index n12
          ! + Brillouin Approximation for when r=r'
          V_H_local(jrcnt)=(integral*dr)+(orbij*2._dp*pi*r_sph**2)

       end do
       !Shut loop over jr

       
       !Gather the local Hartree potential per process into V_H, available to 
       !all process
       call mpi_allgatherv(V_H_local,Ntotal_lcl,MPI_DOUBLE_PRECISION,     & !Local send
                           V_H, rec_size, displ,MPI_DOUBLE_PRECISION,     & !Global receive
                           MPI_COMM_WORLD,MPIierr)

       !Deallocate local Hartree potentials to temporarily save memory
       !deallocate(V_H_local)

       !Time loop over Hartree potential, for all n12 
       call cpu_time(VH_t2)
       VH_time_local=VH_time_local+(VH_t2-VH_t1)



       !**********************************************************************
       ! Section 2
       !
       !**********************************************************************
       !Initialise orbital integrals
       worb_local(:)=0._dp
       worb(:)=0._dp

       !Initialise index for worb_local 
       cnt=0


       !Loop over 2nd composite index n34. Distributed 
       do n34=n34_start,n34_stop
          
          !Iterate w_local index
          !Must be defined before cycle, such that elements of worb_local
          !are consistent with the parallelisation (wastes some memory as a result)
          cnt=cnt+1

          !Return indices in composite index n34 
          call indices(n34, n3_l,n4_l,ka_l,la_l,ka,la,korb,lorb)
          !If ia doesn't have a neighbour la, cycle loop
          if(la==0)cycle

          !Map local NN grid spatial indices to system grid indices using atom ka 
          !on system grid and corresponding atom ka_l on NN grid
          call MapSpatialIndices_NN_sys2(NNatm_pntr,satm_pntr,NNnpts_r,snpts,&
                                         ka,ka_l,jx1,jy1,jz1)


          !*****************************************
          !Integral over r'
          !*****************************************
          integral=0._dp

          !Compute integral over spatial variable r'
          do jro=1,nnz(n3_l,n4_l)

             !Map non-zero product space to NN grid space
             jr_l=orb_grd_pntr(n3_l,n4_l,jro)

             !Unwrap composite index for NN grid space
             call comp_to_ijk(jr_l,NNnpts_r(1),NNnpts_r(2),NNnpts_r(3),&
                              jx_l,jy_l,jz_l)
             !Instead do this by looking up values from index array.
             !See if it affects the result
             !jx_l=NNiptr(jr_l,1)
             !jy_l=NNiptr(jr_l,2)
             !jz_l=NNiptr(jr_l,3)

             !Iterate system grid indices using (jx_l,jy_l,jz_l)
             jx=jx_l+(jx1-1)
             jy=jy_l+(jy1-1)
             jz=jz_l+(jz1-1)

             !Composite spatial counter for system grid
             jr=jz+( (jy+((jx-1)*snpts(2))-1)*snpts(3))

             !Generate (x,y,z) values from (jx,jy,jz) and output for the instance
             !the code crashes
!!$             if(n12==1 .and. n34==1381)then
!!$                x=0._dp
!!$                y=0._dp
!!$                z=0._dp
!!$                call return_grid_xyz(jx,jy,jz,NNrspacing,NNnpts_r, x,y,z)
!!$                write(400,*) x,y,z
!!$             endif

             !Integral
             integral=integral+ (orb_overlap(n3_l,n4_l,jro)*V_H(jr))

          enddo
          !*****************************************
      
          !Store orbital integral (mathematically dr' but dr=dr')
          worb_local(cnt)=integral*dr
          
       enddo
       !Shut loop over n34

       !Gather the local orbital integral  per process into worb
       call mpi_allgatherv(worb_local, n34_stop-n34_start+1 ,MPI_DOUBLE_PRECISION,     & !Local send
                           worb,       rec2, displ2         ,MPI_DOUBLE_PRECISION,     & !Global receive
                           MPI_COMM_WORLD,MPIierr)

       !Convert worb to eV
       worb(:)=worb(:)*unt

       do n34=1,NOI
          !Write worb to file
          if(IO) write(100,*) n12,n34,worb(n34)
       enddo

       !Time construction of worb, for all n12 
       call cpu_time(worb_t2)
       worb_time_local=worb_time_local+(worb_t2-worb_t1)
       
    enddo
    !Shut loop over n12 

    !Deallocate Hartree potentials
    deallocate(V_H_local,V_H)

    !Deallocate orbital integrals
    if(allocated(worb)) deallocate(worb)
    if(allocated(worb_local)) deallocate(worb_local)




    !***************************************************************
    !Routine Timing
    !Important timings are VH and H2p constructions
    !***************************************************************
    !Finish timing subroutine
    call cpu_time(t2_local)

    !Create arrays with only the local run time (_tl) and arrays
    !to contain run times per process (_tg)
    allocate(VH_tl(numprocs),VH_tg(numprocs))
    VH_tg(:)=0._dp
    VH_tl(:)=0._dp

    !Fill local run time into element corresponding to process rank
    VH_tl(rank+1)=VH_time_local

    !Sum local array times such that each process has a global array that contains run times
    !for all processes
    call mpi_allreduce(VH_tl,VH_tg,numprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    !Sum time taken for each section of routine (time in seconds by default)
    !Redundant. Can be done from summing returned arrays above
    call mpi_allreduce(VH_time_local,VH_time,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce(worb_time_local,worb_time,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)
    call mpi_allreduce((t2_local-t1_local),t21,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIierr)


    !Output average times and return in mins
    if(IO)then
       !VH
       write(*,'(A,X,F20.10,A)') 'Avg time taken to construct V_H =',VH_time/(dble(numprocs)*60._dp),' mins'
       write(*,'(A,X,F20.10,A,X,I3)') 'Longest time taken to construct V_H =',maxval(VH_tg)/60._dp,&
                                      ' mins on process',maxloc(VH_tg)-1
       write(*,'(A,X,F20.10,A)')   'The standard deviation for VH construction time =',stnd_dv(VH_tg)/60._dp,'mins'
    endif


    if(IO)write(*,*)'TEST_H2p stops line 4488 of exbnd.f90'
    call mpi_barrier(MPI_COMM_WORLD,MPIierr)
    call stop_mpi()

    !End of the routine. 
  end subroutine TEST_H2p


 !****************************************************************************
  !Construct orb_overlap matrix on NN grid with ALL POINTS
  !****************************************************************************

  subroutine construct_FULL_orb_overlap(orb_overlap)
    
    !Global variables and arrays
    use var,     only: NumNN,orb_bas
    !Orbital functions/profiles
    use orbital, only: radial_function,orbital_profile 
    !Place orbitals on grid
    use smat,    only: place_orbitals 

    implicit none


    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer                 Natm,Ntotal,Nsub    
    integer                 ia,iorb,ja,jorb,n1,n2,ir,iro,ix,iy,iz,i,j
    integer                 x0,y0,z0,xf,yf,zf       !Unused

    !Global Arrays    
    integer,                allocatable :: midpnt(:)
    type(radial_function),  allocatable :: &
                            sorb(:,:,:),pxorb(:,:,:),pyorb(:,:,:),pzorb(:,:,:),sstorb(:,:,:)

    !Modular arrays
    real(dp),     allocatable, intent(out) :: orb_overlap(:,:,:)
    real(dp),     allocatable              :: orbdummy(:,:,:)
    
    !Local Arrays
    complex(dp),            allocatable :: orbital1(:,:,:),orbital2(:,:,:)
    


    !****************
    !Main Routine
    !****************
    !NN grid set up in calling routine and available to the module

    !Generate 3D orbital profiles
    call orbital_profile(NNrspacing,midpnt,sorb,pxorb,pyorb,pzorb,sstorb)

    !Number of atoms in the grid = 1 central atom + 4 NN
    Natm=NumNN

    !Total number of orbital states in NN subsystem
    Nsub=Natm*orb_bas

    !Total number of points in NN grid
    Ntotal=product(NNnpts_r)

    !Orbitals 1 and 2, filled per iteration. Dimensions of the  NN grid
    allocate(orbital1(NNnpts_r(1),NNnpts_r(2),NNnpts_r(3)),orbital2(NNnpts_r(1),NNnpts_r(2),NNnpts_r(3)))

    !orb_overlap has dimensions of the NN grid and contains the overlap between orbitals n1 and n2
    !3rd dimension allocated with upper bound 
    allocate(orb_overlap(Nsub,Nsub,Ntotal))
    orb_overlap(:,:,:)=0._dp


    !*****************************************************************
    !Loop over all atoms in a NN grid to create an array containing all 
    !orbital overlaps/products
    !Order is consistent with loops in orbital integral (important!)
    !*****************************************************************

    do ia=1,Natm
       do iorb=1,orb_bas

          !Composite state index 1
          n1=iorb+((ia-1)*orb_bas)

          do ja=1,Natm
             do jorb=1,orb_bas

                !Composite state index 2
                n2=jorb+((ja-1)*orb_bas)

                !Place appropriate orbitals on grid:orbital1 and orbital2 returned.
                call place_orbitals(orbital1,orbital2,NNatm_pntr,NNiptr,&
                     ia,ja,iorb,jorb,midpnt,sorb,pxorb,pyorb,pzorb,sstorb,&
                     x0,y0,z0,xf,yf,zf,atom_species_NN) 

         
                !Construct an array containing the overlap between two orbitals 
                !Loop over all points in the NNgrid. Can then easily define a composite
                !counter ir and map that to the overlap spatial index iro

                iro=0
                ir=0
  
                do ix=1,NNnpts_r(1)
                   do iy=1,NNnpts_r(2)
                      do iz=1,NNnpts_r(3)

                         !Composite index for NN grid
                         ir=ir+1

                         !Store for all space, not just where the product is non-zero
                         orb_overlap(n1,n2,ir)=real(orbital1(ix,iy,iz))*real(orbital2(ix,iy,iz))

                      enddo
                   enddo
                enddo

             
             enddo
          enddo

       enddo
    enddo


    !Return orb_overlap 
  end subroutine construct_FULL_orb_overlap



  !**********************************************************************
  !Construct the 2-particle kernel in an orthognal, atomic orbital basis
  !**********************************************************************
  subroutine orthogonal_construct_H2p(H2p,rstart,rstop)
    use var, only: dp,Natoms,NumNN,orb_bas,eigen_cnt,r,Ne,Nh,al,evcr,NN_list,ic,atom_type,bond_length,Eg
    implicit none

    !******************
    !Declarations
    !******************

    !Global Variables
    integer,      intent(in)          :: rstart,rstop

    !Global Arrays
    real(dp),     allocatable,        intent(inout) :: H2p(:)
    real(dp),     allocatable         :: U_c(:,:),U_x(:,:)

    !Local Variables
    integer  ::    N,VBT,ia,ja,ia_NN,ia_l,ja_l,iorb,jorb,i,j,i_l,j_l,&
                   Hcnt,row,col,h_row_cnt,h_col_cnt,e_row_cnt,e_col_cnt, cnt
    real(dp) ::    t1_local,t2_local,unt,ieps,Uc,Ux,d_2N,denom,V

    !Local Arrays
    integer,  allocatable :: atom_species(:)
    real(dp), allocatable :: c_vi(:),c_ci(:),c_vj(:),c_cj(:)

 
    !******************
    !Main Routine
    !******************
    call cpu_time(t1_local)
    if(IO) write(*,*)'Constructing H2p in an orthogonal basis'

    !Dimensions of the system (==total number of orbital states)
    N=Natoms*orb_bas

    !VBT
    VBT=eigen_cnt(1,1)

    !Read in tabulated, two-centre orbital integrals (in eV)
    call orbital_integral_energies(U_c,U_x)
    
    !Bulk, high-frequency dielectric constant for CdSe (Bulk static 9.7)
    if(IO .and. ic/=19) write(*,*)'No dielectric constant stored for this material'
    !!if(ic==19)  ieps=(1._dp/6.2)  
    !High-frequency dielectric as a function of NC diameter (enters via band-gap)
    if(ic==19)  ieps=(1._dp/eps_inf_dot(Eg))

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
    H2p(:)=0._dp


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
          !!!!if(ja/=ia) cycle  !On-site only
          ja_l=atom_species(ja)
     
          !NN Coulomb integral taken to be monopole term
          !!!!if(ia/=ja) Uc=1._dp/sqrt(dot_product(r(ia,:)-r(ja,:),r(ia,:)-r(ja,:)))


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
                !!!!if(ia==ja) Uc=U_c(i_l,j_l)
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
    !Beyond NN, only retain terms multiplied by (ii|v|jj)
    do ia=1,Natoms
       do ja=1,Natoms
          !|r-r'|
          denom=sqrt(dot_product(r(ia,:)-r(ja,:),r(ia,:)-r(ja,:)))
          !If atoms ia and ja are NN or same atom, cycle
          if(denom<d_2N) cycle 
          !Inter-atomic integral (ii|v|jj) =~ 1/|r_ia-r_ja|
          V=unt/denom
    
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
                                    -ieps*c_vi(h_row_cnt)*c_cj(e_row_cnt)*c_vi(h_col_cnt)*c_cj(e_col_cnt) )*V

                   enddo
                enddo
                
             enddo
          enddo
       enddo
    enddo

   
    !Output routine runtime
    call cpu_time(t2_local)
    if(IO)write(*,'(A,F20.10,A)') 'Time taken to construct H2p = ',(t2_local-t1_local)/60._dp,' mins'


  end subroutine orthogonal_construct_H2p

  
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
  !Construct the 2-particle kernels in an orthognal, atomic orbital basis
  !H2p=K_(vc)(v'c') and K2=K_(vc)(c'v'), required for the definition of
  !the full BSE Hamiltonian
  !**********************************************************************
  subroutine orthogonal_construct_H2p_K2(H2p,K2,rstart,rstop)
    use var, only: dp,Natoms,NumNN,orb_bas,eigen_cnt,r,Ne,Nh,al,evcr,NN_list,ic,atom_type,bond_length,Eg
    implicit none

    !******************
    !Declarations
    !******************

    !Global Variables
    integer,      intent(in)          :: rstart,rstop

    !Global Arrays
    real(dp),     allocatable,        intent(inout) :: H2p(:),K2(:)
    real(dp),     allocatable         :: U_c(:,:),U_x(:,:)

    !Local Variables
    integer  ::    N,VBT,ia,ja,ia_NN,ia_l,ja_l,iorb,jorb,i,j,i_l,j_l,&
                   Hcnt,row,col,h_row_cnt,h_col_cnt,e_row_cnt,e_col_cnt, cnt
    real(dp) ::    t1_local,t2_local,unt,ieps,Uc,Ux,d_2N,denom,V,c1,c2,c3

    !Local Arrays
    integer,  allocatable :: atom_species(:)
    real(dp), allocatable :: c_vi(:),c_ci(:),c_vj(:),c_cj(:)

 
    !******************
    !Main Routine
    !******************
    call cpu_time(t1_local)
    if(IO) write(*,*)'Constructing H2p and K2 in an orthogonal basis'

    !Dimensions of the system (==total number of orbital states)
    N=Natoms*orb_bas

    !VBT
    VBT=eigen_cnt(1,1)

    !Read in tabulated, two-centre orbital integrals (in eV)
    call orbital_integral_energies(U_c,U_x)
    
    !Bulk, high-frequency dielectric constant for CdSe (Bulk static 9.7)
    if(IO .and. ic/=19) write(*,*)'No dielectric constant stored for this material'
    !!if(ic==19)  ieps=(1._dp/6.2)  
    !High-frequency dielectric as a function of NC diameter (enters via band-gap)
    if(ic==19)  ieps=(1._dp/eps_inf_dot(Eg))

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
    H2p(:)=0._dp
    K2(:) =0._dp

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
          !!!if(ja/=ia) cycle  !On-site only
          ja_l=atom_species(ja)
     
          !NN Coulomb integral taken to be monopole term
          !!!!if(ia/=ja) Uc=1._dp/sqrt(dot_product(r(ia,:)-r(ja,:),r(ia,:)-r(ja,:)))


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
                !!!!if(ia==ja) Uc=U_c(i_l,j_l)
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

                      
                      !Terms * Coulomb integral 
                      !H2p(Hcnt) = H2p(Hcnt) + c1*(2._dp*Uc) + c2*(-ieps*Uc) 
                      !K2(Hcnt)  = K2(Hcnt)  + c1*(2._dp*Uc) + c3*(-ieps*Uc) 
                      !Terms * Exchange integrals 
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
    !Beyond NN, only retain terms multiplied by (ii|v|jj)
    do ia=1,Natoms
       do ja=1,Natoms
          !|r-r'|
          denom=sqrt(dot_product(r(ia,:)-r(ja,:),r(ia,:)-r(ja,:)))
          !If atoms ia and ja are NN or same atom, cycle
          if(denom<d_2N) cycle 
          !Inter-atomic integral (ii|v|jj) =~ 1/|r_ia-r_ja|
          V=unt/denom
    
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

                      !Only retain terms in the kernel multiplied by Uc (denoted as V here)

                      !Packed Two-particle Hamiltonian 
                      !(technically just the kernel elements) K_(vc)(v'c')
                      H2p(Hcnt) = H2p(Hcnt) +  (2._dp*c1-ieps*c2)*V
      
                      !Second Coupling matrix K_(vc)(c'v')
                      K2(Hcnt)  = K2(Hcnt)  +  (2._dp*c1-ieps*c3)*V

                   enddo
                enddo
                
             enddo
          enddo
       enddo
    enddo

   
    !Output routine runtime
    call cpu_time(t2_local)
    if(IO)write(*,'(A,F20.10,A)') 'Time taken to construct H2p = ',(t2_local-t1_local)/60._dp,' mins'


  end subroutine orthogonal_construct_H2p_K2







  !****************************************************************************
  !Compute on-site and NN orbital Coulomb and exchange integrals
  !No point optimising the routine as the integrals can be computed once and
  !tabulated for future use
  !****************************************************************************  
  subroutine orbital_integrals()
    use var,     only: NumNN,orb_bas,s1,s2,s3,ns,bond_length,ic
    use orbital, only: radial_function,orbital_profile 
    use smat,    only: place_orbitals 
    use grids,   only: gen_NN_rspacing
    implicit none


    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer                 Natm
    integer                 ia,iorb,ja,jorb,ir,jr,ix,iy,iz,jx,jy,jz,i,j
    integer                 x0,y0,z0,xf,yf,zf       !Unused but required arguements
    real(dp)                unt,sum,dV,denom,t1,t2,O1,O2
    character               ia_lab*1,ja_lab*1,fname*14

    !Global Arrays    
    integer,                allocatable :: midpnt(:)
    type(radial_function),  allocatable :: &
                            sorb(:,:,:),pxorb(:,:,:),pyorb(:,:,:),pzorb(:,:,:),sstorb(:,:,:)
    complex(dp),            allocatable :: orbital1(:,:,:),orbital2(:,:,:)

    integer,                allocatable ::  NNnpts_r(:),NNiptr(:,:),NNatm_pntr(:),NNxyz_ptr(:,:,:),&
                                            atom_species_NN(:)
    real(dp),               allocatable ::  NNrspacing(:,:),grid(:,:)    


    !****************
    !Main Routine
    !****************
    call cpu_time(t1)
    if(IO) write(*,*) 'Computing on-site and NN Coulomb and exchange integrals'

    !Set up NN grid. Centred on anion (0,0,0) although doesn't matter due to
    !symmetry of the integrals

    !In atom_species, cation=1, anion=2    
    call gen_NN_rspacing(NNnpts_r,NNrspacing,grid,NNiptr,NNatm_pntr,NNxyz_ptr,atom_species_NN)
    deallocate(NNxyz_ptr)

    !Volume element
    dV=NNrspacing(1,1)*NNrspacing(2,2)*NNrspacing(3,3)

    !Generate 3D orbital profiles
    call orbital_profile(NNrspacing,midpnt,sorb,pxorb,pyorb,pzorb,sstorb)

    !Orbitals 1 and 2, filled per iteration. Dimensions of the  NN grid
    allocate(orbital1(NNnpts_r(1),NNnpts_r(2),NNnpts_r(3)),orbital2(NNnpts_r(1),NNnpts_r(2),NNnpts_r(3)))

    !Units. e**2/(4 pi eps_0 a_0) = 1 Ha
    !Unit converts to eV
    unt=14.399644849251

    !Open files
    if(IO)then 
       !open(unit=101,file='coul_int.dat')
       open(unit=101,file='exc_int.dat')
    endif

    !Headers
    if(IO)then 
       !write(101,*)'On-site and NN Coulomb Integrals in Siesta basis (upper triangle)'
       !write(101,*)'Grid Integers:',s1,s2,s3,'. Grid spacing (Ang):',NNrspacing(1,1)
       write(101,*)'On-site and NN Exchange Integrals in Siesta basis (upper triangle)'
       write(101,*)'Grid Integers:',s1,s2,s3,'. Grid spacing (Ang):',NNrspacing(1,1)
    endif


    !Loop over anion and cation pairs (4 combinations) 
    !Produces matrix with same dimensions as Hamiltonian
    !Only need to compute upper triangle due to symmetry 
    do ia=1,ns 
       do ja=1,ns
          do iorb=1,orb_bas
             do jorb=1,orb_bas

                i=iorb+(ia-1)*orb_bas
                j=jorb+(ja-1)*orb_bas
                !Skip lower triangle elements 
                if(i>j) cycle

                !Place appropriate orbitals on grid: orbital1 and orbital2 returned.
                !Don't think x0-zf are utilised (check)
                call place_orbitals(orbital1,orbital2,NNatm_pntr,NNiptr,&
                                    ia,ja,iorb,jorb,midpnt,sorb,pxorb,pyorb,pzorb,sstorb,&
                                    x0,y0,z0,xf,yf,zf,atom_species_NN) 
                
                !Do integral: Loop order consistent with gen_NN_rspacing 
                sum=0._dp
                ir=0
                do ix=1,NNnpts_r(1)
                   do iy=1,NNnpts_r(2)
                      do iz=1,NNnpts_r(3)
                         ir=ir+1
                         
                         !Valid for Coulomb and exchange
                         O1=real(orbital1(ix,iy,iz))
                         if(O1==0._dp) cycle
                         !Valid for exchange
                         O2=real(orbital2(ix,iy,iz))
                         if(O2==0._dp) cycle


                         jr=0
                         do jx=1,NNnpts_r(1)
                            do jy=1,NNnpts_r(2)
                               do jz=1,NNnpts_r(3)
                                  jr=jr+1
                                  
                                  !Denominator 1/|r1-r2|
                                  denom=sqrt(dot_product(grid(:,ir)-grid(:,jr),grid(:,ir)-grid(:,jr)))
                                  !Truncate singularity to sensible value
                                  if(ir==jr) denom=0.001*bond_length(ic)
                                  !sqrt(0.75)*NNrspacing(1,1)

                                  !Sum integral
                                  !Coulomb 
                                  !sum=sum+ O1*O1*&
                                  !         real(orbital2(jx,jy,jz))*real(orbital2(jx,jy,jz))/denom
                                  !Exchange
                                  sum=sum+ O1*O2*&
                                           real(orbital1(jx,jy,jz))*real(orbital2(jx,jy,jz))/denom
                            
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo

                !Output integral values (in eV) 
                if(IO) write(101,*) i, j, unt*sum*dV*dV
                                                       
             enddo
          enddo
       enddo
    enddo

    !Close file
    close(101)

    !Routine time
    call cpu_time(t2)
    if(IO) write(*,*)'Integrals took:',(t2-t1)/60._dp,' mins.'
    
  end subroutine orbital_integrals



  !****************************************************************************
  !Compute 2nd and 3rd NN orbital Coulomb and exchange integrals
  !****************************************************************************  
  subroutine orbital_integrals2()
    use var,     only: dp,d,al,T,bas_atm,bond_length,ic,orb_bas,s1,s2,s3
    use orbital, only: radial_function,orbital_profile 
    use smat,    only: place_orbitals 
    use grids,   only: gen_local_atomic_grid
    use sort,    only: reallocate
    implicit none


    !*****************
    !Declarations
    !*****************

    !Local Variables
    integer                 nx,ny,nz,cnt,Natoms,u1,u2,u3,ia_c,ib,nb,is
    integer                 x0,y0,z0,xf,yf,zf,neighbour,NumNN,iaNN
    integer                 ia,iorb,ja,jorb,ir,jr,ix,iy,iz,jx,jy,jz,i,j,cnt1,cnt2,cnt3,cnt4
    real(dp)                d1,d2,d3,tol,t1,t2,O1,unt,summ,dV,denom,d_val!,Vss,Vsp,Vpp,Vsstp
    character               fname*9,ia_lab*2,ja_lab*2

    !Local Arrays
    integer,                allocatable :: atom_type(:),atom_type_tmp(:),NNindex(:)
    real(dp),               allocatable :: r(:,:),rtmp(:,:),r_centre(:),r1(:),r2(:),Uc(:,:),&
                                           Vss(:),Vsp(:),Vpp(:),Vsstp(:)
    !Global Arrays    
    integer,                allocatable :: midpnt(:)
    type(radial_function),  allocatable :: &
                            sorb(:,:,:),pxorb(:,:,:),pyorb(:,:,:),pzorb(:,:,:),sstorb(:,:,:)
    complex(dp),            allocatable ::  orbital1(:,:,:),orbital2(:,:,:)

    integer,                allocatable ::  npts_r(:),iptr(:,:),atm_pntr(:), atom_species(:)
    real(dp),               allocatable ::  rspacing(:,:)


    !****************
    !Main Routine
    !****************
    call cpu_time(t1)
    if(IO) write(*,*) 'Computing orbital Coulomb integrals'

    !------------------------------------------------------------
    !User-Settings 
    !------------------------------------------------------------
    !Select 1,2, or 3 for 1st to 3rd NN
    neighbour=3

    !Select one for anion-centred grid or 5 for cation-centred grid
    is=5
    !------------------------------------------------------------

    !Guess at a supercell size for central atom and up to 3rd NN
    !(Tests show this supercell to be sufficiently large up 3rd NN interactions)
    nx=6
    ny=6
    nz=6

    !Number of atoms in basis (conventional cell)
    !Note, input script must be using conventional cell
    nb=8

    !Local number of atoms
    Natoms=nx*ny*nz*nb

    !Allocate arrays (local arrays)
    allocate(r(Natoms,d),atom_type(Natoms)) 
    r(:,:)=0._dp
    atom_type(:)=0

    !Generate atomic positions 
    ia=0
    do u3=0,nz-1
       do u2=0,ny-1
          do u1=0,nx-1
             do ib=1,nb
                !Atomic position
                ia=ia+1
                r(ia,1:d) = bas_atm(ib,1:d) +  u1*T(1,1:d) + u2*T(2,1:d) + u3*T(3,1:d)
                !Store atomic type
                if(ib <= int(Nb*0.5)) atom_type(ia)=-1
                if(ib >  int(Nb*0.5)) atom_type(ia)= 1
             enddo
          enddo
       enddo
    enddo


    !Central atom hard-coded. bas_atm(1,:)=anion, bas_atm(2,:)=cation
    allocate(r_centre(d))
    
    !Calculate central atomic coordinates of the system generated.
    !If nx(ny,nz) is odd i.e. 3, ANInt will round 0.5*3. to 2. 
    r_centre(:) = bas_atm(is,:) +  ANInt(0.5*nx)*T(1,:) + ANint(0.5*ny)*T(2,:) + ANInt(0.5*nz)*T(3,:)

    !1st, 2nd and 3rd NN distances
    d1=bond_length(ic)
    d2=al/sqrt(2._dp)
    d3=0.82915619758885*al !i.e. sqrt(11)/4
    !d4=al

    if(IO)then
       write(*,*) 'Separations (Ang): NN,  2nd NN,  3rd NN'
       write(*,*)  d1,d2,d3
    endif

    !1st-3rd NN system
    if(neighbour==1)then
       if(IO .and. is==1)write(*,*)'NN integrals. Anion-Cation'
       if(IO .and. is==5)write(*,*)'NN integrals. Cation-Anion'
       d_val=d1
    elseif(neighbour==2)then
       if(IO .and. is==1)write(*,*)'2nd NN integrals. Anion-Anion'
       if(IO .and. is==5)write(*,*)'2nd NN integrals. Cation-Cation'
       d_val=d2
    elseif(neighbour==3)then
       if(IO .and. is==1)write(*,*)'3rd NN integrals. Anion-Cation '
       if(IO .and. is==5)write(*,*)'3rd NN integrals. Cation-Anion'
       d_val=d3
    endif
   
    !Tolerance
    tol=0.01*bond_length(ic)

    !Temporary arrays with upper size bounds
    allocate(rtmp(Natoms,d),atom_type_tmp(Natoms))
    rtmp(:,:)=0._dp
    atom_type_tmp(:)=0

    !Delete atoms beyond nth NN distance, refining the system size
    !to relevant atoms
    cnt=0
    do ia=1,Natoms
       !If atom is within nth NN distance, keep it
       if(sqrt(dot_product(r(ia,:)-r_centre(:),r(ia,:)-r_centre(:)))<= d_val+tol)then
          cnt=cnt+1
          !Store position, centred w.r.t. central anion
          rtmp(cnt,:)=r(ia,:)-r_centre(:)
          atom_type_tmp(cnt)=atom_type(ia)
          if(sqrt(dot_product(r(ia,:)-r_centre(:),r(ia,:)-r_centre(:)))==0._dp)then
             if(IO)write(*,*)'Central atom type:',atom_type(ia)
             !Define central atom index
             ia_c=cnt
          endif
       endif
    enddo

    !Central atom included in this counts
    if(neighbour==1) write(*,*)'Number of NN',cnt
    if(neighbour==2) write(*,*)'Number of 1st & 2nd NN',cnt
    if(neighbour==3) write(*,*)'Number of 1st, 2nd & 3rd NN',cnt

    !Number of atoms retained
    Natoms=cnt

    !Allocate correct array sizes
    deallocate(r,atom_type)
    allocate(r(Natoms,d),atom_type(Natoms))
    r(1:Natoms,:)=rtmp(1:Natoms,:)
    atom_type(1:Natoms)=atom_type_tmp(1:Natoms)
    deallocate(rtmp,atom_type_tmp)

    !Atom species array required by place_orbitals
    allocate(atom_species(Natoms))
    do ia=1,Natoms
       !cation
       if(atom_type(ia)==1)  atom_species(ia)=1
       !anion
       if(atom_type(ia)==-1) atom_species(ia)=2
    enddo

    !Output atomic positions
    !do ia=1,Natoms
    !   if(IO) write(100,*) r(ia,:)
    !enddo

    !Set up sampling grid based on these atomic positions
    call gen_local_atomic_grid(r, npts_r,rspacing,iptr,atm_pntr)
    if(IO) write(*,*)'Grid Integers:',s1,s2,s3,'. Grid spacing (Ang):',rspacing(1,1)
 

    !--------------------------------------------------------------
    !Set up a NN array containing only NN, 2nd NN or 3rd NN indices
    !--------------------------------------------------------------
    !Allocate with upper bound
    allocate(NNindex(Natoms))
    NNindex(:)=0

    !All neighbours w.r.t. central atom, only 
    if(neighbour==1)then
       cnt=0
       do ia=1,Natoms
          if(sqrt(dot_product(r(ia_c,:)-r(ia,:),r(ia_c,:)-r(ia,:)))<=d1+tol)then
             if(sqrt(dot_product(r(ia_c,:)-r(ia,:),r(ia_c,:)-r(ia,:)))==0._dp)cycle
             cnt=cnt+1
             NNindex(cnt)=ia
          endif
       enddo
    endif

    !Central atom is first value in the array in these two instances
    !as if statements won't find it (centred at (0,0,0))
    if(neighbour==2)then
       cnt=0
       do ia=1,Natoms
          if(sqrt(dot_product(r(ia_c,:)-r(ia,:),r(ia_c,:)-r(ia,:))) > d1+tol .and. &
               sqrt(dot_product(r(ia_c,:)-r(ia,:),r(ia_c,:)-r(ia,:)))<= d2+tol)then
             cnt=cnt+1
             NNindex(cnt)=ia
          endif
       enddo
       if(IO) write(*,*) 'Number of 2nd NN (excluding NN):',cnt
    endif

    if(neighbour==3)then
       cnt=0
       do ia=1,Natoms
          if(sqrt(dot_product(r(ia_c,:)-r(ia,:),r(ia_c,:)-r(ia,:))) > d2+tol .and. &
               sqrt(dot_product(r(ia_c,:)-r(ia,:),r(ia_c,:)-r(ia,:)))<= d3+tol)then
             cnt=cnt+1
             NNindex(cnt)=ia
          endif
       enddo
       if(IO) write(*,*) 'Number of 3rd NN (excluding NN+2ndNN):',cnt
    endif

    !Number of nearest neighbours + central atom
    NumNN=cnt



    !Reallocate with correct size
    call reallocate(NNindex,NumNN)

    !Units. e**2/(4 pi eps_0 a_0) = 1 Ha
    !Unit converts to eV
    unt=14.399644849251


    !Check NNindex works
    !do ia=1,NumNN
    !   if(IO) write(100,*) r(NNindex(ia),:)
    !enddo


    !------------------------------------
    !Compute orbital integrals 
    !------------------------------------

    !Monopole approx
    if(IO)write(*,*) 'Monopole Approximation (eV): ',unt/d_val
    
    !Volume element
    dV=rspacing(1,1)*rspacing(2,2)*rspacing(3,3)

    !Generate 3D orbital profiles
    call orbital_profile(rspacing,midpnt,sorb,pxorb,pyorb,pzorb,sstorb)

    !Orbitals 1 and 2, filled per iteration. Dimensions of the  NN grid
    allocate(orbital1(npts_r(1),npts_r(2),npts_r(3)),orbital2(npts_r(1),npts_r(2),npts_r(3)))

    allocate(r1(d),r2(d),Uc(orb_bas,orb_bas))
    Uc(:,:)=0._dp

    !Average orbital integrals (averaged over NumNN interactions)
    allocate(Vss(NumNN),Vsp(NumNN*6),Vpp(NumNN*9),Vsstp(NumNN*6))
    Vss(:)=0._dp
    Vsp(:)=0._dp
    Vpp(:)=0._dp
    Vsstp(:)=0._dp
    cnt1=0
    cnt2=0
    cnt3=0
    cnt4=0

    !Average values for integral elements
    !Vss=0._dp
    !Vsp=0._dp
    !Vpp=0._dp
    !Vsstp=0._dp

    !Loop over central anion (or cation) and it's set of 1st, 2nd or 3rd NN
    do ia=ia_c,ia_c
       do iaNN=1,NumNN
          ja=NNindex(iaNN)
          
         !if(IO)write(*,'(4(I3,X))') ia,ja,atom_type(ia),atom_type(ja)

          do iorb=1,orb_bas
             do jorb=1,orb_bas
                
                !i=iorb+(ia-1)*orb_bas
                !j=jorb+(ja-1)*orb_bas
                !Skip lower triangle elements 
                !if(i>j) cycle

                !Place appropriate orbitals on grid: orbital1 and orbital2 returned.
                !x0-zf are not utilised 
                call place_orbitals(orbital1,orbital2,atm_pntr,iptr,&
                                    ia,ja,iorb,jorb,midpnt,sorb,pxorb,pyorb,pzorb,sstorb,&
                                    x0,y0,z0,xf,yf,zf,atom_species) 



                !--------------------------
                !Integrate over r1 and r2
                !--------------------------
                !Loop order changed to make array access more efficient
                summ=0._dp      

                !r1 loops
                ir=0
                do iz=1,npts_r(3)
                   do iy=1,npts_r(2)
                      do ix=1,npts_r(1)
                         ir=ir+1
                         r1(:)=(ix-1)*rspacing(1,:) + (iy-1)*rspacing(2,:) +(iz-1)*rspacing(3,:)
                         
                         !If one of the two orbitals is zero, the integral is zero
                         !hence cycle (true more often than not)
                         O1=real(orbital1(ix,iy,iz))
                         if(O1==0._dp) cycle
                         
                         !r2 loops
                         jr=0
                         do jz=1,npts_r(3)
                            do jy=1,npts_r(2)
                               do jx=1,npts_r(1)
                                  jr=jr+1
                                  r2(:)=(jx-1)*rspacing(1,:) + (jy-1)*rspacing(2,:) +(jz-1)*rspacing(3,:)
                                  !Denominator 1/|r1-r2|
                                  denom=sqrt(dot_product(r1(:)-r2(:),r1(:)-r2(:)))
                                  !Truncate singularity to sensible value
                                  if(ir==jr) denom=sqrt(0.75)*rspacing(1,1)
                                  
                                  !Sum integral
                                  !Coulomb 
                                  summ=summ+ O1*O1*&
                                           real(orbital2(jx,jy,jz))*real(orbital2(jx,jy,jz))/denom
                               
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
                
                !Output integral values (in eV) 
                !if(IO) write(*,'(4(I2,X),F12.8)') ia,ja,iorb,jorb, unt*summ*dV*dV
                !Sum averages
                !if(iorb==1 .and. jorb==1)             Vss=Vss+unt*summ*dV*dV
                !if(iorb==1 .and. (jorb>1.and.jorb<5)) Vsp=Vsp+unt*summ*dV*dV     !Technically Vps
                !if(jorb==1 .and. (iorb>1.and.iorb<5)) Vsp=Vsp+unt*summ*dV*dV
                !if(iorb>1  .and. iorb<5 .and. jorb>1 .and. jorb<5) Vpp=Vpp+unt*summ*dV*dV
                !if(iorb==5 .and. (jorb>1.and.jorb<5)) Vsstp=Vsstp+unt*summ*dV*dV  
                !if(jorb==5 .and. (iorb>1.and.iorb<5)) Vsstp=Vsstp+unt*summ*dV*dV !Technically Vpsst
       

                !Store every instance of each interaction in arrays, allowing one to
                !take min, max and avg
                !Vss
                if(iorb==1 .and. jorb==1)then
                   cnt1=cnt1+1
                   Vss(cnt1)=Vss(cnt1)+unt*summ*dV*dV
                endif
                !Technically Vps. Could create a different array if required
                if(iorb==1 .and. (jorb>1.and.jorb<5))then
                   cnt2=cnt2+1
                   Vsp(cnt2)=Vsp(cnt2)+unt*summ*dV*dV     
                endif
                !Vsp
                if(jorb==1 .and. (iorb>1.and.iorb<5))then
                   cnt2=cnt2+1
                   Vsp(cnt2)=Vsp(cnt2)+unt*summ*dV*dV
                endif
                !Vpp
                if(iorb>1  .and. iorb<5 .and. jorb>1 .and. jorb<5)then
                   cnt3=cnt3+1
                   Vpp(cnt3)=Vpp(cnt3)+unt*summ*dV*dV
                endif
                !Vsstp
                if(iorb==5 .and. (jorb>1.and.jorb<5))then
                   cnt4=cnt4+1
                   Vsstp(cnt4)=Vsstp(cnt4)+unt*summ*dV*dV 
                endif
                !Vpsst
                if(jorb==5 .and. (iorb>1.and.iorb<5))then
                   cnt4=cnt4+1
                   Vsstp(cnt4)=Vsstp(cnt4)+unt*summ*dV*dV 
                endif

                !Temporarily store integral values in a local matrix
                Uc(iorb,jorb)=unt*summ*dV*dV

             enddo
          enddo
          
          !Output local matrix with dynamic name, for atoms (ia,ja)
          write(ia_lab,'(I2)') ia 
          write(ja_lab,'(I2)') ja 
          fname=trim(adjustl(ia_lab))//'_'//trim(adjustl(ja_lab))//'.dat'
          if(IO) open(unit=003, file=trim(adjustl(fname)))
          do iorb=1,orb_bas
             if(IO) write(003,'(5(F12.6,X))') Uc(iorb,1:orb_bas)
          enddo
          if(IO)close(003)
          Uc(:,:)=0._dp

       enddo
    enddo

    !Output values for 2-centre integrals between 2nd or 3rd NN
    if(IO)then
       write(*,*) 'Average 2-centre integral elements:'
       write(*,*) 'min,max,avg (eV)'
       write(*,*) 'Vss:',minval(Vss),maxval(Vss),sum(Vss)/dble(size(Vss))
       write(*,*) 'Vsp:',minval(Vsp),maxval(Vsp),sum(Vsp)/dble(size(Vsp))
       write(*,*) 'Vpp:',minval(Vpp),maxval(Vpp),sum(Vpp)/dble(size(Vpp))
       write(*,*) 'Vsstp:',minval(Vsstp),maxval(Vsstp),sum(Vsstp)/dble(size(Vsstp))
    endif

    !Output averages for 2-centre integrals between 2nd or 3rd NN
    !if(IO)then
    !   write(*,*) 'Average 2-centre integral elements:'
    !   write(*,*) 'Vss, Vsp, Vpp, Vsstp'
    !   write(*,*) Vss/dble(NumNN),Vsp/dble(NumNN*6.),Vpp/dble(NumNN*9.),Vsstp/dble(NumNN*6.)
       !Test on fewer terms
    !   !write(*,*) Vss/dble(4.),Vsp/dble(4*6.),Vpp/dble(4*9.),Vsstp/dble(4*6.)
    !endif

    call cpu_time(t2)
    if(IO) write(*,*)'2-centre integral routine took: ',(t2-t1)/60._dp,'mins'
  end subroutine orbital_integrals2



  
  !****************************************************************************
  !Read in and store on-site and NN Coulomb and exchange energies for use
  !in Bethe-Salpeter calculations. All integrals computed in the two-centre
  !approximation using a real basis of single-zeta functions from Siesta
  !****************************************************************************
  subroutine orbital_integral_energies(U_c,U_x)
    use var,     only: dp,ns,orb_bas
    implicit none

    !********************
    !Declarations
    !********************

    !Modular Arrays
    real(dp), allocatable :: U_c(:,:),U_x(:,:)

    !Local Variables
    integer   ::  N,Nt,i,j,k,ierr
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
    N=ns*orb_bas

    !Number of elements in upper triangle of matrix
    Nt=int(0.5*N*(N+1))

    !Allocate arrays. Both are symmetric but easier (and faster?)
    !implementation to fill the whole thing
    allocate(U_c(N,N),U_x(N,N))
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
    call mpi_bcast(U_c,N*N,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIierr)
    call mpi_bcast(U_x,N*N,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIierr)


    !Return U_c and U_x
  end subroutine orbital_integral_energies




  !Same as above, but in serial 
  subroutine serial_orth_construct_H2p(H2p_s)
    use var, only: dp,pi,Natoms,orb_bas,eigen_cnt,r,Ne,Nh,ic,evcr
    implicit none

    !******************
    !Declarations
    !******************

    !Global Arrays
    real(dp),     allocatable,        intent(inout) :: H2p_s(:,:)

    !Local Variables
    integer  :: N,VBT,CBB,ia,ja,iorb,jorb,i,j,Hcnt,row,col,h_row_cnt,h_col_cnt,e_row_cnt,&
                e_col_cnt
    real(dp) :: t1_local,t2_local,V_an,V,r_sph,unt,ieps

    !Local Arrays
    real(dp), allocatable :: c_vi(:),c_ci(:),c_vj(:),c_cj(:)


    !******************
    !Main Routine
    !******************
    if(IO) write(*,*)'Construct_H2p in orthogonal approximation'
    call cpu_time(t1_local)

    !Dimensions of the system (==total number of orbital states)
    N=Natoms*orb_bas

    !VBT
    VBT=eigen_cnt(1,1)
    !CBB
    CBB=VBT+1
    
    !Coulomb potential between atomic orbitals constructed on-the-fly
    !No need to set up a system grid, due to expressions resulting from
    !orthogonality between orbitals located on different atomic sites

    !Units. e**2/(4 pi eps_0 a_0) = 1 Ha
    !Convert atomic positions to Bohr
    !Convert Hartree to eV  
    !unt=(0.5291772109217*27.2113850560)
    unt=14.399644849251

    !Radius of cubic volume element, given as the length from 
    !the centre to the corner (see Brillouin approximation)
    !DUMMY VALUE
    r_sph=1._dp 
    
    !Define the analytic solution to the potential at ia==ja
    !unt set so r_sph can be given in Angstroms and V_an is in eV
    !DUMMY VALUE
    V_an=16._dp
    if(IO) write(*,*)'Using a dummy value for the onsite potential value'

    !NEED TO DIVIDE BY VOL OF UNIT CELL SQUARED
    !V_an=unt*2.*pi*r_sph*r_sph

    !Bulk, static dielectric constant for CdSe
    if(ic==19) ieps=9.7
    if(IO .and. ic/=19) write(*,*)'No dielectric constant stored for this material'

    !1D arrays to hold TB coefficients
    allocate(c_vi(Nh),c_vj(Nh),c_ci(Ne),c_cj(Ne))
    c_vi(:)=0._dp
    c_vj(:)=0._dp
    c_ci(:)=0._dp
    c_cj(:)=0._dp


    do ia=1,Natoms
       do ja=1,Natoms

          !On-site potential
          if(ia==ja) V=V_an
          !Inter-atomic potential 1/|r_ia-r_ja|
          if(ia/=ja) V=unt/(sqrt(dot_product(r(ia,:)-r(ja,:),r(ia,:)-r(ja,:))))

          do iorb=1,orb_bas
             !Global orbital index
             i=iorb+(ia-1)*orb_bas
             !Assign TB coefficients to 1D arrays 
             c_vi(:) = evcr(i, VBT:VBT-Nh+1:-1) 
             c_ci(:) = evcr(i, VBT+1:VBT+Ne)    

             do jorb=1,orb_bas
                j=jorb+(ja-1)*orb_bas
                c_vj(:) = evcr(j, VBT:VBT-Nh+1:-1) 
                c_cj(:) = evcr(j, VBT+1:VBT+Ne)    

                !Row loop distributed => H2p is distributed
                do row=1,Ne*Nh
                   
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
                      H2p_s(row,col) = H2p_s(row,col) + c_vi(h_row_cnt)*c_cj(e_col_cnt)*&
                                       (2._dp*c_ci(e_row_cnt)*c_vj(h_col_cnt) -    &
                                       ieps*c_vi(h_col_cnt)*c_cj(e_row_cnt))*V

                   enddo
                enddo
                !Shut loops over (vc)(v'c')
                
             enddo
          enddo
       enddo
    enddo



    !Output routine runtime
    call cpu_time(t2_local)
    write(*,'(A,F20.10,A)') 'Time taken to construct H2p in serial =  ',(t2_local-t1_local)/60._dp,' mins'

  end subroutine serial_orth_construct_H2p



  !****************************************************************
  !Double-counting Compensation Factor
  !Modified from Hannes Hubener's routine. 
  !Oxford University. 2010-2014
  !****************************************************************
  !Orbital integral matrix has 8-fold symmetry
  !One can therefore loop over 1/8 of the orbital products
  !but needs to account for these permutations in the orbital
  !coefficients which multiply the integral.
  !In some instances, coefficients will be equivalent and a 
  !factor needs to be introduced to compensate for double counting
  !****************************************************************
  !Input:
  !i,j,k,l      Four unique global orbital indices 
  !             (n1,n2,n3,n4 in orbital integral routine)
  !Output:
  !factor     Double-counting factor
  !****************************************************************
 

  !Function for double-counting
  real(dp) function dbfactor(i,j,k,l)
    implicit none

    !***************
    !Declarations
    !***************
    integer, intent(in)  :: i,j,k,l
    real(dp)             :: f1,f2,f3

    !***************
    !Main Routine
    !***************

    !(i/=j/=k/=l), (i=k,j/=l),(i/=k,j=l),(i=l,j/=k),(i/=l,j=k) 
    ! all have 8 unique coefficient products

    !Initialise factors
    f1=1._dp
    f2=1._dp
    f3=1._dp

        
    !Logical statements
    if(i==j)then                       
       f1=0.5                          !i==j,  i==j==k /=l   and  i==j==l /=k  1/2
       if(k==l)then                    
          f2=0.5                       !i==j  k==l    1/4
          if(i==k) f3=0.5              !i==j==k==l    1/8
       endif

    else                               
       if(k==l)then                   
          f2=0.5                       !k==l,   i==k==l /=j    and    j==k==l /=i   1/2
       else                                    
          if(i==k .and. j==l) f3=0.5   !i==k  j==l   1/2
          if(i==l .and. j==k) f3=0.5   !i==l  j==k   1/2
       endif
    endif

    !Double counting factor
    dbfactor=f1*f2*f3

    return
  end function dbfactor




!!$  real(dp) function dbfactor(i,j,k,l)
!!$    use var,only:dp
!!$    implicit none
!!$
!!$    !***************
!!$    !Declarations
!!$    !***************
!!$    integer, intent(in)  :: i,j,k,l
!!$    real(dp)             :: f1,f2,f3
!!$
!!$    !***************
!!$    !Main Routine
!!$    !***************
!!$
!!$    !(i/=j/=k/=l), (i=k,j/=l),(i/=k,j=l),(i=l,j/=k),(i/=l,j=k) 
!!$    ! all have 8 unique coefficient products
!!$
!!$    !Initialise factors
!!$    f1=1._dp
!!$    f2=1._dp
!!$    f3=1._dp
!!$
!!$        
!!$    !Logical statements
!!$
!!$    if(i==j)then                       
!!$       f1=0.5                          !i==j   only
!!$       if(k==l)then                    
!!$          f2=0.5                       !i==j  k==l    k/=i and l/=j  
!!$          if(i==k) f3=0.5              !i==j==k==l
!!$       else
!!$          if(i==k) f3=0.5              !i==j==k /=l   Shouldn't these be total 1/2 not 1/4?
!!$          if(i==l) f3=0.5              !i==j==l /=k   Shouldn't these be total 1/2 not 1/4?
!!$       endif
!!$
!!$    else                               !i/=j   only
!!$       if(k==l)then                   
!!$          f2=0.5                       !i/=j  k==l  
!!$          if(k==i) f3=0.5              !i==k==l /=j  Shouldn't these be total 1/2 not 1/4?
!!$          if(k==j) f3=0.5              !j==k==l /=i  Shouldn't these be total 1/2 not 1/4?
!!$       else                                    
!!$          if(i==k .and. j==l) f3=0.5   !i==k  j==l
!!$          if(i==l .and. j==k) f3=0.5   !i==l  j==k
!!$       endif
!!$    endif
!!$
!!$    !Double counting factor
!!$    dbfactor=f1*f2*f3
!!$
!!$    return
!!$  end function dbfactor



  !End of module
end module exbnd
