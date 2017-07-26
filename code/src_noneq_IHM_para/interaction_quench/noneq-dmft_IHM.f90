!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  date:   February 28, 2013
!
!----------------------------------------------------------------------------
module NONEQ_DMFT_MOD
  
  use CONST_MOD
  use PARM_MOD
  use DMFT_MOD
  use GREEN_MOD
  use IPT_MOD
  use INTEGRAL_MOD
  implicit none
  
  private
  public :: start_noneq_dmft
  
contains
  
  subroutine start_noneq_dmft(parm_,dmft_)
    type(parm), intent(in)        :: parm_
    type(dmft), intent(inout)     :: dmft_
    integer                       :: n ! time step
    integer                       :: i, j
    
    
    do n=2, parm_%N_t+1 ! maximum time is t=t_{n}=(n-1)*dt
       dmft_%converged=.false.
       dmft_%G0_diffA=1.d0
       dmft_%G0_diffB=1.d0
       dmft_%iteration=1
       dmft_%time_step=n
       call initialize_noneq_greenB(parm_,dmft_)
       call initialize_noneq_greenA(parm_,dmft_) !TODO make sure h factor is implimented correctly
       !call initialize_noneq_greenB(parm_,dmft_)
       write (*,'(A,f15.10)') 't =', dble(n-1)*parm_%dt
       
       ! nonequilibrium DMFT self-consistency loop
       do while(.not. dmft_%converged .and. dmft_%iteration<=parm_%N_iter)
          call noneq_impurity_solutionA(parm_,dmft_)
          call noneq_dmft_self_consistencyB(parm_,dmft_)
	  call noneq_impurity_solutionB(parm_,dmft_)
          call noneq_dmft_self_consistencyA(parm_,dmft_)
          dmft_%iteration=dmft_%iteration+1
          if (dmft_%G0_diffA<=parm_%tolerance .and. dmft_%G0_diffB<=parm_%tolerance) then
             dmft_%converged=.true.
             write (*,'(A)') 'Nonequilibrium DMFT is converged'
          end if
          if (dmft_%G0_diffA>parm_%tolerance .and. dmft_%G0_diffB>parm_%tolerance .and. dmft_%iteration>parm_%N_iter) then
             write (*,'(A)') 'Nonequilibrium DMFT is NOT converged !!!!!'
          end if
       end do
       
       ! d/dt G0_{old} <= d/dt G0_{new}
       call deallocate_KB_derivative(dmft_%weiss_green_derA)
       call allocate_KB_derivative(parm_,dmft_%weiss_green_derA,n)
       call deallocate_KB_derivative(dmft_%weiss_green_derB)
       call allocate_KB_derivative(parm_,dmft_%weiss_green_derB,n)

       do i=1, n
          dmft_%weiss_green_derA%retarded(i)=dmft_%weiss_green_der_newA%retarded(i)
          dmft_%weiss_green_derA%lesser(i)=dmft_%weiss_green_der_newA%lesser(i)
	  dmft_%weiss_green_derB%retarded(i)=dmft_%weiss_green_der_newB%retarded(i)
          dmft_%weiss_green_derB%lesser(i)=dmft_%weiss_green_der_newB%lesser(i)
       end do
       do j=1, parm_%N_tau+1
          dmft_%weiss_green_derA%left_mixing(j)=dmft_%weiss_green_der_newA%left_mixing(j)
	  dmft_%weiss_green_derB%left_mixing(j)=dmft_%weiss_green_der_newB%left_mixing(j)
       end do
       call deallocate_KB_derivative(dmft_%weiss_green_der_newA)
       call allocate_KB_derivative(parm_,dmft_%weiss_green_der_newA,n+1)
       call deallocate_KB_derivative(dmft_%weiss_green_der_newB)
       call allocate_KB_derivative(parm_,dmft_%weiss_green_der_newB,n+1)
       
       call measure_density(parm_,dmft_)
       call measure_double_occupancy(parm_,dmft_)
       call measure_kinetic_energy(parm_,dmft_)
       call measure_kinetic_energy2(parm_,dmft_)
       call measure_interaction_energy(parm_,dmft_)
       call measure_total_energy(parm_,dmft_)
       call print_green_reatarded(parm_,dmft_)
    end do
    
  end subroutine start_noneq_dmft
  
!************************ HERE WE HAVE OUR IMPURITY SOLVER ***************************************************    
  subroutine noneq_impurity_solutionA(parm_,dmft_)
    type(parm), intent(in)       :: parm_
    type(dmft), intent(inout)    :: dmft_
    
    if (parm_%solver=='IPT') then
       call noneq_iptA(parm_,dmft_)
    end if
    
  end subroutine noneq_impurity_solutionA

  subroutine noneq_impurity_solutionB(parm_,dmft_)
    type(parm), intent(in)       :: parm_
    type(dmft), intent(inout)    :: dmft_
    
    if (parm_%solver=='IPT') then
       call noneq_iptB(parm_,dmft_)
    end if
    
  end subroutine noneq_impurity_solutionB
  
!****************************** NONEQ SELF CONSISTENCY CALCULATION WAS DONE HERE *********************************************     
  subroutine noneq_dmft_self_consistencyA(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    type(KB)                        :: kernel
    integer                         :: i, j, k, n
    complex(kind(0d0)), allocatable :: h(:)
    
    n=dmft_%time_step
    
    if (parm_%dos=='semicircular') then
       ! solve the lattice Dyson equation: id/dt G0(t,t')-G*G0(t,t')=delta(t,t')
       allocate(h(n))
       h(:)=dmft_%hA
       call Volterra_intdiff(parm_,n,h,dmft_%local_greenB,dmft_%weiss_green_newA,dmft_%weiss_green_derA,dmft_%weiss_green_der_newA)
       deallocate(h)
    end if
    
    ! evaluate |G0_{new}-G0_{old}|
    dmft_%G0_diffA=0.d0
    do j=1, n
       dmft_%G0_diffA=dmft_%G0_diffA+abs(dmft_%weiss_green_newA%retarded(n,j)-dmft_%weiss_greenA%retarded(n,j))
       dmft_%G0_diffA=dmft_%G0_diffA+abs(dmft_%weiss_green_newA%lesser(n,j)-dmft_%weiss_greenA%lesser(n,j))
    end do
    do j=1, parm_%N_tau+1
       dmft_%G0_diffA=dmft_%G0_diffA+abs(dmft_%weiss_green_newA%left_mixing(n,j)-dmft_%weiss_greenA%left_mixing(n,j))
    end do
    do i=1, n-1
       dmft_%G0_diffA=dmft_%G0_diffA+abs(dmft_%weiss_green_newA%lesser(i,n)-dmft_%weiss_greenA%lesser(i,n))
    end do
    write (*,'(A,I2)') '  Iteration #', dmft_%iteration
    write (*,'(A,f15.10)') ' |G0_new-G0_old| =', dmft_%G0_diffA
    
    ! G0_{old} <= G0_{new}
    do j=1, n
       dmft_%weiss_greenA%retarded(n,j)=dmft_%weiss_green_newA%retarded(n,j)
       dmft_%weiss_greenA%lesser(n,j)=dmft_%weiss_green_newA%lesser(n,j)
    end do
    do i=1, n-1
       dmft_%weiss_greenA%lesser(i,n)=dmft_%weiss_green_newA%lesser(i,n)
    end do
    dmft_%weiss_greenA%left_mixing(n,:)=dmft_%weiss_green_newA%left_mixing(n,:)
    
    ! Kadanoff-Baym G0 => contour-ordered G0
    do j=1, n
       dmft_%weiss_green_TA%c12(n,j)=dmft_%weiss_greenA%lesser(n,j)
       dmft_%weiss_green_TA%c21(n,j)=dmft_%weiss_greenA%lesser(n,j)+dmft_%weiss_greenA%retarded(n,j)
    end do
    do i=1, n-1
       dmft_%weiss_green_TA%c12(i,n)=dmft_%weiss_greenA%lesser(i,n)
       dmft_%weiss_green_TA%c21(i,n)=dmft_%weiss_greenA%lesser(i,n)-conjg(dmft_%weiss_greenA%retarded(n,i))
    end do
    do j=1, parm_%N_tau+1
       dmft_%weiss_green_TA%c13(n,j)=dmft_%weiss_greenA%left_mixing(n,j)
    end do
    do i=1, parm_%N_tau+1
       dmft_%weiss_green_TA%c31(i,n)=conjg(dmft_%weiss_greenA%left_mixing(n,parm_%N_tau-i+2))
    end do
    
    ! Hermite conjugate
    do j=1, n
       dmft_%weiss_green_TA%c12(n,j)=0.5d0*(dmft_%weiss_green_TA%c12(n,j)-conjg(dmft_%weiss_green_TA%c12(j,n)))
       dmft_%weiss_green_TA%c21(n,j)=0.5d0*(dmft_%weiss_green_TA%c21(n,j)-conjg(dmft_%weiss_green_TA%c21(j,n)))
    end do
    do i=1, n-1
       dmft_%weiss_green_TA%c12(i,n)=-conjg(dmft_%weiss_green_TA%c12(n,i))
       dmft_%weiss_green_TA%c21(i,n)=-conjg(dmft_%weiss_green_TA%c21(n,i))
    end do
    do j=1, parm_%N_tau+1
       dmft_%weiss_green_TA%c13(n,j)=0.5d0*(dmft_%weiss_green_TA%c13(n,j)+conjg(dmft_%weiss_green_TA%c31(parm_%N_tau-j+2,n)))
    end do
    do i=1, parm_%N_tau+1
       dmft_%weiss_green_TA%c31(i,n)=conjg(dmft_%weiss_green_TA%c13(n,parm_%N_tau-i+2))
    end do
    
  end subroutine noneq_dmft_self_consistencyA

	!////////////////////////////////////////

  subroutine noneq_dmft_self_consistencyB(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    type(KB)                        :: kernel
    integer                         :: i, j, k, n
    complex(kind(0d0)), allocatable :: h(:)
    
    n=dmft_%time_step
    
    if (parm_%dos=='semicircular') then
       ! solve the lattice Dyson equation: id/dt G0(t,t')-G*G0(t,t')=delta(t,t')
       allocate(h(n))
       h(:)= dmft_%hB
       call Volterra_intdiff(parm_,n,h,dmft_%local_greenA,dmft_%weiss_green_newB,dmft_%weiss_green_derB,dmft_%weiss_green_der_newB) !TODO do sothing abt h
       deallocate(h)
    end if
    
    ! evaluate |G0_{new}-G0_{old}|
    dmft_%G0_diffB=0.d0
    do j=1, n
       dmft_%G0_diffB=dmft_%G0_diffB+abs(dmft_%weiss_green_newB%retarded(n,j)-dmft_%weiss_greenB%retarded(n,j))
       dmft_%G0_diffB=dmft_%G0_diffB+abs(dmft_%weiss_green_newB%lesser(n,j)-dmft_%weiss_greenB%lesser(n,j))
    end do
    do j=1, parm_%N_tau+1
       dmft_%G0_diffB=dmft_%G0_diffB+abs(dmft_%weiss_green_newB%left_mixing(n,j)-dmft_%weiss_greenB%left_mixing(n,j))
    end do
    do i=1, n-1
       dmft_%G0_diffB=dmft_%G0_diffB+abs(dmft_%weiss_green_newB%lesser(i,n)-dmft_%weiss_greenB%lesser(i,n))
    end do
    write (*,'(A,I2)') '  Iteration #', dmft_%iteration
    write (*,'(A,f15.10)') ' |G0_new-G0_old| =', dmft_%G0_diffB
    
    ! G0_{old} <= G0_{new}
    do j=1, n
       dmft_%weiss_greenB%retarded(n,j)=dmft_%weiss_green_newB%retarded(n,j)
       dmft_%weiss_greenB%lesser(n,j)=dmft_%weiss_green_newB%lesser(n,j)
    end do
    do i=1, n-1
       dmft_%weiss_greenB%lesser(i,n)=dmft_%weiss_green_newB%lesser(i,n)
    end do
    dmft_%weiss_greenB%left_mixing(n,:)=dmft_%weiss_green_newB%left_mixing(n,:)
    
    ! Kadanoff-Baym G0 => contour-ordered G0
    do j=1, n
       dmft_%weiss_green_TB%c12(n,j)=dmft_%weiss_greenB%lesser(n,j)
       dmft_%weiss_green_TB%c21(n,j)=dmft_%weiss_greenB%lesser(n,j)+dmft_%weiss_greenB%retarded(n,j)
    end do
    do i=1, n-1
       dmft_%weiss_green_TB%c12(i,n)=dmft_%weiss_greenB%lesser(i,n)
       dmft_%weiss_green_TB%c21(i,n)=dmft_%weiss_greenB%lesser(i,n)-conjg(dmft_%weiss_greenB%retarded(n,i))
    end do
    do j=1, parm_%N_tau+1
       dmft_%weiss_green_TB%c13(n,j)=dmft_%weiss_greenB%left_mixing(n,j)

    end do
    do i=1, parm_%N_tau+1
       dmft_%weiss_green_TB%c31(i,n)=conjg(dmft_%weiss_greenB%left_mixing(n,parm_%N_tau-i+2))
    end do
    
    ! Hermite conjugate
    do j=1, n
       dmft_%weiss_green_TB%c12(n,j)=0.5d0*(dmft_%weiss_green_TB%c12(n,j)-conjg(dmft_%weiss_green_TB%c12(j,n)))
       dmft_%weiss_green_TB%c21(n,j)=0.5d0*(dmft_%weiss_green_TB%c21(n,j)-conjg(dmft_%weiss_green_TB%c21(j,n)))
    end do
    do i=1, n-1
       dmft_%weiss_green_TB%c12(i,n)=-conjg(dmft_%weiss_green_TB%c12(n,i))
       dmft_%weiss_green_TB%c21(i,n)=-conjg(dmft_%weiss_green_TB%c21(n,i))
    end do
    do j=1, parm_%N_tau+1
       dmft_%weiss_green_TB%c13(n,j)=0.5d0*(dmft_%weiss_green_TB%c13(n,j)+conjg(dmft_%weiss_green_TB%c31(parm_%N_tau-j+2,n)))
    end do
    do i=1, parm_%N_tau+1
       dmft_%weiss_green_TB%c31(i,n)=conjg(dmft_%weiss_green_TB%c13(n,parm_%N_tau-i+2))
    end do
    
  end subroutine noneq_dmft_self_consistencyB
  
!****************************** INITIALISED NONEQ WEISS GREEN'S FUNCTION *********************************************     
  subroutine initialize_noneq_greenA(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: i, j, k, n
    complex(kind(0d0)), allocatable :: SxG(:)
    
    n=dmft_%time_step
    ! initialize G0(t,t') and G(t,t')
    if (n==2) then
       dmft_%weiss_greenA%retarded(1,1)=-xj
       dmft_%local_greenA%retarded(1,1)=-xj
       dmft_%weiss_green_newA%retarded(1,1)=-xj
       do j=1, parm_%N_tau+1
          dmft_%weiss_greenA%left_mixing(1,j)=-xj*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-j+2)
          dmft_%local_greenA%left_mixing(1,j)=-xj*dmft_%local_greenA%matsubara_t(parm_%N_tau-j+2)
          dmft_%local_greenB%left_mixing(1,j)=-xj*dmft_%local_greenB%matsubara_t(parm_%N_tau-j+2)
          dmft_%weiss_green_newA%left_mixing(1,j)=-xj*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-j+2)
       end do
       dmft_%weiss_greenA%lesser(1,1)=-xj*dmft_%weiss_greenA%matsubara_t(parm_%N_tau+1)
       dmft_%local_greenA%lesser(1,1)=-xj*dmft_%local_greenA%matsubara_t(parm_%N_tau+1)
       dmft_%weiss_green_newA%lesser(1,1)=-xj*dmft_%weiss_greenA%matsubara_t(parm_%N_tau+1)
       
       if (parm_%dos=='semicircular') then
          ! initialize d/dt G0(t,t') = -ihA(t)G0(t,t') -i*(G*G0)(t,t')
          dmft_%weiss_green_derA%retarded(1)= -xj*dmft_%hA*dmft_%weiss_greenA%retarded(1,1)
          dmft_%weiss_green_derA%left_mixing(:)=0.d0
          allocate(SxG(parm_%N_tau+1))
          do j=1, parm_%N_tau+1
             do k=1, j
                SxG(k)=dmft_%local_greenB%left_mixing(1,k)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau+k-j+1)
             end do
             dmft_%weiss_green_derA%left_mixing(j)=dmft_%weiss_green_derA%left_mixing(j)+xj*parm_%dtau*trapezoid_z(SxG,1,j)
             do k=j, parm_%N_tau+1
                SxG(k)=dmft_%local_greenB%left_mixing(1,k)*dmft_%weiss_greenA%matsubara_t(k-j+1)
             end do
             dmft_%weiss_green_derA%left_mixing(j)=dmft_%weiss_green_derA%left_mixing(j)-xj*parm_%dtau*trapezoid_z(SxG,j,parm_%N_tau+1) -xj*dmft_%hA*dmft_%weiss_greenA%left_mixing(1,j)
          end do
          do k=1, parm_%N_tau+1
             SxG(k)=dmft_%local_greenB%left_mixing(1,k)*conjg(dmft_%weiss_greenA%left_mixing(1,parm_%N_tau-k+2)) 
          end do
          dmft_%weiss_green_derA%lesser(1)=-xj*(-xj)*parm_%dtau*trapezoid_z(SxG,1,parm_%N_tau+1) -xj*dmft_%hA*dmft_%weiss_greenA%lesser(1,1)
          deallocate(SxG)
       end if
       
       ! guess G0(t,t') in the next step
       do j=1, n
          dmft_%weiss_greenA%retarded(n,j)=dmft_%weiss_greenA%retarded(1,1)
       end do
       do j=1, parm_%N_tau+1
          dmft_%weiss_greenA%left_mixing(n,j)=dmft_%weiss_greenA%left_mixing(1,j)
       end do
       do j=1, n
          dmft_%weiss_greenA%lesser(n,j)=dmft_%weiss_greenA%lesser(1,1)
       end do
       do i=1, n-1
          dmft_%weiss_greenA%lesser(i,n)=dmft_%weiss_greenA%lesser(1,1)
       end do
       
       ! contour-ordered G0 <= Kadanoff-Baym G0
       do j=1, n
          do i=1, n
             dmft_%weiss_green_TA%c12(i,j)=dmft_%weiss_greenA%lesser(i,j)
             dmft_%weiss_green_TA%c21(i,j)=dmft_%weiss_greenA%lesser(i,j)-xj
          end do
       end do
       do j=1, parm_%N_tau+1
          do i=1, n
             dmft_%weiss_green_TA%c13(i,j)=-xj*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-j+2)
             dmft_%weiss_green_TA%c31(j,i)=xj*dmft_%weiss_greenA%matsubara_t(j)
          end do
       end do
    else if (n>=3) then
       ! guess G0(t,t') in the next step by quadratic extrapolation
       dmft_%weiss_greenA%retarded(n,n)=-xj
       do k=1, n-2
          dmft_%weiss_greenA%retarded(n,k)=2.d0*dmft_%weiss_greenA%retarded(n-1,k)-dmft_%weiss_greenA%retarded(n-2,k)
       end do
       dmft_%weiss_greenA%retarded(n,n-1)=0.5d0*(dmft_%weiss_greenA%retarded(n,n)+dmft_%weiss_greenA%retarded(n,n-2))
       do k=1, parm_%N_tau+1
          dmft_%weiss_greenA%left_mixing(n,k)=2.d0*dmft_%weiss_greenA%left_mixing(n-1,k)-dmft_%weiss_greenA%left_mixing(n-2,k)
       end do
       do k=1, n-1
          dmft_%weiss_greenA%lesser(n,k)=2.d0*dmft_%weiss_greenA%lesser(n-1,k)-dmft_%weiss_greenA%lesser(n-2,k)
          dmft_%weiss_greenA%lesser(k,n)=2.d0*dmft_%weiss_greenA%lesser(k,n-1)-dmft_%weiss_greenA%lesser(k,n-2)
       end do
       dmft_%weiss_greenA%lesser(n,n)=2.d0*dmft_%weiss_greenA%lesser(n-1,n-1)-dmft_%weiss_greenA%lesser(n-2,n-2)
       
       ! contour-ordered G0 <= Kadanoff-Baym G0
       do k=1, n-1
          dmft_%weiss_green_TA%c12(n,k)=dmft_%weiss_greenA%lesser(n,k)
          dmft_%weiss_green_TA%c12(k,n)=dmft_%weiss_greenA%lesser(k,n)
          dmft_%weiss_green_TA%c21(n,k)=dmft_%weiss_greenA%lesser(n,k)+dmft_%weiss_greenA%retarded(n,k)
          dmft_%weiss_green_TA%c21(k,n)=dmft_%weiss_greenA%lesser(k,n)-conjg(dmft_%weiss_greenA%retarded(n,k))
       end do
       dmft_%weiss_green_TA%c12(n,n)=dmft_%weiss_greenA%lesser(n,n)
       dmft_%weiss_green_TA%c21(n,n)=dmft_%weiss_greenA%lesser(n,n)-xj
       do k=1, parm_%N_tau+1
          dmft_%weiss_green_TA%c13(n,k)=dmft_%weiss_greenA%left_mixing(n,k)
          dmft_%weiss_green_TA%c31(k,n)=conjg(dmft_%weiss_greenA%left_mixing(n,parm_%N_tau-k+2))
       end do
    end if
    !TODO need to incorporate h in the initialised weiss_der calculation at t=0
  end subroutine initialize_noneq_greenA
  
!////////////////////////////

  subroutine initialize_noneq_greenB(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: i, j, k, n
    complex(kind(0d0)), allocatable :: SxG(:)
    
    n=dmft_%time_step
    ! initialize G0(t,t') and G(t,t')
    if (n==2) then
       dmft_%weiss_greenB%retarded(1,1)=-xj
       dmft_%local_greenB%retarded(1,1)=-xj
       dmft_%weiss_green_newB%retarded(1,1)=-xj
       do j=1, parm_%N_tau+1
          dmft_%weiss_greenB%left_mixing(1,j)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-j+2)
          dmft_%local_greenB%left_mixing(1,j)=-xj*dmft_%local_greenB%matsubara_t(parm_%N_tau-j+2)
          dmft_%local_greenA%left_mixing(1,j)=-xj*dmft_%local_greenA%matsubara_t(parm_%N_tau-j+2)
          dmft_%weiss_green_newB%left_mixing(1,j)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-j+2)
       end do
       dmft_%weiss_greenB%lesser(1,1)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau+1)
       dmft_%local_greenB%lesser(1,1)=-xj*dmft_%local_greenB%matsubara_t(parm_%N_tau+1)
       dmft_%weiss_green_newB%lesser(1,1)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau+1)
       
       if (parm_%dos=='semicircular') then
          ! initialize d/dt G0(t,t') = -i*(G*G0)(t,t')
          dmft_%weiss_green_derB%retarded(1)=-xj*dmft_%hB*dmft_%weiss_greenB%retarded(1,1)
          dmft_%weiss_green_derB%left_mixing(:)=0.d0
          allocate(SxG(parm_%N_tau+1))
          do j=1, parm_%N_tau+1
             do k=1, j
                SxG(k)=dmft_%local_greenA%left_mixing(1,k)*dmft_%weiss_greenB%matsubara_t(parm_%N_tau+k-j+1)
             end do
             dmft_%weiss_green_derB%left_mixing(j)=dmft_%weiss_green_derB%left_mixing(j)+xj*parm_%dtau*trapezoid_z(SxG,1,j)
             do k=j, parm_%N_tau+1
                SxG(k)=dmft_%local_greenA%left_mixing(1,k)*dmft_%weiss_greenB%matsubara_t(k-j+1)
             end do
             dmft_%weiss_green_derB%left_mixing(j)=dmft_%weiss_green_derB%left_mixing(j)-xj*parm_%dtau*trapezoid_z(SxG,j,parm_%N_tau+1)  -xj*dmft_%hB*dmft_%weiss_greenB%left_mixing(1,j)
          end do
          do k=1, parm_%N_tau+1
             SxG(k)=dmft_%local_greenA%left_mixing(1,k)*conjg(dmft_%weiss_greenB%left_mixing(1,parm_%N_tau-k+2))
          end do
          dmft_%weiss_green_derB%lesser(1)=-xj*(-xj)*parm_%dtau*trapezoid_z(SxG,1,parm_%N_tau+1) -xj*dmft_%hB*dmft_%weiss_greenB%lesser(1,1)
          deallocate(SxG)
       end if
       
       ! guess G0(t,t') in the next step
       do j=1, n
          dmft_%weiss_greenB%retarded(n,j)=dmft_%weiss_greenB%retarded(1,1)
       end do
       do j=1, parm_%N_tau+1
          dmft_%weiss_greenB%left_mixing(n,j)=dmft_%weiss_greenB%left_mixing(1,j)
       end do
       do j=1, n
          dmft_%weiss_greenB%lesser(n,j)=dmft_%weiss_greenB%lesser(1,1)
       end do
       do i=1, n-1
          dmft_%weiss_greenB%lesser(i,n)=dmft_%weiss_greenB%lesser(1,1)
       end do
       
       ! contour-ordered G0 <= Kadanoff-Baym G0
       do j=1, n
          do i=1, n
             dmft_%weiss_green_TB%c12(i,j)=dmft_%weiss_greenB%lesser(i,j)
             dmft_%weiss_green_TB%c21(i,j)=dmft_%weiss_greenB%lesser(i,j)-xj
          end do
       end do
       do j=1, parm_%N_tau+1
          do i=1, n
             dmft_%weiss_green_TB%c13(i,j)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-j+2)
             dmft_%weiss_green_TB%c31(j,i)=xj*dmft_%weiss_greenB%matsubara_t(j)
          end do
       end do
    else if (n>=3) then
       ! guess G0(t,t') in the next step by quadratic extrapolation
       dmft_%weiss_greenB%retarded(n,n)=-xj
       do k=1, n-2
          dmft_%weiss_greenB%retarded(n,k)=2.d0*dmft_%weiss_greenB%retarded(n-1,k)-dmft_%weiss_greenB%retarded(n-2,k)
       end do
       dmft_%weiss_greenB%retarded(n,n-1)=0.5d0*(dmft_%weiss_greenB%retarded(n,n)+dmft_%weiss_greenB%retarded(n,n-2))
       do k=1, parm_%N_tau+1
          dmft_%weiss_greenB%left_mixing(n,k)=2.d0*dmft_%weiss_greenB%left_mixing(n-1,k)-dmft_%weiss_greenB%left_mixing(n-2,k)
       end do

       do k=1, n-1
          dmft_%weiss_greenB%lesser(n,k)=2.d0*dmft_%weiss_greenB%lesser(n-1,k)-dmft_%weiss_greenB%lesser(n-2,k)
          dmft_%weiss_greenB%lesser(k,n)=2.d0*dmft_%weiss_greenB%lesser(k,n-1)-dmft_%weiss_greenB%lesser(k,n-2)
       end do
       dmft_%weiss_greenB%lesser(n,n)=2.d0*dmft_%weiss_greenB%lesser(n-1,n-1)-dmft_%weiss_greenB%lesser(n-2,n-2)
       
       ! contour-ordered G0 <= Kadanoff-Baym G0
       do k=1, n-1
          dmft_%weiss_green_TB%c12(n,k)=dmft_%weiss_greenB%lesser(n,k)
          dmft_%weiss_green_TB%c12(k,n)=dmft_%weiss_greenB%lesser(k,n)
          dmft_%weiss_green_TB%c21(n,k)=dmft_%weiss_greenB%lesser(n,k)+dmft_%weiss_greenB%retarded(n,k)
          dmft_%weiss_green_TB%c21(k,n)=dmft_%weiss_greenB%lesser(k,n)-conjg(dmft_%weiss_greenB%retarded(n,k))
       end do
       dmft_%weiss_green_TB%c12(n,n)=dmft_%weiss_greenB%lesser(n,n)
       dmft_%weiss_green_TB%c21(n,n)=dmft_%weiss_greenB%lesser(n,n)-xj
       do k=1, parm_%N_tau+1
          dmft_%weiss_green_TB%c13(n,k)=dmft_%weiss_greenB%left_mixing(n,k)
          dmft_%weiss_green_TB%c31(k,n)=conjg(dmft_%weiss_greenB%left_mixing(n,parm_%N_tau-k+2))
       end do
    end if
    !TODO need to incorporate h in the initialised weiss_der calculation at t=0
  end subroutine initialize_noneq_greenB

!****************************** MEASURED QUANTITY *********************************************    
  subroutine measure_density(parm_,dmft_)
    !  n(t)=<c^+(t)c(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: n
    double precision                :: t
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    dmft_%densityA(n)=imag(dmft_%local_greenA%lesser(n,n))
    dmft_%densityB(n)=imag(dmft_%local_greenB%lesser(n,n))
    
    open (unit=7, file='densityA', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%densityA(n)
    close (7)

    open (unit=17, file='densityB', status='old', action='write', position='append')
    write (17,'(f15.10,f15.10)') t, dmft_%densityB(n)
    close (17)
  end subroutine measure_density
!//////////////////////////////////////////////////////

  subroutine measure_double_occupancy(parm_,dmft_)
    !  d(t)=<n_up(t)*n_do(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: j, k, n
    double precision                :: t
    complex(kind(0d0)), allocatable :: SxGA(:)
    complex(kind(0d0)), allocatable :: SxGB(:)
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    ! d(t)=n_up(t)*n_do(t)-i/U*[self_energy*G]^{<}(t,t)
    !     =n_up(t)*n_do(t)-i/U*[self_energy^{R}*G^{<}+self_energy^{<}*G^{A}+self_energy^{Left}*G^{Right}](t,t)
    if (U(parm_,t)==0.d0) then
       dmft_%double_occupancyA(n)=dmft_%densityA(n)**2
       dmft_%double_occupancyB(n)=dmft_%densityB(n)**2
    else
       dmft_%double_occupancyA(n)=dmft_%densityA(n)**2
       dmft_%double_occupancyB(n)=dmft_%densityB(n)**2
       allocate(SxGA(max(parm_%N_tau+1,n)))
       allocate(SxGB(max(parm_%N_tau+1,n)))
       do k=1, parm_%N_tau+1
          SxGA(k)=dmft_%self_energyA%left_mixing(n,k)*conjg(dmft_%local_greenA%left_mixing(n,parm_%N_tau-k+2))
	  SxGB(k)=dmft_%self_energyB%left_mixing(n,k)*conjg(dmft_%local_greenB%left_mixing(n,parm_%N_tau-k+2))
       end do
       dmft_%double_occupancyA(n)=dmft_%double_occupancyA(n)+1.d0/U(parm_,t)*parm_%dtau*imag((-xj)*trapezoid_z(SxGA,1,parm_%N_tau+1))
       dmft_%double_occupancyB(n)=dmft_%double_occupancyB(n)+1.d0/U(parm_,t)*parm_%dtau*imag((-xj)*trapezoid_z(SxGB,1,parm_%N_tau+1))
       do k=1, n
          SxGA(k)=dmft_%self_energyA%retarded(n,k)*dmft_%local_greenA%lesser(k,n)
	  SxGB(k)=dmft_%self_energyB%retarded(n,k)*dmft_%local_greenB%lesser(k,n)
       end do
       dmft_%double_occupancyA(n)=dmft_%double_occupancyA(n)+1.d0/U(parm_,t)*parm_%dt*imag(trapezoid_z(SxGA,1,n))
       dmft_%double_occupancyB(n)=dmft_%double_occupancyB(n)+1.d0/U(parm_,t)*parm_%dt*imag(trapezoid_z(SxGB,1,n))
       do k=1, n
          SxGA(k)=dmft_%self_energyA%lesser(n,k)*conjg(dmft_%local_greenA%retarded(n,k))
	  SxGB(k)=dmft_%self_energyB%lesser(n,k)*conjg(dmft_%local_greenB%retarded(n,k))
       end do
       dmft_%double_occupancyA(n)=dmft_%double_occupancyA(n)+1.d0/U(parm_,t)*parm_%dt*imag(trapezoid_z(SxGA,1,n))
       dmft_%double_occupancyB(n)=dmft_%double_occupancyB(n)+1.d0/U(parm_,t)*parm_%dt*imag(trapezoid_z(SxGB,1,n))
       deallocate(SxGA)
       deallocate(SxGB)
    end if
    
    open (unit=7, file='double-occupancyA', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%double_occupancyA(n)
    close (7)
    open (unit=7, file='double-occupancyB', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%double_occupancyB(n)
    close (7)
    
  end subroutine measure_double_occupancy
  
!////////////////////////////////////////////////
!/
  subroutine measure_kinetic_energy(parm_,dmft_)
    !  E_{kin}(t)=2\sum_{k} E(k)<c_{k}^+(t)c_{k}(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: k, n
    double precision                :: t
    complex(kind(0d0)), allocatable :: GxGA(:)
    complex(kind(0d0)), allocatable :: GxGB(:)
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    dmft_%kinetic_energyA(n)=0.d0
    dmft_%kinetic_energyB(n)=0.d0
    
    if (parm_%dos=='semicircular') then
       ! E_{kin}(t)=2*Im[G^{R}*G^{<}+G^{<}*G^{A}+G^{Left}*G^{Right}](t,t)
       allocate(GxGA(max(parm_%N_tau+1,n)))
       allocate(GxGB(max(parm_%N_tau+1,n)))
       do k=1, parm_%N_tau+1
          GxGA(k)=dmft_%local_greenA%left_mixing(n,k)*conjg(dmft_%local_greenA%left_mixing(n,parm_%N_tau-k+2))
          GxGB(k)=dmft_%local_greenB%left_mixing(n,k)*conjg(dmft_%local_greenB%left_mixing(n,parm_%N_tau-k+2))
       end do
       dmft_%kinetic_energyA(n)=dmft_%kinetic_energyA(n)+2.d0*parm_%dtau*imag((-xj)*trapezoid_z(GxGA,1,parm_%N_tau+1))
       dmft_%kinetic_energyB(n)=dmft_%kinetic_energyB(n)+2.d0*parm_%dtau*imag((-xj)*trapezoid_z(GxGB,1,parm_%N_tau+1))
       do k=1, n
          GxGA(k)=dmft_%local_greenA%retarded(n,k)*dmft_%local_greenA%lesser(k,n)
          GxGB(k)=dmft_%local_greenB%retarded(n,k)*dmft_%local_greenB%lesser(k,n)
       end do
       dmft_%kinetic_energyA(n)=dmft_%kinetic_energyA(n)+2.d0*parm_%dt*imag(trapezoid_z(GxGA,1,n))
       dmft_%kinetic_energyB(n)=dmft_%kinetic_energyB(n)+2.d0*parm_%dt*imag(trapezoid_z(GxGB,1,n))
       do k=1, n
          GxGA(k)=dmft_%local_greenA%lesser(n,k)*conjg(dmft_%local_greenA%retarded(n,k))
          GxGB(k)=dmft_%local_greenB%lesser(n,k)*conjg(dmft_%local_greenB%retarded(n,k))
       end do
       dmft_%kinetic_energyA(n)=dmft_%kinetic_energyA(n)+2.d0*parm_%dt*imag(trapezoid_z(GxGA,1,n))
       dmft_%kinetic_energyB(n)=dmft_%kinetic_energyB(n)+2.d0*parm_%dt*imag(trapezoid_z(GxGB,1,n))
       deallocate(GxGA)
       deallocate(GxGB)
    end if
    
    open (unit=7, file='kinetic-energyA', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energyA(n)
    close (7)
   
    open (unit=7, file='kinetic-energyB', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energyB(n)
    close (7)
    
  end subroutine measure_kinetic_energy
!////////////////////////////////////////
  subroutine measure_kinetic_energy2(parm_,dmft_)
    !  E_{kin}(t)=2\sum_{k} E(k)<c_{k}^+(t)c_{k}(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: k, n
    double precision                :: t
    complex(kind(0d0)), allocatable :: GAxGB(:)
    complex(kind(0d0)), allocatable :: GBxGA(:)
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    dmft_%kinetic_energyA(n)=0.d0
    dmft_%kinetic_energyB(n)=0.d0
    
    if (parm_%dos=='semicircular') then
       ! E_{kin}(t)=2*Im[G^{R}*G^{<}+G^{<}*G^{A}+G^{Left}*G^{Right}](t,t)
       allocate(GBxGA(max(parm_%N_tau+1,n)))
       allocate(GAxGB(max(parm_%N_tau+1,n)))
       do k=1, parm_%N_tau+1
          GAxGB(k)=dmft_%local_greenA%left_mixing(n,k)*conjg(dmft_%local_greenB%left_mixing(n,parm_%N_tau-k+2))
          GBxGA(k)=dmft_%local_greenB%left_mixing(n,k)*conjg(dmft_%local_greenA%left_mixing(n,parm_%N_tau-k+2))
       end do
       dmft_%kinetic_energyA(n)=dmft_%kinetic_energyA(n)+2.d0*parm_%dtau*imag((-xj)*trapezoid_z(GBxGA,1,parm_%N_tau+1))
       dmft_%kinetic_energyB(n)=dmft_%kinetic_energyB(n)+2.d0*parm_%dtau*imag((-xj)*trapezoid_z(GAxGB,1,parm_%N_tau+1))
       do k=1, n
          GBxGA(k)=dmft_%local_greenB%retarded(n,k)*dmft_%local_greenA%lesser(k,n)
          GAxGB(k)=dmft_%local_greenA%retarded(n,k)*dmft_%local_greenB%lesser(k,n)
       end do
       dmft_%kinetic_energyA(n)=dmft_%kinetic_energyA(n)+2.d0*parm_%dt*imag(trapezoid_z(GBxGA,1,n))
       dmft_%kinetic_energyB(n)=dmft_%kinetic_energyB(n)+2.d0*parm_%dt*imag(trapezoid_z(GAxGB,1,n))
       do k=1, n
          GBxGA(k)=dmft_%local_greenB%lesser(n,k)*conjg(dmft_%local_greenA%retarded(n,k))
          GAxGB(k)=dmft_%local_greenA%lesser(n,k)*conjg(dmft_%local_greenB%retarded(n,k))
       end do
       dmft_%kinetic_energyA(n)=dmft_%kinetic_energyA(n)+2.d0*parm_%dt*imag(trapezoid_z(GBxGA,1,n))
       dmft_%kinetic_energyB(n)=dmft_%kinetic_energyB(n)+2.d0*parm_%dt*imag(trapezoid_z(GAxGB,1,n))
       deallocate(GBxGA)
       deallocate(GAxGB)
    end if
    
    open (unit=7, file='kinetic-energyA_', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energyA(n)
    close (7)
   
    open (unit=7, file='kinetic-energyB_',  action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energyB(n)
    close (7)
    
  end subroutine measure_kinetic_energy2


!////////////////////////////////////////
  
 subroutine measure_interaction_energy(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: n
    double precision                :: t
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    open (unit=7, file='interaction-energyA', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, U(parm_,t)*(dmft_%double_occupancyA(n)-dmft_%densityA(n)+0.25d0)
    close (7)
    open (unit=7, file='interaction-energyB', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, U(parm_,t)*(dmft_%double_occupancyB(n)-dmft_%densityB(n)+0.25d0)
    close (7)
    
  end subroutine measure_interaction_energy
  
!///////////////////////////////////
  
  subroutine measure_total_energy(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: n
    double precision                :: t
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    open (unit=7, file='total-energyA', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energyA(n)+U(parm_,t)*(dmft_%double_occupancyA(n)-dmft_%densityA(n)+0.25d0)
    close (7)
    open (unit=7, file='total-energyB', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energyB(n)+U(parm_,t)*(dmft_%double_occupancyB(n)-dmft_%densityB(n)+0.25d0)
    close (7)
    
  end subroutine measure_total_energy
  
!****************************** Printing Data *********************************************   
   subroutine print_green_reatarded(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                	    :: j
    double precision                :: w
    
    print*, "printing Green' function************************"
    open (unit=12, file='GfA_reatarded', action='write')
    do j=1, parm_%n_t +1     
       write (12,*) dmft_%local_greenA%retarded(j,:)
    end do
    
    close (12)

    open (unit=13, file='GfB_reatarded',  action='write')
   do j=1, parm_%n_t+1      
       write (13,*) dmft_%local_greenB%retarded(j,:)
    end do
    close (13)
    
  end subroutine print_green_reatarded

  
end module NONEQ_DMFT_MOD
