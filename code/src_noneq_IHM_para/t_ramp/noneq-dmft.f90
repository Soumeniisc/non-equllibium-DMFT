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
       dmft_%G0_diff=1.d0
       dmft_%iteration=1
       dmft_%time_step=n
       call initialize_noneq_green(parm_,dmft_)
       write (*,'(A,f15.10)') 't =', dble(n-1)*parm_%dt
       
       ! nonequilibrium DMFT self-consistency loop
       do while(.not. dmft_%converged .and. dmft_%iteration<=parm_%N_iter)
          call noneq_impurity_solution(parm_,dmft_)
          call noneq_dmft_self_consistency(parm_,dmft_)
          dmft_%iteration=dmft_%iteration+1
          if (dmft_%G0_diff<=parm_%tolerance) then
             dmft_%converged=.true.
             write (*,'(A)') 'Nonequilibrium DMFT is converged'
          end if
          if (dmft_%G0_diff>parm_%tolerance .and. dmft_%iteration>parm_%N_iter) then
             write (*,'(A)') 'Nonequilibrium DMFT is NOT converged !!!!!'
          end if
       end do
       
       ! d/dt G0_{old} <= d/dt G0_{new}
       call deallocate_KB_derivative(dmft_%weiss_green_der)
       call allocate_KB_derivative(parm_,dmft_%weiss_green_der,n)
       do i=1, n
          dmft_%weiss_green_der%retarded(i)=dmft_%weiss_green_der_new%retarded(i)
          dmft_%weiss_green_der%lesser(i)=dmft_%weiss_green_der_new%lesser(i)
       end do
       do j=1, parm_%N_tau+1
          dmft_%weiss_green_der%left_mixing(j)=dmft_%weiss_green_der_new%left_mixing(j)
       end do
       call deallocate_KB_derivative(dmft_%weiss_green_der_new)
       call allocate_KB_derivative(parm_,dmft_%weiss_green_der_new,n+1)
       
       call measure_density(parm_,dmft_)
       call measure_double_occupancy(parm_,dmft_)
       call measure_kinetic_energy(parm_,dmft_)
       call measure_interaction_energy(parm_,dmft_)
       call measure_total_energy(parm_,dmft_)
    end do
    
  end subroutine start_noneq_dmft
  
  
  subroutine noneq_impurity_solution(parm_,dmft_)
    type(parm), intent(in)       :: parm_
    type(dmft), intent(inout)    :: dmft_
    
    if (parm_%solver=='IPT') then
       call noneq_ipt(parm_,dmft_)
    end if
    
  end subroutine noneq_impurity_solution
  
  
  subroutine noneq_dmft_self_consistency(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    type(KB)                        :: kernel
    integer                         :: i, j, k, n
    complex(kind(0d0)), allocatable :: h(:)
    
    n=dmft_%time_step
    
    if (parm_%dos=='semicircular') then
       ! solve the lattice Dyson equation: id/dt G0(t,t')-G*G0(t,t')=delta(t,t')
       allocate(h(n))
       h(:)=0.d0
       call Volterra_intdiff(parm_,n,h,dmft_%local_green,dmft_%weiss_green_new,dmft_%weiss_green_der,dmft_%weiss_green_der_new)
       deallocate(h)
    end if
    
    ! evaluate |G0_{new}-G0_{old}|
    dmft_%G0_diff=0.d0
    do j=1, n
       dmft_%G0_diff=dmft_%G0_diff+abs(dmft_%weiss_green_new%retarded(n,j)-dmft_%weiss_green%retarded(n,j))
       dmft_%G0_diff=dmft_%G0_diff+abs(dmft_%weiss_green_new%lesser(n,j)-dmft_%weiss_green%lesser(n,j))
    end do
    do j=1, parm_%N_tau+1
       dmft_%G0_diff=dmft_%G0_diff+abs(dmft_%weiss_green_new%left_mixing(n,j)-dmft_%weiss_green%left_mixing(n,j))
    end do
    do i=1, n-1
       dmft_%G0_diff=dmft_%G0_diff+abs(dmft_%weiss_green_new%lesser(i,n)-dmft_%weiss_green%lesser(i,n))
    end do
    write (*,'(A,I2)') '  Iteration #', dmft_%iteration
    write (*,'(A,f15.10)') ' |G0_new-G0_old| =', dmft_%G0_diff
    
    ! G0_{old} <= G0_{new}
    do j=1, n
       dmft_%weiss_green%retarded(n,j)=dmft_%weiss_green_new%retarded(n,j)
       dmft_%weiss_green%lesser(n,j)=dmft_%weiss_green_new%lesser(n,j)
    end do
    do i=1, n-1
       dmft_%weiss_green%lesser(i,n)=dmft_%weiss_green_new%lesser(i,n)
    end do
    dmft_%weiss_green%left_mixing(n,:)=dmft_%weiss_green_new%left_mixing(n,:)
    
    ! Kadanoff-Baym G0 => contour-ordered G0
    do j=1, n
       dmft_%weiss_green_T%c12(n,j)=dmft_%weiss_green%lesser(n,j)
       dmft_%weiss_green_T%c21(n,j)=dmft_%weiss_green%lesser(n,j)+dmft_%weiss_green%retarded(n,j)
    end do
    do i=1, n-1
       dmft_%weiss_green_T%c12(i,n)=dmft_%weiss_green%lesser(i,n)
       dmft_%weiss_green_T%c21(i,n)=dmft_%weiss_green%lesser(i,n)-conjg(dmft_%weiss_green%retarded(n,i))
    end do
    do j=1, parm_%N_tau+1
       dmft_%weiss_green_T%c13(n,j)=dmft_%weiss_green%left_mixing(n,j)
    end do
    do i=1, parm_%N_tau+1
       dmft_%weiss_green_T%c31(i,n)=conjg(dmft_%weiss_green%left_mixing(n,parm_%N_tau-i+2))
    end do
    
    ! Hermite conjugate
    do j=1, n
       dmft_%weiss_green_T%c12(n,j)=0.5d0*(dmft_%weiss_green_T%c12(n,j)-conjg(dmft_%weiss_green_T%c12(j,n)))
       dmft_%weiss_green_T%c21(n,j)=0.5d0*(dmft_%weiss_green_T%c21(n,j)-conjg(dmft_%weiss_green_T%c21(j,n)))
    end do
    do i=1, n-1
       dmft_%weiss_green_T%c12(i,n)=-conjg(dmft_%weiss_green_T%c12(n,i))
       dmft_%weiss_green_T%c21(i,n)=-conjg(dmft_%weiss_green_T%c21(n,i))
    end do
    do j=1, parm_%N_tau+1
       dmft_%weiss_green_T%c13(n,j)=0.5d0*(dmft_%weiss_green_T%c13(n,j)+conjg(dmft_%weiss_green_T%c31(parm_%N_tau-j+2,n)))
    end do
    do i=1, parm_%N_tau+1
       dmft_%weiss_green_T%c31(i,n)=conjg(dmft_%weiss_green_T%c13(n,parm_%N_tau-i+2))
    end do
    
  end subroutine noneq_dmft_self_consistency
  
  
  subroutine initialize_noneq_green(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: i, j, k, n
    complex(kind(0d0)), allocatable :: SxG(:)
    
    n=dmft_%time_step
    ! initialize G0(t,t') and G(t,t')
    if (n==2) then
       dmft_%weiss_green%retarded(1,1)=-xj
       dmft_%local_green%retarded(1,1)=-xj
       dmft_%weiss_green_new%retarded(1,1)=-xj
       do j=1, parm_%N_tau+1
          dmft_%weiss_green%left_mixing(1,j)=-xj*dmft_%weiss_green%matsubara_t(parm_%N_tau-j+2)
          dmft_%local_green%left_mixing(1,j)=-xj*dmft_%local_green%matsubara_t(parm_%N_tau-j+2)
          dmft_%weiss_green_new%left_mixing(1,j)=-xj*dmft_%weiss_green%matsubara_t(parm_%N_tau-j+2)
       end do
       dmft_%weiss_green%lesser(1,1)=-xj*dmft_%weiss_green%matsubara_t(parm_%N_tau+1)
       dmft_%local_green%lesser(1,1)=-xj*dmft_%local_green%matsubara_t(parm_%N_tau+1)
       dmft_%weiss_green_new%lesser(1,1)=-xj*dmft_%weiss_green%matsubara_t(parm_%N_tau+1)
       
       if (parm_%dos=='semicircular') then
          ! initialize d/dt G0(t,t') = -i*(G*G0)(t,t')
          dmft_%weiss_green_der%retarded(1)=0.d0
          dmft_%weiss_green_der%left_mixing(:)=0.d0
          allocate(SxG(parm_%N_tau+1))
          do j=1, parm_%N_tau+1
             do k=1, j
                SxG(k)=dmft_%local_green%left_mixing(1,k)*dmft_%weiss_green%matsubara_t(parm_%N_tau+k-j+1)
             end do
             dmft_%weiss_green_der%left_mixing(j)=dmft_%weiss_green_der%left_mixing(j)+xj*parm_%dtau*trapezoid_z(SxG,1,j)
             do k=j, parm_%N_tau+1
                SxG(k)=dmft_%local_green%left_mixing(1,k)*dmft_%weiss_green%matsubara_t(k-j+1)
             end do
             dmft_%weiss_green_der%left_mixing(j)=dmft_%weiss_green_der%left_mixing(j)-xj*parm_%dtau*trapezoid_z(SxG,j,parm_%N_tau+1)
          end do
          do k=1, parm_%N_tau+1
             SxG(k)=dmft_%local_green%left_mixing(1,k)*conjg(dmft_%weiss_green%left_mixing(1,parm_%N_tau-k+2))
          end do
          dmft_%weiss_green_der%lesser(1)=-xj*(-xj)*parm_%dtau*trapezoid_z(SxG,1,parm_%N_tau+1)
          deallocate(SxG)
       end if
       
       ! guess G0(t,t') in the next step
       do j=1, n
          dmft_%weiss_green%retarded(n,j)=dmft_%weiss_green%retarded(1,1)
       end do
       do j=1, parm_%N_tau+1
          dmft_%weiss_green%left_mixing(n,j)=dmft_%weiss_green%left_mixing(1,j)
       end do
       do j=1, n
          dmft_%weiss_green%lesser(n,j)=dmft_%weiss_green%lesser(1,1)
       end do
       do i=1, n-1
          dmft_%weiss_green%lesser(i,n)=dmft_%weiss_green%lesser(1,1)
       end do
       
       ! contour-ordered G0 <= Kadanoff-Baym G0
       do j=1, n
          do i=1, n
             dmft_%weiss_green_T%c12(i,j)=dmft_%weiss_green%lesser(i,j)
             dmft_%weiss_green_T%c21(i,j)=dmft_%weiss_green%lesser(i,j)-xj
          end do
       end do
       do j=1, parm_%N_tau+1
          do i=1, n
             dmft_%weiss_green_T%c13(i,j)=-xj*dmft_%weiss_green%matsubara_t(parm_%N_tau-j+2)
             dmft_%weiss_green_T%c31(j,i)=xj*dmft_%weiss_green%matsubara_t(j)
          end do
       end do
    else if (n>=3) then
       ! guess G0(t,t') in the next step by quadratic extrapolation
       dmft_%weiss_green%retarded(n,n)=-xj
       do k=1, n-2
          dmft_%weiss_green%retarded(n,k)=2.d0*dmft_%weiss_green%retarded(n-1,k)-dmft_%weiss_green%retarded(n-2,k)
       end do
       dmft_%weiss_green%retarded(n,n-1)=0.5d0*(dmft_%weiss_green%retarded(n,n)+dmft_%weiss_green%retarded(n,n-2))
       do k=1, parm_%N_tau+1
          dmft_%weiss_green%left_mixing(n,k)=2.d0*dmft_%weiss_green%left_mixing(n-1,k)-dmft_%weiss_green%left_mixing(n-2,k)
       end do
       do k=1, n-1
          dmft_%weiss_green%lesser(n,k)=2.d0*dmft_%weiss_green%lesser(n-1,k)-dmft_%weiss_green%lesser(n-2,k)
          dmft_%weiss_green%lesser(k,n)=2.d0*dmft_%weiss_green%lesser(k,n-1)-dmft_%weiss_green%lesser(k,n-2)
       end do
       dmft_%weiss_green%lesser(n,n)=2.d0*dmft_%weiss_green%lesser(n-1,n-1)-dmft_%weiss_green%lesser(n-2,n-2)
       
       ! contour-ordered G0 <= Kadanoff-Baym G0
       do k=1, n-1
          dmft_%weiss_green_T%c12(n,k)=dmft_%weiss_green%lesser(n,k)
          dmft_%weiss_green_T%c12(k,n)=dmft_%weiss_green%lesser(k,n)
          dmft_%weiss_green_T%c21(n,k)=dmft_%weiss_green%lesser(n,k)+dmft_%weiss_green%retarded(n,k)
          dmft_%weiss_green_T%c21(k,n)=dmft_%weiss_green%lesser(k,n)-conjg(dmft_%weiss_green%retarded(n,k))
       end do
       dmft_%weiss_green_T%c12(n,n)=dmft_%weiss_green%lesser(n,n)
       dmft_%weiss_green_T%c21(n,n)=dmft_%weiss_green%lesser(n,n)-xj
       do k=1, parm_%N_tau+1
          dmft_%weiss_green_T%c13(n,k)=dmft_%weiss_green%left_mixing(n,k)
          dmft_%weiss_green_T%c31(k,n)=conjg(dmft_%weiss_green%left_mixing(n,parm_%N_tau-k+2))
       end do
    end if
    
  end subroutine initialize_noneq_green
  
  
  subroutine measure_density(parm_,dmft_)
    !  n(t)=<c^+(t)c(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: n
    double precision                :: t
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    dmft_%density(n)=imag(dmft_%local_green%lesser(n,n))
    
    open (unit=7, file='density', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%density(n)
    close (7)
    
  end subroutine measure_density
  
  
  subroutine measure_double_occupancy(parm_,dmft_)
    !  d(t)=<n_up(t)*n_do(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: j, k, n
    double precision                :: t
    complex(kind(0d0)), allocatable :: SxG(:)
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    ! d(t)=n_up(t)*n_do(t)-i/U*[self_energy*G]^{<}(t,t)
    !     =n_up(t)*n_do(t)-i/U*[self_energy^{R}*G^{<}+self_energy^{<}*G^{A}+self_energy^{Left}*G^{Right}](t,t)
    if (U(parm_,t)==0.d0) then
       dmft_%double_occupancy(n)=dmft_%density(n)**2
    else
       dmft_%double_occupancy(n)=dmft_%density(n)**2
       allocate(SxG(max(parm_%N_tau+1,n)))
       do k=1, parm_%N_tau+1
          SxG(k)=dmft_%self_energy%left_mixing(n,k)*conjg(dmft_%local_green%left_mixing(n,parm_%N_tau-k+2))
       end do
       dmft_%double_occupancy(n)=dmft_%double_occupancy(n)+1.d0/U(parm_,t)*parm_%dtau*imag((-xj)*trapezoid_z(SxG,1,parm_%N_tau+1))
       do k=1, n
          SxG(k)=dmft_%self_energy%retarded(n,k)*dmft_%local_green%lesser(k,n)
       end do
       dmft_%double_occupancy(n)=dmft_%double_occupancy(n)+1.d0/U(parm_,t)*parm_%dt*imag(trapezoid_z(SxG,1,n))
       do k=1, n
          SxG(k)=dmft_%self_energy%lesser(n,k)*conjg(dmft_%local_green%retarded(n,k))
       end do
       dmft_%double_occupancy(n)=dmft_%double_occupancy(n)+1.d0/U(parm_,t)*parm_%dt*imag(trapezoid_z(SxG,1,n))
       deallocate(SxG)
    end if
    
    open (unit=7, file='double-occupancy', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%double_occupancy(n)
    close (7)
    
  end subroutine measure_double_occupancy
  
  
  subroutine measure_kinetic_energy(parm_,dmft_)
    !  E_{kin}(t)=2\sum_{k} E(k)<c_{k}^+(t)c_{k}(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: k, n
    double precision                :: t
    complex(kind(0d0)), allocatable :: GxG(:)
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    dmft_%kinetic_energy(n)=0.d0
    
    if (parm_%dos=='semicircular') then
       ! E_{kin}(t)=2*Im[G^{R}*G^{<}+G^{<}*G^{A}+G^{Left}*G^{Right}](t,t)
       allocate(GxG(max(parm_%N_tau+1,n)))
       do k=1, parm_%N_tau+1
          GxG(k)=dmft_%local_green%left_mixing(n,k)*conjg(dmft_%local_green%left_mixing(n,parm_%N_tau-k+2))
       end do
       dmft_%kinetic_energy(n)=dmft_%kinetic_energy(n)+2.d0*parm_%dtau*imag((-xj)*trapezoid_z(GxG,1,parm_%N_tau+1))
       do k=1, n
          GxG(k)=dmft_%local_green%retarded(n,k)*dmft_%local_green%lesser(k,n)
       end do
       dmft_%kinetic_energy(n)=dmft_%kinetic_energy(n)+2.d0*parm_%dt*imag(trapezoid_z(GxG,1,n))
       do k=1, n
          GxG(k)=dmft_%local_green%lesser(n,k)*conjg(dmft_%local_green%retarded(n,k))
       end do
       dmft_%kinetic_energy(n)=dmft_%kinetic_energy(n)+2.d0*parm_%dt*imag(trapezoid_z(GxG,1,n))
       deallocate(GxG)
    end if
    
    open (unit=7, file='kinetic-energy', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energy(n)
    close (7)
    
  end subroutine measure_kinetic_energy
  
  
  subroutine measure_interaction_energy(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: n
    double precision                :: t
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    open (unit=7, file='interaction-energy', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, U(parm_,t)*(dmft_%double_occupancy(n)-dmft_%density(n)+0.25d0)
    close (7)
    
  end subroutine measure_interaction_energy
  
  
  subroutine measure_total_energy(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: n
    double precision                :: t
    
    n=dmft_%time_step
    t=dble(n-1)*parm_%dt
    open (unit=7, file='total-energy', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energy(n)+U(parm_,t)*(dmft_%double_occupancy(n)-dmft_%density(n)+0.25d0)
    close (7)
    
  end subroutine measure_total_energy
  
  
end module NONEQ_DMFT_MOD
