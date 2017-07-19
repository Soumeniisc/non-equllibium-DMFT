!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  data:   February 28, 2013
!
!----------------------------------------------------------------------------
module EQ_DMFT_MOD
  
  use CONST_MOD
  use PARM_MOD
  use DMFT_MOD
  use GREEN_MOD
  use IPT_MOD
  use FFT_MOD
  use INTEGRAL_MOD
  implicit none
  
  private
  public :: start_eq_dmft
  
contains
  
  subroutine start_eq_dmft(parm_,dmft_)
    type(parm), intent(in)        :: parm_
    type(dmft), intent(inout)     :: dmft_
    integer                       :: i, j
    double precision              :: w, tau
    
    dmft_%converged=.false.
    dmft_%G0_diff=1.d0
    dmft_%time_step=1
    dmft_%iteration=1
    call initialize_eq_green(parm_,dmft_)
    ! equilibrium DMFT self-consistency loop
    do while (.not. dmft_%converged .and. dmft_%iteration<=parm_%N_iter)
       call eq_impurity_solution(parm_,dmft_)
       call eq_dmft_self_consistency(parm_,dmft_)
       dmft_%iteration=dmft_%iteration+1
       if (dmft_%G0_diff<=parm_%tolerance) then
          dmft_%converged=.true.
          write (*,'(A)') 'Equilibrium DMFT is converged'
       end if
       if (dmft_%iteration==parm_%N_iter+1 .and. dmft_%G0_diff>parm_%tolerance) then
          write (*,'(A)') 'Equilibrium DMFT is NOT converged !!!!!'
       end if
    end do
    
    ! Fourier transformation: G(w) -> G(t)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_green%matsubara_w(j)=dmft_%local_green%matsubara_w(j)-1.d0/(xj*w)-(parm_%e2+0.25d0*parm_%U_i**2)/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%local_green%matsubara_w,dmft_%local_green%matsubara_t)
    do i=1, parm_%N_tau+1
       tau=dble(i-1)*parm_%dtau
       dmft_%local_green%matsubara_t(i)=dmft_%local_green%matsubara_t(i)-0.5d0+0.25d0*(parm_%e2+0.25d0*parm_%U_i**2)*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_green%matsubara_w(j)=dmft_%local_green%matsubara_w(j)+1.d0/(xj*w)+(parm_%e2+0.25d0*parm_%U_i**2)/(xj*w)/(xj*w)/(xj*w)
    end do

    open (unit=7, file='G_local_tau',  action='write')
    do j=1, parm_%N_tau+1
    write (7,*)  (j-1)*parm_%dtau, dmft_%local_green%matsubara_t(j)
    end do

    call initialize_self_energy(parm_,dmft_)
    
    call measure_density(parm_,dmft_)
    call measure_double_occupancy(parm_,dmft_)
    call measure_kinetic_energy(parm_,dmft_)
    call measure_interaction_energy(parm_,dmft_)
    call measure_total_energy(parm_,dmft_)
    
  end subroutine start_eq_dmft
  
  
  subroutine eq_impurity_solution(parm_,dmft_)
    type(parm), intent(in)       :: parm_
    type(dmft), intent(inout)    :: dmft_
    
    if (parm_%solver=='IPT') then
       call eq_ipt(parm_,dmft_)
    end if
    
  end subroutine eq_impurity_solution
  
  
  subroutine initialize_self_energy(parm_,dmft_)
    type(parm), intent(in)       :: parm_
    type(dmft), intent(inout)    :: dmft_
    
    if (parm_%solver=='IPT') then
       call initialize_self_energy_ipt(parm_,dmft_)
    end if
    
  end subroutine initialize_self_energy
  
  
  subroutine eq_dmft_self_consistency(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: i, j, k
    double precision                :: w, tau
    
    if (parm_%dos=='semicircular') then
       ! solve the lattice Dyson equation: G0(w)=1/[i*w-G(w)]
       do j=1, parm_%N_tau
          w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
          dmft_%weiss_green%matsubara_w(j)=1.d0/(xj*w-dmft_%local_green%matsubara_w(j))
       end do
    end if
    
    ! Fourier transformation: G0(w) -> G0(t)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)-1.d0/(xj*w)-parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%weiss_green%matsubara_w,dmft_%weiss_green_new%matsubara_t)
    do i=1, parm_%N_tau+1
       tau=dble(i-1)*parm_%dtau
       dmft_%weiss_green_new%matsubara_t(i)=dmft_%weiss_green_new%matsubara_t(i)-0.5d0+0.25d0*parm_%e2*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+1.d0/(xj*w)+parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    
    ! evaluate |G0_{new}-G0_{old}|
    dmft_%G0_diff=0.d0
    do i=1, parm_%N_tau+1
       dmft_%G0_diff=dmft_%G0_diff+abs(dmft_%weiss_green_new%matsubara_t(i)-dmft_%weiss_green%matsubara_t(i))
    end do
    write (*,'(A,I3)') '  Iteration #', dmft_%iteration
    write (*,'(A,f15.10)') ' |G0_new-G0_old| =', dmft_%G0_diff
    ! G0_{old} <= G0_{new}
    dmft_%weiss_green%matsubara_t(:)=dmft_%weiss_green_new%matsubara_t(:)
    
  end subroutine eq_dmft_self_consistency
  
  
  subroutine initialize_eq_green(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    double precision                :: ek, e_max, de, w, tau
    integer                         :: i, j, k
    
    ! initialize G0(w)
    if (parm_%dos=='semicircular') then
       e_max=2.d0
       de=e_max/dble(parm_%N_e)
       dmft_%weiss_green%matsubara_w(j)=0.d0
       do j=1, parm_%N_tau
          w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
          do k=2, 2*parm_%N_e
             ek=dble(k-parm_%N_e-1)*de
             dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+de*sqrt(4.d0-ek**2)/(2.d0*pi)/(xj*w-ek)
          end do
       end do
    end if
    
    ! Fourier transformation: G0(w) -> G0(t)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)-1.d0/(xj*w)-parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%weiss_green%matsubara_w,dmft_%weiss_green%matsubara_t)
    do i=1, parm_%N_tau+1
       tau=dble(i-1)*parm_%dtau
       dmft_%weiss_green%matsubara_t(i)=dmft_%weiss_green%matsubara_t(i)-0.5d0+0.25d0*parm_%e2*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+1.d0/(xj*w)+parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    
  end subroutine initialize_eq_green
  
  
  subroutine measure_density(parm_,dmft_)
    !  n(t)=<c^+(t)c(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    double precision                :: t
    
    t=0.d0
    dmft_%density(1)=-dmft_%local_green%matsubara_t(parm_%N_tau+1)
    
    open (unit=7, file='density', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%density(1)
    close (7)
    
  end subroutine measure_density
  
  
  subroutine measure_double_occupancy(parm_,dmft_)
    !  d(t)=<n_up(t)*n_do(t)>
    type(parm), intent(in)        :: parm_
    type(dmft), intent(inout)     :: dmft_
    integer                       :: i
    double precision              :: t
    double precision, allocatable :: SxG(:)
    
    t=0.d0
    ! d(t)=n_up(t)*n_do(t)-1/U*\int dtau self_energy^{M}(beta-tau)*G^{M}(tau)
    if (parm_%U_i==0.d0) then
       dmft_%double_occupancy(1)=dmft_%density(1)**2
    else
       allocate(SxG(parm_%N_tau+1))
       do i=1, parm_%N_tau+1
          SxG(i)=dmft_%self_energy%matsubara_t(parm_%N_tau-i+2)*dmft_%local_green%matsubara_t(i)
       end do
       dmft_%double_occupancy(1)=dmft_%density(1)**2-1.d0/parm_%U_i*parm_%dtau*trapezoid_d(SxG,1,parm_%N_tau+1)
       deallocate(SxG)
    end if
    
    open (unit=7, file='double-occupancy', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%double_occupancy(1)
    close (7)
    
  end subroutine measure_double_occupancy
  
  
  subroutine measure_kinetic_energy(parm_,dmft_)
    !  E_{kin}(t)=2\sum_{k} E(k)<c_{k}^+(t)c_{k}(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: i
    double precision                :: t
    double precision, allocatable   :: GxG(:)
    
    t=0.d0
    if (parm_%dos=='semicircular') then
       ! E_{kin}(0)=-2*\int_0^{beta} dtau G^M(tau)*G^M(beta-tau)
       allocate(GxG(parm_%N_tau+1))
       do i=1, parm_%N_tau+1
          GxG(i)=dmft_%local_green%matsubara_t(i)*dmft_%local_green%matsubara_t(parm_%N_tau-i+2)
       end do
       dmft_%kinetic_energy(1)=-2.d0*parm_%dtau*trapezoid_d(GxG,1,parm_%N_tau+1)
       deallocate(GxG)
    end if
    
    open (unit=7, file='kinetic-energy', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energy(1)
    close (7)
    
  end subroutine measure_kinetic_energy
  
  
  subroutine measure_interaction_energy(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    double precision                :: t
    
    t=0.d0
    open (unit=7, file='interaction-energy', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, parm_%U_i*(dmft_%double_occupancy(1)-dmft_%density(1)+0.25d0)
    close (7)
    
  end subroutine measure_interaction_energy
  
  
  subroutine measure_total_energy(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    double precision                :: t
    
    t=0.d0
    open (unit=7, file='total-energy', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10)') t, dmft_%kinetic_energy(1)+parm_%U_i*(dmft_%double_occupancy(1)-dmft_%density(1)+0.25d0)
    close (7)
    
  end subroutine measure_total_energy
  
  
end module EQ_DMFT_MOD
