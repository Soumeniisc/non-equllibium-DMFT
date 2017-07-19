!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  date:   February 28, 2013
!
!----------------------------------------------------------------------------
module IPT_MOD
  
  use CONST_MOD
  use PARM_MOD
  use DMFT_MOD
  use GREEN_MOD
  use FFT_MOD
  implicit none
  
  private
  public :: eq_ipt,                     &
            initialize_self_energy_ipt, &
            noneq_ipt
  
contains
  
  subroutine eq_ipt(parm_,dmft_)
    type(parm), intent(in)       :: parm_
    type(dmft), intent(inout)    :: dmft_
    integer                      :: i, j
    double precision             :: w, tau
    
    ! second-order perturbation in U at half filling
    ! self_energy(t) = -U(0-)*U(0-)*G0(t)*G0(-t)*G0(t)
    do i=1, parm_%N_tau+1
       dmft_%self_energy%matsubara_t(i)=parm_%U_i**2*dmft_%weiss_green%matsubara_t(i)*dmft_%weiss_green%matsubara_t(parm_%N_tau-i+2)*dmft_%weiss_green%matsubara_t(i)
    end do
    
    ! Fourier transformation
    ! self_energy(t) -> self_energy(w)
    call fft_t2w(parm_,dmft_%self_energy%matsubara_t,dmft_%self_energy%matsubara_w)
    
    ! solve the impurity Dyson equation
    ! G(w)=[G0(w)^{-1}-self_energy(w)]^{-1}
    do j=1, parm_%N_tau
       dmft_%local_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)/(1.d0-dmft_%weiss_green%matsubara_w(j)*dmft_%self_energy%matsubara_w(j))
    end do
    
  end subroutine eq_ipt
  
  
  subroutine initialize_self_energy_ipt(parm_,dmft_)
    type(parm), intent(in)     :: parm_
    type(dmft), intent(inout)  :: dmft_
    integer                    :: j

    ! initialize the self-energy for real-time evolution
    ! self_energy_{12}(0,0) = i^3*U(0+)*U(0+)*G0(0-)*G0(0+)*G0(0-)
    dmft_%self_energy_T%c12(1,1)=U(parm_,0.d0)*U(parm_,0.d0)*(-xj)*dmft_%weiss_green%matsubara_t(parm_%N_tau+1) &
                                 *xj*dmft_%weiss_green%matsubara_t(1)*(-xj)*dmft_%weiss_green%matsubara_t(parm_%N_tau+1)
    ! self_energy_{21}(0,0) = i^3*U(0+)*U(0+)*G0(0+)*G0(0-)*G0(0+)
    dmft_%self_energy_T%c21(1,1)=U(parm_,0.d0)*U(parm_,0.d0)*xj*dmft_%weiss_green%matsubara_t(1) &
                                 *(-xj)*dmft_%weiss_green%matsubara_t(parm_%N_tau+1)*xj*dmft_%weiss_green%matsubara_t(1)
    do j=1, parm_%N_tau+1
       ! self_energy_{13}(0,t) = i^3*U(0+)*U(0-)*G0(-t)*G0(t)*G0(-t)
       dmft_%self_energy_T%c13(1,j)=U(parm_,0.d0)*parm_%U_i*(-xj)*dmft_%weiss_green%matsubara_t(parm_%N_tau-j+2) &
                                    *xj*dmft_%weiss_green%matsubara_t(j)*(-xj)*dmft_%weiss_green%matsubara_t(parm_%N_tau-j+2)
       ! self_energy_{31}(t,0) = i^3*U(0-)*U(0+)*G0(t)*G0(-t)*G0(t)
       dmft_%self_energy_T%c31(j,1)=parm_%U_i*U(parm_,0.d0)*xj*dmft_%weiss_green%matsubara_t(j) &
                                    *(-xj)*dmft_%weiss_green%matsubara_t(parm_%N_tau-j+2)*xj*dmft_%weiss_green%matsubara_t(j)
    end do
    dmft_%self_energy%retarded(1,1)=dmft_%self_energy_T%c21(1,1)-dmft_%self_energy_T%c12(1,1)
    dmft_%self_energy%lesser(1,1)=dmft_%self_energy_T%c12(1,1)
    dmft_%self_energy%left_mixing(1,:)=dmft_%self_energy_T%c13(1,:)
    
  end subroutine initialize_self_energy_ipt
  
  
  subroutine noneq_ipt(parm_,dmft_)
    type(parm), intent(in)     :: parm_
    type(dmft), intent(inout)  :: dmft_
    type(KB)                   :: kernel
    double precision           :: t1, t2
    integer                    :: i, j, n
    
    ! second-order perturbation in U at half filling
    ! self_energy(t,t') = U(t)*U(t')*G0(t,t')*G0(t',t)*G0(t,t')
    n=dmft_%time_step
    do j=1, n
       t2=dble(j-1)*parm_%dt
       t1=dble(n-1)*parm_%dt
       dmft_%self_energy_T%c12(n,j)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_T%c12(n,j)*dmft_%weiss_green_T%c21(j,n)*dmft_%weiss_green_T%c12(n,j)
       dmft_%self_energy_T%c21(n,j)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_T%c21(n,j)*dmft_%weiss_green_T%c12(j,n)*dmft_%weiss_green_T%c21(n,j)
    end do
    t2=dble(n-1)*parm_%dt
    do i=1, n-1
       t1=dble(i-1)*parm_%dt
       dmft_%self_energy_T%c12(i,n)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_T%c12(i,n)*dmft_%weiss_green_T%c21(n,i)*dmft_%weiss_green_T%c12(i,n)
       dmft_%self_energy_T%c21(i,n)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_T%c21(i,n)*dmft_%weiss_green_T%c12(n,i)*dmft_%weiss_green_T%c21(i,n)
    end do
    t1=dble(n-1)*parm_%dt
    dmft_%self_energy_T%c13(n,:)=U(parm_,t1)*parm_%U_i*dmft_%weiss_green_T%c13(n,:)*dmft_%weiss_green_T%c31(:,n)*dmft_%weiss_green_T%c13(n,:)
    t2=dble(n-1)*parm_%dt
    dmft_%self_energy_T%c31(:,n)=parm_%U_i*U(parm_,t2)*dmft_%weiss_green_T%c31(:,n)*dmft_%weiss_green_T%c13(n,:)*dmft_%weiss_green_T%c31(:,n)
    
    ! Kadanoff-Baym self-energy <= contour-ordered self-energy
    do j=1, n
       dmft_%self_energy%retarded(n,j)=dmft_%self_energy_T%c21(n,j)-dmft_%self_energy_T%c12(n,j)
       dmft_%self_energy%lesser(n,j)=dmft_%self_energy_T%c12(n,j)
    end do
    do i=1, n-1
       dmft_%self_energy%lesser(i,n)=dmft_%self_energy_T%c12(i,n)
    end do
    dmft_%self_energy%left_mixing(n,:)=dmft_%self_energy_T%c13(n,:)
    
    ! solve the Impurity Dyson equation: G = G0+G0*self_energy*G
    ! K=G0*self_energy
    call allocate_KB(parm_,kernel)
    call convolute_KB(parm_,n,dmft_%weiss_green,dmft_%self_energy,kernel)
    
    ! G=[1-K]^{-1}*G0
    call Volterra_int(parm_,n,dmft_%weiss_green,kernel,dmft_%local_green)
    call deallocate_KB(kernel)
    
  end subroutine noneq_ipt
  
  
end module IPT_MOD
