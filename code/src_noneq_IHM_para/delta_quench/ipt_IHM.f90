!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  date:   February 28, 2013
!
!----------------------------------------------------------------------------
!last modified by soume on 12/04/17 modied over equllibium IHM  ipt_IHM.f90 file in src directory. i have initialize_self_energy_iptA(B) and noneq_iptA(B) subroutine added.
module IPT_MOD
  
  use CONST_MOD
  use PARM_MOD
  use DMFT_MOD
  use GREEN_MOD
  use FFT_MOD
  implicit none
  
  private
  public :: eq_iptA,                     &
            eq_iptB,			 &
	    initialize_self_energy_iptA, &
	    initialize_self_energy_iptB, &
            noneq_iptA,			 &	
	    noneq_iptB
  
contains
  
  subroutine eq_iptA(parm_,dmft_)
    type(parm), intent(in)       :: parm_
    type(dmft), intent(inout)    :: dmft_
    integer                      :: i, j
    double precision             :: w, tau
    
    ! second-order perturbation in U at half filling
    ! self_energy(t) = -U(0-)*U(0-)*G0(t)*G0(-t)*G0(t)
    do i=1, parm_%N_tau+1
       dmft_%self_energyA%matsubara_t(i)=parm_%U_i**2*dmft_%weiss_greenA%matsubara_t(i)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-i+2)*dmft_%weiss_greenA%matsubara_t(i)  
    end do
    
    ! Fourier transformation
    ! self_energy(t) -> self_energy(w)
    call fft_t2w(parm_,dmft_%self_energyA%matsubara_t,dmft_%self_energyA%matsubara_w)

   
    ! solve the impurity Dyson equation
    ! G(w)=[G0(w)^{-1}-self_energy(w)]^{-1}
    do j=1, parm_%N_tau
       dmft_%local_greenA%matsubara_w(j)=dmft_%weiss_greenA%matsubara_w(j)/(1.d0 - dmft_%weiss_greenA%matsubara_w(j)*dmft_%self_energyA%matsubara_w(j))
     
    end do
     do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_greenA%matsubara_w(j)=dmft_%local_greenA%matsubara_w(j)-1.d0/(xj*w)-(parm_%e2)/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%local_greenA%matsubara_w,dmft_%local_greenA%matsubara_t)
    do i=1, parm_%N_tau+1
       tau=dble(i-1)*parm_%dtau
       dmft_%local_greenA%matsubara_t(i)=dmft_%local_greenA%matsubara_t(i)-0.5d0+0.25d0*(parm_%e2)*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_greenA%matsubara_w(j)=dmft_%local_greenA%matsubara_w(j)+1.d0/(xj*w)+(parm_%e2)/(xj*w)/(xj*w)/(xj*w)
    end do
    print*, "before solver A nA", dmft_%nA
    dmft_%nA = -dmft_%local_greenA%matsubara_t(parm_%N_tau+1)
    print*, "after solver A nA", dmft_%nA
    ! this thing will be neded to calculate intrigal diiferential equationin non eq problem and its needed in initialisation too
    dmft_%hA = parm_%delta_i + parm_%U_i*dmft_%nA - parm_%U_i/2.0
  end subroutine eq_iptA

  subroutine eq_iptB(parm_,dmft_)
    type(parm), intent(in)       :: parm_
    type(dmft), intent(inout)    :: dmft_
    integer                      :: i, j
    double precision             :: w, tau
    
    ! second-order perturbation in U at half filling
    ! self_energy(t) = -U(0-)*U(0-)*G0(t)*G0(-t)*G0(t)
    do i=1, parm_%N_tau+1
       dmft_%self_energyB%matsubara_t(i)=parm_%U_i**2*dmft_%weiss_greenB%matsubara_t(i)*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-i+2)*dmft_%weiss_greenB%matsubara_t(i)  ! U*n_B
    end do
    
    ! Fourier transformation
    ! self_energy(t) -> self_energy(w)
    call fft_t2w(parm_,dmft_%self_energyB%matsubara_t,dmft_%self_energyB%matsubara_w)

   

    ! solve the impurity Dyson equation
    ! G(w)=[G0(w)^{-1}-self_energy(w)]^{-1}
    do j=1, parm_%N_tau
       dmft_%local_greenB%matsubara_w(j)=dmft_%weiss_greenB%matsubara_w(j)/(1.d0  - dmft_%weiss_greenB%matsubara_w(j)*dmft_%self_energyB%matsubara_w(j))
    end do
    
     ! doing fourier transform
     do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_greenB%matsubara_w(j)=dmft_%local_greenB%matsubara_w(j)-1.d0/(xj*w)-(parm_%e2)/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%local_greenB%matsubara_w,dmft_%local_greenB%matsubara_t)
    do i=1, parm_%N_tau+1
       tau=dble(i-1)*parm_%dtau
       dmft_%local_greenB%matsubara_t(i)=dmft_%local_greenB%matsubara_t(i)-0.5d0+0.25d0*(parm_%e2)*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_greenB%matsubara_w(j)=dmft_%local_greenB%matsubara_w(j)+1.d0/(xj*w)+(parm_%e2)/(xj*w)/(xj*w)/(xj*w)
    end do
    dmft_%nB = -dmft_%local_greenB%matsubara_t(parm_%N_tau+1)
    dmft_%hB = -parm_%delta_i + parm_%U_i*dmft_%nB - parm_%U_i/2.0
  end subroutine eq_iptB

!*********************************************************************

  subroutine initialize_self_energy_iptA(parm_,dmft_) !TODO carefully chek the intialization for A sublattice
    type(parm), intent(in)     :: parm_
    type(dmft), intent(inout)  :: dmft_
    integer                    :: j

    ! initialize the self-energy for real-time evolution
    ! self_energy_{12}(0,0) = i^3*U(0+)*U(0+)*G0(0-)*G0(0+)*G0(0-)
    dmft_%self_energy_TA%c12(1,1)=U(parm_,0.d0)*U(parm_,0.d0)*(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau+1) &
                                 *xj*dmft_%weiss_greenA%matsubara_t(1)*(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau+1)
    ! self_energy_{21}(0,0) = i^3*U(0+)*U(0+)*G0(0+)*G0(0-)*G0(0+)
    dmft_%self_energy_TA%c21(1,1)=U(parm_,0.d0)*U(parm_,0.d0)*xj*dmft_%weiss_greenA%matsubara_t(1) &
                                 *(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau+1)*xj*dmft_%weiss_greenA%matsubara_t(1)
    do j=1, parm_%N_tau+1
       ! self_energy_{13}(0,t) = i^3*U(0+)*U(0-)*G0(-t)*G0(t)*G0(-t)
       dmft_%self_energy_TA%c13(1,j)=U(parm_,0.d0)*parm_%U_i*(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-j+2) &
                                    *xj*dmft_%weiss_greenA%matsubara_t(j)*(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-j+2)
       ! self_energy_{31}(t,0) = i^3*U(0-)*U(0+)*G0(t)*G0(-t)*G0(t)
       dmft_%self_energy_TA%c31(j,1)=parm_%U_i*U(parm_,0.d0)*xj*dmft_%weiss_greenA%matsubara_t(j) &
                                    *(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-j+2)*xj*dmft_%weiss_greenA%matsubara_t(j)
    end do
    dmft_%self_energyA%retarded(1,1)=dmft_%self_energy_TA%c21(1,1)-dmft_%self_energy_TA%c12(1,1)
    dmft_%self_energyA%lesser(1,1)=dmft_%self_energy_TA%c12(1,1)
    dmft_%self_energyA%left_mixing(1,:)=dmft_%self_energy_TA%c13(1,:)
    
  end subroutine initialize_self_energy_iptA
  
   subroutine initialize_self_energy_iptB(parm_,dmft_)!TODO carefully chek the intialization for B sublattice
    type(parm), intent(in)     :: parm_
    type(dmft), intent(inout)  :: dmft_
    integer                    :: j

    ! initialize the self-energy for real-time evolution
    ! self_energy_{12}(0,0) = i^3*U(0+)*U(0+)*G0(0-)*G0(0+)*G0(0-)
    dmft_%self_energy_TB%c12(1,1)=U(parm_,0.d0)*U(parm_,0.d0)*(-xj)*dmft_%weiss_greenB%matsubara_t(parm_%N_tau+1) &
                                 *xj*dmft_%weiss_greenB%matsubara_t(1)*(-xj)*dmft_%weiss_greenB%matsubara_t(parm_%N_tau+1)
    ! self_energy_{21}(0,0) = i^3*U(0+)*U(0+)*G0(0+)*G0(0-)*G0(0+)
    dmft_%self_energy_TB%c21(1,1)=U(parm_,0.d0)*U(parm_,0.d0)*xj*dmft_%weiss_greenB%matsubara_t(1) &
                                 *(-xj)*dmft_%weiss_greenB%matsubara_t(parm_%N_tau+1)*xj*dmft_%weiss_greenB%matsubara_t(1)
    do j=1, parm_%N_tau+1
       ! self_energy_{13}(0,t) = i^3*U(0+)*U(0-)*G0(-t)*G0(t)*G0(-t)
       dmft_%self_energy_TB%c13(1,j)=U(parm_,0.d0)*parm_%U_i*(-xj)*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-j+2) &
                                    *xj*dmft_%weiss_greenB%matsubara_t(j)*(-xj)*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-j+2)
       ! self_energy_{31}(t,0) = i^3*U(0-)*U(0+)*G0(t)*G0(-t)*G0(t)
       dmft_%self_energy_TB%c31(j,1)=parm_%U_i*U(parm_,0.d0)*xj*dmft_%weiss_greenB%matsubara_t(j) &
                                    *(-xj)*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-j+2)*xj*dmft_%weiss_greenB%matsubara_t(j)
    end do
    dmft_%self_energyB%retarded(1,1)=dmft_%self_energy_TB%c21(1,1)-dmft_%self_energy_TB%c12(1,1)
    dmft_%self_energyB%lesser(1,1)=dmft_%self_energy_TB%c12(1,1)
    dmft_%self_energyB%left_mixing(1,:)=dmft_%self_energy_TB%c13(1,:)
    
  end subroutine initialize_self_energy_iptB
  
!****************************************************************************************

  subroutine noneq_iptA(parm_,dmft_) !TODO carefully chek the self energy for A sublattice and local GF sroe in A sublattice
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
       dmft_%self_energy_TA%c12(n,j)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TA%c12(n,j)*dmft_%weiss_green_TA%c21(j,n)*dmft_%weiss_green_TA%c12(n,j)
       dmft_%self_energy_TA%c21(n,j)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TA%c21(n,j)*dmft_%weiss_green_TA%c12(j,n)*dmft_%weiss_green_TA%c21(n,j)
    end do
    t2=dble(n-1)*parm_%dt
    do i=1, n-1
       t1=dble(i-1)*parm_%dt
       dmft_%self_energy_TA%c12(i,n)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TA%c12(i,n)*dmft_%weiss_green_TA%c21(n,i)*dmft_%weiss_green_TA%c12(i,n)
       dmft_%self_energy_TA%c21(i,n)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TA%c21(i,n)*dmft_%weiss_green_TA%c12(n,i)*dmft_%weiss_green_TA%c21(i,n)
    end do
    t1=dble(n-1)*parm_%dt
    dmft_%self_energy_TA%c13(n,:)=U(parm_,t1)*parm_%U_i*dmft_%weiss_green_TA%c13(n,:)*dmft_%weiss_green_TA%c31(:,n)*dmft_%weiss_green_TA%c13(n,:)
    t2=dble(n-1)*parm_%dt
    dmft_%self_energy_TA%c31(:,n)=parm_%U_i*U(parm_,t2)*dmft_%weiss_green_TA%c31(:,n)*dmft_%weiss_green_TA%c13(n,:)*dmft_%weiss_green_TA%c31(:,n)
    
    ! Kadanoff-Baym self-energy <= contour-ordered self-energy
    do j=1, n
       dmft_%self_energyA%retarded(n,j)=dmft_%self_energy_TA%c21(n,j)-dmft_%self_energy_TA%c12(n,j)
       dmft_%self_energyA%lesser(n,j)=dmft_%self_energy_TA%c12(n,j)
    end do
    do i=1, n-1
       dmft_%self_energyA%lesser(i,n)=dmft_%self_energy_TA%c12(i,n)
    end do
    dmft_%self_energyA%left_mixing(n,:)=dmft_%self_energy_TA%c13(n,:)
    
    ! solve the Impurity Dyson equation: G = G0+G0*self_energy*G
    ! K=G0*self_energy
    call allocate_KB(parm_,kernel)
    call convolute_KB(parm_,n,dmft_%weiss_greenA,dmft_%self_energyA,kernel)
    
    ! G=[1-K]^{-1}*G0
    call Volterra_int(parm_,n,dmft_%weiss_greenA,kernel,dmft_%local_greenA)
    call deallocate_KB(kernel)

    dmft_%nA = imag(dmft_%local_greenA%lesser(n,n))
    !print*, "after solver A nA", dmft_%nA
    ! this thing will be neded to calculate intrigal diiferential equationin non eq problem and its needed in initialisation too
    t2=dble(n-1)*parm_%dt
    dmft_%hA = delta_(parm_,t2) + U(parm_,t2)*dmft_%nA - U(parm_,t2)/2.0
  end subroutine noneq_iptA

  subroutine noneq_iptB(parm_,dmft_)!TODO carefully chek the self energy for B sublattice and local GF sroe in B sublattice
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
       dmft_%self_energy_TB%c12(n,j)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TB%c12(n,j)*dmft_%weiss_green_TB%c21(j,n)*dmft_%weiss_green_TB%c12(n,j)
       dmft_%self_energy_TB%c21(n,j)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TB%c21(n,j)*dmft_%weiss_green_TB%c12(j,n)*dmft_%weiss_green_TB%c21(n,j)
    end do
    t2=dble(n-1)*parm_%dt
    do i=1, n-1
       t1=dble(i-1)*parm_%dt
       dmft_%self_energy_TB%c12(i,n)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TB%c12(i,n)*dmft_%weiss_green_TB%c21(n,i)*dmft_%weiss_green_TB%c12(i,n)
       dmft_%self_energy_TB%c21(i,n)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TB%c21(i,n)*dmft_%weiss_green_TB%c12(n,i)*dmft_%weiss_green_TB%c21(i,n)
    end do
    t1=dble(n-1)*parm_%dt
    dmft_%self_energy_TB%c13(n,:)=U(parm_,t1)*parm_%U_i*dmft_%weiss_green_TB%c13(n,:)*dmft_%weiss_green_TB%c31(:,n)*dmft_%weiss_green_TB%c13(n,:)
    t2=dble(n-1)*parm_%dt
    dmft_%self_energy_TB%c31(:,n)=parm_%U_i*U(parm_,t2)*dmft_%weiss_green_TB%c31(:,n)*dmft_%weiss_green_TB%c13(n,:)*dmft_%weiss_green_TB%c31(:,n)
    
    ! Kadanoff-Baym self-energy <= contour-ordered self-energy
    do j=1, n
       dmft_%self_energyB%retarded(n,j)=dmft_%self_energy_TB%c21(n,j)-dmft_%self_energy_TB%c12(n,j)
       dmft_%self_energyB%lesser(n,j)=dmft_%self_energy_TB%c12(n,j)
    end do
    do i=1, n-1
       dmft_%self_energyB%lesser(i,n)=dmft_%self_energy_TB%c12(i,n)
    end do
    dmft_%self_energyB%left_mixing(n,:)=dmft_%self_energy_TB%c13(n,:)
    
    ! solve the Impurity Dyson equation: G = G0+G0*self_energy*G
    ! K=G0*self_energy
    call allocate_KB(parm_,kernel)
    call convolute_KB(parm_,n,dmft_%weiss_greenB,dmft_%self_energyB,kernel)
    
    ! G=[1-K]^{-1}*G0
    call Volterra_int(parm_,n,dmft_%weiss_greenB,kernel,dmft_%local_greenB)
    call deallocate_KB(kernel)
    
    dmft_%nB = imag(dmft_%local_greenB%lesser(n,n))
    ! this thing will be neded to calculate intrigal diiferential equationin non eq problem and its needed in initialisation too
    t2=dble(n-1)*parm_%dt
    dmft_%hB = -delta_(parm_,t2) + U(parm_,t2)*dmft_%nB - U(parm_,t2)/2.0
  end subroutine noneq_iptB

  
end module IPT_MOD
