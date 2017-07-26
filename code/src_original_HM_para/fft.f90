!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  date:   February 28, 2013
!
!----------------------------------------------------------------------------
module FFT_MOD
  
  use CONST_MOD
  use PARM_MOD
  implicit none

  private
  public :: fft_t2w, &
            fft_w2t
  
contains
  
  subroutine fft_t2w(parm_,y,z)
    !----------------------------------------------------------------------------
    !  This subroutine computes the Fourier transformation
    !
    !     z(k) = sum_{j=1}^{n+1} w_{n+1,j}*y(j)*exp(i*w_{k}*tau_{j})*dtau.
    !     (k = 1, 2,..., n)
    !
    !     tau_{j} = (j-1)*dtau     (dtau = beta/n)
    !
    !     w_{k} = (2*k-n-1)*pi/beta
    !
    !     w_{n+1,j} = 1/2  for j=1, n+1
    !               = 1    for 2<=j<=n
    !----------------------------------------------------------------------------
    type(parm), intent(in)            :: parm_
    double precision, intent(in)      :: y(:)
    complex(kind(0d0)), intent(out)   :: z(:)
    complex(kind(0d0)), allocatable   :: w(:)
    integer                           :: plan(8), j, k
    
    include 'fftw3.f'
    allocate(w(parm_%N_tau))
    
    w(1)=0.5d0*y(1)
    do j=2, parm_%N_tau
       w(j)=y(j)*exp(-xj*dble(parm_%N_tau-1)*dble(j-1)*pi/dble(parm_%N_tau))
    end do
    call dfftw_plan_dft_1d(plan,parm_%N_tau,w,z,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)
    
    do k=1, parm_%N_tau
       z(k)=(z(k)+0.5d0*y(parm_%N_tau+1)*exp(-xj*dble(parm_%N_tau-1)*pi))*parm_%beta/dble(parm_%N_tau)
    end do
    deallocate(w)
    
  end subroutine fft_t2w
  
  
  subroutine fft_w2t(parm_,z,y)
    !----------------------------------------------------------------------------
    !  This subroutine computes the Fourier transformation
    !
    !     y(j) = 1/beta*sum_{k=1}^{n} z(k)*exp(-i*w_{k}*tau_{j}).
    !     (j = 1, 2,..., n+1)
    !
    !     tau_{j} = (j-1)*dtau     (dtau = beta/n)
    !
    !     w_{k} = (2*k-n-1)*pi/beta
    !----------------------------------------------------------------------------
    type(parm), intent(in)          :: parm_
    complex(kind(0d0)), intent(in)  :: z(:)
    double precision, intent(out)   :: y(:)
    complex(kind(0d0)), allocatable :: w(:)
    integer                         :: plan(8), j, k
    
    include 'fftw3.f'
    allocate(w(parm_%N_tau))
    
    call dfftw_plan_dft_1d(plan,parm_%N_tau,z,w,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)
    
    do j=1, parm_%N_tau
       y(j)=dble(w(j)*exp(xj*dble(parm_%N_tau-1)*dble(j-1)*pi/dble(parm_%N_tau)))/parm_%beta
    end do
    y(parm_%N_tau+1)=0.d0
    do k=1, parm_%N_tau
       y(parm_%N_tau+1)=y(parm_%N_tau+1)+dble(z(k)*exp(xj*dble(parm_%N_tau+1)*pi))/parm_%beta
    end do
    deallocate(w)
    
  end subroutine fft_w2t
  
  
end module FFT_MOD
