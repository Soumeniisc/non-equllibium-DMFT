!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  date:   February 28, 2013
!
!----------------------------------------------------------------------------
module GREEN_MOD
  
  use CONST_MOD
  use PARM_MOD
  use INTEGRAL_MOD
  implicit none
  
  private
  public :: allocate_KB,              &
            deallocate_KB,            &
            allocate_KB_derivative,   &
            deallocate_KB_derivative, &
            allocate_contour,         &
            deallocate_contour,       &
            convolute_KB,             &
            Volterra_int,             &
            Volterra_intdiff
  
  ! Kadanoff-Baym Green function: G(t,t')
  type, public :: KB
     complex(kind(0d0)), pointer :: matsubara_w(:)   ! Matsubara Green function G^{M}(w)
     double precision, pointer   :: matsubara_t(:)   ! Matsubara Green function G^{M}(tau)
     complex(kind(0d0)), pointer :: retarded(:,:)    ! retarded Green function G^{R}(t,t')
     complex(kind(0d0)), pointer :: left_mixing(:,:) ! left-mixing Green function G^{Left}(t,tau')
     complex(kind(0d0)), pointer :: lesser(:,:)      ! lesser Green function G^{<}(t,t')
  end type KB
  
  ! derivative of Kadanoff-Baym Green function: d/dt G(t,t')
  type, public :: KB_derivative
     complex(kind(0d0)), pointer :: retarded(:)    ! derivative of retarded Green function d/dt G^{R}(t_{n},t_{j})
     complex(kind(0d0)), pointer :: left_mixing(:) ! derivative of left-mixing Green function d/dt G^{Left}(t_{n},tau_{j})
     complex(kind(0d0)), pointer :: lesser(:)      ! derivative of lesser Green function d/dt G^{<}(t_{n},t_{j})
  end type KB_derivative
  
  ! contour-ordered Green function: G(t,t')
  type, public :: contour
     complex(kind(0d0)), pointer :: c12(:,:) ! G^{12}(t,t')
     complex(kind(0d0)), pointer :: c21(:,:) ! G^{21}(t,t')
     complex(kind(0d0)), pointer :: c13(:,:) ! G^{13}(t,tau')
     complex(kind(0d0)), pointer :: c31(:,:) ! G^{31}(tau,t')
  end type contour
  
contains
  
  subroutine allocate_KB(parm_,KB_)
    type(parm), intent(in) :: parm_
    type(KB), intent(out)  :: KB_
    
    allocate(KB_%matsubara_w(parm_%N_tau))
    allocate(KB_%matsubara_t(parm_%N_tau+1))
    allocate(KB_%retarded(parm_%N_t+1,parm_%N_t+1))
    allocate(KB_%left_mixing(parm_%N_t+1,parm_%N_tau+1))
    allocate(KB_%lesser(parm_%N_t+1,parm_%N_t+1))
    
  end subroutine allocate_KB
  
  
  subroutine deallocate_KB(KB_)
    type(KB), intent(out) :: KB_
    
    deallocate(KB_%matsubara_w)
    deallocate(KB_%matsubara_t)
    deallocate(KB_%retarded)
    deallocate(KB_%left_mixing)
    deallocate(KB_%lesser)
    
  end subroutine deallocate_KB
  
  
  subroutine allocate_KB_derivative(parm_,KB_der_,n)
    type(parm), intent(in)           :: parm_
    type(KB_derivative), intent(out) :: KB_der_
    integer                          :: n
    
    allocate(KB_der_%retarded(n))
    allocate(KB_der_%left_mixing(parm_%N_tau+1))
    allocate(KB_der_%lesser(n))
    
    KB_der_%retarded=0.d0
    KB_der_%left_mixing=0.d0
    KB_der_%lesser=0.d0
    
  end subroutine allocate_KB_derivative
  
  
  subroutine deallocate_KB_derivative(KB_der_)
    type(KB_derivative), intent(out) :: KB_der_
    
    deallocate(KB_der_%retarded)
    deallocate(KB_der_%left_mixing)
    deallocate(KB_der_%lesser)
    
  end subroutine deallocate_KB_derivative
  
  
  subroutine allocate_contour(parm_,contour_)
    type(parm), intent(in)     :: parm_
    type(contour), intent(out) :: contour_
    
    allocate(contour_%c12(parm_%N_t+1,parm_%N_t+1))
    allocate(contour_%c21(parm_%N_t+1,parm_%N_t+1))
    allocate(contour_%c13(parm_%N_t+1,parm_%N_tau+1))
    allocate(contour_%c31(parm_%N_tau+1,parm_%N_t+1))
    
    contour_%c12=0.d0
    contour_%c21=0.d0
    contour_%c13=0.d0
    contour_%c31=0.d0
    
  end subroutine allocate_contour
  
  
  subroutine deallocate_contour(contour_)
    type(contour), intent(out) :: contour_
    
    deallocate(contour_%c12)
    deallocate(contour_%c21)
    deallocate(contour_%c13)
    deallocate(contour_%c31)
    
  end subroutine deallocate_contour
  
  
  subroutine convolute_KB(parm_,n,A,B,C)
    !----------------------------------------------------------------------------
    !  This subroutine calculates a convolution of two Kadanoff-Baym function A 
    !  and B,
    !              C(t,t')=(A*B)(t,t'),
    !
    !  for t=n*dt or t'=n*dt. Integrals are evaluated by the trapezoidal rule.
    !----------------------------------------------------------------------------
    type(parm), intent(in)             :: parm_
    integer, intent(in)                :: n       ! time step
    type(KB), intent(in)               :: A, B
    type(KB), intent(inout)            :: C
    integer                            :: i, j, k
    complex(kind(0d0)), allocatable    :: AxB(:)
    
    ! retarded component
    ! C^{R}(t,t')=\int_{t'}^{t} ds A^{R}(t,s)*B^{R}(s,t')
    allocate(AxB(max(parm_%N_tau+1,n)))
    do j=1, n
       C%retarded(n,j)=0.d0
       do k=j, n
          AxB(k)=A%retarded(n,k)*B%retarded(k,j)
       end do
       C%retarded(n,j)=C%retarded(n,j)+parm_%dt*trapezoid_z(AxB,j,n)
    end do
    
    ! left-mixing component
    ! C^{Left}(t,tau')=\int_0^{beta} dtau A^{Left}(t,tau)*B^{M}(tau,tau')
    !                  +\int_0^{t} ds A^{R}(t,s)*B^{Left}(s,tau')
    do j=1, parm_%N_tau+1
       C%left_mixing(n,j)=0.d0
       do k=1, j
          AxB(k)=A%left_mixing(n,k)*B%matsubara_t(parm_%N_tau+k-j+1)
       end do
       C%left_mixing(n,j)=C%left_mixing(n,j)-parm_%dtau*trapezoid_z(AxB,1,j)
       do k=j, parm_%N_tau+1
          AxB(k)=A%left_mixing(n,k)*B%matsubara_t(k-j+1)
       end do
       C%left_mixing(n,j)=C%left_mixing(n,j)+parm_%dtau*trapezoid_z(AxB,j,parm_%N_tau+1)
       do k=1, n
          AxB(k)=A%retarded(n,k)*B%left_mixing(k,j)
       end do
       C%left_mixing(n,j)=C%left_mixing(n,j)+parm_%dt*trapezoid_z(AxB,1,n)
    end do
    
    ! lesser component
    ! C^{<}(t,t')=-i\int_0^{beta} dtau A^{Left}(t,tau)*B^{Right}(tau,t')
    !             +\int_0^{t'} ds A^{<}(t,s)*B^{A}(s,t')
    !             +\int_0^{t} ds A^{R}(t,s)*B^{<}(s,t')
    ! (i,j)=(n,j)
    do j=1, n-1
       C%lesser(n,j)=0.d0
       do k=1, parm_%N_tau+1
          AxB(k)=A%left_mixing(n,k)*conjg(B%left_mixing(j,parm_%N_tau-k+2))
       end do
       C%lesser(n,j)=C%lesser(n,j)-xj*parm_%dtau*trapezoid_z(AxB,1,parm_%N_tau+1)
       do k=1, j
          AxB(k)=A%lesser(n,k)*conjg(B%retarded(j,k))
       end do
       C%lesser(n,j)=C%lesser(n,j)+parm_%dt*trapezoid_z(AxB,1,j)
       do k=1, n
          AxB(k)=A%retarded(n,k)*B%lesser(k,j)
       end do
       C%lesser(n,j)=C%lesser(n,j)+parm_%dt*trapezoid_z(AxB,1,n)
    end do
    
    ! (i,j)=(i,n)
    do i=1, n
       C%lesser(i,n)=0.d0
       do k=1, parm_%N_tau+1
          AxB(k)=A%left_mixing(i,k)*conjg(B%left_mixing(n,parm_%N_tau-k+2))
       end do
       C%lesser(i,n)=C%lesser(i,n)-xj*parm_%dtau*trapezoid_z(AxB,1,parm_%N_tau+1)
       do k=1, n
          AxB(k)=A%lesser(i,k)*conjg(B%retarded(n,k))
       end do
       C%lesser(i,n)=C%lesser(i,n)+parm_%dt*trapezoid_z(AxB,1,n)
       do k=1, i
          AxB(k)=A%retarded(i,k)*B%lesser(k,n)
       end do
       C%lesser(i,n)=C%lesser(i,n)+parm_%dt*trapezoid_z(AxB,1,i)
    end do
    deallocate(AxB)
    
  end subroutine convolute_KB
  
  
  subroutine Volterra_int(parm_,n,G0,K,G)
    !----------------------------------------------------------------------------
    !  This subroutine solves a Volterra integral equation of the second kind,
    ! 
    !              G(t,t')=G0(t,t')+(K*G)(t,t'),
    !
    !  for t=n*dt or t'=n*dt. The integral equation is solved by the second-order 
    !  implicit Runge-Kutta method.
    !----------------------------------------------------------------------------
    type(parm), intent(in)             :: parm_
    integer, intent(in)                :: n       ! time step
    type(KB), intent(in)               :: G0, K
    type(KB), intent(out)              :: G
    integer                            :: i, j, l
    complex(kind(0d0)), allocatable    :: KxG(:)
    
    allocate(KxG(max(parm_%N_tau+1,n)))
    
    ! retarded Green function
    ! G^{R}(t,t')-\int_{t'}^{t} ds K^{R}(t,s)*G^{R}(s,t') = G0^{R}(t,t')
    G%retarded(n,n)=G0%retarded(n,n)
    do j=1, n-1
       G%retarded(n,j)=G0%retarded(n,j)
       do l=j, n-1
          KxG(l)=K%retarded(n,l)*G%retarded(l,j)
       end do
       G%retarded(n,j)=G%retarded(n,j)+parm_%dt*trapezoid_half_edge(KxG,j,n-1)
       G%retarded(n,j)=G%retarded(n,j)/(1.d0-0.5d0*parm_%dt*K%retarded(n,n))
    end do
    
    ! left-mixing Green function
    ! G^{Left}(t,tau')-\int_0^{t} ds K^{R}(t,s)*G^{Left}(s,tau')
    !    = G0^{Left}(t,tau')+\int_0^{beta} dtau K^{Left}(t,tau)*G^{M}(tau,tau')
    do j=1, parm_%N_tau+1
       G%left_mixing(n,j)=G0%left_mixing(n,j)
       do l=1, j
          KxG(l)=K%left_mixing(n,l)*G%matsubara_t(parm_%N_tau+l-j+1)
       end do
       G%left_mixing(n,j)=G%left_mixing(n,j)-parm_%dtau*trapezoid_z(KxG,1,j)
       do l=j, parm_%N_tau+1
          KxG(l)=K%left_mixing(n,l)*G%matsubara_t(l-j+1)
       end do
       G%left_mixing(n,j)=G%left_mixing(n,j)+parm_%dtau*trapezoid_z(KxG,j,parm_%N_tau+1)
       do l=1, n-1
          KxG(l)=K%retarded(n,l)*G%left_mixing(l,j)
       end do
       G%left_mixing(n,j)=G%left_mixing(n,j)+parm_%dt*trapezoid_half_edge(KxG,1,n-1)
       G%left_mixing(n,j)=G%left_mixing(n,j)/(1.d0-0.5d0*parm_%dt*K%retarded(n,n))
    end do
    
    ! lesser Green function
    ! G^{<}(t,t')-\int_0^{t} ds K^{R}(t,s)*G^{<}(s,t')
    !    = G0^{<}(t,t')-i\int_0^{beta} dtau K^{Left}(t,tau)*G^{Right}(tau,t')
    !      +\int_0^{t'} ds K^{<}(t,s)*G^{A}(s,t')
    ! G^{<}(t_{n},t_{j})
    do j=1, n-1
       G%lesser(n,j)=G0%lesser(n,j)
       do l=1, parm_%N_tau+1
          KxG(l)=K%left_mixing(n,l)*conjg(G%left_mixing(j,parm_%N_tau-l+2))
       end do
       G%lesser(n,j)=G%lesser(n,j)-xj*parm_%dtau*trapezoid_z(KxG,1,parm_%N_tau+1)
       do l=1, j
          KxG(l)=K%lesser(n,l)*conjg(G%retarded(j,l))
       end do
       G%lesser(n,j)=G%lesser(n,j)+parm_%dt*trapezoid_z(KxG,1,j)
       do l=1, n-1
          KxG(l)=K%retarded(n,l)*G%lesser(l,j)
       end do
       G%lesser(n,j)=G%lesser(n,j)+parm_%dt*trapezoid_half_edge(KxG,1,n-1)
       G%lesser(n,j)=G%lesser(n,j)/(1.d0-0.5d0*parm_%dt*K%retarded(n,n))
    end do
    ! Hermite conjugate
    ! G^{<}(t_{i},t_{n})
    do i=1, n-1
       G%lesser(i,n)=-conjg(G%lesser(n,i))
    end do
    ! G^{<}(t_{n},t_{n})
    G%lesser(n,n)=G0%lesser(n,n)
    do l=1, parm_%N_tau+1
       KxG(l)=K%left_mixing(n,l)*conjg(G%left_mixing(n,parm_%N_tau-l+2))
    end do
    G%lesser(n,n)=G%lesser(n,n)-xj*parm_%dtau*trapezoid_z(KxG,1,parm_%N_tau+1)
    do l=1, n
       KxG(l)=K%lesser(n,l)*conjg(G%retarded(n,l))
    end do
    G%lesser(n,n)=G%lesser(n,n)+parm_%dt*trapezoid_z(KxG,1,n)
    do l=1, n-1
       KxG(l)=K%retarded(n,l)*G%lesser(l,n)
    end do
    G%lesser(n,n)=G%lesser(n,n)+parm_%dt*trapezoid_half_edge(KxG,1,n-1)
    G%lesser(n,n)=G%lesser(n,n)/(1.d0-0.5d0*parm_%dt*K%retarded(n,n))
    
    deallocate(KxG)
    
  end subroutine Volterra_int
  
  
  subroutine Volterra_intdiff(parm_,n,h,K,G,G_der,G_der_new)
    !----------------------------------------------------------------------------
    !  This subroutine solves a Volterra integrodifferential equation of 
    !  the second kind,
    ! 
    !              [i*d/dt-h(t)]G(t,t')=delta(t,t')+(K*G)(t,t'),
    !
    !  for t=n*dt or t'=n*dt. The integrodifferential equation is solved by 
    !  the second-order implicit Runge-Kutta method.
    !----------------------------------------------------------------------------
    type(parm), intent(in)                 :: parm_
    integer, intent(in)                    :: n       ! time step
    complex(kind(0d0))                     :: h(:)
    type(KB), intent(in)                   :: K
    type(KB), intent(out)                  :: G
    type(KB_derivative), intent(in)        :: G_der
    type(KB_derivative), intent(out)       :: G_der_new
    integer                                :: i, j, l
    complex(kind(0d0)), allocatable        :: KxG(:)
    complex(kind(0d0))                     :: G_der_lesser ! d/dt G^{<}(t_{n-1},t_{n})
    
    allocate(KxG(max(parm_%N_tau+1,n)))
    
    ! retarded Green function
    ! d/dt G^{R}(t,t') = -i*delta(t,t')-i*h(t)*G^{R}(t,t')
    !                    -i\int_{t'}^{t} ds K^{R}(t,s)*G^{R}(s,t')
    G%retarded(n,n)=-xj
    G_der_new%retarded(n)=-xj*h(n)*G%retarded(n,n)
    do j=1, n-1
       G%retarded(n,j)=G%retarded(n-1,j)+0.5d0*parm_%dt*G_der%retarded(j)
       do l=j, n-1
          KxG(l)=K%retarded(n,l)*G%retarded(l,j)
       end do
       G_der_new%retarded(j)=-xj*parm_%dt*trapezoid_half_edge(KxG,j,n-1)
       G%retarded(n,j)=G%retarded(n,j)+0.5d0*parm_%dt*G_der_new%retarded(j)
       G%retarded(n,j)=G%retarded(n,j)/(1.d0+0.5d0*xj*parm_%dt*h(n)+0.25d0*xj*parm_%dt**2*K%retarded(n,n))
       G_der_new%retarded(j)=G_der_new%retarded(j)-xj*h(n)*G%retarded(n,j)-0.5d0*xj*parm_%dt*K%retarded(n,n)*G%retarded(n,j)
    end do
    
    ! left-mixing Green function
    ! d/dt G^{Left}(t,tau') = -i*h(t)*G^{Left}(t,tau')
    !                         -i\int_0^{beta} dtau K^{Left}(t,tau)*G^{M}(tau,tau')
    do j=1, parm_%N_tau+1
       G%left_mixing(n,j)=G%left_mixing(n-1,j)+0.5d0*parm_%dt*G_der%left_mixing(j)
       do l=1, j
          KxG(l)=K%left_mixing(n,l)*G%matsubara_t(parm_%N_tau+l-j+1)
       end do
       G_der_new%left_mixing(j)=xj*parm_%dtau*trapezoid_z(KxG,1,j)
       do l=j, parm_%N_tau+1
          KxG(l)=K%left_mixing(n,l)*G%matsubara_t(l-j+1)
       end do
       G_der_new%left_mixing(j)=G_der_new%left_mixing(j)-xj*parm_%dtau*trapezoid_z(KxG,j,parm_%N_tau+1)
       do l=1, n-1
          KxG(l)=K%retarded(n,l)*G%left_mixing(l,j)
       end do
       G_der_new%left_mixing(j)=G_der_new%left_mixing(j)-xj*parm_%dt*trapezoid_half_edge(KxG,1,n-1)
       G%left_mixing(n,j)=G%left_mixing(n,j)+0.5d0*parm_%dt*G_der_new%left_mixing(j)
       G%left_mixing(n,j)=G%left_mixing(n,j)/(1.d0+0.5d0*xj*parm_%dt*h(n)+0.25d0*xj*parm_%dt**2*K%retarded(n,n))
       G_der_new%left_mixing(j)=G_der_new%left_mixing(j)-xj*h(n)*G%left_mixing(n,j)-0.5d0*xj*parm_%dt*K%retarded(n,n)*G%left_mixing(n,j)
    end do
    
    ! lesser Green function
    ! d/dt G^{<}(t,t') = -i*h(t)*G^{<}(t,t')
    !                    -i*(-i)*\int_0^{beta} dtau K^{Left}(t,tau)*G^{Right}(tau,t')
    !                    -i*\int_0^{t'} ds K^{<}(t,s)*G^{A}(s,t')
    !                    -i*\int_0^{t} ds K^{R}(t,s)*G^{<}(s,t')
    !
    ! G^{<}(t_{n},t_{j}), d/dt G^{<}(t_{n},t_{j})
    do j=1, n-1
       G%lesser(n,j)=G%lesser(n-1,j)+0.5d0*parm_%dt*G_der%lesser(j)
       do l=1, parm_%N_tau+1
          KxG(l)=K%left_mixing(n,l)*conjg(G%left_mixing(j,parm_%N_tau-l+2))
       end do
       G_der_new%lesser(j)=-xj*(-xj)*parm_%dtau*trapezoid_z(KxG,1,parm_%N_tau+1)
       do l=1, j
          KxG(l)=K%lesser(n,l)*conjg(G%retarded(j,l))
       end do
       G_der_new%lesser(j)=G_der_new%lesser(j)-xj*parm_%dt*trapezoid_z(KxG,1,j)
       do l=1, n-1
          KxG(l)=K%retarded(n,l)*G%lesser(l,j)
       end do
       G_der_new%lesser(j)=G_der_new%lesser(j)-xj*parm_%dt*trapezoid_half_edge(KxG,1,n-1)
       G%lesser(n,j)=G%lesser(n,j)+0.5d0*parm_%dt*G_der_new%lesser(j)
       G%lesser(n,j)=G%lesser(n,j)/(1.d0+0.5d0*xj*parm_%dt*h(n)+0.25d0*xj*parm_%dt**2*K%retarded(n,n))
       G_der_new%lesser(j)=G_der_new%lesser(j)-xj*h(n)*G%lesser(n,j)-0.5d0*xj*parm_%dt*K%retarded(n,n)*G%lesser(n,j)
    end do
    ! Hermite conjugate
    do i=1, n-1
       G%lesser(i,n)=-conjg(G%lesser(n,i))
    end do
    ! d/dt G^{<}(t_{n-1},t_{n})
    G_der_lesser=-xj*h(n-1)*G%lesser(n-1,n)
    do l=1, parm_%N_tau+1
       KxG(l)=K%left_mixing(n-1,l)*conjg(G%left_mixing(n,parm_%N_tau-l+2))
    end do
    G_der_lesser=G_der_lesser-xj*(-xj)*parm_%dtau*trapezoid_z(KxG,1,parm_%N_tau+1)
    do l=1, n
       KxG(l)=K%lesser(n-1,l)*conjg(G%retarded(n,l))
    end do
    G_der_lesser=G_der_lesser-xj*parm_%dt*trapezoid_z(KxG,1,n)
    do l=1, n-1
       KxG(l)=K%retarded(n-1,l)*G%lesser(l,n)
    end do
    G_der_lesser=G_der_lesser-xj*parm_%dt*trapezoid_z(KxG,1,n-1)
    ! G^{<}(t_{n},t_{n}), d/dt G^{<}(t_{n},t_{n})
    G%lesser(n,n)=G%lesser(n-1,n)+0.5d0*parm_%dt*G_der_lesser
    do l=1, parm_%N_tau+1
       KxG(l)=K%left_mixing(n,l)*conjg(G%left_mixing(n,parm_%N_tau-l+2))
    end do
    G_der_new%lesser(n)=-xj*(-xj)*parm_%dtau*trapezoid_z(KxG,1,parm_%N_tau+1)
    do l=1, n
       KxG(l)=K%lesser(n,l)*conjg(G%retarded(n,l))
    end do
    G_der_new%lesser(n)=G_der_new%lesser(n)-xj*parm_%dt*trapezoid_z(KxG,1,n)
    do l=1, n-1
       KxG(l)=K%retarded(n,l)*G%lesser(l,n)
    end do
    G_der_new%lesser(n)=G_der_new%lesser(n)-xj*parm_%dt*trapezoid_half_edge(KxG,1,n-1)
    G%lesser(n,n)=G%lesser(n,n)+0.5d0*parm_%dt*G_der_new%lesser(n)
    G%lesser(n,n)=G%lesser(n,n)/(1.d0+0.5d0*xj*parm_%dt*h(n)+0.25d0*xj*parm_%dt**2*K%retarded(n,n))
    G_der_new%lesser(n)=G_der_new%lesser(n)-xj*h(n)*G%lesser(n,n)-0.5d0*xj*parm_%dt*K%retarded(n,n)*G%lesser(n,n)
    
    deallocate(KxG)
    
  end subroutine Volterra_intdiff
  
  
end module GREEN_MOD
