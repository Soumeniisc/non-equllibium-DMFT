!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  date:   February 28, 2013
!
!----------------------------------------------------------------------------
module PARM_MOD
  
  use CONST_MOD
  implicit none
  
  private
  public :: set_parm, U
  
  type, public :: parm
     character(100)     :: dos              ! density of states
     double precision   :: e2               ! <e^2>=\int de D(e)*e^2
     double precision   :: temperature      ! temperature of the initial thermal state
     double precision   :: beta             ! inverse temperature (beta=1/temperature)
     double precision   :: t_max            ! maximum time up to which the system evolves
     double precision   :: U_i              ! initial value of the interaction parameter U
     double precision   :: U_f              ! interaction quench from U=U_i to U_f
     integer            :: N_tau            ! # of imaginary-time steps
     integer            :: N_t              ! # of real-time steps
     integer            :: N_e              ! # of energy steps
     integer            :: N_iter           ! # of maximum iterations of DMFT self-consistency
     double precision   :: dtau             ! dtau=beta/N_tau
     double precision   :: dt               ! dt=t_max/N_t
     double precision   :: tolerance        ! tolerance for DMFT convergence
     character(100)     :: solver           ! impurity solver
     integer            :: init             ! # 0 when G0 is non-interacting GF. 1 when G0 is interacting GF of previous U value
     double precision   :: delta            ! # ionic potetial
     double precision   :: mix              ! # G = (1-mix)*G_new + G except first iteration

  end type parm
  
contains
  
  subroutine set_parm(parm_)
    type(parm), intent(inout) :: parm_
    character(100)            :: arg
    integer                   :: n_arg
    
    call getarg(1,parm_%dos)
    if (parm_%dos=='semicircular') then
       parm_%e2=0.d0
    end if
    call getarg(2,arg)
    read(arg,*) parm_%beta
    parm_%temperature=1.d0/parm_%beta
    call getarg(3,arg)
    read(arg,*) parm_%t_max
    call getarg(4,arg)
    read(arg,*) parm_%U_i
    call getarg(5,arg)
    read(arg,*) parm_%U_f
    call getarg(6,arg)
    read(arg,*) parm_%N_tau
    parm_%dtau=parm_%beta/dble(parm_%N_tau)
    call getarg(7,arg)
    read(arg,*) parm_%N_t
    parm_%dt=parm_%t_max/dble(parm_%N_t)
    call getarg(8,arg)
    read(arg,*) parm_%N_e
    call getarg(9,arg)
    read(arg,*) parm_%N_iter
    call getarg(10,arg)
    read(arg,*) parm_%tolerance
    call getarg(11,parm_%solver)
    print*, parm_%solver,"------------------------"
    call getarg(12,arg)
    read(arg,*) parm_%init
    print*, parm_%init,"------init------------------"
    call getarg(13,arg)
    read(arg,*) parm_%mix
    !read(arg,*) parm_%init
    !parm_%init = 1
    call getarg(14,arg)
    read(arg,*) parm_%delta
    print*,"mixing--------------",parm_%mix

    
  end subroutine set_parm
  
  
  double precision function U(parm_,t)
    type(parm), intent(in)       :: parm_
    double precision, intent(in) :: t
    
    U=parm_%U_f
    
  end function U

  double precision function U1(parm_,t)
    type(parm), intent(in)       :: parm_
    double precision, intent(in) :: t
    
    if (t<3.0) then
	U1 = (parm_%U_f - parm_%U_i)*t/3.0
    else
    	U1=parm_%U_f
    end if   
  end function U1
  
  
end module PARM_MOD
