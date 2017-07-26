!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  data:   February 28, 2013
!
!----------------------------------------------------------------------------
! last modified by soume on 12/04/17 modied over equllibium IHM dmft dmft.f90 file in src directory. we defined all the dmft green function and self energy and measuring quentity and file for saving data diffent for A and B sunlatice. veribla name ended with B is for B soblatice.  one end with A is for A sublattice. 
module DMFT_MOD
  
  use CONST_MOD
  use GREEN_MOD
  use PARM_MOD
  use FFT_MOD
  use INTEGRAL_MOD
  implicit none
  
  private
  public :: construct_dmft,           &
            destruct_dmft,            &
            initialize_output_file
            
  type, public :: dmft
     type(KB)                     :: local_greenA         ! local Green function G(t,t')
     type(KB)                     :: weiss_greenA         ! Weiss Green function G0(t,t')
     type(KB)                     :: self_energyA         ! local self-energy Sigma(t,t')
     type(KB)                     :: weiss_green_newA     ! updated Weiss Green function G0(t,t')
     type(KB)                     :: local_greenB         ! local Green function G(t,t')
     type(KB)                     :: weiss_greenB         ! Weiss Green function G0(t,t')
     type(KB)                     :: self_energyB         ! local self-energy Sigma(t,t')
     type(KB)                     :: weiss_green_newB     ! updated Weiss Green function G0(t,t')

     type(KB_derivative)          :: weiss_green_derA     ! derivative of Weiss Green function d/dt G0(t,t')
     type(KB_derivative)          :: weiss_green_der_newA ! updated derivative of Weiss Green function d/dt G0(t,t')
     type(KB_derivative)          :: weiss_green_derB     ! derivative of Weiss Green function d/dt G0(t,t')
     type(KB_derivative)          :: weiss_green_der_newB ! updated derivative of Weiss Green function d/dt G0(t,t')

     type(contour)                :: weiss_green_TA       ! contour-ordered Weiss Green function G0(t,t')
     type(contour)                :: self_energy_TA       ! contour-ordered self-energy Sigma(t,t')
     type(contour)                :: weiss_green_TB       ! contour-ordered Weiss Green function G0(t,t')
     type(contour)                :: self_energy_TB       ! contour-ordered self-energy Sigma(t,t')

     double precision, pointer    :: densityA(:)          ! n(t)=<c^+(t)c(t)>
     double precision, pointer    :: double_occupancyA(:) ! d(t)=<n_up(t)*n_do(t)>
     double precision, pointer    :: kinetic_energyA(:)   ! E_{kin}(t)=2\sum_{k} E(k)*<c_{k}^+(t)c_{k}(t)>
     double precision, pointer    :: densityB(:)          ! n(t)=<c^+(t)c(t)>
     double precision, pointer    :: double_occupancyB(:) ! d(t)=<n_up(t)*n_do(t)>
     double precision, pointer    :: kinetic_energyB(:)   ! E_{kin}(t)=2\sum_{k} E(k)*<c_{k}^+(t)c_{k}(t)>


     double precision             :: G0_diffA            ! convergence measure |G0_{new}-G0_{old}|
     double precision             :: G0_diffB            ! convergence measure |G0_{new}-G0_{old}|

     integer                      :: time_step           ! n: t_n=(n-1)*dt
     integer                      :: iteration           ! # of DMFT iteration
     logical                      :: converged           ! DMFT self-consistency is converged or not

     double precision   :: nA               ! # occupency of A sublattice
     double precision   :: nB               ! # occupency of B sublattice
     double precision   :: hA               ! parm_%delta - U/2 +UnA ! chek whic U should be there. its initialised in start_noneq_dmft
     double precision   :: hB               ! -parm_%delta -U/2 +UnB ! chek whic U should be there.

  end type dmft
  
contains
  
  subroutine construct_dmft(parm_,dmft_)
    type(parm), intent(in)    :: parm_
    type(dmft), intent(inout) :: dmft_
    
    call allocate_KB(parm_,dmft_%local_greenA)
    call allocate_KB(parm_,dmft_%weiss_greenA)
    call allocate_KB(parm_,dmft_%self_energyA)
    call allocate_KB(parm_,dmft_%weiss_green_newA)
    call allocate_KB(parm_,dmft_%local_greenB)
    call allocate_KB(parm_,dmft_%weiss_greenB)
    call allocate_KB(parm_,dmft_%self_energyB)
    call allocate_KB(parm_,dmft_%weiss_green_newB)

    call allocate_KB_derivative(parm_,dmft_%weiss_green_derA,1)
    call allocate_KB_derivative(parm_,dmft_%weiss_green_der_newA,2)
    call allocate_KB_derivative(parm_,dmft_%weiss_green_derB,1)
    call allocate_KB_derivative(parm_,dmft_%weiss_green_der_newB,2)

    call allocate_contour(parm_,dmft_%weiss_green_TA)
    call allocate_contour(parm_,dmft_%self_energy_TA)
    call allocate_contour(parm_,dmft_%weiss_green_TB)
    call allocate_contour(parm_,dmft_%self_energy_TB)

    allocate(dmft_%densityA(parm_%n_t+1))
    allocate(dmft_%double_occupancyA(parm_%n_t+1))
    allocate(dmft_%kinetic_energyA(parm_%n_t+1))
    allocate(dmft_%densityB(parm_%n_t+1))
    allocate(dmft_%double_occupancyB(parm_%n_t+1))
    allocate(dmft_%kinetic_energyB(parm_%n_t+1))
    
  end subroutine construct_dmft
  
  
  subroutine destruct_dmft(parm_,dmft_)
    type(parm), intent(in)    :: parm_
    type(dmft), intent(inout) :: dmft_
    integer                   :: k
    
    call deallocate_KB(dmft_%local_greenA)
    call deallocate_KB(dmft_%weiss_greenA)
    call deallocate_KB(dmft_%self_energyA)
    call deallocate_KB(dmft_%weiss_green_newA)
    call deallocate_KB(dmft_%local_greenB)
    call deallocate_KB(dmft_%weiss_greenB)
    call deallocate_KB(dmft_%self_energyB)
    call deallocate_KB(dmft_%weiss_green_newB)

    call deallocate_KB_derivative(dmft_%weiss_green_derA)
    call deallocate_KB_derivative(dmft_%weiss_green_der_newA)
    call deallocate_KB_derivative(dmft_%weiss_green_derB)
    call deallocate_KB_derivative(dmft_%weiss_green_der_newB)

    call deallocate_contour(dmft_%weiss_green_TA)
    call deallocate_contour(dmft_%self_energy_TA)
    call deallocate_contour(dmft_%weiss_green_TB)
    call deallocate_contour(dmft_%self_energy_TB)

    deallocate(dmft_%densityA)
    deallocate(dmft_%double_occupancyA)
    deallocate(dmft_%kinetic_energyA)
    deallocate(dmft_%densityB)
    deallocate(dmft_%double_occupancyB)
    deallocate(dmft_%kinetic_energyB)
    
  end subroutine destruct_dmft
  
  
  subroutine initialize_output_file()
    
    open (unit=7, file="densityA", status="replace", action="write")
    open (unit=8, file="double-occupancyA", status="replace", action="write")
    open (unit=9, file="kinetic-energyA", status="replace", action="write")
    open (unit=10,file="interaction-energyA", status="replace", action="write")
    open (unit=11,file="total-energyA", status="replace", action="write")
    open (unit=71, file="densityB", status="replace", action="write")
    open (unit=81, file="double-occupancyB", status="replace", action="write")
    open (unit=91, file="kinetic-energyB", status="replace", action="write")
    open (unit=101,file="interaction-energyB", status="replace", action="write")
    open (unit=111,file="total-energyB", status="replace", action="write")

    open (unit=12,file="GfA", status="replace", action="write")
    open (unit=13,file="GfB", status="replace", action="write")

    close (7)
    close (8)
    close (9)
    close (10)
    close (11)
    close (71)
    close (81)
    close (91)
    close (101)
    close (111)

    close (12)
    close (13)
  end subroutine initialize_output_file
  
  
end module DMFT_MOD
