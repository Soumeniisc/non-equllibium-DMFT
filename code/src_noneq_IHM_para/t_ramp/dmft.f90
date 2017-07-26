!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  data:   February 28, 2013
!
!----------------------------------------------------------------------------
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
     type(KB)                     :: local_green         ! local Green function G(t,t')
     type(KB)                     :: weiss_green         ! Weiss Green function G0(t,t')
     type(KB)                     :: self_energy         ! local self-energy Sigma(t,t')
     type(KB)                     :: weiss_green_new     ! updated Weiss Green function G0(t,t')
     type(KB_derivative)          :: weiss_green_der     ! derivative of Weiss Green function d/dt G0(t,t')
     type(KB_derivative)          :: weiss_green_der_new ! updated derivative of Weiss Green function d/dt G0(t,t')
     type(contour)                :: weiss_green_T       ! contour-ordered Weiss Green function G0(t,t')
     type(contour)                :: self_energy_T       ! contour-ordered self-energy Sigma(t,t')
     double precision, pointer    :: density(:)          ! n(t)=<c^+(t)c(t)>
     double precision, pointer    :: double_occupancy(:) ! d(t)=<n_up(t)*n_do(t)>
     double precision, pointer    :: kinetic_energy(:)   ! E_{kin}(t)=2\sum_{k} E(k)*<c_{k}^+(t)c_{k}(t)>
     double precision             :: G0_diff             ! convergence measure |G0_{new}-G0_{old}|
     integer                      :: time_step           ! n: t_n=(n-1)*dt
     integer                      :: iteration           ! # of DMFT iteration
     logical                      :: converged           ! DMFT self-consistency is converged or not
  end type dmft
  
contains
  
  subroutine construct_dmft(parm_,dmft_)
    type(parm), intent(in)    :: parm_
    type(dmft), intent(inout) :: dmft_
    
    call allocate_KB(parm_,dmft_%local_green)
    call allocate_KB(parm_,dmft_%weiss_green)
    call allocate_KB(parm_,dmft_%self_energy)
    call allocate_KB(parm_,dmft_%weiss_green_new)
    call allocate_KB_derivative(parm_,dmft_%weiss_green_der,1)
    call allocate_KB_derivative(parm_,dmft_%weiss_green_der_new,2)
    call allocate_contour(parm_,dmft_%weiss_green_T)
    call allocate_contour(parm_,dmft_%self_energy_T)
    allocate(dmft_%density(parm_%n_t+1))
    allocate(dmft_%double_occupancy(parm_%n_t+1))
    allocate(dmft_%kinetic_energy(parm_%n_t+1))
    
  end subroutine construct_dmft
  
  
  subroutine destruct_dmft(parm_,dmft_)
    type(parm), intent(in)    :: parm_
    type(dmft), intent(inout) :: dmft_
    integer                   :: k
    
    call deallocate_KB(dmft_%local_green)
    call deallocate_KB(dmft_%weiss_green)
    call deallocate_KB(dmft_%self_energy)
    call deallocate_KB(dmft_%weiss_green_new)
    call deallocate_KB_derivative(dmft_%weiss_green_der)
    call deallocate_KB_derivative(dmft_%weiss_green_der_new)
    call deallocate_contour(dmft_%weiss_green_T)
    call deallocate_contour(dmft_%self_energy_T)
    deallocate(dmft_%density)
    deallocate(dmft_%double_occupancy)
    deallocate(dmft_%kinetic_energy)
    
  end subroutine destruct_dmft
  
  
  subroutine initialize_output_file()
    
    open (unit=7, file="density", status="replace", action="write")
    open (unit=8, file="double-occupancy", status="replace", action="write")
    open (unit=9, file="kinetic-energy", status="replace", action="write")
    open (unit=10,file="interaction-energy", status="replace", action="write")
    open (unit=11,file="total-energy", status="replace", action="write")
    open (unit=12,file="Gf", status="replace", action="write")
    close (7)
    close (8)
    close (9)
    close (10)
    close (11)
    close (12)
    
  end subroutine initialize_output_file
  
  
end module DMFT_MOD
