!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  date:   February 28, 2013
!
!----------------------------------------------------------------------------
program NONEQ_DMFT_HUBBARD_IPT
  
  use PARM_MOD
  use DMFT_MOD
  use EQ_DMFT_MOD
  use NONEQ_DMFT_MOD
  implicit none
  
  type(parm) :: parm_
  type(dmft) :: dmft_
  
  call set_parm(parm_)
  call construct_dmft(parm_,dmft_)
  call initialize_output_file()
  call start_eq_dmft(parm_,dmft_)
  call start_noneq_dmft(parm_,dmft_)
  call destruct_dmft(parm_,dmft_)
  
end program NONEQ_DMFT_HUBBARD_IPT

