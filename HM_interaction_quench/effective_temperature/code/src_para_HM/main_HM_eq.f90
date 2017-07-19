program NONEQ_DMFT_HUBBARD_IPT
  
  use PARM_MOD
  use DMFT_MOD
  use EQ_DMFT_MOD
  implicit none
  type(parm) :: parm_
  type(dmft) :: dmft_
  print*,"main is started2"
  call set_parm(parm_)
  print*, "parameter is set"
  call construct_dmft(parm_,dmft_)
  print*, "dmft is constructed"
  call initialize_output_file()
  print*, "output files are intialised"
  call start_eq_dmft(parm_,dmft_)
  !call destruct_dmft(parm_,dmft_)
  
end program NONEQ_DMFT_HUBBARD_IPT

