!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  data:   February 28, 2013
!
!----------------------------------------------------------------------------
module CONST_MOD
  
  implicit none
  
  public :: pi, xj
  
  double precision, parameter   :: pi=3.14159265358979d0
  complex(kind(0d0)), parameter :: xj=(0.d0,1.d0)
  
end module CONST_MOD
