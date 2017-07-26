!----------------------------------------------------------------------------
!
!  author: Naoto Tsuji <tsuji@cms.phys.s.u-tokyo.ac.jp>
!
!          Department of Physics, University of Tokyo
!
!  date:   February 28, 2013
!
!----------------------------------------------------------------------------
module INTEGRAL_MOD
  
  use CONST_MOD
  implicit none
  
  private
  public :: trapezoid_d, trapezoid_z, trapezoid_half_edge
  
contains
  
  double precision function trapezoid_d(f,i,j)
    !----------------------------------------------------------------------------
    !  This function calculates the sum
    !
    !    \sum_{k=i}^{j} w_{k}^{i,j} f_{k}.
    !
    !    w_{k}^{i,j} = 1/2  for k=i,j
    !                = 1    for i<k<j
    !----------------------------------------------------------------------------
    double precision, intent(in) :: f(:)
    integer, intent(in)          :: i, j
    integer                      :: k
    double precision             :: integral
    
    if (j==i) then
       trapezoid_d=0.d0
    else
       integral=0.d0
       integral=integral+0.5d0*f(i)
       do k=i+1, j-1
          integral=integral+f(k)
       end do
       integral=integral+0.5d0*f(j)
       trapezoid_d=integral
    end if
  end function trapezoid_d
  
  
  complex(kind(0d0)) function trapezoid_z(f,i,j)
    !----------------------------------------------------------------------------
    !  This function calculates the sum
    !
    !    \sum_{k=i}^{j} w_{k}^{i,j} f_{k}.
    !
    !    w_{k}^{i,j} = 1/2  for k=i,j
    !                = 1    for i<k<j
    !----------------------------------------------------------------------------
    complex(kind(0d0)), intent(in) :: f(:)
    integer, intent(in)            :: i, j
    integer                        :: k
    complex(kind(0d0))             :: integral
    
    if (j==i) then
       trapezoid_z=0.d0
    else
       integral=0.d0
       integral=integral+0.5d0*f(i)
       do k=i+1, j-1
          integral=integral+f(k)
       end do
       integral=integral+0.5d0*f(j)
       trapezoid_z=integral
    end if
  end function trapezoid_z
  
  
  complex(kind(0d0)) function trapezoid_half_edge(f,i,j)
    !----------------------------------------------------------------------------
    !  This function calculates the sum
    !
    !    \sum_{k=i}^{j} w_{k}^{i,j} f_{k}.
    !
    !    w_{k}^{i,j} = 1/2  for k=i
    !                = 1    for i<k<=j
    !----------------------------------------------------------------------------
    complex(kind(0d0)), intent(in) :: f(:)
    integer, intent(in)            :: i, j
    integer                        :: k
    complex(kind(0d0))             :: integral
    
    integral=0.5d0*f(i)
    do k=i+1, j
       integral=integral+f(k)
    end do
    trapezoid_half_edge=integral
    
  end function trapezoid_half_edge
  
  
end module INTEGRAL_MOD
