subroutine initialize_noneq_greenB(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: i, j, k, n
    complex(kind(0d0)), allocatable :: SxG(:)
    
    n=dmft_%time_step
    ! initialize G0(t,t') and G(t,t')
    if (n==2) then
       dmft_%weiss_greenB%retarded(1,1)=-xj
       dmft_%local_greenB%retarded(1,1)=-xj
       dmft_%weiss_green_newB%retarded(1,1)=-xj
       do j=1, parm_%N_tau+1
          dmft_%weiss_greenB%left_mixing(1,j)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-j+2)
          dmft_%local_greenB%left_mixing(1,j)=-xj*dmft_%local_greenB%matsubara_t(parm_%N_tau-j+2)
          dmft_%weiss_green_newB%left_mixing(1,j)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-j+2)
       end do
       dmft_%weiss_greenB%lesser(1,1)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau+1)
       dmft_%local_greenB%lesser(1,1)=-xj*dmft_%local_greenB%matsubara_t(parm_%N_tau+1)
       dmft_%weiss_green_newB%lesser(1,1)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau+1)
       
       if (parm_%dos=='semicircular') then
          ! initialize d/dt G0(t,t') = -i*(G*G0)(t,t')
          dmft_%weiss_green_derB%retarded(1)=0.d0
          dmft_%weiss_green_derB%left_mixing(:)=0.d0
          allocate(SxG(parm_%N_tau+1))
          do j=1, parm_%N_tau+1
             do k=1, j
                SxG(k)=dmft_%local_greenB%left_mixing(1,k)*dmft_%weiss_greenB%matsubara_t(parm_%N_tau+k-j+1)
             end do
             dmft_%weiss_green_derB%left_mixing(j)=dmft_%weiss_green_derB%left_mixing(j)+xj*parm_%dtau*trapezoid_z(SxG,1,j)
             do k=j, parm_%N_tau+1
                SxG(k)=dmft_%local_greenB%left_mixing(1,k)*dmft_%weiss_greenB%matsubara_t(k-j+1)
             end do
             dmft_%weiss_green_derB%left_mixing(j)=dmft_%weiss_green_derB%left_mixing(j)-xj*parm_%dtau*trapezoid_z(SxG,j,parm_%N_tau+1)
          end do
          do k=1, parm_%N_tau+1
             SxG(k)=dmft_%local_greenB%left_mixing(1,k)*conjg(dmft_%weiss_greenB%left_mixing(1,parm_%N_tau-k+2))
          end do
          dmft_%weiss_green_derB%lesser(1)=-xj*(-xj)*parm_%dtau*trapezoid_z(SxG,1,parm_%N_tau+1)
          deallocate(SxG)
       end if
       
       ! guess G0(t,t') in the next step
       do j=1, n
          dmft_%weiss_greenB%retarded(n,j)=dmft_%weiss_greenB%retarded(1,1)
       end do
       do j=1, parm_%N_tau+1
          dmft_%weiss_greenB%left_mixing(n,j)=dmft_%weiss_greenB%left_mixing(1,j)
       end do
       do j=1, n
          dmft_%weiss_greenB%lesser(n,j)=dmft_%weiss_greenB%lesser(1,1)
       end do
       do i=1, n-1
          dmft_%weiss_greenB%lesser(i,n)=dmft_%weiss_greenB%lesser(1,1)
       end do
       
       ! contour-ordered G0 <= Kadanoff-Baym G0
       do j=1, n
          do i=1, n
             dmft_%weiss_green_TB%c12(i,j)=dmft_%weiss_greenB%lesser(i,j)
             dmft_%weiss_green_TB%c21(i,j)=dmft_%weiss_greenB%lesser(i,j)-xj
          end do
       end do
       do j=1, parm_%N_tau+1
          do i=1, n
             dmft_%weiss_green_TB%c13(i,j)=-xj*dmft_%weiss_greenB%matsubara_t(parm_%N_tau-j+2)
             dmft_%weiss_green_TB%c31(j,i)=xj*dmft_%weiss_greenB%matsubara_t(j)
          end do
       end do
    else if (n>=3) then
       ! guess G0(t,t') in the next step by quadratic extrapolation
       dmft_%weiss_greenB%retarded(n,n)=-xj
       do k=1, n-2
          dmft_%weiss_greenB%retarded(n,k)=2.d0*dmft_%weiss_greenB%retarded(n-1,k)-dmft_%weiss_greenB%retarded(n-2,k)
       end do
       dmft_%weiss_greenB%retarded(n,n-1)=0.5d0*(dmft_%weiss_greenB%retarded(n,n)+dmft_%weiss_greenB%retarded(n,n-2))
       do k=1, parm_%N_tau+1
          dmft_%weiss_greenB%left_mixing(n,k)=2.d0*dmft_%weiss_greenB%left_mixing(n-1,k)-dmft_%weiss_greenB%left_mixing(n-2,k)
       end do
       do k=1, n-1
          dmft_%weiss_greenB%lesser(n,k)=2.d0*dmft_%weiss_greenB%lesser(n-1,k)-dmft_%weiss_greenB%lesser(n-2,k)
          dmft_%weiss_greenB%lesser(k,n)=2.d0*dmft_%weiss_greenB%lesser(k,n-1)-dmft_%weiss_greenB%lesser(k,n-2)
       end do
       dmft_%weiss_greenB%lesser(n,n)=2.d0*dmft_%weiss_greenB%lesser(n-1,n-1)-dmft_%weiss_greenB%lesser(n-2,n-2)
       
       ! contour-ordered G0 <= Kadanoff-Baym G0
       do k=1, n-1
          dmft_%weiss_green_TB%c12(n,k)=dmft_%weiss_greenB%lesser(n,k)
          dmft_%weiss_green_TB%c12(k,n)=dmft_%weiss_greenB%lesser(k,n)
          dmft_%weiss_green_TB%c21(n,k)=dmft_%weiss_greenB%lesser(n,k)+dmft_%weiss_greenB%retarded(n,k)
          dmft_%weiss_green_TB%c21(k,n)=dmft_%weiss_greenB%lesser(k,n)-conjg(dmft_%weiss_greenB%retarded(n,k))
       end do
       dmft_%weiss_green_TB%c12(n,n)=dmft_%weiss_greenB%lesser(n,n)
       dmft_%weiss_green_TB%c21(n,n)=dmft_%weiss_greenB%lesser(n,n)-xj
       do k=1, parm_%N_tau+1
          dmft_%weiss_green_TB%c13(n,k)=dmft_%weiss_greenB%left_mixing(n,k)
          dmft_%weiss_green_TB%c31(k,n)=conjg(dmft_%weiss_greenB%left_mixing(n,parm_%N_tau-k+2))
       end do
    end if
    !TODO need to incorporate h in the initialised weiss_der calculation at t=0
  end subroutine initialize_noneq_greenB

