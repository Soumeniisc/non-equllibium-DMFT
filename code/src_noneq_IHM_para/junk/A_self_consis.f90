  subroutine noneq_dmft_self_consistencyA(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    type(KB)                        :: kernel
    integer                         :: i, j, k, n
    complex(kind(0d0)), allocatable :: h(:)
    
    n=dmft_%time_step
    
    if (parm_%dos=='semicircular') then
       ! solve the lattice Dyson equation: id/dt G0(t,t')-G*G0(t,t')=delta(t,t')
       allocate(h(n))
       h(:)=0.d0
       call Volterra_intdiff(parm_,n,h,dmft_%local_greenB,dmft_%weiss_green_newA,dmft_%weiss_green_derA,dmft_%weiss_green_der_newA)
       deallocate(h)
    end if
    
    ! evaluate |G0_{new}-G0_{old}|
    dmft_%G0_diffA=0.d0
    do j=1, n
       dmft_%G0_diffA=dmft_%G0_diffA+abs(dmft_%weiss_green_newA%retarded(n,j)-dmft_%weiss_greenA%retarded(n,j))
       dmft_%G0_diffA=dmft_%G0_diffA+abs(dmft_%weiss_green_newA%lesser(n,j)-dmft_%weiss_greenA%lesser(n,j))
    end do
    do j=1, parm_%N_tau+1
       dmft_%G0_diffA=dmft_%G0_diffA+abs(dmft_%weiss_green_newA%left_mixing(n,j)-dmft_%weiss_greenA%left_mixing(n,j))
    end do
    do i=1, n-1
       dmft_%G0_diffA=dmft_%G0_diffA+abs(dmft_%weiss_green_newA%lesser(i,n)-dmft_%weiss_greenA%lesser(i,n))
    end do
    write (*,'(A,I2)') '  Iteration #', dmft_%iteration
    write (*,'(A,f15.10)') ' |G0_new-G0_old| =', dmft_%G0_diffA
    
    ! G0_{old} <= G0_{new}
    do j=1, n
       dmft_%weiss_greenA%retarded(n,j)=dmft_%weiss_green_newA%retarded(n,j)
       dmft_%weiss_greenA%lesser(n,j)=dmft_%weiss_green_newA%lesser(n,j)
    end do
    do i=1, n-1
       dmft_%weiss_greenA%lesser(i,n)=dmft_%weiss_green_newA%lesser(i,n)
    end do
    dmft_%weiss_greenA%left_mixing(n,:)=dmft_%weiss_green_newA%left_mixing(n,:)
    
    ! Kadanoff-Baym G0 => contour-ordered G0
    do j=1, n
       dmft_%weiss_green_TA%c12(n,j)=dmft_%weiss_greenA%lesser(n,j)
       dmft_%weiss_green_TA%c21(n,j)=dmft_%weiss_greenA%lesser(n,j)+dmft_%weiss_greenA%retarded(n,j)
    end do
    do i=1, n-1
       dmft_%weiss_green_TA%c12(i,n)=dmft_%weiss_greenA%lesser(i,n)
       dmft_%weiss_green_TA%c21(i,n)=dmft_%weiss_greenA%lesser(i,n)-conjg(dmft_%weiss_greenA%retarded(n,i))
    end do
    do j=1, parm_%N_tau+1
       dmft_%weiss_green_TA%c13(n,j)=dmft_%weiss_greenA%left_mixing(n,j)
    end do
    do i=1, parm_%N_tau+1
       dmft_%weiss_green_TA%c31(i,n)=conjg(dmft_%weiss_greenA%left_mixing(n,parm_%N_tau-i+2))
    end do
    
    ! Hermite conjugate
    do j=1, n
       dmft_%weiss_green_TA%c12(n,j)=0.5d0*(dmft_%weiss_green_TA%c12(n,j)-conjg(dmft_%weiss_green_TA%c12(j,n)))
       dmft_%weiss_green_TA%c21(n,j)=0.5d0*(dmft_%weiss_green_TA%c21(n,j)-conjg(dmft_%weiss_green_TA%c21(j,n)))
    end do
    do i=1, n-1
       dmft_%weiss_green_TA%c12(i,n)=-conjg(dmft_%weiss_green_TA%c12(n,i))
       dmft_%weiss_green_TA%c21(i,n)=-conjg(dmft_%weiss_green_TA%c21(n,i))
    end do
    do j=1, parm_%N_tau+1
       dmft_%weiss_green_TA%c13(n,j)=0.5d0*(dmft_%weiss_green_TA%c13(n,j)+conjg(dmft_%weiss_green_TA%c31(parm_%N_tau-j+2,n)))
    end do
    do i=1, parm_%N_tau+1
       dmft_%weiss_green_TA%c31(i,n)=conjg(dmft_%weiss_green_TA%c13(n,parm_%N_tau-i+2))
    end do
    
  end subroutine noneq_dmft_self_consistencyA
