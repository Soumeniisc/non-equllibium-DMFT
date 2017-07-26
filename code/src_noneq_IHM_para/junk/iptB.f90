subroutine noneq_iptB(parm_,dmft_)!TODO carefully chek the self energy for B sublattice and local GF sroe in B sublattice
    type(parm), intent(in)     :: parm_
    type(dmft), intent(inout)  :: dmft_
    type(KB)                   :: kernel
    double precision           :: t1, t2
    integer                    :: i, j, n
    
    ! second-order perturbation in U at half filling
    ! self_energy(t,t') = U(t)*U(t')*G0(t,t')*G0(t',t)*G0(t,t')
    n=dmft_%time_step
    do j=1, n
       t2=dble(j-1)*parm_%dt
       t1=dble(n-1)*parm_%dt
       dmft_%self_energy_TB%c12(n,j)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TB%c12(n,j)*dmft_%weiss_green_TB%c21(j,n)*dmft_%weiss_green_TB%c12(n,j)
       dmft_%self_energy_TB%c21(n,j)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TB%c21(n,j)*dmft_%weiss_green_TB%c12(j,n)*dmft_%weiss_green_TB%c21(n,j)
    end do
    t2=dble(n-1)*parm_%dt
    do i=1, n-1
       t1=dble(i-1)*parm_%dt
       dmft_%self_energy_TB%c12(i,n)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TB%c12(i,n)*dmft_%weiss_green_TB%c21(n,i)*dmft_%weiss_green_TB%c12(i,n)
       dmft_%self_energy_TB%c21(i,n)=U(parm_,t1)*U(parm_,t2)*dmft_%weiss_green_TB%c21(i,n)*dmft_%weiss_green_TB%c12(n,i)*dmft_%weiss_green_TB%c21(i,n)
    end do
    t1=dble(n-1)*parm_%dt
    dmft_%self_energy_TB%c13(n,:)=U(parm_,t1)*parm_%U_i*dmft_%weiss_green_TB%c13(n,:)*dmft_%weiss_green_TB%c31(:,n)*dmft_%weiss_green_TB%c13(n,:)
    t2=dble(n-1)*parm_%dt
    dmft_%self_energy_TB%c31(:,n)=parm_%U_i*U(parm_,t2)*dmft_%weiss_green_TB%c31(:,n)*dmft_%weiss_green_TB%c13(n,:)*dmft_%weiss_green_TB%c31(:,n)
    
    ! Kadanoff-Baym self-energy <= contour-ordered self-energy
    do j=1, n
       dmft_%self_energyB%retarded(n,j)=dmft_%self_energy_TB%c21(n,j)-dmft_%self_energy_TB%c12(n,j)
       dmft_%self_energyB%lesser(n,j)=dmft_%self_energy_TB%c12(n,j)
    end do
    do i=1, n-1
       dmft_%self_energyB%lesser(i,n)=dmft_%self_energy_TB%c12(i,n)
    end do
    dmft_%self_energyB%left_mixing(n,:)=dmft_%self_energy_TB%c13(n,:)
    
    ! solve the Impurity Dyson equation: G = G0+G0*self_energy*G
    ! K=G0*self_energy
    call allocate_KB(parm_,kernel)
    call convolute_KB(parm_,n,dmft_%weiss_greenB,dmft_%self_energyB,kernel)
    
    ! G=[1-K]^{-1}*G0
    call Volterra_int(parm_,n,dmft_%weiss_greenB,kernel,dmft_%local_greenB)
    call deallocate_KB(kernel)
    
  end subroutine noneq_iptB

