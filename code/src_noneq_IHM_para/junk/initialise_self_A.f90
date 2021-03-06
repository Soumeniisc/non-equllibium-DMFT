subroutine initialize_self_energy_iptA(parm_,dmft_) !TODO carefully chek the intialization for A sublattice
    type(parm), intent(in)     :: parm_
    type(dmft), intent(inout)  :: dmft_
    integer                    :: j

    ! initialize the self-energy for real-time evolution
    ! self_energy_{12}(0,0) = i^3*U(0+)*U(0+)*G0(0-)*G0(0+)*G0(0-)
    dmft_%self_energy_TA%c12(1,1)=U(parm_,0.d0)*U(parm_,0.d0)*(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau+1) &
                                 *xj*dmft_%weiss_greenA%matsubara_t(1)*(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau+1)
    ! self_energy_{21}(0,0) = i^3*U(0+)*U(0+)*G0(0+)*G0(0-)*G0(0+)
    dmft_%self_energy_TA%c21(1,1)=U(parm_,0.d0)*U(parm_,0.d0)*xj*dmft_%weiss_greenA%matsubara_t(1) &
                                 *(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau+1)*xj*dmft_%weiss_greenA%matsubara_t(1)
    do j=1, parm_%N_tau+1
       ! self_energy_{13}(0,t) = i^3*U(0+)*U(0-)*G0(-t)*G0(t)*G0(-t)
       dmft_%self_energy_TA%c13(1,j)=U(parm_,0.d0)*parm_%U_i*(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-j+2) &
                                    *xj*dmft_%weiss_greenA%matsubara_t(j)*(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-j+2)
       ! self_energy_{31}(t,0) = i^3*U(0-)*U(0+)*G0(t)*G0(-t)*G0(t)
       dmft_%self_energy_TA%c31(j,1)=parm_%U_i*U(parm_,0.d0)*xj*dmft_%weiss_greenA%matsubara_t(j) &
                                    *(-xj)*dmft_%weiss_greenA%matsubara_t(parm_%N_tau-j+2)*xj*dmft_%weiss_greenA%matsubara_t(j)
    end do
    dmft_%self_energyA%retarded(1,1)=dmft_%self_energy_TA%c21(1,1)-dmft_%self_energy_TA%c12(1,1)
    dmft_%self_energyA%lesser(1,1)=dmft_%self_energy_TA%c12(1,1)
    dmft_%self_energyA%left_mixing(1,:)=dmft_%self_energy_TA%c13(1,:)
    
  end subroutine initialize_self_energy_iptA
