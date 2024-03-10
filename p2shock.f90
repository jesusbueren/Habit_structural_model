!subroutine p2shock(sigma_y,sigma_hs,sigma_cg,rho_ey,rand_ey)
!use nrtype;use state_space_dim
!implicit none
!real(DP),intent(in)::sigma_hs,sigma_cg,sigma_y
!real(DP),intent(in)::rho_ey
!real(DP),dimension(G_df,G_educ,G_types),intent(out)::rand_ey
!real(DP)::eps_pro,eps_hs,eps_col
!interface
!    double precision function c4_normal_01( )
!                implicit none
!    end function c4_normal_01
!end interface
!
!    rand_ey=0.0d0
!    !Random normal
!    eps_pro=c4_normal_01( )
!    print*,eps_pro*sqrt(sigma_y) 
!    rand_ey(:,:,1)=eps_pro*sqrt(sigma_y) 
!    !conditional distribution given eps_pro
!    eps_hs=c4_normal_01( )
!    print*,eps_hs*sqrt(sigma_hs)
!    rand_ey(:,2,:)=rand_ey(:,2,:)+eps_hs*sqrt(sigma_hs)
!    eps_col=c4_normal_01( )
!    print*,eps_col*sqrt(sigma_cg)
!    rand_ey(:,3,:)=rand_ey(:,3,:)+eps_col*sqrt(sigma_cg)
!
!
!end subroutine
    
subroutine p2shock(mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey,rand_ey)
use nrtype;use state_space_dim
implicit none
real(DP),intent(in)::mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey
real(DP),dimension(G_df,G_educ,G_types),intent(out)::rand_ey
real(DP)::eps_y,eps_e
integer::i_l,j_l,e_l,y_l
interface
    double precision function c4_normal_01( )
                implicit none
    end function c4_normal_01
end interface

    !call RANDOM_NUMBER(eps_e)
    !eps_e=(eps_e-0.5d0)*2.0d0*sqrt(var_e)+mu_e
    !call RANDOM_NUMBER(eps_y)
    !eps_y=(eps_y-0.5d0)*2.0d0*sqrt(var_y)+mu_y
    eps_e=c4_normal_01( )
    eps_y=rho_ey*eps_e+sqrt(1.0d0-rho_ey**2.0d0)*c4_normal_01( )
    eps_e=eps_e*sqrt(var_e)+mu_e
    !call sample_cauchy(var_y,eps_y)
    !call gumbel_shock((/0.0d0,var_y/),eps_y)
    !call sample_gamma(1.3d0, var_y, eps_y)
    !eps_y=eps_y+mu_y
    eps_y=eps_y*sqrt(var_y)+mu_y !exp(eps_y*sqrt(var_y))+ mu_y


    
    rand_ey=0.0d0
    do e_l=1,G_educ; do y_l=1,G_types
        if (e_l==2) then
            rand_ey(1,e_l,y_l)=rand_ey(1,e_l,y_l)+eps_e
        elseif (e_l==3) then
            rand_ey(1,e_l,y_l)=rand_ey(1,e_l,y_l)+eps_e*tau_cg
        end if
        if (y_l==1) then
            rand_ey(1,e_l,y_l)=rand_ey(1,e_l,y_l)+eps_y
        end if
        rand_ey(1,e_l,y_l)=rand_ey(1,e_l,y_l)!+c4_normal_01( )*sqrt(var_nu)
    end do; end do

end subroutine
    
SUBROUTINE gumbel_shock(gumbel,u_sample)
use nrtype
implicit none
real(DP),dimension(2),intent(in):: gumbel
real(DP),intent(out):: u_sample
real(DP):: u

call RANDOM_NUMBER(u)

u_sample = gumbel(1) - gumbel(2) * log(-log(u))


end subroutine
    
subroutine sample_exponential(lambda, sample)
use nrtype
implicit none
  real(DP), intent(in) :: lambda
  real(DP), intent(out) :: sample

  ! Generate a random number from a uniform distribution [0, 1)
  real(DP) :: u
  
  call random_number(u)

  ! Inverse transform sampling for exponential distribution
  sample = -log(u) / lambda

    end subroutine sample_exponential    
    
subroutine sample_cauchy(scale, sample)
use nrtype
    implicit none
    real(DP), intent(in) :: scale
    real(DP), intent(out) :: sample
    real(DP) :: u
    integer :: i

    call random_number(u)
    sample = scale * tan((u - 0.50d0) * atan(1.0d0) * 4.0d0)

    end subroutine sample_cauchy
    
    
subroutine sample_gamma(shape, scale, sample)
use nrtype
    implicit none
    real(DP), intent(in) :: shape, scale
    real(DP), intent(out) :: sample
    real(DP) :: u1, u2, gamma_value, delta, temp

    call random_number(u1)
    call random_number(u2)

    ! Inverse transform sampling for gamma distribution
    gamma_value = shape * scale
    delta = sqrt(2.0d0 * shape - 1.0d0)
    temp = delta * sqrt(-2.0d0 * log(u1))
    sample = gamma_value + temp * cos(2.0d0 * atan(1.0d0) * u2)
end subroutine sample_gamma

    