subroutine p2shock(sigma_y,sigma_hs,sigma_cg,rho_ey,rand_ey)
use nrtype;use state_space_dim
implicit none
real(DP),intent(in)::sigma_hs,sigma_cg,sigma_y
real(DP),intent(in)::rho_ey
real(DP),dimension(G_df,G_educ,G_types),intent(out)::rand_ey
real(DP)::eps_pro,eps_hs,eps_col
interface
    double precision function c4_normal_01( )
                implicit none
    end function c4_normal_01
end interface

    rand_ey=0.0d0
    !Random normal
    eps_pro=c4_normal_01( )
    rand_ey(:,:,1)=eps_pro*sqrt(sigma_y) 
    !conditional distribution given eps_pro
    eps_hs=c4_normal_01( )
    rand_ey(:,2,:)=rand_ey(:,2,:)+eps_hs*sqrt(sigma_hs)
    eps_col=c4_normal_01( )
    rand_ey(:,3,:)=rand_ey(:,3,:)+eps_col*sqrt(sigma_cg)


end subroutine