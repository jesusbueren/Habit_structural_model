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
    
subroutine p2shock(gumbel_y,gumbel_hs,gumbel_cg,rand_ey)
use nrtype;use state_space_dim
implicit none
real(DP),dimension(2),intent(in)::gumbel_y,gumbel_hs,gumbel_cg
real(DP),dimension(G_df,G_educ,G_types),intent(out)::rand_ey
real(DP)::eps_pro,eps_hs,eps_cg


integer::i_l,j_l,e_l,y_l
interface
    double precision function c4_normal_01( )
                implicit none
    end function c4_normal_01
end interface

    

    call gumbel_shock(gumbel_hs,eps_hs)
    call gumbel_shock(gumbel_cg,eps_cg)
    call gumbel_shock(gumbel_y,eps_pro)
    
    rand_ey=0.0d0
    do e_l=1,G_educ; do y_l=1,G_types
        if (e_l==2) then
            rand_ey(1,e_l,y_l)=rand_ey(1,e_l,y_l)+eps_hs
        elseif (e_l==3) then
            rand_ey(1,e_l,y_l)=rand_ey(1,e_l,y_l)+eps_cg
        end if
        if (y_l==1) then
            rand_ey(1,e_l,y_l)=rand_ey(1,e_l,y_l)+eps_pro
        end if
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