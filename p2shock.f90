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
    
subroutine p2shock(sigma_y,sigma_hs,sigma_cg,rho_ey,rand_ey)
use nrtype;use state_space_dim
implicit none
real(DP),intent(in)::sigma_hs,sigma_cg,sigma_y
real(DP),dimension(2),intent(in)::rho_ey
real(DP),dimension(G_df,G_educ,G_types),intent(out)::rand_ey
real(DP)::eps_pro,eps_hs,eps_col
real(DP),dimension(3,3)::var
real(DP),dimension(3)::L
real(DP),dimension(3,1)::u
integer::i_l,j_l
interface
    double precision function c4_normal_01( )
                implicit none
    end function c4_normal_01
end interface

    var=0.0d0
    var(1,1)=sigma_y
    var(2,2)=sigma_hs
    var(3,3)=sigma_cg
    var(1,2)=rho_ey(1)*sqrt(sigma_hs)*sqrt(sigma_y)
    var(1,3)=rho_ey(2)*sqrt(sigma_cg)*sqrt(sigma_y)
    var(2,1)=var(1,2)
    var(3,1)=var(1,3)
    
    call choldc(var,L)
    do i_l=1,3; do j_l=1,3
        if (i_l==j_l) then
            var(i_l,j_l)=L(i_l)
        elseif (j_l>i_l) then
            var(i_l,j_l)=0.0d0
        end if
    end do; end do
    
    do i_l=1,3
        u(i_l,1)=c4_normal_01( )
    end do
    
    u=matmul(var,u)
    
    rand_ey=0.0d0
    rand_ey(:,:,1)=u(1,1)
    rand_ey(:,2,:)=rand_ey(:,2,:)+u(2,1)
    rand_ey(:,3,:)=rand_ey(:,3,:)+u(3,1)
    
end subroutine
    
SUBROUTINE choldc(a,p)
USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
IMPLICIT NONE
REAL(DP), DIMENSION(3,3), INTENT(INOUT) :: a
REAL(DP), DIMENSION(3), INTENT(OUT) :: p
INTEGER(I4B) :: i,n
REAL(SP) :: summ
n=assert_eq(size(a,1),size(a,2),size(p),'choldc')
do i=1,n
	summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
	if (summ <= 0.0) call nrerror('choldc failed')
	p(i)=sqrt(summ)
	a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
end do
END SUBROUTINE choldc