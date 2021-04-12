subroutine  discretize_shocks(rho_av,s2_nu_av,nzz,mag_p,pr0_p,pr_p)
    use nrtype
    implicit none
    integer,intent(in)::nzz
    real(DP),intent(in)::rho_av,s2_nu_av
    real(DP),dimension(nzz,1),intent(out)::mag_p,pr0_p
    real(DP),dimension(nzz,nzz),intent(out)::pr_p
    real(DP),dimension(nzz)::x
    real(DP),dimension(nzz-1,1)::dz
    integer::i_l,j_l
    real(DP)::x1,x2,pr_1,pr_2
    character::error
    
    if (nzz==2) then
        x=(/-0.7071067811,0.7071067811/)
    elseif (nzz==3) then
        x=(/-1.22474487139158,0.0,1.22474487139158/)
    elseif (nzz==4) then
        x=(/-1.650680123,-0.5246476232,0.5246476232,1.650680123/)
    elseif (nzz==5) then
        x=(/-2.02018287,-0.9585724646,0.0,0.9585724646,2.02018287/) 
    else
        print*,'error: introduce more points'
        read*,error
    end if
    
    !Persistent shock
    do i_l=1,nzz
        mag_p(i_l,1)=sqrt(2.0d0)*sqrt(s2_nu_av/(1.0d0-rho_av**2))*x(i_l)
        if (i_l>1) then
            dz(i_l-1,1)=mag_p(i_l,1)-mag_p(i_l-1,1)
        end if
    end do
    
    i_l=1
    x1=(mag_p(i_l,1)+dz(i_l,1)/2.0d0)/sqrt(s2_nu_av/(1.0d0-rho_av**2))
    pr0_p(i_l,1)=0.5d0+0.5d0*erf(x1/sqrt(2.0d0))
    
    i_l=nzz
    x1=(mag_p(i_l,1)-dz(i_l-1,1)/2.0d0)/sqrt(s2_nu_av/(1.0d0-rho_av**2))
    pr0_p(i_l,1)=1.0d0-(0.5d0+0.5d0*erf(x1/sqrt(2.0d0)))
    
    do i_l=2,nzz-1
        x1=(mag_p(i_l,1)-dz(i_l-1,1)/2d0)/sqrt(s2_nu_av/(1.0d0-rho_av**2))
        pr_1=0.5d0+0.5d0*erf(x1/sqrt(2.0d0))
        x2=(mag_p(i_l,1)+dz(i_l,1)/2d0)/sqrt(s2_nu_av/(1.0d0-rho_av**2))
        pr_2=0.5d0+0.5d0*erf(x2/sqrt(2.0d0))
        pr0_p(i_l,1)=pr_2-pr_1
    end do   
    
    do i_l=1,nzz
        do j_l=2,nzz-1
            x1=(mag_p(j_l,1)-dz(j_l-1,1)/2.0d0-rho_av*mag_p(i_l,1))/sqrt(s2_nu_av)
            pr_1=0.5d0+0.5d0*erf(x1/sqrt(2.0d0))
            x2=(mag_p(j_l,1)+dz(j_l,1)/2.0d0-rho_av*mag_p(i_l,1))/sqrt(s2_nu_av)
            pr_2=0.5d0+0.5d0*erf(x2/sqrt(2.0d0))
            pr_p(i_l,j_l)=pr_2-pr_1
        end do

        j_l=1
        x1=(mag_p(j_l,1)+dz(j_l,1)/2.0d0-rho_av*mag_p(i_l,1))/sqrt(s2_nu_av)
        pr_p(i_l,j_l)=0.5d0+0.5d0*erf(x1/sqrt(2.0d0)) 
        
        j_l=nzz
        x1=(mag_p(j_l,1)-dz(j_l-1,1)/2.0d0-rho_av*mag_p(i_l,1))/sqrt(s2_nu_av)
        pr_p(i_l,j_l)=1.0d0-(0.5d0+0.5d0*erf(x1/sqrt(2.0d0))) 
        
        pr_p(i_l,:)=pr_p(i_l,:)/sum(pr_p(i_l,:))
        print*,sum(pr_p(i_l,:))
        print*,''
    end do

end subroutine