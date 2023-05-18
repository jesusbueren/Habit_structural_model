subroutine  discretize_iid_shocks(s2_nu_av,nzz,mag_p,pr0_p)
    implicit none
    integer,intent(in)::nzz
    double precision,intent(in)::s2_nu_av
    double precision,dimension(nzz,1),intent(out)::mag_p
    double precision,dimension(nzz,1),intent(out)::pr0_p
    double precision,dimension(nzz)::x
    double precision,dimension(nzz-1,1)::dz
    integer::i_l,j_l
    double precision::x1,x2,pr_1,pr_2
    character::error
    
    x(1)=-2.5d0 !-2.33d0
    x(nzz)=2.5d0 !=2.33d0
    
    do i_l=2,nzz-1
        x(i_l)=x(i_l-1)+(x(nzz)-x(1))/dble(nzz-1)
    end do

    
    !Persistent shock   
    mag_p=-9.0d0
    pr0_p=-9.0d0
    
    do i_l=1,nzz
        mag_p(i_l,1)=sqrt(2.0d0)*sqrt(s2_nu_av)*x(i_l)
        if (i_l>1) then
            dz(i_l-1,1)=mag_p(i_l,1)-mag_p(i_l-1,1)
        end if
    end do
    
    i_l=1
    x1=(mag_p(i_l,1)+dz(i_l,1)/2.0d0)/sqrt(s2_nu_av)
    pr0_p(i_l,1)=0.5d0+0.5d0*erf(x1/sqrt(2.0d0))
    
    i_l=nzz
    x1=(mag_p(i_l,1)-dz(i_l-1,1)/2.0d0)/sqrt(s2_nu_av)
    pr0_p(i_l,1)=1.0d0-(0.5d0+0.5d0*erf(x1/sqrt(2.0d0)))
    
    do i_l=2,nzz-1
        x1=(mag_p(i_l,1)-dz(i_l-1,1)/2d0)/sqrt(s2_nu_av)
        pr_1=0.5d0+0.5d0*erf(x1/sqrt(2.0d0))
        x2=(mag_p(i_l,1)+dz(i_l,1)/2d0)/sqrt(s2_nu_av)
        pr_2=0.5d0+0.5d0*erf(x2/sqrt(2.0d0))
        pr0_p(i_l,1)=pr_2-pr_1
    end do   

    end subroutine
    
subroutine  transition_pr(rho_av,s2_nu_av,nzz,mag_1,mag_2,pr_p_new)
    implicit none
    integer,intent(in)::nzz
    double precision,intent(in)::rho_av,s2_nu_av
    double precision,dimension(nzz,1),intent(in)::mag_1,mag_2
    double precision,dimension(nzz,nzz),intent(out)::pr_p_new
    double precision,dimension(nzz)::x
    double precision,dimension(nzz-1,1)::dz
    integer::i_l,j_l
    double precision::x1,x2,pr_1,pr_2
    character::error
    
    do i_l=1,nzz
        if (i_l>1) then
            dz(i_l-1,1)=mag_2(i_l,1)-mag_2(i_l-1,1)
        end if
    end do
    
    
    do i_l=1,nzz
        do j_l=2,nzz-1
            x1=(mag_2(j_l,1)-dz(j_l-1,1)/2.0d0-rho_av*mag_1(i_l,1))/sqrt(s2_nu_av)
            pr_1=0.5d0+0.5d0*erf(x1/sqrt(2.0d0))
            x2=(mag_2(j_l,1)+dz(j_l,1)/2.0d0-rho_av*mag_1(i_l,1))/sqrt(s2_nu_av)
            pr_2=0.5d0+0.5d0*erf(x2/sqrt(2.0d0))
            pr_p_new(i_l,j_l)=pr_2-pr_1
        end do
    
        j_l=1
        x1=(mag_2(j_l,1)+dz(j_l,1)/2.0d0-rho_av*mag_1(i_l,1))/sqrt(s2_nu_av)
        pr_p_new(i_l,j_l)=0.5d0+0.5d0*erf(x1/sqrt(2.0d0)) 
        
        j_l=nzz
        x1=(mag_2(j_l,1)-dz(j_l-1,1)/2.0d0-rho_av*mag_1(i_l,1))/sqrt(s2_nu_av)
        pr_p_new(i_l,j_l)=1.0d0-(0.5d0+0.5d0*erf(x1/sqrt(2.0d0))) 
        
        pr_p_new(i_l,:)=pr_p_new(i_l,:)/sum(pr_p_new(i_l,:))

    end do

end subroutine
    
subroutine  discretize_shocks(rho_av,s2_nu_av,nzz,mag_p,pr0_p,pr_p_new)
    implicit none
    integer,intent(in)::nzz
    double precision,intent(in)::rho_av,s2_nu_av
    double precision,dimension(nzz,1),intent(out)::mag_p
    double precision,dimension(nzz,1),intent(out)::pr0_p
    double precision,dimension(nzz,nzz),intent(out)::pr_p_new
    double precision,dimension(nzz)::x
    double precision,dimension(nzz-1,1)::dz
    integer::i_l,j_l
    double precision::x1,x2,pr_1,pr_2
    character::error
    
    x(1)=-3.0d0
    x(nzz)=3.0d0
    
    do i_l=2,nzz-1
        x(i_l)=x(i_l-1)+(x(nzz)-x(1))/dble(nzz-1)
    end do

    
    !Persistent shock   
    mag_p=-9.0d0
    pr0_p=-9.0d0
    pr_p_new=0.0d0
    
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
            pr_p_new(i_l,j_l)=pr_2-pr_1
        end do
    
        j_l=1
        x1=(mag_p(j_l,1)+dz(j_l,1)/2.0d0-rho_av*mag_p(i_l,1))/sqrt(s2_nu_av)
        pr_p_new(i_l,j_l)=0.5d0+0.5d0*erf(x1/sqrt(2.0d0)) 
        
        j_l=nzz
        x1=(mag_p(j_l,1)-dz(j_l-1,1)/2.0d0-rho_av*mag_p(i_l,1))/sqrt(s2_nu_av)
        pr_p_new(i_l,j_l)=1.0d0-(0.5d0+0.5d0*erf(x1/sqrt(2.0d0))) 
        
        pr_p_new(i_l,:)=pr_p_new(i_l,:)/sum(pr_p_new(i_l,:))!pr_p_new(3,:)

    end do

end subroutine
        
    