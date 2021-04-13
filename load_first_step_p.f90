subroutine  load_first_step_p()
    use var_first_step
    implicit none
    
    !Load health transition pbs of the retirees
    call load_transitions_retirees()
    
    !Load medical expenditures of the retirees.
    call load_medical_retirees()
    
end subroutine
    
    
subroutine load_transitions_retirees()
    use global_var; use nrtype;use var_first_step; use state_space_dim
    implicit none
        integer,parameter::iterations=58640
        real(dp),dimension(clusters**2*covariates,iterations)::c_tr_all
        real(dp),dimension(covariates_habits*habits*types,iterations)::c_habits_all
        real(dp),dimension(clusters**2*covariates,1)::c_tr
        real(dp),dimension(covariates,clusters,clusters+1)::beta
        real(DP),dimension(clusters+1,clusters+1,T,types,L_gender,L_educ)::H
        integer::t_l,h_l
        real(DP),dimension(G_h+1,G_types,G_educ)::pr_sr

        open(unit=9,file=path_s//'c_tr.txt')
                read(9,'(F20.8)') c_tr_all
        close(9)
    
        c_tr(:,1)=sum(c_tr_all,2)/dble(iterations)
    
        beta=0.0d0
        beta(:,:,1:clusters)=reshape(c_tr,(/covariates,clusters,clusters/))
    
        call extrapolated_transitions(beta,H)
    
        H_sm(:,:,1:T,:,:)=H(:,:,1:T,:,1,:)
        !For individuals younger than fifty I assume that they don't die
        do t_l=1,T_50-1
            H_sm(1:G_h,G_h+1,T_50,:,:)=0.0d0
            pr_sr=sum(H_sm(:,:,T_50,:,:),2)
            do h_l=1,G_h
                H_sm(:,h_l,t_l,:,:)=H_sm(:,h_l,T_50,:,:)/pr_sr
            end do
        end do
!H_sm(1,1,:,2,1)
end subroutine
    
subroutine load_medical_retirees()
    use nrtype;use var_first_step;use state_space_dim
    implicit none
    integer,parameter::K=13
    real(DP),dimension(K+2)::data_csv
    real(DP),dimension(K,1)::X,beta
    real(DP)::rho,sigma2_u,age,h,z,d1,d2,d3
    real(DP),dimension(G_nzz,1)::mag_p,pr0_p
    integer::t_l,e_l,h_l,z_l
    
    open(unit=9,file=path//'data\medical_exp.csv')
            read(9,*) data_csv
    close(9)
    
    beta(:,1)=data_csv(1:K)
    rho=data_csv(K+1)
    sigma2_u=data_csv(K+2)**2
    
    call discretize_shocks(rho,sigma2_u,G_nzz,mag_p,pr0_p,Pi_m)
    
    do t_l=1,T;do e_l=1,G_educ;do h_l=1,G_h;do z_l=1,G_nzz
        age=dble(first_age_sm+(t_l-1)*2)
        d1=0.0d0
        d2=0.0d0
        d3=0.0d0
        if (e_l==1) then
            d1=1.0d0
        elseif (e_l==2) then
            d2=1.0d0
        elseif (e_l==3) then
            d3=1.0d0
        end if
        h=dble(h_l-1)
        z=mag_p(z_l,1)
        X(:,1)=(/d1,d2,d3,age,age**2.0d0,h,d1*h,d2*h,d3*h,d1*age,d2*age,d3*age,1.0d0/)
        m_grid(e_l,t_l,h_l,z_l)=exp(sum(X(:,1)*beta(:,1))+z)
    end do;end do;end do;end do
       

end subroutine
