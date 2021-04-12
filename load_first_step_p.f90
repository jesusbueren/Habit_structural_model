subroutine  load_first_step_p()
    use var_first_step
    implicit none
    
    !Load health transition pbs of the retirees
    call load_transitions_retirees()
    
    !Load medical expenditures of the retirees.
    call load_medical_retirees()
    
end subroutine
    
    
subroutine load_transitions_retirees()
    use global_var; use nrtype;use var_first_step
    implicit none
        integer,parameter::iterations=58640
        real(dp),dimension(clusters**2*covariates,iterations)::c_tr_all
        real(dp),dimension(covariates_habits*habits*types,iterations)::c_habits_all
        real(dp),dimension(clusters**2*covariates,1)::c_tr
        real(dp),dimension(covariates,clusters,clusters+1)::beta
        real(DP),dimension(clusters,types,L_gender,L_educ)::init_cond
        real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H
        real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE

        open(unit=9,file=path_s//'c_tr.txt')
                read(9,'(F20.8)') c_tr_all
        close(9)
    
        c_tr(:,1)=sum(c_tr_all,2)/dble(iterations)
    
        beta=0.0d0
        beta(:,:,1:clusters)=reshape(c_tr,(/covariates,clusters,clusters/))
    
        init_cond=0.0d0
        init_cond(1,:,:,:)=1.0d0
    
        call transitions(beta,init_cond,H,LE)
    
        H_fs(:,:,T_50:T,:,:)=H(:,:,:,:,1,:)
        

end subroutine
    
subroutine load_medical_retirees()
    use nrtype;use var_first_step;use state_space_dim
    implicit none
    integer,parameter::K=8
    real(DP),dimension(K+2)::data_csv
    real(DP),dimension(K,1)::X,beta
    real(DP)::rho,sigma2_u,age,pi_income,h,z
    real(DP),dimension(G_nzz,1)::mag_p,pr0_p
    integer::t_l,pi_l,h_l,z_l
    
    open(unit=9,file=path//'data\medical_exp.csv')
            read(9,*) data_csv
    close(9)
    
    beta(:,1)=data_csv(1:K)
    rho=data_csv(K+1)
    sigma2_u=data_csv(K+2)**2
    
    call discretize_shocks(rho,sigma2_u,G_nzz,mag_p,pr0_p,Pi_m)
    
    do t_l=T_R,T;do pi_l=1,G_PI;do h_l=1,G_h;do z_l=1,G_nzz
        age=dble(first_age+(t_l-1)*2)
        pi_income=PI_grid(pi_l)
        h=dble(h_l-1)
        z=mag_p(z_l,1)
        X(:,1)=(/pi_income, pi_income**2.0d0, age, age**2.0d0,h,h*pi_income,pi_income*age,1.0d0/)
        m_grid(pi_l,t_l,h_l,z_l)=exp(sum(X(:,1)*beta(:,1))+z)
    end do;end do;end do;end do
        

end subroutine
