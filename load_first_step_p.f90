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