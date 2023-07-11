subroutine  load_first_step_p()
    implicit none
    
    !Create grid of capital
    call create_grids()
    
    !Load health transition pbs of the retirees
    call load_transitions_retirees()
    
    !Load medical expenditures of the retirees.
    call load_medical_expenses()
    
    !Load income risk of workers
    call load_income_risk()

end subroutine
    
subroutine create_grids()
    use nrtype; use var_first_step
    implicit none
    integer::i_l

    do i_l=1,G_nkk
        if (i_l==1) then
            a_grid(i_l)=a_min
        else
            a_grid(i_l)=a_grid(i_l-1)+step
        end if
    end do

end subroutine
    
    
subroutine load_transitions_retirees()
    use global_var; use nrtype;use var_first_step; use state_space_dim
    implicit none
        integer,parameter::iterations=2728
        real(dp),dimension(covariates,clusters,L_gender,L_educ,iterations)::c_tr_all
        real(dp),dimension(covariates,clusters,L_gender,L_educ,iterations)::c_tr_d_all
        real(dp),dimension(covariates_habits*habits*types,iterations)::c_habits_all

        real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_h
        real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_d
        real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H
        integer::t_l,h_l,e_l
        real(DP),dimension(G_h+1,G_types,G_educ)::pr_sr
        real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE
        real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts)::joint_yh
        real(DP),dimension(generations,L_gender,L_educ,types,cohorts,iterations)::fraction_t
        real(DP),dimension(generations,G_h,L_gender,L_educ,types,iterations)::fraction_h

        open(unit=9,file=path_s//'c_tr.txt')
                read(9,'(F20.8)') c_tr_all
        close(9)
        open(unit=9,file=path_s//'c_tr_d.txt')
                read(9,'(F20.8)') c_tr_d_all
        close(9)
        
        open(unit=9,file=path_s//'fraction_t.txt')
                read(9,'(F7.4)') fraction_t
        close(9)
        
        open(unit=9,file=path_s//'fraction_h.txt')
                read(9,'(F7.4)') fraction_h
        close(9)
        
        fraction_types=sum(fraction_t(1,1,:,:,reference_cohort,:),3)/dble(iterations) 
        do e_l=1,G_educ
            fraction_types(e_l,:)=fraction_types(e_l,:)/sum(fraction_types(e_l,:)) 
        end do
        fraction_h_ey=sum(fraction_h(1,:,1,:,:,:),4)/dble(iterations) 
        
        if (reference_cohort==3) then
            fraction_e=(/0.089d0,0.514d0,0.397d0/)
        elseif (reference_cohort==4) then
            fraction_e=(/0.085d0,0.550d0,0.365d0/)
        end if
    
    
        beta_h=sum(c_tr_all,5)/dble(iterations)
        beta_d=sum(c_tr_d_all,5)/dble(iterations)
        
        joint_yh=1.0d0
        call transitions(beta_h,beta_d,H,LE,joint_yh)
        
        H_sm=H(:,:,2:generations,:,1,:) 
        
        !Dead for sure the last period in the model H_sm(1,3,:,3,3)
        H_sm(1:clusters,clusters+1,T-1:T,:,:)=1.0d0
        H_sm(1:clusters,1:clusters,T-1:T,:,:)=0.0d0
        
        !Modify vfct for bequest
        !print*,'trick to change vfct'
        !H_sm(:,clusters+1,:,:,:)=1.0d0



end subroutine
    
subroutine load_medical_expenses()
    use nrtype;use var_first_step;use state_space_dim; use global_var
    implicit none
    integer,parameter::K=14
    real(DP),dimension(14*2)::data_csv
    real(DP),dimension(K,1)::X,beta,beta_var
    real(DP)::age,h,z,d1,d2,d3
    real(DP),dimension(G_nzz,1)::mag_p
    real(DP),dimension(K)::beta_insurance_u65,beta_insurance_o65
    real(DP),dimension(G_educ,T,G_h)::fraction_covered!Medical grid
    integer::t_l,e_l,h_l,z_l
    real(DP),dimension(2)::h_d
    real(DP),dimension(3)::e_d
    real(DP),dimension(G_educ,T,G_h)::sigma2_u
    character::pause_k
    
    open(unit=9,file=path//'data\medical_exp.csv')
            read(9,*) data_csv
    close(9)

    
    beta(:,1)=data_csv(1:K)
    beta_var(:,1)=data_csv(K+1:2*K)
    

    !Mean expenses from stata regression
    m_grid=0.0d0
    do e_l=1,G_educ;do h_l=1,G_h;do t_l=T_50,T
        !Set covariates
        age=dble(first_age_sm+(t_l-1)*2)
        e_d=0.0d0
        e_d(e_l)=1.0d0
        h_d=0.0d0
        h_d(h_l)=1.0d0
        X(:,1)=(/age,age**2.0d0,age**3.0d0,h_d,h_d*age,e_d*age,e_d,1.0d0/)
        !set variance
        sigma2_u(e_l,t_l,h_l)=sum(X(:,1)*beta_var(:,1))
        call discretize_iid_shocks(sigma2_u(e_l,t_l,h_l),G_nzz,mag_p,pr0_p)
        do z_l=1,G_nzz
            z=mag_p(z_l,1)
            m_grid(e_l,t_l,h_l,z_l)=exp(sum(X(:,1)*beta(:,1))+z)
        end do
    end do;end do;end do

end subroutine
    
    
subroutine load_income_risk()
    use nrtype;use var_first_step;use state_space_dim; use global_var
    implicit none
    integer,parameter::K=6
    real(DP),dimension(K,1)::X
    real(DP),dimension(K,G_educ)::beta
    real(DP),dimension(G_educ,G_cohorts):: s2_nu,s2_w,rho,s2_0
    real(DP)::sigma2_u,age,h,z,d1,d2,d3,z2,var_aux
    real(DP),dimension(G_PI,generations,G_educ)::mag_p
    real(DP),dimension(G_PI,G_educ)::mag_tr
    real(DP),dimension(G_PI,G_PI)::Pi_nothing
    real(DP),dimension(G_types)::y_d
    real(DP),dimension(G_cohorts)::cohort_d
    integer::t_l,e_l,h_l,z_l,z_l2,y_l
    character::pause_k
    real(DP),dimension(5,generations,clusters,L_educ,types)::labor_force_participation
    real(DP),dimension(5,generations,clusters,L_educ,types,2)::labor_force_participation_dynamic
    real(DP)::gross_annual_income
    
    !Labor force participation risk at age 25
    open(unit=9,file=path//'metric_model\Results\labor_force_participation.txt')
            read(9,*) labor_force_participation 
    close(9)
    print*,'cond on health, edu and age, labor force participation does not vary by type'
    do t_l=1,T;do h_l=1,G_h;do e_l=1,G_educ
        Pi_p_0(G_PI,t_l,h_l,e_l)=1.0d0-labor_force_participation(5,t_l,h_l,e_l,1)
    end do;end do;end do
    
    !Labor force participation conditional previous labor force participation
    open(unit=9,file=path//'metric_model\Results\labor_force_participation_dynamic.txt')
            read(9,*) labor_force_participation_dynamic
    close(9)
    print*,'cond on health, edu and age, previous labor force participation does not vary by type'
    Pi_p=-9.0d0
    do t_l=1,T;do h_l=1,G_h;do e_l=1,G_educ
        Pi_p(1:G_PI-1,G_PI,t_l,h_l,e_l)=1.0d0-labor_force_participation_dynamic(5,t_l,h_l,e_l,1,2) ! transition in the labor force -> out of the labor force 
        Pi_p(G_PI,G_PI,t_l,h_l,e_l)=1.0d0-labor_force_participation_dynamic(5,t_l,h_l,e_l,1,1) ! transition out of the labor force -> out of the labor force 
    end do;end do;end do

    !Income risk for those being employed
    open(unit=9,file=path//'metric_model\Results\parameters_income.txt')
            read(9,*) beta,s2_nu,s2_w,rho,s2_0
    close(9)
    
   !Compute unconditional probability and magnitude of the transitory shock
    do e_l=1,G_educ
        !Discretize the transitory shock
        call discretize_iid_shocks(s2_w(e_l,1),G_PI,mag_tr(:,e_l),Pi_t(:,e_l))  
    end do
    
    !Compute unconditional probability and magnitude of the persistent shock
    do e_l=1,G_educ;do h_l=1,G_h
        var_aux=s2_0(e_l,1)
        t_l=1
        call discretize_iid_shocks(var_aux,G_PI-1,mag_p(1:G_PI-1,t_l,e_l),Pi_p_0(1:G_PI-1,t_l,h_l,e_l))  !mag_p(1:9,t_l,e_l)
        !correct for the pr of being unemployed
        Pi_p_0(1:G_PI-1,t_l,h_l,e_l)=Pi_p_0(1:G_PI-1,t_l,h_l,e_l)/sum(Pi_p_0(1:G_PI-1,t_l,h_l,e_l))*(1.0d0-Pi_p_0(G_PI,t_l,h_l,e_l)) 
        if (isnan(sum(Pi_p_0))) then
            print*,'pb 1'
        end if 
        do t_l=2,T
            var_aux=rho(e_l,1)**2*var_aux+s2_nu(e_l,1)
1            call discretize_iid_shocks(var_aux,G_PI-1,mag_p(1:G_PI-1,t_l,e_l),Pi_p_0(1:G_PI-1,t_l,h_l,e_l))  
            !correct for the pr of being unemployed
            Pi_p_0(1:G_PI-1,t_l,h_l,e_l)=Pi_p_0(1:G_PI-1,t_l,h_l,e_l)/sum(Pi_p_0(1:G_PI-1,t_l,h_l,e_l))*(1.0d0-Pi_p_0(G_PI,t_l,h_l,e_l))
            if (sum(Pi_p_0(1:G_PI,t_l,h_l,e_l))>1.01d0) then
                print*,sum(Pi_p_0(1:G_PI,t_l,h_l,e_l))
                go to 1
            end if
            if (isnan(sum(Pi_p_0))) then
                print*,'pb 2'
            end if 
        end do
         
    end do;end do


    !Compute the transition probability of the persistent component cond on working
    do h_l=1,G_h;do e_l=1,G_educ;do t_l=1,T
        call transition_pr(rho(e_l,1),s2_nu(e_l,1),G_PI-1,mag_p(1:G_PI-1,t_l,e_l),mag_p(1:G_PI-1,t_l+1,e_l),Pi_p(1:G_PI-1,1:G_PI-1,t_l,h_l,e_l)) 
    end do;end do;end do

    
    !I need to adjust the pr of each persistent shock conditional on being in labor force so that from 1 to G_PI they sum to 1
    !I need to include the pr of moving to employment from unemployment
    do t_l=1,T-1;do h_l=1,G_h;do e_l=1,G_educ;do z_l=1,G_PI
        if (z_l<G_PI) then
            Pi_p(z_l,1:G_PI-1,t_l,h_l,e_l)=Pi_p(z_l,1:G_PI-1,t_l,h_l,e_l)/sum(Pi_p(z_l,1:G_PI-1,t_l,h_l,e_l))*(1.0d0-Pi_p(z_l,G_PI,t_l,h_l,e_l))
        else
            !the unemployment to employement: the productivity shock is equal to the uncond mean
            Pi_p(z_l,1:G_PI-1,t_l,h_l,e_l)=Pi_p_0(1:G_PI-1,t_l,h_l,e_l)/sum(Pi_p_0(1:G_PI-1,t_l,h_l,e_l))*(1.0d0-Pi_p(z_l,G_PI,t_l,h_l,e_l)) 
            !the unemployment to employement: the productivity shock is the lowest possible
            !Pi_p(z_l,1:G_PI-1,t_l,h_l,e_l)=0.0d0
            !Pi_p(z_l,1,t_l,h_l,e_l)=(1.0d0-Pi_p(z_l,G_PI,t_l,h_l,e_l)) 
        end if  
        if (isnan(sum(Pi_p(z_l,:,t_l,h_l,e_l)))) then
            print*,'pb 3',Pi_p(z_l,:,t_l,h_l,e_l)
            print*,z_l,t_l,h_l,e_l
            pause
        end if
        if (sum(Pi_p(z_l,:,t_l,h_l,e_l))<0.99d0) then
            print*,'pb 3',Pi_p(z_l,:,t_l,h_l,e_l)
            print*,z_l,t_l,h_l,e_l
            pause
        end if
    end do;end do;end do;end do
        

    income_grid=0.0d0
    do t_l=1,T_R;do y_l=1,G_types;do e_l=1,G_educ;do h_l=1,G_h;do z_l=1,G_PI;do z_l2=1,G_PI
        if (z_l<G_PI) then
            age=first_age_sm+(t_l-1)*2-70 
            z=mag_p(z_l,t_l,e_l)
            z2=mag_tr(z_l2,e_l)
            y_d=0.0d0
            y_d(y_l)=1.0d0
            !reference_cohort
            cohort_d(4)=1.0d0
            X(:,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(h_l-1),dble(h_l-1)*dble(age)/)
            gross_annual_income=exp(sum(X(:,1)*beta(:,e_l))+z+z2)/1000.0d0 !beta(:,3)
            gross_annual_income2(y_l,e_l,t_l,h_l,z_l,z_l2)=exp(sum(X(:,1)*beta(:,e_l))+z+z2)/1000.0d0
            taxable_income(y_l,e_l,t_l,h_l,z_l,z_l2)=min(gross_annual_income,ss_bar)*2.0d0
            income_grid(y_l,e_l,t_l,h_l,z_l,z_l2)=(av_income*lambda*(gross_annual_income/av_income)**(1.0d0-tau) - tau_mcr*gross_annual_income-tau_ss*min(gross_annual_income,ss_bar))*2.0d0
            income_tax(y_l,e_l,t_l,h_l,z_l,z_l2)=gross_annual_income*2.0d0-av_income*lambda*(gross_annual_income/av_income)**(1.0d0-tau)*2.0d0
        else
            income_grid(y_l,e_l,t_l,h_l,z_l,z_l2)=0.0d0
            gross_annual_income2(y_l,e_l,t_l,h_l,z_l,z_l2)=0.0d0
            taxable_income(y_l,e_l,t_l,h_l,z_l,z_l2)=0.0d0
        end if
        if (income_grid(y_l,e_l,t_l,h_l,z_l,z_l2)<0.0d0) then
            print*,'negative income'           
        end if
    end do;end do;end do;end do;end do;end do
    
    !print*,'chg this!!' gross_annual_income2(3,3,10,:,5,5)
    !income_grid(:,3,:,:,:,:)=income_grid(:,3,:,:,:,:)*1.5d0

    call compute_pension()
    
!income_grid(1,3,:,1,3,3)
end subroutine load_income_risk
    
subroutine compute_pension()
    use nrtype;use state_space_dim;use var_first_step
        implicit none
        integer,parameter::indv_sim_p=100000
        double precision::u
        integer::h,t_l,ind,i_l,y,e,h2,pi_l,df,ts_l,pi_l2,pi_old
        character::continue_k
        real(DP),dimension(T_R)::gross_earnings
        integer,dimension(G_PI)::counter
        real(DP),dimension(indv_sim_p,G_PI)::pension
        real(DP),dimension(indv_sim_p,T)::panel_income
        integer,dimension(T)::counter_income,counter_bh
        integer,dimension(T_R)::counter_un
        real(DP),dimension(T,G_educ,G_types)::av_income_panel,std_income
        real(DP),dimension(T_R,G_educ,G_types)::pr_unemployed,pr_bh
        real(DP)::aime,pia,coeff
        
        open(unit=9,file='wages.txt')
        
        pr_unemployed=0.0d0
        df=G_DF
        do e=1,G_educ; do y=1,G_types
            counter=0
            counter_income=0
            counter_un=0
            counter_bh=0
            panel_income=-9
            do i_l=1,indv_sim_p;gross_earnings=-9.0d0;aime=-9.0d0
                pi_old=-9
                do t_l=1,T_R
                    if (t_l==1) then
                        call RANDOM_NUMBER(u)
                        ind=1
                        h=-9
                        do while (h==-9)
                            if (u<sum(fraction_h_ey(1:ind,e,y)))then
                                h=ind
                            else
                                ind=ind+1
                            end if
                        end do
                        !Sample persisent income shock from uncond distribution
                        call RANDOM_NUMBER(u)
                        ind=1
                        pi_l=-9
                        do while (pi_l==-9)
                            if (u<sum(Pi_p_0(1:ind,t_l,h,e))) then !Pi_p_0(:,t_l,h,e)
                                pi_l=ind
                            else
                                ind=ind+1
                            end if
                        end do

                        !Sample transitory income shock from uncond distribution
                        call RANDOM_NUMBER(u)
                        ind=1
                        ts_l=-9
                        do while (ts_l==-9)
                            if (ind>G_PI) then
                                print*,'error simulating initial persistent shock 1'
                                read*,continue_k
                            end if
                            if (u<sum(Pi_t(1:ind,e))) then
                                ts_l=ind
                            else
                                ind=ind+1
                            end if
                        end do

                    else
                        !Sample persisent income shock from cond distribution
                        if (t_l<=T_R) then
                            call RANDOM_NUMBER(u)
                            ind=1
                            pi_old=pi_l
                            pi_l2=-9
                            do while (pi_l2==-9)
                                if (ind>G_PI) then
                                    print*,'error simulating initial persistent shock 2'
                                    read*,continue_k
                                end if
                                if (u<sum(Pi_p(pi_l,1:ind,t_l,h,e))) then 
                                    pi_l2=ind
                                    pi_l=pi_l2
                                else
                                    ind=ind+1
                                end if
                            end do
            
                            !Sample transitory income shock from uncond distribution
                            call RANDOM_NUMBER(u)
                            ind=1
                            ts_l=-9
                            do while (ts_l==-9)
                                if (ind>G_PI) then
                                    print*,'error simulating initial persistent shock 3'
                                    read*,continue_k
                                end if
                                if (u<sum(Pi_t(1:ind,e))) then
                                    ts_l=ind
                                else
                                    ind=ind+1
                                end if
                            end do
                        end if
            
                        !Sample health status
                        call RANDOM_NUMBER(u)
                        ind=1
                        h2=-9
                        do while (h2==-9)
                            if (ind>G_h+1) then
                                print*,'error simulating initial persistent shock 4'
                                read*,continue_k
                            end if
                            if (u<sum(H_sm(h,1:ind,t_l,y,e))) then !H_sm(h,:,t_l,y,e)
                                h2=ind
                                h=h2
                            else
                                ind=ind+1
                            end if
                        end do
                    end if
                    if (h==2) then
                        counter_bh(t_l)=counter_bh(t_l)+1
                    end if
            
                    gross_earnings(t_l)=taxable_income(y,e,t_l,min(h,G_h),pi_l,ts_l)   
            
                    !If dead go to next individual
                    if (h==G_h+1) then
                        exit
                    end if
                    

                    counter_income(t_l)=counter_income(t_l)+1
                    panel_income(counter_income(t_l),t_l)=income_grid(y,e,t_l,min(h,G_h),pi_l,ts_l) !gross_earnings(t_l) 

            
                    if (t_l<=T_R) then
                        if (h==1 .and. pi_old<G_PI .and. pi_old>0) then
                            counter_un(t_l)=counter_un(t_l)+1
                            if (pi_l==G_PI) then 
                                pr_unemployed(t_l,e,y)=pr_unemployed(t_l,e,y)+1.0d0
                            end if
                        end if
                    end if
            
                end do;
                if (h<G_h+1) then
                    coeff=1.0d0
                    if (gross_earnings(t_R)==0.0) then
                        coeff=0.7d0
                    end if
                    if (sum(gross_earnings(t_R-1:t_R))==0.0) then
                        coeff=0.0d0
                    end if 
                    call sort(gross_earnings,T_R)
                    aime=sum(gross_earnings(T_R-16:T_R))/17.0d0/2.0d0/12.0d0 
                    pia=0.9d0*min(aime,0.895d0)
                    if (aime>0.895) then
                        pia=pia+(min(aime,5.397d0)-0.895d0)*0.32d0
                    end if
                    if (aime>5.397d0) then
                        pia=pia+(aime-5.397d0)*0.15d0
                    end if
                    pia=coeff*pia
                    counter(pi_l)=counter(pi_l)+1
            
                    pension(counter(pi_l),pi_l)=pia*12.0d0*2.0d0    
                    counter_income(T_R+1:T)=counter_income(T_R+1:T)+1
                    panel_income(counter_income(T_R+1:T),T_R+1:T)=pension(counter(pi_l),pi_l)
                end if
            end do
            do pi_l=1,G_PI
                PI_grid(pi_l,e,y)=sum(pension(1:counter(pi_l),pi_l))/dble(counter(pi_l))
            end do
            do t_l=1,T
                av_income_panel(t_l,e,y)=sum(panel_income(1:counter_income(t_l),t_l))/dble(counter_income(t_l)) 
                std_income(t_l,e,y)=sqrt(sum((panel_income(1:counter_income(t_l),t_l)-av_income_panel(t_l,e,y))**2.0d0)/dble(counter_income(t_l)-1))
                if (t_l<=T_R) then
                    pr_unemployed(t_l,e,y)=pr_unemployed(t_l,e,y)/dble(counter_un(t_l)) 
                    pr_bh(t_l,e,y)=counter_bh(t_l)/dble(counter_income(t_l))
                end if             
            end do
        end do;end do
!av_income_panel(:,1,1)
end subroutine  

    

