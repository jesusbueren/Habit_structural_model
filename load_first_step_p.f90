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
        integer,parameter::iterations=4996
        real(dp),dimension(covariates,clusters,L_gender,L_educ,iterations)::c_tr_all
        real(dp),dimension(covariates,clusters,L_gender,L_educ,iterations)::c_tr_d_all
        real(dp),dimension(covariates_habits*habits*types,iterations)::c_habits_all

        real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_h
        real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_d
        real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H
        integer::t_l,h_l
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
        fraction_h_ey=sum(fraction_h(1,:,1,:,:,:),4)/dble(iterations) 
        
        if (reference_cohort==3) then
            fraction_e=(/0.089d0,0.514d0,0.397d0/)!(/0.16d0,0.59d0,0.24d0/)
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



end subroutine
    
subroutine load_medical_expenses()
    use nrtype;use var_first_step;use state_space_dim; use global_var
    implicit none
    integer,parameter::K=15
    real(DP),dimension(K+2)::data_csv
    real(DP),dimension(K,1)::X,beta
    real(DP)::rho,sigma2_u,age,h,z,d1,d2,d3
    real(DP),dimension(G_nzz,1)::mag_p
    real(DP),dimension(K)::beta_insurance_u65,beta_insurance_o65
    real(DP),dimension(G_educ,T,G_h)::fraction_covered!Medical grid
    integer::t_l,e_l,h_l,z_l
    character::pause_k
    
    open(unit=9,file=path//'data\medical_exp.csv')
            read(9,*) data_csv
    close(9)
    open(unit=9,file=path//'data\insurance_u65.csv')
            read(9,*) beta_insurance_u65
    close(9)
    open(unit=9,file=path//'data\insurance_o65.csv')
            read(9,*) beta_insurance_o65
    close(9)
    
    beta(:,1)=data_csv(1:K)
    rho=data_csv(K+1)
    sigma2_u=(data_csv(K+2)**2)*(1-rho**2)
    
    call discretize_shocks(rho,sigma2_u,G_nzz,mag_p,pr0_p,Pi_m)
    
    !Mean expenses from stata regression
    do e_l=1,G_educ;do h_l=1,G_h;do z_l=1,G_nzz;do t_l=1,T
        age=min(dble(first_age_sm+(t_l-1)*2),80.0d0)
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
        X(:,1)=(/d1,d2,d3,age,age**2.0d0,age**3.0d0,h,d1*h,d2*h,d3*h,d1*age,d2*age,d3*age,h*age,1.0d0/)
        m_grid(e_l,t_l,h_l,z_l)=exp(sum(X(:,1)*beta(:,1))+z)
        if (age<65) then
            fraction_covered(e_l,t_l,h_l)=1-sum(X(:,1)*beta_insurance_u65)
        else
            fraction_covered(e_l,t_l,h_l)=1-sum(X(:,1)*beta_insurance_o65)
        end if
        m_grid(e_l,t_l,h_l,z_l)=m_grid(e_l,t_l,h_l,z_l)*(1-fraction_covered(e_l,t_l,h_l))
    end do;end do;end do;end do

end subroutine
    
    
subroutine load_income_risk()
    use nrtype;use var_first_step;use state_space_dim; use global_var
    implicit none
    integer,parameter::K=9
    real(DP),dimension(K,1)::X
    real(DP),dimension(K,G_educ)::beta
    real(DP),dimension(G_educ,G_cohorts):: s2_nu,s2_w,rho,s2_0
    real(DP)::sigma2_u,age,h,z,d1,d2,d3,z2,var_aux
    real(DP),dimension(G_PI,1)::mag_p,mag_tr
    real(DP),dimension(G_PI,G_PI)::Pi_nothing
    real(DP),dimension(G_types)::y_d
    real(DP),dimension(G_cohorts)::cohort_d
    integer::t_l,e_l,h_l,z_l,z_l2,y_l
    character::pause_k
    real(DP),dimension(5,generations,clusters,L_educ,G_cohorts)::labor_force_participation
    real(DP),dimension(5,generations,clusters,L_educ,G_cohorts,2)::labor_force_participation_dynamic
    real(DP)::gross_annual_income
    
    !Labor force participation risk at age 25
    open(unit=9,file=path//'metric_model\Results\labor_force_participation.txt')
            read(9,*) labor_force_participation 
    close(9)
    do t_l=1,T;do h_l=1,G_h;do e_l=1,G_educ
        Pi_p_0(G_PI,t_l,h_l,e_l)=1.0d0-labor_force_participation(5,t_l,h_l,e_l,reference_cohort) !Pi_p_0(10,:,1,1)
    end do;end do;end do
    
    !Labor force participation conditional previous labor force participation
    open(unit=9,file=path//'metric_model\Results\labor_force_participation_dynamic.txt')
            read(9,*) labor_force_participation_dynamic
    close(9)
    do t_l=1,T;do h_l=1,G_h;do e_l=1,G_educ
        Pi_p(1:G_PI-1,G_PI,t_l,h_l,e_l)=1.0d0-labor_force_participation_dynamic(5,t_l,h_l,e_l,reference_cohort,2) ! transition in the labor force -> out of the labor force 
        Pi_p(G_PI,G_PI,t_l,h_l,e_l)=1.0d0-labor_force_participation_dynamic(5,t_l,h_l,e_l,reference_cohort,1) ! transition out of the labor force -> out of the labor force 
    end do;end do;end do
    !Pi_p(10,10,:,1,1)
    !Income risk for those being employed
    open(unit=9,file=path//'metric_model\Results\parameters_income.txt')
            read(9,*) beta,s2_nu,s2_w,rho,s2_0
    close(9)

    !Discretize the transitory shock
    call discretize_shocks(0.0d0,s2_w(1,1),G_PI,mag_tr,Pi_t,Pi_nothing) 
    
    !Compute transition matrix of the persistent component and the grid of the persistent shock
    call discretize_shocks(rho(1,1),s2_nu(1,1),G_PI-1,mag_p(1:G_PI-1,1),Pi_p_0(1:G_PI-1,1,1,1),Pi_p(1:G_PI-1,1:G_PI-1,1,1,1)) 
    do t_l=1,T;do h_l=1,G_h;do e_l=1,G_educ;do z_l=1,G_PI
        if (z_l<G_PI) then
            Pi_p(z_l,1:G_PI-1,t_l,h_l,e_l)=Pi_p(z_l,1:G_PI-1,1,1,1)/sum(Pi_p(z_l,1:G_PI-1,1,1,1))*(1.0d0-Pi_p(z_l,G_PI,t_l,h_l,e_l))            
        end if            
    end do;end do;end do;end do
    
    !Compute the unconditional probability of the persistent shock
    var_aux=s2_0(1,1)
    call pr_grid_shock(G_PI-1,var_aux,mag_p(1:G_PI-1,1),Pi_p_0(1:G_PI-1,1,1,1))
    do t_l=2,T
        var_aux=rho(1,1)**2*var_aux+s2_nu(1,1)
        call pr_grid_shock(G_PI-1,var_aux,mag_p(1:G_PI-1,1),Pi_p_0(1:G_PI-1,t_l,1,1))
    end do
    
    do t_l=1,T;do h_l=1,G_h;do e_l=1,G_educ;do z_l=1,G_PI
        if (z_l<G_PI) then
            Pi_p_0(z_l,t_l,h_l,e_l)=Pi_p_0(z_l,t_l,1,1)/sum(Pi_p_0(1:G_PI-1,t_l,1,1))*(1.0d0-Pi_p_0(G_PI,t_l,h_l,e_l))           
        else
            Pi_p(z_l,1:G_PI-1,t_l,h_l,e_l)=Pi_p_0(1:G_PI-1,t_l,1,1)/sum(Pi_p_0(1:G_PI-1,t_l,1,1))*(1.0d0-Pi_p(z_l,G_PI,t_l,h_l,e_l))
        end if            
    end do;end do;end do;end do
    

    income_grid=0.0d0
    do t_l=1,T_R;do y_l=1,G_types;do e_l=1,G_educ;do h_l=1,G_h;do z_l=1,G_PI;do z_l2=1,G_PI
        if (z_l<G_PI) then
            age=first_age_sm+(t_l-1)*2-70 
            z=mag_p(z_l,1)
            z2=mag_tr(z_l2,1)
            y_d=0.0d0
            y_d(y_l)=1.0d0
            !reference_cohort
            cohort_d(4)=1.0d0
            X(:,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(h_l-1),y_d(2:types),cohort_d(4:5)/) 
            gross_annual_income=exp(sum(X(:,1)*beta(:,e_l))+z+z2)/1000.0d0
            taxable_income(y_l,e_l,t_l,h_l,z_l,z_l2)=min(gross_annual_income,ss_bar)*2.0d0
            income_grid(y_l,e_l,t_l,h_l,z_l,z_l2)=av_income*lambda*(gross_annual_income/av_income)**(1.0d0-tau)*2.0d0 - tau_mcr*gross_annual_income*2.0d0-tau_ss*min(gross_annual_income,ss_bar)*2.0d0       
        else
            income_grid(y_l,e_l,t_l,h_l,z_l,z_l2)=0.0d0
            taxable_income(y_l,e_l,t_l,h_l,z_l,z_l2)=0.0d0
        end if
        if (income_grid(y_l,e_l,t_l,h_l,z_l,z_l2)<0.0d0) then
            print*,''           
        end if
    end do;end do;end do;end do;end do;end do

    call compute_pension()
    
!income_grid(1,3,:,1,3,3)
end subroutine load_income_risk
    
subroutine compute_pension()
    use nrtype;use state_space_dim;use var_first_step
        implicit none
        integer,parameter::indv_sim=50000
        double precision::u
        integer::h,t_l,ind,i_l,y,e,h2,pi_l,df,ts_l,pi_l2
        character::continue_k
        real(DP),dimension(T_R)::gross_earnings
        integer,dimension(G_PI)::counter
        real(DP),dimension(indv_sim,G_PI)::pension
        real(DP),dimension(indv_sim,T)::panel_income
        integer,dimension(T)::counter_income
        integer,dimension(T_R)::counter_un
        real(DP),dimension(T,G_educ,G_types)::av_income_panel,std_income
        real(DP),dimension(T_R,G_educ,G_types)::pr_unemployed
        real(DP)::aime,pia,coeff
        
        open(unit=9,file='wages.txt')
        
        pr_unemployed=0.0d0
        df=G_DF
        do e=1,G_educ; do y=1,G_types
            counter=0
            counter_income=0
            counter_un=0
            panel_income=-9
            do i_l=1,indv_sim;gross_earnings=-9.0d0;aime=-9.0d0
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
                            if (u<sum(Pi_p_0(1:ind,t_l,h,e))) then 
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
                                print*,'error simulating initial persistent shock '
                                read*,continue_k
                            end if
                            if (u<sum(Pi_t(1:ind,1))) then
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
                            pi_l2=-9
                            do while (pi_l2==-9)
                                if (ind>G_PI) then
                                    print*,'error simulating initial persistent shock '
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
                                    print*,'error simulating initial persistent shock '
                                    read*,continue_k
                                end if
                                if (u<sum(Pi_t(1:ind,1))) then
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
                                print*,'error simulating initial persistent shock '
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
            
                    gross_earnings(t_l)=taxable_income(y,e,t_l,min(h,G_h),pi_l,ts_l)   
            
                    !If dead go to next individual
                    if (h==G_h+1) then
                        exit
                    end if
            
                    counter_income(t_l)=counter_income(t_l)+1
                    panel_income(counter_income(t_l),t_l)=income_grid(y,e,t_l,min(h,G_h),pi_l,ts_l)
            
                    if (t_l<=T_R) then
                        if (h==2) then
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
                    pia=0.9*min(aime,0.895)
                    if (aime>0.895) then
                        pia=pia+(min(aime,5.397)-0.895)*0.32
                    end if
                    if (aime>5.397) then
                        pia=pia+(aime-5.397)*0.15
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
                end if
            end do
        end do;end do
!pr_unemployed(:,1,3)
end subroutine  

    

