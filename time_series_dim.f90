subroutine time_series_dim()
    use nrtype; use preference_p; use global_var; use var_first_step; use second_step
    implicit none
    real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF)::a_policy
    real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF)::VSL
    real(DP),dimension(T,G_educ,G_types,4)::asset_distribution
    real(DP),dimension(G_educ)::av_VSL
    real(DP),dimension(G_DF,G_educ,G_types)::av_V_new
    integer::e,y,e_ref,y_ref,df_l,c_l,y_l,e_l,c_l2
    real(DP),dimension(PAR_2)::p_fs,p_fs_new
    real(DP),dimension(G_h+1,G_h+1,T,G_types,G_educ)::H_store
    real(DP):: smthg
    real(DP),dimension(G_types,G_educ,T_R,G_h,G_PI,G_PI)::income_store
    real(DP),dimension(G_PI,G_educ,G_types)::PI_store
    integer,dimension(2)::cohort_v=(/2,4/),born=(/1930,1970/),unit_c=(/10,11/)
    real(DP),dimension(G_df,G_educ,G_types,G_cohorts)::joint_pr_bmk,joint_pr_bmk_de
    real(DP),dimension(G_df,G_educ,G_types,G_cohorts,3)::joint_pr_c,joint_pr_c_de
    real(DP),dimension(G_educ,G_h,T)::p50_delta
    INTERFACE
        subroutine simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta,lambda_ref)
            use nrtype;use state_space_dim;use var_first_step;use preference_p
            implicit none
            real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::a_policy
            real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::VSL
            real(DP),dimension(G_educ,G_types),optional::lambda_ref
            real(DP),dimension(T,G_educ,G_types,4),intent(out)::asset_distribution
            real(DP),dimension(G_DF,G_educ,G_types),intent(out)::av_V_ini
            real(DP),dimension(G_educ),intent(out)::av_VSL
            real(DP),dimension(G_educ,G_h,T),intent(out)::p50_delta
        END subroutine simulate_model
        function obj_function_costs(parameters)
            use nrtype; use preference_p
            implicit none
            real(DP),dimension(:),intent(in)::parameters
            real(DP)::obj_function_costs
        end function
        subroutine direct_effect(p_in,V_in,V_out,str)
            use nrtype; use preference_p
            implicit none
            real(DP),dimension(:),intent(in)::p_in
            real(DP),dimension(G_df,G_educ,G_types),intent(in)::V_in,V_out
            character(len=*), intent(in) :: str
        end subroutine
    end interface
    
    
    !Load structural parameters
    open(unit=9,file='parameter.txt')
        read(9,*) c_floor,beq_cur,beq_mu,betas
    close (9)
    
    open(unit=9,file='b_bar_cost.txt')
        read(9,*)  b_bar
    close(9)
    
    !! 0. Benchmark
    open(unit=9,file='parameters_first_stage.txt')
        read(9,*) p_fs
    close(9)
    
    counterfactual=-1
    !do c_l=2,G_cohorts-1
    !    reference_cohort=c_l
    !    call load_income_risk() !college premium
    !    call load_medical_expenses() !adjust tuition
    !    call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini,p50_delta)
    !    av_V_ini_all(:,:,:,c_l)=av_V_ini
    !end do
    open(unit=9,file='v_ini.txt')
        read(9,*) av_V_ini_all 
    close(9)

    smthg=obj_function_costs(p_fs)
    joint_pr_bmk=joint_pr_model 
    
    !compute direct effect
    joint_pr_bmk_de=joint_pr_bmk 
    call direct_effect(p_fs,av_V_ini_all(:,:,:,2),av_V_ini_all(:,:,:,4),'c3')
    joint_pr_bmk_de(:,:,:,4)=joint_pr_model(:,:,:,4) 
    
    do c_l=1,3
    
        if (c_l==1) then
            !! 1. only changes in wages
            !a.adjust correlation to 1930
            open(unit=9,file='parameters_first_stage.txt')
                read(9,*) p_fs_new
            close(9) 
    
            !b.adjust wages 1990 
            reference_cohort=4
            call load_income_risk() 
    
            !c.adjust tuition 1930 
            reference_cohort=2
            call load_medical_expenses()  
            
        elseif (c_l==2) then
            !! 2. only changes in tuition
            !a.adjust correlation to 1930
            open(unit=9,file='parameters_first_stage.txt')
                read(9,*) p_fs_new
            close(9) 
    
            !b.adjust wages 1930 
            reference_cohort=2
            call load_income_risk() 
    
            !c.adjust tuition 1970 
            reference_cohort=4
            call load_medical_expenses()  
            
        elseif (c_l==3) then
            !! 3. only changes in initial conditions
            !a.adjust correlation to 1970
            open(unit=9,file='parameters_first_stage.txt')
                read(9,*) p_fs_new
            close(9)  
    
            !b.adjust wages 1930 
            reference_cohort=2
            call load_income_risk() 
        
            !c.adjust tuition 1930 
            reference_cohort=2
            call load_medical_expenses()  
        end if
        
    
    !Solve model
    reference_cohort=4
    call solve_model(a_policy,VSL)
    call simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta) 
    av_V_ini_all(:,:,:,reference_cohort)=av_V_ini
    
    !Total effect: direct + selection
    smthg=obj_function_costs(p_fs_new)
    pause
    joint_pr_c(:,:,:,:,c_l)=joint_pr_model
    
    !direct effect
    joint_pr_c_de(:,:,:,:,c_l)=joint_pr_c(:,:,:,:,c_l) 
    call direct_effect(p_fs,av_V_ini_all(:,:,:,2),av_V_ini_all(:,:,:,4),'c4')

    joint_pr_c_de(:,:,:,4,c_l)=joint_pr_model(:,:,:,4)
    

    if (c_l==1) then
        open(unit=10,file='C:\Users\jbueren\Dropbox\habits\Slides\v2\tables\counterfactual_ts_w.txt')
         write(10,'(A125)'),'& \multicolumn{3}{c}{Full model} & \multicolumn{3}{c}{Wage increase} \\'
    elseif (c_l==2) then
        open(unit=10,file='C:\Users\jbueren\Dropbox\habits\Slides\v2\tables\counterfactual_ts_tf.txt')
         write(10,'(A125)'),'& \multicolumn{3}{c}{Full model} & \multicolumn{3}{c}{Tuition increase} \\'
    elseif (c_l==3) then
        open(unit=10,file='C:\Users\jbueren\Dropbox\habits\Slides\v2\tables\counterfactual_ts_ic.txt')
         write(10,'(A125)'),'& \multicolumn{3}{c}{Full model} & \multicolumn{3}{c}{Correlation increase} \\'
    end if
        
    
   
    write(10,'(A125)'),'\cline{2-4} \cline{5-7}'
    write(10,'(A125)'),'& \textsc{cg} & \textsc{hsd} & $\Delta$ & \textsc{cg} & \textsc{hsd} & $\Delta$ \\'
    write(10,'(A125)'),'\cline{2-7}'
    write(10,'(A125)'),'\\[0ex]'
    !Delta Pr(e)
    write(10,'(A125,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A10)'),'$\Delta\text{Pr}\left( \text{e} \right)$ &',(sum(joint_pr_bmk(:,3,:,4))-sum(joint_pr_bmk(:,3,:,2)))*100.0d0, &
                                                                                        '&',(sum(joint_pr_bmk(:,1,:,4))-sum(joint_pr_bmk(:,1,:,2)))*100.0d0, &
                                                                                        '&',(sum(joint_pr_bmk(:,3,:,4))-sum(joint_pr_bmk(:,3,:,2))-(sum(joint_pr_bmk(:,1,:,4))-sum(joint_pr_bmk(:,1,:,2))))*100.0d0, &
                                                                                        '&',(sum(joint_pr_c(:,3,:,4,c_l))-sum(joint_pr_c(:,3,:,2,c_l)))*100.0d0, &
                                                                                        '&',(sum(joint_pr_c(:,1,:,4,c_l))-sum(joint_pr_c(:,1,:,2,c_l)))*100.0d0, &
                                                                                        '&',(sum(joint_pr_c(:,3,:,4,c_l))-sum(joint_pr_c(:,3,:,2,c_l))-(sum(joint_pr_c(:,1,:,4,c_l))-sum(joint_pr_c(:,1,:,2,c_l))))*100.0d0,'\\[2ex]'
    !Delta Pr(y=pro|e)
    write(10,'(A125,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A10)'),'$\Delta\text{Pr}\left( y=\textsc{pro}\, |\, \text{e} \right)$  &',(joint_pr_bmk(1,3,1,4)/sum(joint_pr_bmk(1,3,:,4))-joint_pr_bmk(1,3,1,2)/sum(joint_pr_bmk(1,3,:,2)))*100.0d0, &
                                                                                        '&',(joint_pr_bmk(1,1,1,4)/sum(joint_pr_bmk(1,1,:,4))-joint_pr_bmk(1,1,1,2)/sum(joint_pr_bmk(1,1,:,2)))*100.0d0, &
                                                                                        '&',((joint_pr_bmk(1,3,1,4)/sum(joint_pr_bmk(1,3,:,4))-joint_pr_bmk(1,3,1,2)/sum(joint_pr_bmk(1,3,:,2)))*100.0d0)-((joint_pr_bmk(1,1,1,4)/sum(joint_pr_bmk(1,1,:,4))-joint_pr_bmk(1,1,1,2)/sum(joint_pr_bmk(1,1,:,2)))*100.0d0), &
                                                                                        '&',(joint_pr_c(1,3,1,4,c_l)/sum(joint_pr_c(1,3,:,4,c_l))-joint_pr_c(1,3,1,2,c_l)/sum(joint_pr_c(1,3,:,2,c_l)))*100.0d0, &
                                                                                        '&',(joint_pr_c(1,1,1,4,c_l)/sum(joint_pr_c(1,1,:,4,c_l))-joint_pr_c(1,1,1,2,c_l)/sum(joint_pr_c(1,1,:,2,c_l)))*100.0d0, &
                                                                                        '&',((joint_pr_c(1,3,1,4,c_l)/sum(joint_pr_c(1,3,:,4,c_l))-joint_pr_c(1,3,1,2,c_l)/sum(joint_pr_c(1,3,:,2,c_l)))*100.0d0)-((joint_pr_c(1,1,1,4,c_l)/sum(joint_pr_c(1,1,:,4,c_l))-joint_pr_c(1,1,1,2,c_l)/sum(joint_pr_c(1,1,:,2,c_l)))*100.0d0),'\\'
    !Delta Pr(y=pro|e): direct effect
    write(10,'(A125,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A10)'),'\ \ \ \ direct effect  &',(joint_pr_bmk_de(1,3,1,4)/sum(joint_pr_bmk_de(1,3,:,4))-joint_pr_bmk_de(1,3,1,2)/sum(joint_pr_bmk_de(1,3,:,2)))*100.0d0, &
                                                                                        '&',(joint_pr_bmk_de(1,1,1,4)/sum(joint_pr_bmk_de(1,1,:,4))-joint_pr_bmk_de(1,1,1,2)/sum(joint_pr_bmk_de(1,1,:,2)))*100.0d0, &
                                                                                        '&',((joint_pr_bmk_de(1,3,1,4)/sum(joint_pr_bmk_de(1,3,:,4))-joint_pr_bmk_de(1,3,1,2)/sum(joint_pr_bmk_de(1,3,:,2)))*100.0d0)-((joint_pr_bmk_de(1,1,1,4)/sum(joint_pr_bmk_de(1,1,:,4))-joint_pr_bmk_de(1,1,1,2)/sum(joint_pr_bmk_de(1,1,:,2)))*100.0d0), &
                                                                                        '&',(joint_pr_c_de(1,3,1,4,c_l)/sum(joint_pr_c_de(1,3,:,4,c_l))-joint_pr_c_de(1,3,1,2,c_l)/sum(joint_pr_c_de(1,3,:,2,c_l)))*100.0d0, &
                                                                                        '&',(joint_pr_c_de(1,1,1,4,c_l)/sum(joint_pr_c_de(1,1,:,4,c_l))-joint_pr_c_de(1,1,1,2,c_l)/sum(joint_pr_c_de(1,1,:,2,c_l)))*100.0d0, &
                                                                                        '&',((joint_pr_c_de(1,3,1,4,c_l)/sum(joint_pr_c_de(1,3,:,4,c_l))-joint_pr_c_de(1,3,1,2,c_l)/sum(joint_pr_c_de(1,3,:,2,c_l)))*100.0d0)-((joint_pr_c_de(1,1,1,4,c_l)/sum(joint_pr_c_de(1,1,:,4,c_l))-joint_pr_c_de(1,1,1,2,c_l)/sum(joint_pr_c_de(1,1,:,2,c_l)))*100.0d0),'\\[2ex]'
    
    !Delta LE(e)
    write(10,'(A125,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A10)'),'$\Delta$LE$(\text{e})$   &',(sum(joint_pr_bmk(1,3,:,4)/sum(joint_pr_bmk(1,3,:,4))*LE_sm(:,3,3))-sum(joint_pr_bmk(1,3,:,2)/sum(joint_pr_bmk(1,3,:,2))*LE_sm(:,3,3))), &
                                                                                        '&',(sum(joint_pr_bmk(1,1,:,4)/sum(joint_pr_bmk(1,1,:,4))*LE_sm(:,1,3))-sum(joint_pr_bmk(1,1,:,2)/sum(joint_pr_bmk(1,1,:,2))*LE_sm(:,1,3))), &
                                                                                        '&',((sum(joint_pr_bmk(1,3,:,4)/sum(joint_pr_bmk(1,3,:,4))*LE_sm(:,3,3))-sum(joint_pr_bmk(1,3,:,2)/sum(joint_pr_bmk(1,3,:,2))*LE_sm(:,3,3))))-((sum(joint_pr_bmk(1,1,:,4)/sum(joint_pr_bmk(1,1,:,4))*LE_sm(:,1,3))-sum(joint_pr_bmk(1,1,:,2)/sum(joint_pr_bmk(1,1,:,2))*LE_sm(:,1,3)))), &
                                                                                        '&',(sum(joint_pr_c(1,3,:,4,c_l)/sum(joint_pr_c(1,3,:,4,c_l))*LE_sm(:,3,3))-sum(joint_pr_c(1,3,:,2,c_l)/sum(joint_pr_c(1,3,:,2,c_l))*LE_sm(:,3,3))), &
                                                                                        '&',(sum(joint_pr_c(1,1,:,4,c_l)/sum(joint_pr_c(1,1,:,4,c_l))*LE_sm(:,1,3))-sum(joint_pr_c(1,1,:,2,c_l)/sum(joint_pr_c(1,1,:,2,c_l))*LE_sm(:,1,3))), &
                                                                                        '&',((sum(joint_pr_c(1,3,:,4,c_l)/sum(joint_pr_c(1,3,:,4,c_l))*LE_sm(:,3,3))-sum(joint_pr_c(1,3,:,2,c_l)/sum(joint_pr_c(1,3,:,2,c_l))*LE_sm(:,3,3))))-((sum(joint_pr_c(1,1,:,4,c_l)/sum(joint_pr_c(1,1,:,4,c_l))*LE_sm(:,1,3))-sum(joint_pr_c(1,1,:,2,c_l)/sum(joint_pr_c(1,1,:,2,c_l))*LE_sm(:,1,3)))),'\\'
    
    !Delta LE(e): direct effect
    write(10,'(A125,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A2,F25.1,A10)'),'\ \ \ \ direct effect   &',(sum(joint_pr_bmk_de(1,3,:,4)/sum(joint_pr_bmk_de(1,3,:,4))*LE_sm(:,3,3))-sum(joint_pr_bmk_de(1,3,:,2)/sum(joint_pr_bmk_de(1,3,:,2))*LE_sm(:,3,3))), &
                                                                                        '&',(sum(joint_pr_bmk_de(1,1,:,4)/sum(joint_pr_bmk_de(1,1,:,4))*LE_sm(:,1,3))-sum(joint_pr_bmk_de(1,1,:,2)/sum(joint_pr_bmk_de(1,1,:,2))*LE_sm(:,1,3))), &
                                                                                        '&',((sum(joint_pr_bmk_de(1,3,:,4)/sum(joint_pr_bmk_de(1,3,:,4))*LE_sm(:,3,3))-sum(joint_pr_bmk_de(1,3,:,2)/sum(joint_pr_bmk_de(1,3,:,2))*LE_sm(:,3,3))))-((sum(joint_pr_bmk_de(1,1,:,4)/sum(joint_pr_bmk_de(1,1,:,4))*LE_sm(:,1,3))-sum(joint_pr_bmk_de(1,1,:,2)/sum(joint_pr_bmk_de(1,1,:,2))*LE_sm(:,1,3)))), &
                                                                                        '&',(sum(joint_pr_c_de(1,3,:,4,c_l)/sum(joint_pr_c_de(1,3,:,4,c_l))*LE_sm(:,3,3))-sum(joint_pr_c_de(1,3,:,2,c_l)/sum(joint_pr_c_de(1,3,:,2,c_l))*LE_sm(:,3,3))), &
                                                                                        '&',(sum(joint_pr_c_de(1,1,:,4,c_l)/sum(joint_pr_c_de(1,1,:,4,c_l))*LE_sm(:,1,3))-sum(joint_pr_c_de(1,1,:,2,c_l)/sum(joint_pr_c_de(1,1,:,2,c_l))*LE_sm(:,1,3))), &
                                                                                        '&',((sum(joint_pr_c_de(1,3,:,4,c_l)/sum(joint_pr_c_de(1,3,:,4,c_l))*LE_sm(:,3,3))-sum(joint_pr_c_de(1,3,:,2,c_l)/sum(joint_pr_c_de(1,3,:,2,c_l))*LE_sm(:,3,3))))-((sum(joint_pr_c_de(1,1,:,4,c_l)/sum(joint_pr_c_de(1,1,:,4,c_l))*LE_sm(:,1,3))-sum(joint_pr_c_de(1,1,:,2,c_l)/sum(joint_pr_c_de(1,1,:,2,c_l))*LE_sm(:,1,3)))),'\\[2ex]'
        
    close(10)
    end do
   
    
    end subroutine
    
    
