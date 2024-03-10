subroutine cross_sectional_dim()
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
    
    open(unit=9,file='parameters_first_stage.txt')
        read(9,*) p_fs
    close(9)
    
    !Read value function in the benchmark model
    open(unit=9,file='v_ini.txt')
        read(9,*) av_V_ini_all 
    close(9)
    
    !Solve the value function if everyone had the same income as college
    reference_cohort=2
    call load_income_risk() !college premium
    call load_medical_expenses() !adjust tuition
    income_store=income_grid
    income_grid(:,1,:,:,:,:)=income_grid(:,3,:,:,:,:)
    income_grid(:,1,:,:,:,:)=income_grid(:,3,:,:,:,:)
    PI_grid(:,1,:)=PI_grid(:,3,:)
    PI_grid(:,2,:)=PI_grid(:,3,:)
    Pi_p(:,:,:,:,1)=Pi_p(:,:,:,:,3)
    Pi_p(:,:,:,:,2)=Pi_p(:,:,:,:,3)
    Pi_p_0(:,:,:,1)=Pi_p_0(:,:,:,3)
    Pi_p_0(:,:,:,2)=Pi_p_0(:,:,:,3)
    Pi_t(:,1)=Pi_t(:,3)
    Pi_t(:,2)=Pi_t(:,3)
    open(unit=9,file='v_out1.txt')
        write(9,*) av_V_ini
    close(9)
    call solve_model(a_policy,VSL)
    call simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta) 
    call direct_effect(p_fs,av_V_ini_all(:,:,:,2),av_V_ini,'c1')
    
    !Solve the value function if everyone had the same health transitions as college
    H_store=H_sm
    H_sm(:,:,:,:,1)=H_sm(:,:,:,:,3)
    H_sm(:,:,:,:,2)=H_sm(:,:,:,:,3)
    reference_cohort=2
    call load_income_risk() !college premium
    call load_medical_expenses() !adjust tuition
    call solve_model(a_policy,VSL)
    call simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta)
    
    open(unit=9,file='v_out2.txt')
        write(9,*) av_V_ini
    close(9)
    call direct_effect(p_fs,av_V_ini_all(:,:,:,2),av_V_ini,'c2')
    
    pause

end subroutine
    
    
    
!subroutine cross_sectional_dim()
!    use nrtype; use preference_p; use global_var; use var_first_step; use second_step
!    implicit none
!    real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF)::a_policy
!    real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF)::VSL
!    real(DP),dimension(T,G_educ,G_types,4)::asset_distribution
!    real(DP),dimension(G_educ)::av_VSL
!    real(DP),dimension(G_DF,G_educ,G_types,G_cohorts)::av_V_bk
!    integer::e,y,e_ref,y_ref,df_l,c_l,y_l,e_l,c_l2
!    real(DP)::lambda_max,lambda_min
!    real(DP),dimension(G_educ,G_types)::lambda_c
!    real(DP),dimension(indv_sim,G_h,T)::panel_delta_assets
!    integer,dimension(G_h,T)::counter_h
!    real(DP),dimension(G_educ,G_h,T)::p50_delta
!    real(DP),dimension(PAR_2)::p_fs,p_bk
!    real(DP),dimension(G_h+1,G_h+1,T,G_types,G_educ)::H_store
!    real(DP):: smthg
!    real(DP),dimension(G_types,G_educ,T_R,G_h,G_PI,G_PI)::income_store
!    real(DP),dimension(G_PI,G_educ,G_types)::PI_store
!    integer,dimension(2)::cohort_v=(/2,4/),born=(/1930,1970/),unit_c=(/10,11/)
!    INTERFACE
!        subroutine simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta,lambda_ref)
!            use nrtype;use state_space_dim;use var_first_step;use preference_p
!            implicit none
!            real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::a_policy
!            real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::VSL
!            real(DP),dimension(G_educ,G_types),optional::lambda_ref
!            real(DP),dimension(T,G_educ,G_types,4),intent(out)::asset_distribution
!            real(DP),dimension(G_DF,G_educ,G_types),intent(out)::av_V_ini
!            real(DP),dimension(G_educ),intent(out)::av_VSL
!            real(DP),dimension(G_educ,G_h,T),intent(out)::p50_delta
!        END subroutine simulate_model
!        function obj_function_costs(parameters)
!            use nrtype; use preference_p
!            implicit none
!            real(DP),dimension(:),intent(in)::parameters
!            real(DP)::obj_function_costs
!        end function
!        subroutine direct_effect(parameters_in,parameters_out,V_in,V_out)
!            use nrtype; use preference_p; use global_var;use second_step
!            implicit none
!            real(DP),dimension(:),intent(in)::parameters_in,parameters_out
!            real(DP),dimension(G_df,G_educ,G_types),intent(in)::V_in,V_out
!        end subroutine
!    end interface
!    
!    
!   open(unit=9,file='parameter.txt')
!        read(9,*) c_floor,beq_cur,beq_mu
!    close (9)
!    
!    open(unit=9,file='b_bar_cost.txt')
!        read(9,*)  b_bar
!    close(9)
!    
!    
!    print*,'----------------------'
!    print*,"Parameter"
!    print('(A20,<G_DF>F5.2)'),"beta",betas
!    print('(A20,F10.2)'),"beq cur",beq_cur
!    print('(A20,F10.2)'),"c floor",c_floor
!    print('(A20,F10.2)'),"beq mu",beq_mu
!    print('(A20,F10.2)'),"b_bar",b_bar(1)
!    
!    !!! 0. Benchmark
!    open(unit=9,file='parameters_first_stage.txt')
!        read(9,*) p_fs
!    close(9)
!    p_bk=p_fs
!    
!    open(unit=9,file='v_ini.txt')
!        read(9,*) av_V_ini_all 
!    close(9)
!    av_V_bk=av_V_ini_all
!    counterfactual=-1
!    smthg=obj_function_costs(p_fs)
!    !pause
!    do c_l2=1,2
!        if (c_l2==1) then
!            open(unit=unit_c(c_l2),file='C:\Users\jbueren\Dropbox\habits\Slides\v2\tables\counterfactual_c1.txt')
!        else
!            open(unit=unit_c(c_l2),file='C:\Users\jbueren\Dropbox\habits\Slides\v2\tables\counterfactual_c2.txt')
!        end if
!        c_l=cohort_v(c_l2)
!        write(unit_c(c_l2),'(A2,A50,A3)'),'&','\multicolumn{2}{c}{$\text{Pr}(y=\textsc{pro}|e)$}',' \\'
!        write(unit_c(c_l2),'(A19,I4,A2,A40,A40,A40)'),'\textbf{Born in ',born(c_l2),'}&','$e=\textsc{hsd}$ &','$e=\textsc{cg}$ &','$\Delta$LE  \\'
!        print'(A25,A25,A25,A25)','','Dropout Detrimental','College Detrimental','LE Gradient'
!        print'(A25,A25,A25,A25)','','(%)','(%)','(years)'
!        write(unit_c(c_l2),'(A25,F25.3,A2,F25.3,A2,F25.3,A2)'),'\hline'
!        print'(A25,F25.3,A2,F25.3,A2,F25.3,A2)','Benchmark &',joint_pr_model(1,1,2,c_l)/sum(joint_pr_model(1,1,:,c_l)),'&',joint_pr_model(1,3,2,c_l)/sum(joint_pr_model(1,3,:,c_l)),'&',sum(joint_pr_model(1,3,:,c_l)/sum(joint_pr_model(1,3,:,c_l))*LE_sm(:,3,3))-sum(joint_pr_model(1,1,:,c_l)/sum(joint_pr_model(1,1,:,c_l))*LE_sm(:,1,3)),'\\'
!        write(unit_c(c_l2),'(A25,F23.1,A2,F23.1,A2,F23.1,A10)'),'Benchmark &',joint_pr_model(1,1,2,c_l)/sum(joint_pr_model(1,1,:,c_l))*100.0d0,'&',joint_pr_model(1,3,2,c_l)/sum(joint_pr_model(1,3,:,c_l))*100.0d0,'&',sum(joint_pr_model(1,3,:,c_l)/sum(joint_pr_model(1,3,:,c_l))*LE_sm(:,3,3))-sum(joint_pr_model(1,1,:,c_l)/sum(joint_pr_model(1,1,:,c_l))*LE_sm(:,1,3)),'\\[5pt]'
!    end do
!    
!    !!! 1. No correlation 
!    !counterfactual=1
!    !open(unit=9,file='parameters_first_stage_prem.txt')
!    !    read(9,*) p_fs
!    !close(9)
!    !open(unit=9,file='v_ini.txt')
!    !    read(9,*) av_V_ini_all 
!    !close(9)
!    !p_fs(7:10)=-1.0/0.0d0
!    !
!    !smthg=obj_function_costs(p_fs)
!    !do c_l2=1,2
!    !    c_l=cohort_v(c_l2)
!    !    print'(A25,F25.3,F25.3,F25.3)','No correlation',joint_pr_model(1,1,2,c_l)/sum(joint_pr_model(1,1,:,c_l)),joint_pr_model(1,3,2,c_l)/sum(joint_pr_model(1,3,:,c_l)),sum(joint_pr_model(1,3,:,c_l)/sum(joint_pr_model(1,3,:,c_l))*LE_sm(:,3,3))-sum(joint_pr_model(1,1,:,c_l)/sum(joint_pr_model(1,1,:,c_l))*LE_sm(:,1,3))
!    !    write(unit_c(c_l2),'(A50,F23.1,A2,F23.1,A2,F23.1,A10)'),'\bluetext{$\rho^c_{\textsc{pro},\text{e}}=0$} &',joint_pr_model(1,1,2,c_l)/sum(joint_pr_model(1,1,:,c_l))*100.0d0,'&',joint_pr_model(1,3,2,c_l)/sum(joint_pr_model(1,3,:,c_l))*100.0d0,'&',sum(joint_pr_model(1,3,:,c_l)/sum(joint_pr_model(1,3,:,c_l))*LE_sm(:,3,3))-sum(joint_pr_model(1,1,:,c_l)/sum(joint_pr_model(1,1,:,c_l))*LE_sm(:,1,3)),'\\[5pt]'
!    !end do
!    
!    
!    !! 2. Same income
!    !counterfactual=2
!    !reference_cohort=2
!    !call load_income_risk() !college premium
!    !call load_medical_expenses() !adjust tuition
!    !income_store=income_grid
!    !income_grid(:,3,:,:,:,:)=income_grid(:,1,:,:,:,:)
!    !income_grid(:,2,:,:,:,:)=income_grid(:,1,:,:,:,:)
!    !PI_grid(:,3,:)=PI_grid(:,1,:)
!    !PI_grid(:,2,:)=PI_grid(:,1,:)
!    !Pi_p(:,:,:,:,3)=Pi_p(:,:,:,:,1)
!    !Pi_p(:,:,:,:,2)=Pi_p(:,:,:,:,1)
!    !Pi_p_0(:,:,:,3)=Pi_p_0(:,:,:,1)
!    !Pi_p_0(:,:,:,2)=Pi_p_0(:,:,:,1)
!    !Pi_t(:,3)=Pi_t(:,1)
!    !Pi_t(:,2)=Pi_t(:,1)
!    !call solve_model(a_policy,VSL)
!    !call simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta) 
!    !av_V_ini_all(:,:,:,reference_cohort)=av_V_ini
!    !
!    !reference_cohort=4
!    !call load_income_risk() !college premium
!    !call load_medical_expenses() !adjust tuition
!    !income_store=income_grid
!    !income_grid(:,3,:,:,:,:)=income_grid(:,1,:,:,:,:)
!    !income_grid(:,2,:,:,:,:)=income_grid(:,1,:,:,:,:)
!    !PI_grid(:,3,:)=PI_grid(:,1,:)
!    !PI_grid(:,2,:)=PI_grid(:,1,:)
!    !Pi_p(:,:,:,:,3)=Pi_p(:,:,:,:,1)
!    !Pi_p(:,:,:,:,2)=Pi_p(:,:,:,:,1)
!    !Pi_p_0(:,:,:,3)=Pi_p_0(:,:,:,1)
!    !Pi_p_0(:,:,:,2)=Pi_p_0(:,:,:,1)
!    !Pi_t(:,3)=Pi_t(:,1)
!    !Pi_t(:,2)=Pi_t(:,1)
!    !call solve_model(a_policy,VSL)
!    !call simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta) 
!    !av_V_ini_all(:,:,:,reference_cohort)=av_V_ini
!
!    open(unit=9,file='v_I.txt')
!        read(9,*) av_V_ini_all 
!    close(9)
!    !pause
!    
!    call direct_effect(p_bk(1:8),p_bk(1:8),av_V_bk(:,:,:,2),av_V_ini_all(:,:,:,2))    
!    
!    
!    
!    
!    
!    !pause
!    open(unit=9,file='parameters_first_stage.txt')
!        read(9,*) p_fs
!    close(9)
!
!    smthg=obj_function_costs(p_fs)
!    do c_l2=1,2
!        c_l=cohort_v(c_l2)
!        print'(A25,F25.3,F25.3,F25.3)','Equal Income',joint_pr_model(1,1,2,c_l)/sum(joint_pr_model(1,1,:,c_l)),joint_pr_model(1,3,2,c_l)/sum(joint_pr_model(1,3,:,c_l)),sum(joint_pr_model(1,3,:,c_l)/sum(joint_pr_model(1,3,:,c_l))*LE_sm(:,3,3))-sum(joint_pr_model(1,1,:,c_l)/sum(joint_pr_model(1,1,:,c_l))*LE_sm(:,1,3))
!        write(unit_c(c_l2),'(A140,F23.1,A2,F23.1,A2,F23.1,A10)'),'\bluetext{$w^{\textsc{cg}}_t=w^{\textsc{hsd}}_t$} &',joint_pr_model(1,1,2,c_l)/sum(joint_pr_model(1,1,:,c_l))*100.0d0,'&',joint_pr_model(1,3,2,c_l)/sum(joint_pr_model(1,3,:,c_l))*100.0d0,'&',sum(joint_pr_model(1,3,:,c_l)/sum(joint_pr_model(1,3,:,c_l))*LE_sm(:,3,3))-sum(joint_pr_model(1,1,:,c_l)/sum(joint_pr_model(1,1,:,c_l))*LE_sm(:,1,3)),'\\[5pt]'
!    end do
!
!    !pause
!    !! 3. Same health transitions
!    !counterfactual=3
!    !H_store=H_sm
!    !H_sm(:,:,:,:,3)=H_sm(:,:,:,:,1)
!    !H_sm(:,:,:,:,2)=H_sm(:,:,:,:,1)
!    !reference_cohort=2
!    !call load_income_risk() !college premium
!    !call load_medical_expenses() !adjust tuition
!    !call solve_model(a_policy,VSL)
!    !call simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta) 
!    !av_V_ini_all(:,:,:,reference_cohort)=av_V_ini
!    !
!    !reference_cohort=4
!    !call load_income_risk() !college premium
!    !call load_medical_expenses() !adjust tuition
!    !call solve_model(a_policy,VSL)
!    !call simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta) 
!    !av_V_ini_all(:,:,:,reference_cohort)=av_V_ini
!    
!     open(unit=9,file='v_H.txt')
!        read(9,*) av_V_ini_all 
!    close(9)
!    
!    open(unit=9,file='parameters_first_stage.txt')
!        read(9,*) p_fs
!    close(9)
!    
!    call direct_effect(p_bk(1:8),p_bk(1:8),av_V_bk(:,:,:,2),av_V_ini_all(:,:,:,2))  
!
!    smthg=obj_function_costs(p_fs)
!    
!    do c_l2=1,2
!        c_l=cohort_v(c_l2)
!        print'(A25,F25.3,F25.3,F25.3)','Equal Health ',joint_pr_model(1,1,2,c_l)/sum(joint_pr_model(1,1,:,c_l)),joint_pr_model(1,3,2,c_l)/sum(joint_pr_model(1,3,:,c_l)),sum(joint_pr_model(1,3,:,c_l)/sum(joint_pr_model(1,3,:,c_l))*LE_sm(:,1,3))-sum(joint_pr_model(1,1,:,c_l)/sum(joint_pr_model(1,1,:,c_l))*LE_sm(:,1,3))
!        write(unit_c(c_l2),'(A140,F23.1,A2,F23.1,A2,F23.1,A10)'),'\bluetext{$\Gamma^{\textsc{cg}}_t\left( h^\prime|h \right)=\Gamma^{\textsc{hsd}}_t\left( h^\prime|h \right)$} &',joint_pr_model(1,1,2,c_l)/sum(joint_pr_model(1,1,:,c_l))*100.0d0,'&',joint_pr_model(1,3,2,c_l)/sum(joint_pr_model(1,3,:,c_l))*100.0d0,'&',sum(joint_pr_model(1,3,:,c_l)/sum(joint_pr_model(1,3,:,c_l))*LE_sm(:,1,3))-sum(joint_pr_model(1,1,:,c_l)/sum(joint_pr_model(1,1,:,c_l))*LE_sm(:,1,3)),'\\'
!        write(unit_c(c_l2),'(A25,F23.3,A2,F23.3,A2,F23.3,A2)'),'\hline'
!    end do
!    H_sm=H_store
!
!    
!    
!    
!    
!    
!    
!  
!    
!    !!0
!    !lambda_c=0.0d0
!    !do df_l=1,G_df
!    !    do e=1,G_educ;do y=1,G_types
!    !        e_ref=e;y_ref=2
!    !        if (e/=e_ref .or. y/=y_ref ) then !reference category
!    !            lambda_max=1.0
!    !            lambda_min=0.0d0
!    !
!    !        1   lambda_c(e,y)=(lambda_max+lambda_min)/2
!    !        
!    !            call simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_new,p50_delta,lambda_c)
!    !            if (abs(av_V_new(df_l,e,y)-av_V_ini(df_l,e_ref,y_ref))/av_V_ini(df_l,e_ref,y_ref)>1d-5) then
!    !                if (av_V_new(df_l,e,y)>av_V_ini(df_l,e_ref,y_ref)) then
!    !                    lambda_min=lambda_c(e,y)
!    !                else
!    !                    lambda_max=lambda_c(e,y)
!    !                end if
!    !                !print*,'diff',av_V_new(df_l,e,y),av_V_ini2(df_l,e_ref,y_ref),lambda_c(e,y)
!    !                go to 1
!    !            end if
!    !        end if
!    !        print '(A3,A3,I3,A3,I3,A3,I3,F10.3)','CE','df',df_l,'e',e,'y',y,lambda_c(e,y)
!    !    end do;end do
!    !end do
!    !
!    
!    
!    
!    
!    
!    
!    end subroutine
!    
!    
!    