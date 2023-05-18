subroutine welfare_analysis()
    use nrtype; use preference_p; use global_var; use var_first_step
    implicit none
    real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF)::a_policy
    real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF)::VSL
    real(DP),dimension(T,G_educ,G_types,4)::asset_distribution
    real(DP),dimension(G_DF,G_educ,G_types)::av_VSL,av_V_ini,av_V_new
    integer::e,y,e_ref,y_ref,df_l
    real(DP)::lambda_max,lambda_min
    real(DP),dimension(G_educ,G_types)::lambda_c
    INTERFACE
        subroutine simulate_model(e,y,a_policy,VSL,asset_distribution,av_VSL,av_V_ini,lambda_ref)
            use nrtype;use state_space_dim;use var_first_step;use preference_p
            implicit none
            integer,intent(in)::e,y
            real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::a_policy
            real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::VSL
            real(DP),optional::lambda_ref
            real(DP),dimension(T,4),intent(out)::asset_distribution
            real(DP),dimension(G_DF),intent(out)::av_VSL,av_V_ini
        END subroutine simulate_model
    end interface
    
    
    open(unit=9,file='parameter.txt')
        read(9,*) betas,beq_cur,c_floor,beq_mu,pr_betas
    close (9)
    
    !Define b_bar
    b_bar=0.05d0
 
    
    print*,'----------------------'
    print*,"Parameter"
    print('(A20,<G_DF>F5.2)'),"beta",betas
    print('(A20,F10.2)'),"beq cur",beq_cur
    print('(A20,F10.2)'),"c floor",c_floor
    print('(A20,F10.2)'),"beq mu",beq_mu
    print('(A20,F10.2)'),"b_bar",b_bar
    
    
    !call solve_model(a_policy,VSL)
    
    open(unit=9,file='C:\Users\jbueren\OneDrive - Istituto Universitario Europeo\endo_health\sol.txt')
        read(9,*) a_policy,VSL
    close (9)
    
    
    do e=1,G_educ;do y=1,G_types
        call simulate_model(e,y,a_policy,VSL,asset_distribution(:,e,y,:),av_VSL(:,e,y),av_V_ini(:,e,y)) !av_V_ini(1,3,:)
    end do; end do
    
    
    !0
    lambda_c=0.0d0
    do df_l=1,G_df
        do e=1,G_educ;do y=1,G_types
            e_ref=e;y_ref=3
            if (e/=e_ref .or. y/=y_ref ) then !reference category
                lambda_max=1.0
                lambda_min=0.0d0
    
            1   lambda_c(e,y)=(lambda_max+lambda_min)/2
            
                call simulate_model(e,y,a_policy,VSL,asset_distribution(:,e,y,:),av_VSL(:,e,y),av_V_new(:,e,y),lambda_c(e,y))
                if (abs(av_V_new(df_l,e,y)-av_V_ini(df_l,e_ref,y_ref))/av_V_ini(df_l,e_ref,y_ref)>1d-5) then
                    if (av_V_new(df_l,e,y)>av_V_ini(df_l,e_ref,y_ref)) then
                        lambda_min=lambda_c(e,y)
                    else
                        lambda_max=lambda_c(e,y)
                    end if
                    !print*,'diff',av_V_new(df_l,e,y),av_V_ini(df_l,e_ref,y_ref),lambda_c(e,y)
                    go to 1
                end if
            end if
            print '(A3,A3,I3,A3,I3,A3,I3,F10.3)','CE','df',df_l,'e',e,'y',y,lambda_c(e,y)
        end do;end do
    end do
    
    
    
    
    
    
    
end subroutine