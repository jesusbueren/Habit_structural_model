subroutine solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini,p50_delta)
    use nrtype;use state_space_dim
    implicit none
    real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF)::a_policy
    real(DP),dimension(T,G_educ,G_types,4),intent(out)::asset_distribution
    real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF)::VSL
    real(DP),dimension(G_DF,G_educ,G_types),intent(out)::av_V_ini
    real(DP),dimension(G_educ),intent(out)::av_VSL
    integer::e,y
    real(DP),dimension(indv_sim,G_h,T)::panel_delta_assets
    integer,dimension(G_h,T)::counter_h
    real(DP),dimension(G_educ,G_h,T),intent(out)::p50_delta
    INTERFACE
        subroutine simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta,lambda_ref)
            use nrtype;use state_space_dim;use var_first_step;use preference_p
            implicit none
            real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::a_policy
            real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::VSL
            real(DP),dimension(G_educ,G_types),optional::lambda_ref
            real(DP),dimension(T,G_educ,G_types,4),intent(out)::asset_distribution
            real(DP),dimension(G_educ),intent(out)::av_VSL
            real(DP),dimension(G_DF,G_educ,G_types),intent(out)::av_V_ini
            real(DP),dimension(G_educ,G_h,T),intent(out)::p50_delta
        END subroutine simulate_model
    end interface
    
    call solve_model(a_policy,VSL)
    
    !open(unit=9,file='C:\Users\jbueren\OneDrive - Istituto Universitario Europeo\endo_health\sol.txt')
    !    read(9,*) a_policy,VSL
    !close (9)
    
        

    call simulate_model(a_policy,VSL,asset_distribution,av_VSL,av_V_ini,p50_delta)

    
    
    
    
end subroutine
    