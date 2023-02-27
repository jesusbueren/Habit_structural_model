subroutine solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini)
    use nrtype;use state_space_dim
    implicit none
    real(DP),dimension(G_nkk,T,G_nzz,G_h,G_PI,G_educ,G_types,G_DF)::a_policy
    real(DP),dimension(G_nkk,G_nzz,G_h,G_PI,G_educ,G_types,G_DF)::V_ini
    real(DP),dimension(T,G_educ,G_types,4),intent(out)::asset_distribution
    real(DP),dimension(G_nkk,G_nzz,G_h,G_PI,G_educ,G_types,G_DF)::VSL
    real(DP),dimension(G_educ,G_types),intent(out)::av_VSL,av_V_ini
    
    call solve_model(a_policy,V_ini,VSL)
    

    call simulate_model(a_policy,VSL,V_ini,asset_distribution,av_VSL,av_V_ini)
    
    
    
    
end subroutine
    