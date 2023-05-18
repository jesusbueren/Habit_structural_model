subroutine solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini)
    use nrtype;use state_space_dim
    implicit none
    real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF)::a_policy
    real(DP),dimension(T,G_educ,G_types,4),intent(out)::asset_distribution
    real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF)::VSL
    real(DP),dimension(G_DF,G_educ,G_types),intent(out)::av_VSL,av_V_ini
    integer::e,y
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
    
    call solve_model(a_policy,VSL)
    
    !open(unit=9,file='C:\Users\jbueren\OneDrive - Istituto Universitario Europeo\endo_health\sol.txt')
    !    write(9,*) a_policy,VSL
    !close (9)
    
        
    do e=1,G_educ; do y=1,G_types
        call simulate_model(e,y,a_policy,VSL,asset_distribution(:,e,y,:),av_VSL(:,e,y),av_V_ini(:,e,y))
    end do;end do
    
    
    
    
end subroutine
    