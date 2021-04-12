module state_space_dim
    use nrtype
    implicit none
    integer,parameter::first_age=25, last_age=98 !First and last age in the model
    integer,parameter::T=(last_age-first_age+1)/2, T_R=(65-first_age+1)/2, T_50=(50-first_age+1)/2 !T: Model periods, T_R: time at retirement, T_50: period at age 50
    integer,parameter::G_educ=3,G_PI=5,G_types=2,G_DF=1,G_h=2 !G_educ: education levels,G_PI: pension income,G_types: health behavior types,G_DF: discount factors,G_h: health types
    integer,parameter::G_nzz=5,G_nkk=50  !L_nzz: grid for the persistent shock 
    character(LEN=42)::path="C:\Users\jbueren\Google Drive\endo_health\"  
end module
    
module var_first_step
    use nrtype;use state_space_dim
    implicit none
    real(DP),dimension(G_h+1,G_h+1,T,G_types,G_educ)::H_fs=-9.0d0 !Health transition probabilities in the model
    real(DP),dimension(G_PI)::PI_grid=(/4.3d0,15.1d0,20.4d0,25.7d0,47.5d0/) ! Median of permanent income for each quantile
    real(DP),dimension(G_nzz,G_nzz)::Pi_m=-9.0d0 !Transition probabilities of the medical shock
    real(DP),dimension(G_PI,T,G_h,G_nzz)::m_grid=-9.0d0 !Medical grid
end module