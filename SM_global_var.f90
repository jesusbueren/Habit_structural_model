module state_space_dim
    use nrtype
    implicit none
    integer,parameter::first_age_sm=26, last_age_sm=98,age_svl=27 !First and last age in the model
    integer,parameter::T=(last_age_sm-first_age_sm+1)/2, T_R=(67-first_age_sm+1)/2, T_50=(50-first_age_sm+1)/2,T_svl=(age_svl-first_age_sm+1)/2 !T: Model periods, T_R: time at retirement, T_50: period at age 50
    integer,parameter::G_educ=3,G_PI=10,G_types=3,G_DF=2,G_h=2,G_cohorts=5 !G_educ: education levels,G_PI: pension income,G_types: health behavior types,G_DF: discount factors,G_h: health types
    integer,parameter::G_nzz=10,G_nkk=100 !G_nzz: grid for the persistent medical shock ; G_nkk: grid for assets  
    integer,parameter::indv_sim=20000
    real(DP)::tol=1.0d-11
end module
    
module var_first_step
    use nrtype;use state_space_dim
    implicit none
    integer,parameter::reference_cohort=3
    real(DP),parameter::a_min=0.0d0,a_max=5000.0d0,r=(1.0d0+0.02d0)**2.0d0-1.0d0
    real(DP)::lambda=0.940772d0,tau=0.158466d0,tau_mcr=0.029d0,tau_ss=0.124d0,ss_bar=128.400d0,av_income=51.4d0
    real(DP)::step=(a_max-a_min)/dble(G_nkk-1)  
    real(DP),dimension(G_nkk)::a_grid
    real(DP),dimension(G_h+1,G_h+1,T,G_types,G_educ)::H_sm=-9.0d0 !Health transition probabilities in the model
    real(DP),dimension(G_PI,G_educ,G_types)::PI_grid,PI_grid2
    real(DP),dimension(G_nzz,1)::pr0_p
    real(DP),dimension(G_educ,T,G_h,G_nzz)::m_grid=-9.0d0 !Medical grid
    real(DP),dimension(G_PI,G_PI,T,G_h,G_educ)::Pi_p !Transition probabilities of the persistent income shock
    real(DP),dimension(G_PI,T,G_h,G_educ)::Pi_p_0 !Stationary distribution of the persistent income shock
    real(DP),dimension(G_PI,G_educ)::Pi_t !Transition probabilities of the transitory income shock
    real(DP),dimension(G_types,G_educ,T_R,G_h,G_PI,G_PI)::income_grid=-9.0d0,taxable_income=-9.0d0,income_tax=-9.0d0,gross_annual_income2=-9.0d0 !income grid
    real(DP),dimension(T)::n_bar=(/2.25d0,2.50807d0,2.81898d0,3.08435d0,3.25147d0,3.36645d0,3.43471d0,3.41466d0,3.34425d0,3.2399d0,3.07015d0,2.84822d0,2.71829d0,2.53466d0,2.36906d0,2.27883d0,2.14743d0,2.10747d0,2.06412d0,2.0274d0,&
                                1.98d0,1.93d0,1.9d0,1.84d0,1.8d0,1.73d0,1.7d0,1.64d0,1.61d0,1.62d0,1.57d0,1.52d0,1.52d0,1.52d0,1.52d0,1.52d0/)
    real(DP),dimension(G_educ,G_types)::fraction_types
    real(DP),dimension(G_educ)::fraction_e
    real(DP),dimension(G_h,G_educ,G_types)::fraction_h_ey
end module
    
    
module preference_p
    use nrtype; use state_space_dim
    implicit none
    !initial guess of parameters
    real(DP),dimension(G_DF)::betas=-9.0d0,pr_betas_u
    real(DP),dimension(G_types,G_educ)::pr_betas=1.0d0
    real(DP),parameter::RRA=2.0d0
    real(DP)::beq_cur=100.0d0,c_floor=-9.0d0,beq_mu=0.0d0,best_obj_fct=1.0d0/0.0d0,b_bar=3.55d0,delta_h=0.0d0,c_bar=0.0d0
    real(DP)::RRA_beq=RRA
    integer,parameter::PAR=14
    real(DP)::beta_max=1.1d0, beta_min=0.6d0
end module
    
module initial_p
    use nrtype; use state_space_dim
    implicit none
    !initial guess of parameters
    real(DP),dimension(G_DF)::betas_ini=(/0.67d0,0.89d0/)
    real(DP)::beq_cur_ini=345.3d0 ,c_floor_ini=3.15d0,beq_mu_ini=2415.0d0
    real(DP),dimension(G_types,G_educ)::pr_betas_ini=reshape((/0.452,0.58,0.79,0.30,0.54,0.7,0.21,0.68,0.88/),shape(pr_betas_ini))
end module
    

    
    