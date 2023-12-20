module state_space_dim
    use nrtype
    implicit none
    integer,parameter::first_age_sm=26, last_age_sm=98,age_svl=36 !First and last age in the model
    integer,parameter::T=(last_age_sm-first_age_sm+1)/2, T_R=(67-first_age_sm+1)/2, T_50=(50-first_age_sm+1)/2,T_svl=(age_svl-first_age_sm+1)/2 !T: Model periods, T_R: time at retirement, T_50: period at age 50
    integer,parameter::G_educ=3,G_PI=10,G_types=2,G_DF=1,G_h=2,G_cohorts=5 !G_educ: education levels,G_PI: pension income,G_types: health behavior types,G_DF: discount factors,G_h: health types
    integer,parameter::G_nzz=10,G_nkk=200 !G_nzz: grid for the persistent medical shock ; G_nkk: grid for assets  
    integer,parameter::indv_sim=500000
    real(DP)::tol=1.0d-11
end module
    
module var_first_step
    use nrtype;use state_space_dim
    implicit none
    integer::reference_cohort=2
    real(DP),parameter::a_min=0.0d0,a_max=10000.0d0
    real(DP)::r=(1.0d0+0.02d0)**2.0d0-1.0d0
    real(DP)::tau_mcr=0.029d0,tau_ss=0.124d0,av_income=51.4d0,ss_bar=128.400d0,lambda=0.873964d0,tau=0.108002d0
    real(DP)::step=(a_max-a_min)/dble(G_nkk-1)  
    real(DP),dimension(G_nkk)::a_grid
    real(DP),dimension(G_h+1,G_h+1,T,G_types,G_educ)::H_sm=-9.0d0 !Health transition probabilities in the model
    real(DP),dimension(G_types,G_educ,G_h+1)::LE_sm=-9.0d0
    real(DP),dimension(G_PI,G_educ,G_types)::PI_grid,PI_grid2
    real(DP),dimension(G_nzz,1)::pr0_p
    real(DP),dimension(G_educ,T,G_h,G_nzz)::m_grid=-9.0d0 !Medical grid
    real(DP),dimension(G_PI,G_PI,T,G_h,G_educ)::Pi_p !Transition probabilities of the persistent income shock
    real(DP),dimension(G_PI,T,G_h,G_educ)::Pi_p_0 !Stationary distribution of the persistent income shock
    real(DP),dimension(G_PI,G_educ)::Pi_t !Transition probabilities of the transitory income shock
    real(DP),dimension(G_types,G_educ,T_R,G_h,G_PI,G_PI)::income_grid=-9.0d0,taxable_income=-9.0d0,income_tax=-9.0d0,gross_annual_income2=-9.0d0 !income grid
    real(DP),dimension(T)::n_bar=(/2.25d0,2.50807d0,2.81898d0,3.08435d0,3.25147d0,3.36645d0,3.43471d0,3.41466d0,3.34425d0,3.2399d0,3.07015d0,2.84822d0,2.71829d0,2.53466d0,2.36906d0,2.27883d0,2.14743d0,2.10747d0,2.06412d0,2.0274d0,&
                                1.98d0,1.93d0,1.9d0,1.84d0,1.8d0,1.73d0,1.7d0,1.64d0,1.61d0,1.62d0,1.57d0,1.52d0,1.52d0,1.52d0,1.52d0,1.52d0/)
    real(DP),dimension(G_educ,G_types,G_cohorts)::fraction_types
    real(DP),dimension(G_educ,G_cohorts)::fraction_e
    real(DP),dimension(G_h,G_educ,G_types)::fraction_h_ey
    real(DP),dimension(G_cohorts)::born_cohort=(/10,30,50,70,80/),tuition=(/2.826d0,2.197d0,6.673d0,10.001d0,-9.0d0/)
end module
    
    
module preference_p
    use nrtype; use state_space_dim
    implicit none
    !initial guess of parameters
    real(DP),dimension(G_DF)::betas=0.98d0,pr_betas_u
    real(DP),dimension(G_types,G_educ)::pr_betas=1.0d0
    real(DP),parameter::RRA=1.0d0
    real(DP)::beq_cur=100.0d0,c_floor=-9.0d0,beq_mu=0.0d0,best_obj_fct=1.0d0/0.0d0,delta_h=0.0d0,c_bar=0.0d0
    real(DP),dimension(G_h)::b_bar=0.0d0
    real(DP)::RRA_beq=RRA
    integer,parameter::PAR=4
    real(DP)::beta_max=1.0d0, beta_min=0.95d0
end module
    
module initial_p
    use nrtype; use state_space_dim
    implicit none
    !initial guess of parameters
    real(DP),dimension(G_DF)::betas_ini=0.98d0 !(/0.89d0,0.937d0/)
    real(DP)::beq_cur_ini=100.76d0 ,c_floor_ini=5.8d0,beq_mu_ini=15.35d0,delta_h_ini=0.2d0
    real(DP),dimension(G_types,G_educ)::pr_betas_ini=1.0d0 !reshape((/0.511d0,0.575d0,0.518d0,0.385d0,0.552d0,0.566d0,0.263d0,0.597d0,0.614d0/),shape(pr_betas_ini))
end module
    
    
module second_step
    use nrtype; use state_space_dim
    implicit none
    integer::PAR_2=6,counterfactual=0,DGP_sim=0
    real(dp),dimension(G_df)::pr_beta_un
    real(DP),dimension(G_df,G_educ,G_types,G_cohorts)::joint_pr
    real(DP),dimension(G_df,G_educ,G_types)::av_V_ini  
    real(DP),dimension(G_df,G_educ,G_types,G_cohorts)::av_V_ini_all     
    real(DP),dimension(G_df,G_educ,G_types,G_cohorts)::joint_pr_model
end module  
    
    


    
    