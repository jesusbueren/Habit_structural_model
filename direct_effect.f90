subroutine direct_effect(parameters_in,parameters_out,V_in,V_out)
    use nrtype; use preference_p; use global_var;use second_step
    implicit none
    real(DP),dimension(:),intent(in)::parameters_in,parameters_out
    real(DP),dimension(G_df,G_educ,G_types),intent(in)::V_in,V_out
    real(DP),dimension(G_df,G_educ,G_types)::rand_ey,U
    integer,dimension(3)::maxloc_U
    real(DP)::sigma_hs,sigma_cg,sigma_y
    real(DP)::eps_pro,eps_hs,eps_col,eps_01_hs,eps_01_col
    real(DP)::rho_ey
    integer:: e_l,y_l,df_l,i_l
    integer(8),dimension(1)::seed=321
    real(DP),dimension(G_educ)::cost_e    
    real(DP),dimension(G_types)::cost_y
    real(DP),dimension(G_df,G_educ,G_types)::cost_ey
    real(DP),dimension(G_df,G_educ,G_types)::joint_pr_in,joint_pr_out
    interface
        subroutine input2p(parameters,cost_ey,sigma_y,sigma_hs,sigma_cg,rho_ey)
            use nrtype; use state_space_dim
            implicit none
            real(DP),dimension(:),intent(in)::parameters    
            real(DP),dimension(G_df,G_educ,G_types),intent(out)::cost_ey
            real(DP),intent(out)::sigma_hs,sigma_cg,sigma_y,rho_ey
        end subroutine
    end interface

    
    !print*,parameters
    
    rho_ey=0.0d0
    call random_seed(PUT=seed)

    joint_pr_in=0.0d0
    joint_pr_out=0.0d0
    
    do i_l=1,500000
        !Best option for the individual with parameters_in and V_in 
        call input2p(parameters_in(1:6),cost_ey,sigma_y,sigma_hs,sigma_cg,rho_ey)
        call random_seed(GET=seed)
        call p2shock(sigma_y,sigma_hs,sigma_cg,rho_ey,rand_ey)
        U=V_in-cost_ey+rand_ey  
        maxloc_U=maxloc(U)
        e_l=maxloc_U(2) !education decision made
        y_l=maxloc_U(3) !lifestyle decision made
        joint_pr_in(maxloc_U(1),maxloc_U(2),maxloc_U(3))=joint_pr_in(maxloc_U(1),maxloc_U(2),maxloc_U(3))+1.0d0
        
        !Best option for the individual with parameters_out and V_out if education choice constrained to choice made with parameters_in and V_in 
        call input2p(parameters_out(1:6),cost_ey,sigma_y,sigma_hs,sigma_cg,rho_ey)
        call random_seed(PUT=seed)
        call p2shock(sigma_y,sigma_hs,sigma_cg,rho_ey,rand_ey)
        U=-1.0d0/0.0d0
        U(:,e_l,:)=V_out(:,e_l,:)-cost_ey(:,e_l,:)+rand_ey(:,e_l,:) 
        maxloc_U=maxloc(U)
        joint_pr_out(maxloc_U(1),maxloc_U(2),maxloc_U(3))=joint_pr_out(maxloc_U(1),maxloc_U(2),maxloc_U(3))+1.0d0

    end do


    joint_pr_in=joint_pr_in/sum(joint_pr_in)   
    joint_pr_out=joint_pr_out/sum(joint_pr_out)   

    

    
end subroutine