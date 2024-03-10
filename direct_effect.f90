subroutine direct_effect(parameters_in,V_in,V_out,str)
    use nrtype; use preference_p; use global_var;use second_step
    implicit none
    real(DP),dimension(:),intent(in)::parameters_in
    real(DP),dimension(G_df,G_educ,G_types),intent(in)::V_in,V_out
    real(DP),dimension(G_df,G_educ,G_types)::rand_ey,U,V_s
    integer,dimension(3)::maxloc_U
    real(DP)::sigma_hs,sigma_cg,sigma_nu,rho_hs_cg
    real(DP)::eps_pro,eps_hs,eps_col,eps_01_hs,eps_01_col
    real(DP)::mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey
    integer:: e_l,y_l,df_l,i_l,e_l_te,y_l_de,y_l_te
    integer(8),dimension(1)::seed=321
    real(DP),dimension(G_educ)::cost_e    
    real(DP),dimension(3)::cost_y,sigma_y
    real(DP),dimension(G_df,G_educ,G_types)::cost_ey
    real(DP),dimension(G_df,G_educ,G_types)::joint_pr_in,joint_pr_out,joint_pr_out2
    character(len=*), intent(in) :: str
    interface
        subroutine input2p(parameters,mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey)
            use nrtype; use state_space_dim
            implicit none
            real(DP),dimension(:),intent(in)::parameters
            real(DP),intent(out)::mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey
        end subroutine
    end interface
    

    call random_seed(PUT=seed)

    joint_pr_in=0.0d0
    joint_pr_out=0.0d0
    open(unit=9,file='eps_'//str//'.txt')
    do i_l=1,200000
        !Best option for the individual with parameters_in and V_in 
        call input2p(parameters_in,mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey)
        call p2shock(mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey,rand_ey)
        U=V_in-rand_ey  
        maxloc_U=maxloc(U)
        e_l=maxloc_U(2) !education decision made
        y_l=maxloc_U(3) !lifestyle decision made
        joint_pr_in(maxloc_U(1),maxloc_U(2),maxloc_U(3))=joint_pr_in(maxloc_U(1),maxloc_U(2),maxloc_U(3))+1.0d0
        
        !Total effect best option for the individual with V_out if education choice unconstrained 
        U=-1.0d0/0.0d0
        U=V_out-rand_ey   
        maxloc_U=maxloc(U)
        e_l_te=maxloc_U(2) 
        y_l_te=maxloc_U(3) 
        
        !Direct effect: best option for the individual with V_out if education choice constrained to choice made with parameters_in and V_in 
        U=-1.0d0/0.0d0
        U(:,e_l,:)=V_out(:,e_l,:)-rand_ey(:,e_l,:) 
        maxloc_U=maxloc(U)
        joint_pr_out(maxloc_U(1),maxloc_U(2),maxloc_U(3))=joint_pr_out(maxloc_U(1),maxloc_U(2),maxloc_U(3))+1.0d0
        y_l_de=maxloc_U(3) 
        
        write(9,'(<5>I3,F10.3)'),e_l,e_l_te,y_l,y_l_de,y_l_te,rand_ey(1,1,1) 
        
        !Selection effect: decision made in terms of y by CG & HSG if they had to same returns as HSD on their health behavior
        U=-1.0d0/0.0d0
        V_s(:,e_l,:)=V_in(:,1,:)
        U(:,e_l,:)=V_s(:,e_l,:)-rand_ey(:,e_l,:) 
        maxloc_U=maxloc(U)
        joint_pr_out2(maxloc_U(1),maxloc_U(2),maxloc_U(3))=joint_pr_out2(maxloc_U(1),maxloc_U(2),maxloc_U(3))+1.0d0

    end do
    close(9)

    joint_pr_in=joint_pr_in/sum(joint_pr_in)   
    joint_pr_out=joint_pr_out/sum(joint_pr_out)   
    joint_pr_out2=joint_pr_out2/sum(joint_pr_out2)   

print*,'direct effect computed'
pause
    
end subroutine