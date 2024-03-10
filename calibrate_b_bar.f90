subroutine calibrate_b_bar()
    use nrtype; use preference_p; use global_var; use var_first_step; use second_step
    implicit none
    real(DP)::obj_function,b_bar_max,b_bar_min,VSL_data
    real(DP),dimension(T,G_h,G_educ,G_types,4)::asset_distribution
    real(DP),dimension(G_educ)::av_VSL

    real(DP),dimension(G_df)::k
    real(DP),dimension(8,generations,types,L_educ)::dist_assets_data
    real(DP),dimension(G_educ,G_h,T)::p50_delta
    real(DP),dimension(PAR_2+1,PAR_2)::p_g
    real(DP),dimension(PAR_2,PAR_2)::xi
    real(DP),dimension(PAR_2+1)::y
    real(DP)::ftol,u,fret=1.0d-5,fret_old=999.0d0
    integer::y_l,t_l,e_l,it,df_l,p_l,iter,c_l
    character::end_k
    interface
        function obj_function_costs(parameters)
            use nrtype; use preference_p
            implicit none
            real(DP),dimension(:),intent(in)::parameters
            real(DP)::obj_function_costs
        end function  
        subroutine amoeba(p,y,ftol,func,iter)
                use nrtype
                implicit none
                INTEGER(I4B), INTENT(OUT) :: iter
	            REAL(DP), INTENT(IN) :: ftol
	            REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
	            REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p
	            INTERFACE
		            FUNCTION func(x)
		            USE nrtype
		            IMPLICIT NONE
		            REAL(DP), DIMENSION(:), INTENT(IN) :: x
		            REAL(DP) :: func
		            END FUNCTION func
	            END INTERFACE
        end subroutine
        SUBROUTINE powell(p,xi,ftol,iter,fret)
	        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	        USE nr, ONLY : linmin
	        IMPLICIT NONE
	        REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
	        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: xi
	        INTEGER(I4B), INTENT(OUT) :: iter
	        REAL(DP), INTENT(IN) :: ftol
	        REAL(DP), INTENT(OUT) :: fret
        end subroutine
        function inverse_logistic_mapping(x, beta_min, beta_max) result(inverse_mapped_value)
    use nrtype
    implicit none
        real(DP), intent(in) ::beta_min, beta_max
        real(DP), intent(in) :: x
        real(DP) :: inverse_mapped_value
    end function inverse_logistic_mapping
    end interface
    
    open(unit=9,file='parameter.txt')
        read(9,*) c_floor,beq_cur,beq_mu
    close (9)

    print*,'----------------------'
    print*,"Parameter"
    print('(A20,<2>F5.2)'),"beta",betas
    print('(A20,F10.2)'),"beq cur",beq_cur
    print('(A20,F10.2)'),"c floor",c_floor
    print('(A20,F10.2)'),"beq mu",beq_mu
    print('(A20,<9>F6.3)'),"pr low beta",pr_betas
    


    pr_beta_un(1)=1.0d0
    
    
    if (G_DF==1) then
        pr_betas=1.0d0
    end if
    
    b_bar_min=-2.0d0
    b_bar_max=2.0d0
    it=0
    b_bar=-9.0d0
    VSL_data=1000.0d0
    do while (abs(av_VSL(2)-VSL_data)>50.0d0 .and. it<10) 
        b_bar=(b_bar_max+b_bar_min)/2.0d0
        !Solve model
        call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini,p50_delta)
        if (av_VSL(2)>VSL_data) then
            b_bar_max=b_bar(1)
        else
            b_bar_min=b_bar(1)
        end if
        it=it+1
        print*,it,b_bar(1),av_VSL(1),av_VSL(2)   
    end do
    
    open(unit=9,file='b_bar_cost.txt')
        write(9,*)  b_bar
    close(9)

    
    !Solve model for different cohorts
    av_V_ini_all=0.0d0
    av_V_ini_all(:,:,:,reference_cohort)=av_V_ini
    do c_l=2,G_cohorts-1,2
        reference_cohort=c_l
        call load_income_risk() !college premium
        call load_medical_expenses() !adjust tuition
        call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini,p50_delta)
        av_V_ini_all(:,:,:,c_l)=av_V_ini
        print*,it,b_bar(1),av_VSL(1),av_VSL(2) 
    end do
    reference_cohort=2
    
    open(unit=9,file='v_ini.txt')
        write(9,*) av_V_ini_all 
    close(9)
    
    do c_l=1,G_cohorts;do df_l=1,G_DF;do e_l=1,G_educ;do y_l=1,types
        if (df_l==1) then
            joint_pr(df_l,e_l,y_l,c_l)=fraction_types(e_l,y_l,c_l)*fraction_e(e_l,c_l)*pr_betas(y_l,e_l)
        else
            joint_pr(df_l,e_l,y_l,c_l)=fraction_types(e_l,y_l,c_l)*fraction_e(e_l,c_l)*(1.0d0-pr_betas(y_l,e_l))
        end if
    end do;end do;end do;end do
    
    !!Generate a DGP to check identification (comment if don't want to check for identification)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !p_g(1,1:3)=(/11.8d0,9.0d0,35.2d0/) !av cost of protective, High-school, College   
    !p_g(1,4:6)=(/11.2d0,2.0d0,15.0d0/) !Variance of shock of protective, High-school, College 
    !p_g(1,4:6)=log(p_g(1,4:6))
    !p_g(1,7:9)=(/0.1d0,0.3d0,0.4d0/) !Variance-cov matrix of correlations
    !p_g(1,7:9)=log((p_g(1,7:9)+1.0)/(2.0d0-p_g(1,7:9)-1.0d0))
    !
    !DGP_sim=1
    !print*,obj_function_costs(p_g(1,:))
    


    !Estimate costs:
    !!!!!!!!!!!!!!!!    
    !(df_l,e_l,y_l,c_l)
    p_g(1,1:2)=(/4.53d0,6.55d0 /)!mean taste shock for protective and education 
    p_g(1,3)=log(5.17d0) !multiplicative factor to the cost of college 
    p_g(1,4:5)=(/log(3.3d0),log(30.69d0)/) !variance of taste shock of protective and education
    !p_g(1,6)=inverse_logistic_mapping(0.1d0, -1.0d0, 1.0d0) 
    it=0

1    do p_l=1,PAR_2+1
        if (p_l>1) then
            p_g(p_l,:)=p_g(1,:)
            call random_number(u)
            p_g(p_l,p_l-1)=p_g(1,p_l-1)*0.9d0 
        end if
        y(p_l)=obj_function_costs(p_g(p_l,:))
        pause
    end do

    print*,y 
    ftol=1d-8
    call amoeba(p_g,y,ftol,obj_function_costs,iter) 
    open(unit=9,file='parameters_first_stage.txt')
        write(9,*) p_g(1,:)
    close(9)
    !pause
    print '(A20,<PAR_2>F10.5)','parameters A',p_g(1,:)
    print*,'obj fct',y(1)
    !pause
    xi=0.0d0
    do p_l=1,PAR_2
        xi(p_l,p_l)=1.0d0
    end do

    ftol=1d-3
    call powell(p_g(1,:),xi,ftol,iter,fret)
    open(unit=9,file='parameters_first_stage.txt')
        write(9,*) p_g(1,:)
    close(9)
    it=it+1
    print '(A20,<PAR_2>F10.5)','parameters P',p_g(1,:)
    print*,'obj fct',fret
    print*,'*************************************************'
    print*,'crit',abs(fret_old-fret)/fret
    if (abs(fret_old-fret)/fret>1.0d-2) then
        fret_old=fret
        print*,'it',it
        !pause
        go to 1
    end if 
    
    open(unit=9,file='parameters_first_stage.txt')
        write(9,*) p_g(1,:)
    close(9)

  
   print*,'End optimization'
pause
end subroutine  
    
    
    
function obj_function_costs(parameters)
    use nrtype; use preference_p; use global_var;use second_step
    implicit none
    real(DP),dimension(:),intent(in)::parameters
    real(DP),dimension(G_df,G_educ,G_types)::rand_ey,U
    integer,dimension(3)::maxloc_U
    real(DP)::obj_function_costs
    real(DP)::mu_y,mu_e,var_y,var_e,tau_cg,var_nu
    real(DP)::rho_ey
    real(DP),dimension(G_df,G_educ,G_types)::cost_ey
    integer:: e_l,y_l,df_l,i_l,c_l
    integer(8),dimension(1)::seed=321
    interface
        subroutine input2p(parameters,mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey)
            use nrtype; use state_space_dim
            implicit none
            real(DP),dimension(:),intent(in)::parameters
            real(DP),intent(out)::mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey
        end subroutine
        subroutine direct_effect(p_in,V_in,V_out)
            use nrtype; use preference_p
            implicit none
            real(DP),dimension(:),intent(in)::p_in
            real(DP),dimension(G_df,G_educ,G_types),intent(in)::V_in,V_out
        end subroutine
    end interface
    
    !print*,parameters


    call random_seed(PUT=seed)
    joint_pr_model=0.0d0
    !open(unit=9,file='eps.txt')
    call input2p(parameters,mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey)
    print '(A20,<10>F10.2)','parameters P',mu_y,mu_e,tau_cg,var_y,var_e,var_nu,rho_ey
    do c_l=2,G_cohorts-1,2
        do i_l=1,200000
        call p2shock(mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey,rand_ey)
        U=av_V_ini_all(:,:,:,c_l)-rand_ey 
        maxloc_U=maxloc(U)
        joint_pr_model(maxloc_U(1),maxloc_U(2),maxloc_U(3),c_l)=joint_pr_model(maxloc_U(1),maxloc_U(2),maxloc_U(3),c_l)+1.0d0 
        !if (maxloc_U(2)==3) then
        !    print*,''
        !end if
        !write(9,'(I3,I3,F10.3)'),c_l,maxloc_U(2),maxloc_U(3),rand_ey(1,1,1),log(rand_ey(1,2,2)),log(rand_ey(1,3,2))
    end do;end do
    !close(9)
    do c_l=2,G_cohorts-1,2
        joint_pr_model(:,:,:,c_l)=joint_pr_model(:,:,:,c_l)/sum(joint_pr_model(:,:,:,c_l))      
    end do
    
    if (counterfactual==0) then
        open(unit=9,file='costs_e_y.txt')
    end if
    
    if (DGP_sim==1) then
        joint_pr=joint_pr_model
        DGP_sim=0
    end if
    
    do df_l=1,G_df;do c_l=2,G_cohorts-1,2;do e_l=1,G_educ;do y_l=1,G_types
        write(9,'(A4,I4,A4,I4,A4,I4,A4,I4,F20.5,F20.5,F20.5,F20.5)'),'df_l',df_l,'c_l',c_l,'e_l',e_l,'y_l',y_l,cost_ey(df_l,e_l,y_l),av_V_ini_all(df_l,e_l,y_l,c_l),joint_pr(df_l,e_l,y_l,c_l),joint_pr_model(df_l,e_l,y_l,c_l)
    end do;end do;end do;end do
    close(9)
       
    
    obj_function_costs=0.0d0
    
    !Match pr(y=1|e)
    do df_l=1,G_df;do e_l=1,G_educ;do c_l=2,G_cohorts-1,2
        !if (e_l==3 .or. e_l==1) then
        !    obj_function_costs=obj_function_costs+((joint_pr_model(df_l,e_l,1,c_l)/sum(joint_pr_model(df_l,e_l,:,c_l))-joint_pr(df_l,e_l,1,c_l)/sum(joint_pr(df_l,e_l,:,c_l)))/(joint_pr(df_l,e_l,1,c_l)/sum(joint_pr(df_l,e_l,:,c_l))))**2.0d0*100.0d0
        !else
            obj_function_costs=obj_function_costs+((joint_pr_model(df_l,e_l,1,c_l)/sum(joint_pr_model(df_l,e_l,:,c_l))-joint_pr(df_l,e_l,1,c_l)/sum(joint_pr(df_l,e_l,:,c_l)))/(joint_pr(df_l,e_l,1,c_l)/sum(joint_pr(df_l,e_l,:,c_l))))**2.0d0
        !end if
    end do;end do;end do
    
    !Match pr(e)
    do df_l=1,G_df;do e_l=2,G_educ;do c_l=2,G_cohorts-1,2
        obj_function_costs=obj_function_costs+((sum(joint_pr_model(df_l,e_l,:,c_l))-sum(joint_pr(df_l,e_l,:,c_l)))/sum(joint_pr(df_l,e_l,:,c_l)))**2.0d0
    end do;end do;end do
    
    !Match pr(y)
    !do df_l=1,G_df;do y_l=2,G_types;do c_l=2,G_cohorts-1,2
    !    obj_function_costs=obj_function_costs+(sum(joint_pr_model(df_l,:,y_l,c_l))-sum(joint_pr(df_l,:,y_l,c_l)))**2.0d0
    !end do;end do;end do
    

    
    if (isnan(obj_function_costs))then
        obj_function_costs=10000.0d0
    end if
    
    
    print*,obj_function_costs
    
    open(unit=9,file='parameters_first_stage_prem.txt')
        write(9,*) parameters
    close(9)
    
    !call direct_effect(parameters,av_V_ini_all(:,:,:,2),av_V_ini_all(:,:,:,4))

    
    end function
    
    subroutine input2p(parameters,mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey)
        use nrtype; use state_space_dim
        implicit none
        real(DP),dimension(:),intent(in)::parameters
        real(DP),intent(out)::mu_y,mu_e,var_y,var_e,tau_cg,var_nu,rho_ey
        interface
            function logistic_mapping(beta, beta_min, beta_max) result(mapped_value)
            use nrtype
            implicit none
                real(DP), intent(in) :: beta
                real(DP), intent(in) ::beta_min, beta_max
                real(DP) :: mapped_value
                real(DP):: sigmoid_beta
            end function logistic_mapping
        end interface
    
        mu_y=parameters(1)
        mu_e=parameters(2)
        tau_cg=exp(parameters(3))
    
        !Location
        var_y=exp(parameters(4))
        var_e=exp(parameters(5))
        var_nu=0.0d0
        rho_ey=0.0d0!logistic_mapping(parameters(6), -1.0d0, 1.0d0)
        
    
    end subroutine
    
    

    
    
    
    double precision function c4_normal_01 (  )
!------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  double precision, parameter :: r4_pi=3.14159265358979323846264338327950288419716939937510
  double precision:: v1
  double precision:: v2
  double precision:: x_c
  double precision:: x_r 
  call random_number(v1)
  call random_number(v2)
  x_r =sqrt(-2.0d0*log(v1))*cos(2.0d0*r4_pi*v2)
  x_c =sqrt(-2.0d0*log(v1))*sin(2.0d0*r4_pi*v2)
  c4_normal_01=x_r
  return
end