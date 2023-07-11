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
    real(DP),dimension(PAR_2+1)::y
    real(DP)::ftol=1d-12
    integer::y_l,t_l,e_l,it,df_l,p_l,iter
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

    end interface
    
    open(unit=9,file='parameter.txt')
        read(9,*) betas,c_floor,beq_cur,beq_mu,pr_betas,obj_function
    close (9)
 
    print*,'----------------------'
    print*,"Parameter"
    print('(A20,<2>F5.2)'),"beta",betas
    print('(A20,F10.2)'),"beq cur",beq_cur
    print('(A20,F10.2)'),"c floor",c_floor
    print('(A20,F10.2)'),"beq mu",beq_mu
    print('(A20,<9>F6.3)'),"pr low beta",pr_betas

    df_l=1
    pr_beta_un(df_l)=0.0d0
    do e_l=1,G_educ;do y_l=1,G_types
        pr_beta_un(df_l)=pr_beta_un(df_l)+pr_betas(y_l,e_l)*fraction_types(e_l,y_l)*fraction_e(e_l)
    end do;end do
    pr_beta_un(G_df)=1.0d0-pr_beta_un(1)

    if (G_DF==1) then
        pr_betas=1.0d0
    end if
    
    b_bar_min=0.0d0
    b_bar_max=100.0d0
    it=0
    
    VSL_data=6000.0d0

    
    !do while (abs(av_VSL(1)-VSL_data)>100.0d0 .and. it<10) 
        b_bar=10.0!(b_bar_max+b_bar_min)/2.0d0
        !Solve model
        call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini,p50_delta)
        if (av_VSL(1)>VSL_data) then
            b_bar_max=b_bar(1)
        else
            b_bar_min=b_bar(1)
        end if
        it=it+1
        print*,it,b_bar,av_VSL(1)   
    !end do
        
    open(unit=9,file='v_ini.txt')
        write(9,*) av_V_ini
    close(9)
    pause
    !Estimate costs:
    !!!!!!!!!!!!!!!!    
    
    !Specification I. 9 cost parameter, var=1.
    do df_l=1,G_DF;
        do e_l=1,G_educ;do y_l=1,types
            if (df_l==1) then
                joint_pr(df_l,e_l,y_l)=fraction_types(e_l,y_l)*fraction_e(e_l)*pr_betas(y_l,e_l)
            else
                joint_pr(df_l,e_l,y_l)=fraction_types(e_l,y_l)*fraction_e(e_l)*(1.0d0-pr_betas(y_l,e_l))
            end if
        end do;end do
        joint_pr(df_l,:,:)=joint_pr(df_l,:,:)/sum(joint_pr(df_l,:,:))
        !print*,'cov',sum((log(joint_pr(df_l,:,:))-sum(log(joint_pr(df_l,:,:)))/dble(G_educ*G_types))*(av_V_ini(df_l,:,:)-sum(av_V_ini(df_l,:,:))/dble(G_educ*G_types)))/dble(G_educ*G_types)
        !print*,'var x',sum((av_V_ini(df_l,:,:)-sum(av_V_ini(df_l,:,:))/dble(G_educ*G_types))**2.0d0)/dble(G_educ*G_types)
        !k(df_l)=sum((log(joint_pr(df_l,:,:))-sum(log(joint_pr(df_l,:,:)))/dble(G_educ*G_types))*(av_V_ini(df_l,:,:)-sum(av_V_ini(df_l,:,:))/dble(G_educ*G_types)))/sum((av_V_ini(df_l,:,:)-sum(av_V_ini(df_l,:,:))/dble(G_educ*G_types))**2.0d0)
        !k(df_l)=1.0d0/k(df_l)
        !print*,'k=',k(df_l)
        !print*,'cost: e x y'
        !
        !do e_l=1,L_educ;do y_l=1,types
        !    cost_ey(df_l,e_l,y_l)=av_V_ini(df_l,e_l,y_l)-av_V_ini(df_l,1,3)-k(df_l)*(log(joint_pr(df_l,e_l,y_l))-log(joint_pr(df_l,1,3)))
        !    print'(A4,I4,A4,I4,A4,I4,A4,F20.5,A5,F20.5)','df_l',df_l,'e_l',e_l,'y_l',y_l,'%',joint_pr(df_l,e_l,y_l),'cost',cost_ey(df_l,e_l,y_l)
        !    write(9,'(A4,I4,A4,I4,A4,I4,F20.5,F20.5,F20.5,F20.5)'),'df_l',df_l,'e_l',e_l,'y_l',y_l,cost_ey(df_l,e_l,y_l),av_V_ini(df_l,e_l,y_l),joint_pr(df_l,e_l,y_l),exp(1.0d0/k(df_l)*av_V_ini(df_l,e_l,y_l))/sum(exp(1.0d0/k(df_l)*av_V_ini(df_l,:,:)))
        !end do;end do
    end do

    
    p_g(1,1:PAR_2-1)=(/5.0d0,5.0d0,4.0d0,15.0d0/)
    p_g(1,PAR_2)=8.0d0
    do p_l=1,PAR_2+1
        if (p_l>1) then
            p_g(p_l,:)=p_g(1,:)
            p_g(p_l,p_l-1)=p_g(1,p_l-1)*0.95d0
        end if
        y(p_l)=obj_function_costs(p_g(p_l,:))
        !pause
    end do
    
    print*,size(p_g,2),size(p_g,1)-1,size(y)-1

    
    call amoeba(p_g,y,ftol,obj_function_costs,iter)
    
    
    
    open(unit=9,file='b_bar_costs.txt')
        write(9,*) b_bar,cost_ey,k
    close(9)
    

    

   
    
    
    end subroutine  
    
    
    
function obj_function_costs(parameters)
    use nrtype; use preference_p; use global_var;use second_step
    implicit none
    real(DP),dimension(:),intent(in)::parameters
    real(DP),dimension(G_df,G_educ,G_types)::joint_pr_model
    real(DP)::obj_function_costs
    real(DP)::sigma
    integer:: e_l,y_l,df_l
    
    print*,parameters
    cost_y(3)=0.0d0
    cost_e(1)=0.0d0
    cost_y(1:2)=parameters(1:2)
    cost_e(2:3)=parameters(3:4)
    sigma=parameters(5)
    
    do df_l=1,G_df;do e_l=1,G_educ;do y_l=1,G_types
        cost_ey(df_l,e_l,y_l)=cost_e(e_l)+cost_y(y_l)
        !exp(1.0d0/sigma*av_V_ini(df_l,e_l,y_l))/sum(exp(1.0d0/sigma*av_V_ini(df_l,:,:)))
    end do;end do; end do
    
    open(unit=9,file='costs_e_y.txt')
    do df_l=1,G_df+1;do e_l=1,G_educ;do y_l=1,G_types
        if (df_l<G_df+1) then
            joint_pr_model(df_l,e_l,y_l)=exp(1.0d0/sigma*(av_V_ini(df_l,e_l,y_l)-cost_ey(df_l,e_l,y_l)))/sum(exp(1.0d0/sigma*(av_V_ini(df_l,:,:)-cost_ey(df_l,:,:))))
            write(9,'(A4,I4,A4,I4,A4,I4,F20.5,F20.5,F20.5,F20.5)'),'df_l',df_l,'e_l',e_l,'y_l',y_l,cost_ey(df_l,e_l,y_l),av_V_ini(df_l,e_l,y_l),joint_pr(df_l,e_l,y_l),joint_pr_model(df_l,e_l,y_l)
        else
            write(9,'(A4,I4,A4,I4,A4,I4,F20.5,F20.5,F20.5,F20.5)'),'df_l',df_l,'e_l',e_l,'y_l',y_l,sum(cost_ey(:,e_l,y_l)*pr_beta_un),sum(av_V_ini(:,e_l,y_l)*pr_beta_un),sum(joint_pr(:,e_l,y_l)*pr_beta_un),sum(joint_pr_model(:,e_l,y_l)*pr_beta_un)
        end if
    end do;end do;end do
    close(9)

            

    
    obj_function_costs=sum((joint_pr_model-joint_pr)**2.0d0)
    
    print*,obj_function_costs

    
    
end function
    