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
    real(DP)::ftol=1d-15,u
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
        read(9,*) c_floor,beq_cur,beq_mu,obj_function
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
    
    b_bar_min=-10.0d0
    b_bar_max=50.0d0
    it=0
    
    VSL_data=1000.0d0

    
    !do while (abs(av_VSL(1)-VSL_data)>100.0d0 .and. it<10) 
    !    b_bar=(b_bar_max+b_bar_min)/2.0d0
    !    !Solve model
    !    call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini,p50_delta)
    !    if (av_VSL(1)>VSL_data) then
    !        b_bar_max=b_bar(1)
    !    else
    !        b_bar_min=b_bar(1)
    !    end if
    !    it=it+1
    !    print*,it,b_bar(1),av_VSL(1)   
    !end do
 
    open(unit=9,file='v_ini.txt')
        read(9,*) av_V_ini
    close(9)
    
    open(unit=9,file='b_bar_cost.txt')
        read(9,*)  b_bar
    close(9)

    
    do df_l=1,G_DF;do e_l=1,G_educ;do y_l=1,types
        if (df_l==1) then
            joint_pr(df_l,e_l,y_l)=fraction_types(e_l,y_l)*fraction_e(e_l)*pr_betas(y_l,e_l)
        else
            joint_pr(df_l,e_l,y_l)=fraction_types(e_l,y_l)*fraction_e(e_l)*(1.0d0-pr_betas(y_l,e_l))
        end if
    end do;end do;end do


    !Estimate costs:
    !!!!!!!!!!!!!!!!    
    
    p_g(1,1:PAR_2-2)=(/50.7d0,50.9d0,50.7d0/)
    p_g(1,PAR_2-1)=log(20.0d0)
    p_g(1,PAR_2)=log(20.2d0)
    it=0

1    do p_l=1,PAR_2+1
        if (p_l>1) then
            p_g(p_l,:)=p_g(1,:)
            call random_number(u)
            p_g(p_l,p_l-1)=p_g(1,p_l-1)*0.95d0
        end if
        y(p_l)=obj_function_costs(p_g(p_l,:))

    end do

    print*,y
    pause
    call amoeba(p_g,y,ftol,obj_function_costs,iter) 
    it=it+1
    print('(A20,<PAR_2>F10.5)'),'parameters',p_g(1,:)
    print*,'obj fct',y(1)
    print*,'*************************************************'
    if (it<10) then
        go to 1
    end if
    
     open(unit=9,file='b_bar_cost.txt')
        write(9,*) b_bar,cost_ey,exp(p_g(1,PAR_2))
    close(9)   
   
    

    end subroutine  
    
    
    
function obj_function_costs(parameters)
    use nrtype; use preference_p; use global_var;use second_step
    implicit none
    real(DP),dimension(:),intent(in)::parameters
    real(DP),dimension(G_df,G_educ,G_types)::joint_pr_model,rand_ey,rand_ndo,U
    integer,dimension(3)::maxloc_U
    real(DP)::obj_function_costs
    real(DP)::sigma,sigma_neqdo
    integer:: e_l,y_l,df_l,i_l
    integer(8),dimension(1)::seed=321
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    

    !print*,parameters
    cost_y(G_types)=0.0d0
    cost_y(1:G_types-1)=parameters(1:G_types-1)
    cost_e(1)=0.0d0
    cost_e(2:3)=parameters(G_types-1+1:PAR_2-2)
    
    sigma_neqdo=exp(parameters(PAR_2-1))
    sigma=exp(parameters(PAR_2))
    !print*,'cost_y',cost_y
    !print*,'cost_e',cost_e
    
     do df_l=1,G_df;do e_l=1,G_educ;do y_l=1,G_types
        cost_ey(df_l,e_l,y_l)=cost_e(e_l)+cost_y(y_l)
     end do;end do; end do
    
     call random_seed(PUT=seed)
     
    joint_pr_model=0.0d0
    do i_l=1,100000
        rand_ndo=0.0d0
        rand_ndo(:,2:G_educ,:)=c4_normal_01( )*sqrt(sigma_neqdo)
        do df_l=1,G_df;do e_l=1,G_educ;do y_l=1,G_types
            rand_ey(df_l,e_l,y_l)=c4_normal_01( )*sqrt(sigma)
        end do;end do;end do
        U=av_V_ini-cost_ey+rand_ndo+rand_ey
        maxloc_U=maxloc(U)
        joint_pr_model(maxloc_U(1),maxloc_U(2),maxloc_U(3))=joint_pr_model(maxloc_U(1),maxloc_U(2),maxloc_U(3))+1.0d0
    end do
        
    joint_pr_model=joint_pr_model/sum(joint_pr_model)    

    
   
    
    open(unit=9,file='costs_e_y.txt')
    do df_l=1,G_df+1;do e_l=1,G_educ;do y_l=1,G_types
        if (df_l<G_df+1) then
            write(9,'(A4,I4,A4,I4,A4,I4,F20.5,F20.5,F20.5,F20.5)'),'df_l',df_l,'e_l',e_l,'y_l',y_l,cost_ey(df_l,e_l,y_l),av_V_ini(df_l,e_l,y_l),joint_pr(df_l,e_l,y_l),joint_pr_model(df_l,e_l,y_l)
        else
            write(9,'(A4,I4,A4,I4,A4,I4,F20.5,F20.5,F20.5,F20.5)'),'df_l',df_l,'e_l',e_l,'y_l',y_l,sum(cost_ey(:,e_l,y_l)*pr_beta_un),sum(av_V_ini(:,e_l,y_l)*pr_beta_un),sum(joint_pr(:,e_l,y_l)*pr_beta_un),sum(joint_pr_model(:,e_l,y_l)*pr_beta_un)
        end if
    end do;end do;end do
    close(9)
       

    
    obj_function_costs=sum((joint_pr_model-joint_pr)**2.0d0)
    
    if (isnan(obj_function_costs))then
        print*,''
    end if
    
    
    !print*,obj_function_costs

    !pause
    
end function
    
    
    
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