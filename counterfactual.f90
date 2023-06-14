subroutine counterfactuals()
    use nrtype; use preference_p; use global_var; use var_first_step
    implicit none
    real(DP)::obj_function,av_VSL_do,b_bar_max,b_bar_min,VSL_data
    real(DP),dimension(G_df)::k
    real(DP),dimension(T,G_h,G_educ,G_types,4)::asset_distribution
    real(DP),dimension(G_DF,G_educ,G_types)::av_VSL,av_V_ini,joint_pr,cost_ey,cost_zero
    real(DP),dimension(8,generations,types,L_educ)::dist_assets_data
    integer::y_l,t_l,e_l,it,df_l
    real(DP)::G,G_new,lambda_max,lambda_min,G_benchmark
    real(DP),dimension(3,6)::counterfactual_table
    integer,parameter::K2=(G_educ-1)*(G_types-1)+1
    real(DP),dimension(K2,G_educ,G_types)::X_T
    real(DP),dimension(K2,K2)::X_TX
    real(DP),dimension(G_types)::y_d
    real(DP),dimension(G_educ)::e_d
    real(DP),dimension(K2,G_educ*G_types)::X_T_new
    character::end_k
    real(DP),dimension(K2,1)::beta
    real(DP),dimension(G_educ,G_h,T)::p50_delta
    
    !load preferences
    open(unit=9,file='parameter.txt')
        read(9,*) betas,c_floor,beq_cur,beq_mu,pr_betas
    close (9)
    open(unit=9,file='b_bar_costs.txt')
        read(9,*)  b_bar,cost_ey,k
    close(9)
    
    cost_zero=0.0d0
    
    !Benchmark economy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Recover value function from benchmark model
    open(unit=9,file='v_ini.txt')
        read(9,*) av_V_ini
    close(9)
    
    !compute shares of (y,e,df)
    pr_betas_u=0.0d0
    do e_l=1,L_educ;do y_l=1,G_types;do df_l=1,G_DF
        joint_pr(df_l,e_l,y_l)=exp((av_V_ini(df_l,e_l,y_l)-cost_ey(df_l,e_l,y_l))/k(df_l))/sum(exp((av_V_ini(df_l,:,:)-cost_ey(df_l,:,:))/k(df_l))) 
        print*,'df_l',df_l,'e_l',e_l,'y_l',y_l,joint_pr(df_l,e_l,y_l)
        if (df_l==1) then
            pr_betas_u(df_l)=pr_betas_u(df_l)+pr_betas(y_l,e_l)*fraction_types(e_l,y_l)*fraction_e(e_l)
        end if
    end do;end do;end do
    pr_betas_u(2)=1.0d0-pr_betas_u(1) 

    
    !simulate model to compute tax revenue
    call compute_labor_income_revenue(joint_pr,G_benchmark,av_V_ini,cost_ey,counterfactual_table(1,:))

    
    !Counterfactual 1: no complementarities in costs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Recover value function from benchmark model
    open(unit=9,file='v_ini.txt')
        read(9,*) av_V_ini
    close(9)
    
    !compute shares of (y,e,df)
    do df_l=1,G_DF
        do e_l=1,L_educ;do y_l=1,G_types   
            y_d=0.0d0
            e_d=0.0d0
            y_d(y_l)=1.0d0
            e_d(e_l)=1.0d0
            X_t(:,e_l,y_l)=(/1.0d0, y_d(2:G_types), e_d(2:L_educ)/)
        end do;end do
        X_TX=matmul(reshape(X_t,(/K2,G_educ*G_types/)),transpose(reshape(X_t,(/K2,G_educ*G_types/))))
        call inverse(X_TX,X_TX,K2)
        beta=matmul(X_TX,matmul(reshape(X_t,(/K2,G_educ*G_types/)),reshape(cost_ey(df_l,:,:),(/G_educ*G_types,1/))))
        do e_l=1,L_educ;do y_l=1,G_types   
            cost_zero(df_l,e_l,y_l)=sum(X_t(:,e_l,y_l)*beta(:,1))
        end do;end do
    end do
    
    do e_l=1,L_educ;do y_l=1,G_types;do df_l=1,G_DF
        joint_pr(df_l,e_l,y_l)=exp((av_V_ini(df_l,e_l,y_l)-cost_zero(df_l,e_l,y_l))/k(df_l))/sum(exp((av_V_ini(df_l,:,:)-cost_zero(df_l,:,:))/k(df_l)))  
    end do;end do;end do

    !simulate model to compute tax revenue
    call compute_labor_income_revenue(joint_pr,G,av_V_ini,cost_zero,counterfactual_table(2,:))

    
    
    !Counterfactual 2: US -> DNK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tau=0.305d0
    lambda=1.015d0
    lambda_max=1.2d0
    lambda_min=0.8d0
    G_new=0.9d0
    do while (abs(G_new-G_benchmark)>1.0d0)
        
        call load_first_step_p()
        call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini,p50_delta)    

         !compute new shares of (y,e,df)
        do e_l=1,L_educ;do y_l=1,types;do df_l=1,G_DF
            joint_pr(df_l,e_l,y_l)=exp((av_V_ini(df_l,e_l,y_l)-cost_ey(df_l,e_l,y_l))/k(df_l))/sum(exp((av_V_ini(df_l,:,:)-cost_ey(df_l,:,:))/k(df_l)))
            print*,'e_l',e_l,'y_l',y_l,joint_pr(df_l,e_l,y_l)
        end do;end do;end do
    
        call load_first_step_p()
        call compute_labor_income_revenue(joint_pr,G_new,av_V_ini,cost_ey,counterfactual_table(3,:))
        
        
        if (G_new>G_benchmark) then
            lambda_min=lambda
        else
            lambda_max=lambda
        end if
        lambda=(lambda_min+lambda_max)/2.0d0
        print*,lambda,lambda_min,lambda_max
    end do
    
    open(unit=9,file='counterfactual.txt')
    do df_l=1,G_df;do e_l=1,L_educ;do y_l=1,types
        write(9,'(A4,I4,A4,I4,A4,I4,F20.3,F10.3,F10.3)'),'df_l',df_l,'e_l',e_l,'y_l',y_l,cost_ey(df_l,e_l,y_l),av_V_ini(df_l,e_l,y_l),joint_pr(df_l,e_l,y_l)
    end do; end do;end do
    close(9)
    
    print '(A4,A4,<6>A10)','C','N','Col (%)','Pro (%)', 'LE_50','Delta LE','mean(V)','var(V)'
    do e_l=1,3
        print '(A4,I4,<6>F10.3)','c_l',e_l,counterfactual_table(e_l,:)
    end do
    
    
    print*,'pause in counterfactuals'
    pause
    
end subroutine
    
subroutine compute_labor_income_revenue(joint_pr,G,av_V_ini,cost_ey,counterfactuals)
    use nrtype;use state_space_dim;use var_first_step; use preference_p
        implicit none
        real(DP),dimension(G_df,G_educ,G_types),intent(in)::joint_pr,cost_ey
        real(DP),dimension(G_DF,G_educ,G_types),intent(in)::av_V_ini
        real(DP),dimension(6)::counterfactuals
        real(DP),intent(out)::G
        integer,parameter::indv_sim2=100000
        double precision::u
        integer::h,t_l,ind,i_l,y,e,h2,pi_l,df,ts_l,pi_l2,pi_old,counter
        character::continue_k
        real(DP),dimension(indv_sim2*T_R)::panel_tax,panel_income
        integer,dimension(T)::counter_income
        real(DP),dimension(T,G_educ,G_types)::av_income_panel
        real(DP)::aime,pia,coeff
        integer,dimension(G_educ,G_types,G_h,T)::joint_pr_t
        real(DP),dimension(indv_sim2)::V_i
        real(DP),dimension(G_educ)::life_exp,counter_e
        
        counter_income=0
        counter=0
        counter_e=0

        
            panel_tax=-9;panel_income=-9.0d0
            do i_l=1,indv_sim2;aime=-9.0d0
                !Sample discount factor
                call RANDOM_NUMBER(u)
                if (u<pr_betas_u(1)) then
                    df=1
                else
                    df=2
                end if
                !Sample type
                call RANDOM_NUMBER(u)
                ind=1
                e=-9
                do while (e==-9)
                    if (u<sum(joint_pr(df,1:ind,:)))then
                        e=ind
                    else
                        ind=ind+1
                    end if
                end do
                call RANDOM_NUMBER(u)
                ind=1
                y=-9
                do while (y==-9)
                    if (u<sum(joint_pr(df,e,1:ind))/sum(joint_pr(df,e,:)))then
                        y=ind
                    else
                        ind=ind+1
                    end if
                end do

                V_i(i_l)=av_V_ini(1,e,y)!-cost_ey(e,y)
                do t_l=1,T
                    if (t_l==1) then
                        call RANDOM_NUMBER(u)
                        ind=1
                        h=-9
                        do while (h==-9)
                            if (u<sum(fraction_h_ey(1:ind,e,y)))then
                                h=ind
                            else
                                ind=ind+1
                            end if
                        end do
                        !Sample persisent income shock from uncond distribution
                        call RANDOM_NUMBER(u)
                        ind=1
                        pi_l=-9
                        do while (pi_l==-9)
                            if (u<sum(Pi_p_0(1:ind,t_l,h,e))) then 
                                pi_l=ind
                            else
                                ind=ind+1
                            end if
                        end do

                        !Sample transitory income shock from uncond distribution
                        call RANDOM_NUMBER(u)
                        ind=1
                        ts_l=-9
                        do while (ts_l==-9)
                            if (ind>G_PI) then
                                print*,'error simulating initial persistent shock 1'
                                read*,continue_k
                            end if
                            if (u<sum(Pi_t(1:ind,e))) then
                                ts_l=ind
                            else
                                ind=ind+1
                            end if
                        end do

                    else
                        !Sample persisent income shock from cond distribution
                        if (t_l<=T_R) then
                            call RANDOM_NUMBER(u)
                            ind=1
                            pi_old=pi_l
                            pi_l2=-9
                            do while (pi_l2==-9)
                                if (ind>G_PI) then
                                    print*,'error simulating initial persistent shock 2',pi_l,t_l,h,e,Pi_p(pi_l,1:ind,t_l,h,e)
                                    read*,continue_k
                                end if
                                if (u<sum(Pi_p(pi_l,1:ind,t_l,h,e))) then 
                                    pi_l2=ind
                                    pi_l=pi_l2
                                else
                                    ind=ind+1
                                end if
                            end do
            
                            !Sample transitory income shock from uncond distribution
                            call RANDOM_NUMBER(u)
                            ind=1
                            ts_l=-9
                            do while (ts_l==-9)
                                if (ind>G_PI) then
                                    print*,'error simulating initial persistent shock 3'
                                    read*,continue_k
                                end if
                                if (u<sum(Pi_t(1:ind,e))) then
                                    ts_l=ind
                                else
                                    ind=ind+1
                                end if
                            end do
                        end if
            
                        !Sample health status
                        call RANDOM_NUMBER(u)
                        ind=1
                        h2=-9
                        do while (h2==-9)
                            if (ind>G_h+1) then
                                print*,'error simulating initial persistent shock 4'
                                read*,continue_k
                            end if
                            if (u<sum(H_sm(h,1:ind,t_l,y,e))) then !H_sm(h,:,t_l,y,e)
                                h2=ind
                                h=h2
                            else
                                ind=ind+1
                            end if
                        end do
                    end if
                    
            
                    !If dead go to next individual
                    if (h==G_h+1 .and. t_l>=(52-first_age_sm)/2+1 ) then   
                        counter_e(e)=counter_e(e)+1.0d0
                        !print'(I5,F8.2)',t_l,dble((t_l-1)*2+first_age_sm-50)
                        life_exp(e)=1.0d0/counter_e(e)*(min((t_l-1)*2+first_age_sm+2,95))+(counter_e(e)-1.0d0)/counter_e(e)*life_exp(e) !dble((t_l-1)*2+first_age_sm-50)
                        exit
                    elseif (h==G_h+1) then
                        exit
                    end if
                           
                    counter_income(t_l)=counter_income(t_l)+1
                    joint_pr_t(e,y,h,t_l)=joint_pr_t(e,y,h,t_l)+1.0d0
                    if (t_l<=T_R) then
                        counter=counter+1     
                        panel_tax(counter)=income_tax(y,e,t_l,min(h,G_h),pi_l,ts_l)
                        panel_income(counter)=gross_annual_income2(y,e,t_l,min(h,G_h),pi_l,ts_l)
                        if (panel_income(counter)==-9.0d0) then
                            print*,''
                        end if
                    end if
                end do;
            end do

            G=sum(panel_tax)/sum(counter_income)
            
            !Share college
            counterfactuals(1)=sum(joint_pr(1,3,:)*pr_betas_u(1))+ sum(joint_pr(2,3,:)*pr_betas_u(2))
            !Share protective
            counterfactuals(2)=sum(joint_pr(1,:,1)*pr_betas_u(1))+sum(joint_pr(2,:,1)*pr_betas_u(2))
            !Av LE
            counterfactuals(3)=sum((sum(joint_pr(1,:,:),2)*pr_betas_u(1)+sum(joint_pr(2,:,:),2)*pr_betas_u(2))*life_exp)-50.0d0 
            !LE gradient
            counterfactuals(4)=life_exp(3)-life_exp(1)
            !Mean life-time utility
            counterfactuals(5)=sum(V_i)/dble(indv_sim2)
            !Variance of life-time utility
            counterfactuals(6)=sum((V_i-sum(V_i)/dble(indv_sim2))**2.0d0)/dble(indv_sim2)
            
            print*,'G raised',G
            
            call sort(panel_income,counter)
            
            print*,'marginal tax rates'
            do e=1,10
                ind=counter/10*e-counter/10/2
                print*,e,1-lambda*(1-tau)*(panel_income(ind)/av_income)**(-tau)
            end do
            
    end subroutine    
    
    
subroutine inverse(a,c,n)
use nrtype
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real(DP) :: a(n,n), c(n,n)
real(DP):: L(n,n), U(n,n), b(n), d(n), x(n)
real(DP) :: coeff
integer :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0d0
U=0.0d0
b=0.0d0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0d0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0d0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0d0
end do
end subroutine inverse