subroutine solve_model(a_policy,VSL)
    use nrtype;use state_space_dim;use var_first_step; use preference_p
    implicit none
    real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF),intent(out)::a_policy
    real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF),intent(out)::VSL
    real(DP),dimension(G_nkk)::aux_a,aux_c,aux_v
    integer::x_l,t_l,h_l,pi_l,e_l,y_l,df_l,x_l2,a_l2,loc,ind,a_l_max,a_l2_2,loc2,a_l1
    real(DP)::jumps,alpha,savings_at_floor,alpha_left,alpha_right,beta_left,beta_right,x_policy_left,V_x2_left,x_policy_right,V_x2_right,eps,x_star,x_star2,alpha_right2,alpha_left2
    real(DP),dimension(G_nkk,T,G_h+1,G_PI,G_educ,G_types,G_DF)::V
    real(DP),dimension(G_nkk)::x_policy,dV_x2,V_x2,x_policy_s,a_grid_s,x_policy_new,u_md
    integer,dimension(G_nkk)::discard_points,index_sort
    INTERFACE
        FUNCTION beq_fct (a)
            use nrtype
            implicit none
             REAL(DP):: beq_fct
             REAL(DP), INTENT(IN) :: a
        END FUNCTION beq_fct
        FUNCTION u_fct(a,n,h,y)
            use nrtype
            implicit none
             REAL(DP):: u_fct
             REAL(DP), INTENT(IN) :: a,n
             integer,intent(in):: h,y
        END FUNCTION u_fct
        FUNCTION ExpConVal (a,FV,t_l,h_l,pi_l,e_l,y_l)
            use nrtype; use state_space_dim
            implicit none
             REAL(DP):: ExpConVal
             REAL(DP), INTENT(IN) :: a
             REAL(DP),dimension(G_nkk,G_h+1,G_PI), INTENT(IN) ::FV
             integer, INTENT(IN) :: t_l,h_l,pi_l,e_l,y_l
        END FUNCTION ExpConVal
        FUNCTION ExpMarginalUtility (a,FSavings,t_l,h_l,pi_l,e_l,y_l)
            use nrtype; use preference_p;use var_first_step
            implicit none
                REAL(DP):: ExpMarginalUtility
                REAL(DP), INTENT(IN) :: a
                REAL(DP),dimension(G_nkk,G_h,G_PI), INTENT(IN) ::FSavings
                integer, INTENT(IN) :: t_l,h_l,pi_l,e_l,y_l
        end function ExpMarginalUtility 
        FUNCTION value_of_stat_life (a_today,a_prime,FV,t_l,h_l,pi_l,e_l,y_l,df_l)
        use nrtype; use preference_p;use var_first_step
        implicit none
            REAL(DP):: value_of_stat_life
            REAL(DP), INTENT(IN) :: a_today,a_prime
            REAL(DP),dimension(G_nkk,G_h+1,G_PI), INTENT(IN) ::FV
            integer, INTENT(IN) :: t_l,h_l,pi_l,e_l,y_l,df_l
        end function value_of_stat_life
        FUNCTION marginal_u_fct (a,n,h)
        use nrtype; use preference_p
        implicit none
            REAL(DP):: marginal_u_fct
            REAL(DP), INTENT(IN) :: a,n
            integer,intent(in):: h
        END FUNCTION marginal_u_fct
        
    end interface
      
    
    a_policy=-9.0d0
    V=-9.0d0
    
    print*,'Compute value of bequest'
    do x_l=1,G_nkk
        V(x_l,:,G_h+1,:,:,:,:)=beq_fct((1.0d0+r)*a_grid(x_l))
        a_policy(:,T,:,:,:,:,:)=0.0d0
    end do
    
    print*,'Started solving model'
    
    !$omp parallel default(shared) private(aux_a,aux_c,aux_v,x_l,t_l,h_l,pi_l,e_l,y_l,df_l,x_l2,a_l2,loc,loc2,ind,jumps,alpha,x_policy,u_md,dV_x2,V_x2,a_l_max,discard_points,x_star2,a_l1,alpha_right2,alpha_left2,x_policy_new,x_policy_s,a_grid_s,index_sort)
    !$omp  do collapse(2)
    do e_l=1,G_educ;do y_l=1,G_types;do df_l=1,G_DF
    do t_l=T-1,1,-1
        do h_l=1,G_h; do pi_l=1,G_PI
        
        !!DC-EGM 
        !!!!!!!!
        !Value of being in Medicaid
        !If no beq consume everything in last period
        if (beq_mu==0.0d0 .and. t_l==T-1) then
            do x_l=1,G_nkk
                if (a_grid(x_l)<=c_floor)then
                    V(x_l,t_l,h_l,pi_l,:,:,:)=u_fct(c_floor,n_bar(t_l),h_l,y_l)
                else
                    V(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)=u_fct(a_grid(x_l),n_bar(t_l),h_l,y_l)
                    a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)=0.0d0
                end if
            end do 
        else
2           jumps=0
            x_policy=-9.0d0
            V_x2=-9.0d0
            dV_x2=-9.0d0
            discard_points=0
            !Std EGM
            a_l1=1
            do a_l2=1,G_nkk
                dV_x2(a_l2)=betas(df_l)*(1.0d0+r)*ExpMarginalUtility(a_grid(a_l2),a_policy(:,t_l+1,:,:,e_l,y_l,df_l),t_l,h_l,pi_l,e_l,y_l)
                x_policy(a_l2)=(dV_x2(a_l2)*n_bar(t_l))**(-1.0d0/RRA)*n_bar(t_l)+a_grid(a_l2)
                if (isnan(x_policy(a_l2))) then
                    print*,'isnan solve model'
                else
                    V_x2(a_l2)=u_fct(x_policy(a_l2)-a_grid(a_l2),n_bar(t_l),h_l,y_l)+betas(df_l)*ExpConVal(a_grid(a_l2),V(:,t_l+1,:,:,e_l,y_l,df_l),t_l,h_l,pi_l,e_l,y_l)!(1.0d0-betas(y_l))*
                end if
                u_md(a_l2)=u_fct(x_policy(a_l2),n_bar(t_l),h_l,y_l)+betas(df_l)*ExpConVal(0.0d0,V(:,t_l+1,:,:,e_l,y_l,df_l),t_l,h_l,pi_l,e_l,y_l) !(1.0d0-betas(y_l))*
                if (u_md(a_l2)>V_x2(a_l2) .or. dV_x2(a_l2)==0.0d0) then
                    discard_points(a_l2)=1
                    a_l1=a_l2+1
                    x_policy(1)=x_policy(a_l2)
                end if
                if (a_l2>1) then; if (x_policy(a_l2)<x_policy(a_l2-1) .and. discard_points(a_l2)==0 ) then
                    jumps=jumps+1
                end if; end if
            end do
            x_policy_new=x_policy
            
            
            !Discard non-optimal points
            if (jumps>0) then
                a_l2=a_l1+1  
                do while (a_l2<G_nkk+1) 
                    if (x_policy(a_l2)<x_policy(a_l1)) then
                        !Construct firt-order taylor expansion around the left and the right of the value function
                        alpha_left2=marginal_u_fct(x_policy(a_l1)-a_grid(a_l1),n_bar(t_l),h_l) 
                        alpha_right2=marginal_u_fct(x_policy(a_l2)-a_grid(a_l2),n_bar(t_l),h_l) 

                        x_star2=(V_x2(a_l1)-alpha_left2*x_policy(a_l1)-(V_x2(a_l2)-alpha_right2*x_policy(a_l2)))/(alpha_right2-alpha_left2)  
                        !Discard point to the right of x_star
                        ind=0
                        !do while (ind==0)
                            if (x_policy(a_l1)>x_star2) then
                                discard_points(a_l1)=1
                                x_policy_new(a_l1)=1.0d0/0.0d0
                                do while (discard_points(a_l1)==1)
                                    a_l1=a_l1-1
                                    if (a_l1==0) then
                                        exit                                       
                                    end if
                                end do
                                if (a_l1==0) then
                                    ind=1
                                end if
                            else
                                ind=1
                            end if
                        !end do
                        
                        !Discard points to the left of x_star
                        if (x_policy(a_l2)<x_star2) then
                            discard_points(a_l2)=1
                            x_policy_new(a_l2)=1.0d0/0.0d0
                            a_l2=a_l2+1
                        end if
                        
                        if (a_l1==0) then
                            a_l1=a_l2
                            a_l2=a_l2+1
                        end if
                    else
                        a_l1=a_l1+1
                            do while (discard_points(a_l1)==1 .and. a_l1<G_nkk)
                                a_l1=a_l1+1
                            end do
                        a_l2=a_l2+1
                    end if               
                end do
                !Sort the endogenous grid
                ind=0
                do a_l2=1,G_nkk
                    if (discard_points(a_l2)==0) then
                        ind=ind+1
                        x_policy_s(ind)=x_policy_new(a_l2)
                        a_grid_s(ind)=a_grid(a_l2)
                    end if
                end do
            else
                if (sum(discard_points)>1) then
                    ind=0
                    do a_l2=1,G_nkk
                        if (discard_points(a_l2)==0) then
                            ind=ind+1
                            x_policy_s(ind)=x_policy_new(a_l2)
                            a_grid_s(ind)=a_grid(a_l2)
                        end if
                    end do               
                else
                    x_policy_s=x_policy
                    a_grid_s=a_grid
                end if
            end if
            
            if (sum(discard_points)>150)then
                go to 2
                print*,'suspicious'
            end if
            

            !Recover savings policy function in the original grid
1           loc=1
            loc2=2
            do x_l=1,G_nkk
                if (a_grid(x_l)<=x_policy_s(1)) then
                    a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)=0.0d0
                    V(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)=u_fct(max(a_grid(x_l),c_floor),n_bar(t_l),h_l,y_l)+betas(df_l)*ExpConVal(0.0d0,V(:,t_l+1,:,:,e_l,y_l,df_l),t_l,h_l,pi_l,e_l,y_l) !(1.0d0-betas(y_l))*
                else
                    a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)=-9.0d0 
                    do while (a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)==-9.0d0)
                        if (loc2==G_nkk-sum(discard_points)) then
                            alpha=(a_grid(x_l)-x_policy_s(loc-1))/(x_policy_s(loc)-x_policy_s(loc-1))
                            a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)=(a_grid_s(loc)-a_grid_s(loc-1))*alpha+a_grid_s(loc-1)                           
                        elseif (a_grid(x_l)>=x_policy_s(loc) .and. a_grid(x_l)<=x_policy_s(loc2)) then
                            alpha=(a_grid(x_l)-x_policy_s(loc))/(x_policy_s(loc2)-x_policy_s(loc))
                            a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)=(a_grid_s(loc2)-a_grid_s(loc))*alpha+a_grid_s(loc)
                            V(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)=u_fct(a_grid(x_l)-a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l),n_bar(t_l),h_l,y_l)+betas(df_l)*ExpConVal(a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l),V(:,t_l+1,:,:,e_l,y_l,df_l),t_l,h_l,pi_l,e_l,y_l) !(1.0d0-betas(y_l))*
                        else
                            loc=loc2
                            loc2=loc2+1
                        end if
                    end do
                end if
                aux_a(x_l)=a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)
                aux_v(x_l)=v(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)

                if (t_l==T_svl) then
                    VSL(x_l,h_l,pi_l,e_l,y_l,df_l)=value_of_stat_life(a_grid(x_l),a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l),V(:,t_l+1,:,:,e_l,y_l,df_l),t_l,h_l,pi_l,e_l,y_l,df_l)
                end if
            end do
            
            do x_l=1,G_nkk
                if (a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)>a_grid(x_l) .and. a_grid(x_l)>c_floor) then
                    print*,'error'
                end if
                if (x_l>1) then
                    if (a_policy(x_l-1,t_l,h_l,pi_l,e_l,y_l,df_l)>a_policy(x_l,t_l,h_l,pi_l,e_l,y_l,df_l)) then
                        print*,'error'
                        go to 2
                    end if
                end if
            end do
              
        end if
    
    end do;end do
    end do; 
    end do;end do;end do
    !$omp end do
    !$omp end parallel 
    
    print*,'Finished solving model'

end subroutine
    
FUNCTION beq_fct (a)
        use nrtype; use preference_p
        implicit none
            REAL(DP):: beq_fct
            REAL(DP), INTENT(IN) :: a
            if (RRA/=1.0d0) then
                beq_fct=beq_mu**(-RRA)*(a+beq_cur)**(1.0d0-RRA)/(1.0d0-RRA)
            else
                beq_fct=beq_mu**(-RRA)*log(a+beq_cur)
            end if
END FUNCTION beq_fct
FUNCTION u_fct(c,n,h,y)
        use nrtype; use preference_p
        implicit none
            REAL(DP):: u_fct
            REAL(DP), INTENT(IN) :: c,n
            integer,intent(in):: h,y
            if (RRA/=1.0d0) then
                u_fct=(1.0d0-delta_h*dble(h-1))*(c/n)**(1.0d0-RRA)/(1.0d0-RRA)+b_bar !-1.0d0/(1.0d0-RRA)
            else
                u_fct=(1.0d0-delta_h*dble(h-1))*log(c/n)+b_bar
            end if
            
END FUNCTION u_fct
FUNCTION ExpConVal (a,FV,t_l,h_l,pi_l,e_l,y_l)
            use nrtype; use preference_p;use var_first_step
            implicit none
             REAL(DP):: ExpConVal
             REAL(DP), INTENT(IN) :: a
             REAL(DP),dimension(G_nkk,G_h+1,G_PI), INTENT(IN) ::FV
             integer, INTENT(IN) :: t_l,h_l,pi_l,e_l,y_l
             REAL(DP)::cash_on_hand,alpha,sum
             integer:: xi_l2, h_l2,loc_coh,pi_l2,ts_l2
             step=(a_max-a_min)/dble(G_nkk-1)
             ExpConVal=0.0d0
             sum=0.0d0
             if (t_l>=T_R) then
                 do xi_l2=1,G_nzz; do h_l2=1,G_h+1
                     cash_on_hand=max((1.0d0+r)*a+PI_grid(pi_l,e_l,y_l)-m_grid(e_l,t_l,h_l,xi_l2),c_floor)
                     !locate cash on hand in capital grid
                     loc_coh=max(min(int(cash_on_hand/step)+1,G_nkk),1)
                     if (loc_coh<G_nkk) then
                        alpha=(cash_on_hand-a_grid(loc_coh))/step
                        ExpConVal=ExpConVal+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*((FV(loc_coh+1,h_l2,pi_l)-FV(loc_coh,h_l2,pi_l))*alpha+FV(loc_coh,h_l2,pi_l))
                    else
                        alpha=(cash_on_hand-a_grid(G_nkk-1))/step
                        ExpConVal=ExpConVal+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*((FV(G_nkk,h_l2,pi_l)-FV(G_nkk-1,h_l2,pi_l))*alpha+FV(G_nkk-1,h_l2,pi_l))
                     end if  
                 end do;end do
            else
                do xi_l2=1,G_nzz; do h_l2=1,G_h+1; do pi_l2=1,G_PI; do ts_l2=1,G_PI
                    cash_on_hand=max((1.0d0+r)*a+income_grid(y_l,e_l,t_l,min(h_l2,G_h),pi_l2,ts_l2)-m_grid(e_l,t_l,min(h_l2,G_h),xi_l2),c_floor) 
                     !locate cash on hand in capital grid
                     loc_coh=max(min(int(cash_on_hand/step)+1,G_nkk),1)
                     if (loc_coh<G_nkk) then !Pi_p(pi_l,:,t_l,h_l,e_l)
                        alpha=(cash_on_hand-a_grid(loc_coh))/step
                        ExpConVal=ExpConVal+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*((FV(loc_coh+1,h_l2,pi_l2)-FV(loc_coh,h_l2,pi_l2))*alpha+FV(loc_coh,h_l2,pi_l2))
                        sum=sum+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)
                    else
                        alpha=(cash_on_hand-a_grid(G_nkk-1))/step
                        ExpConVal=ExpConVal+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*FV(G_nkk,h_l2,pi_l2) !((FV(G_nkk,h_l2,pi_l2)-FV(G_nkk-1,h_l2,pi_l2))*alpha+FV(G_nkk-1,h_l2,pi_l2)) 
                    end if  
                 end do;end do;end do;end do
            end if
    
        
END FUNCTION ExpConVal
    
FUNCTION marginal_u_fct (c,n,h)
    use nrtype; use preference_p
    implicit none
        REAL(DP):: marginal_u_fct
        REAL(DP), INTENT(IN) :: c,n
        integer,intent(in):: h
        marginal_u_fct=(1.0d0-delta_h*dble(h-1))*(1.0d0/n)*(c/n)**(-RRA)
END FUNCTION marginal_u_fct
    
FUNCTION ExpMarginalUtility (a,FSavings,t_l,h_l,pi_l,e_l,y_l)
            use nrtype; use preference_p;use var_first_step
            implicit none
             REAL(DP):: ExpMarginalUtility,a2
             REAL(DP), INTENT(IN) :: a
             REAL(DP),dimension(G_nkk,G_h,G_PI), INTENT(IN) ::FSavings
             integer, INTENT(IN) :: t_l,h_l,pi_l,e_l,y_l
             REAL(DP)::cash_on_hand,alpha,const,Pr_zero_mu,non_zero_mu
             integer:: xi_l2, h_l2,loc_coh,pi_l2,ts_l2
             interface
                FUNCTION marginal_beq_fct (a)
                    use nrtype; use preference_p
                    implicit none
                        REAL(DP):: marginal_beq_fct
                        REAL(DP), INTENT(IN) :: a
                END FUNCTION marginal_beq_fct
                FUNCTION marginal_u_fct (a,n,h)
                    use nrtype; use preference_p
                    implicit none
                        REAL(DP):: marginal_u_fct
                        REAL(DP), INTENT(IN) :: a,n
                        integer,intent(in):: h
                END FUNCTION marginal_u_fct
             end interface

             step=(a_max-a_min)/dble(G_nkk-1)
             ExpMarginalUtility=0.0d0
             Pr_zero_mu=0.0d0
             non_zero_mu=0.0d0

             if (t_l>=T_R) then
                    do xi_l2=1,G_nzz; do h_l2=1,G_h+1
                        cash_on_hand=max((1.0d0+r)*a+PI_grid(pi_l,e_l,y_l)-m_grid(e_l,t_l,min(h_l2,G_h),xi_l2),c_floor)
                        if (h_l2<G_h+1) then
                            !locate cash on hand in capital grid                         
                            loc_coh=max(min(int(cash_on_hand/step)+1,G_nkk),1)
                            if (cash_on_hand==c_floor) then
                                ExpMarginalUtility=ExpMarginalUtility+0.0d0 
                                Pr_zero_mu=Pr_zero_mu+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)
                                if (ExpMarginalUtility<0.0d0) then
                                    print*,cash_on_hand
                                    pause
                                end if
                            elseif (loc_coh<G_nkk) then
                                if (loc_coh==1) then
                                    a2=FSavings(loc_coh,h_l2,pi_l)
                                    ExpMarginalUtility=ExpMarginalUtility+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2)
                                else
                                    alpha=(cash_on_hand-a_grid(loc_coh))/step
                                    a2=(FSavings(loc_coh+1,h_l2,pi_l)-FSavings(loc_coh,h_l2,pi_l))*alpha+FSavings(loc_coh,h_l2,pi_l) 
                                    ExpMarginalUtility=ExpMarginalUtility+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2)
                                end if                                    
                                if (ExpMarginalUtility<0.0d0) then
                                    print*,marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2),cash_on_hand,cash_on_hand-a2
                                    pause
                                end if
                                non_zero_mu=non_zero_mu+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2)
                            else
                                alpha=(FSavings(G_nkk,h_l2,pi_l)-FSavings(G_nkk-1,h_l2,pi_l))/(a_grid(G_nkk)-a_grid(G_nkk-1))
                                const=FSavings(G_nkk-1,h_l2,pi_l)-alpha*a_grid(G_nkk-1)
                                a2=alpha*cash_on_hand+const
                                ExpMarginalUtility=ExpMarginalUtility+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2) 
                                if (ExpMarginalUtility<0.0d0) then
                                    print*,marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2),cash_on_hand,cash_on_hand-a2
                                    pause
                                end if
                                non_zero_mu=non_zero_mu+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2)
                            end if 
                        else
                            cash_on_hand=max((1.0d0+r)*a+PI_grid(pi_l,e_l,y_l)-m_grid(e_l,t_l,min(h_l2,G_h),xi_l2),0.0d0)
                            ExpMarginalUtility=ExpMarginalUtility+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*marginal_beq_fct(cash_on_hand)
                        end if
                    end do;end do
            else
                   do xi_l2=1,G_nzz; do h_l2=1,G_h+1; do pi_l2=1,G_PI; do ts_l2=1,G_PI
                        cash_on_hand=max((1.0d0+r)*a+income_grid(y_l,e_l,t_l,min(h_l2,G_h),pi_l2,ts_l2)-m_grid(e_l,t_l,min(h_l2,G_h),xi_l2),c_floor)
                        if (h_l2<G_h+1) then
                            !locate cash on hand in capital grid                         
                            loc_coh=max(min(int(cash_on_hand/step)+1,G_nkk),1)
                            if (cash_on_hand==c_floor) then
                                ExpMarginalUtility=ExpMarginalUtility+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*0.0d0 
                                Pr_zero_mu=Pr_zero_mu+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)
                                if (ExpMarginalUtility<0.0d0) then
                                    print*,marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2),cash_on_hand,cash_on_hand-a2
                                    pause
                                end if
                            elseif (loc_coh<G_nkk) then
                                if (loc_coh==1) then
                                    a2=FSavings(loc_coh,h_l2,pi_l2)
                                    ExpMarginalUtility=ExpMarginalUtility+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2)
                                else
                                    alpha=(cash_on_hand-a_grid(loc_coh))/step 
                                    a2=(FSavings(loc_coh+1,h_l2,pi_l2)-FSavings(loc_coh,h_l2,pi_l2))*alpha+FSavings(loc_coh,h_l2,pi_l2)
                                    ExpMarginalUtility=ExpMarginalUtility+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2)
                                end if
                                if (ExpMarginalUtility<0.0d0) then
                                    print*,marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2),cash_on_hand,cash_on_hand-a2
                                    pause
                                end if
                                non_zero_mu=non_zero_mu+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2)
                            else
                                alpha=(FSavings(G_nkk,h_l2,pi_l2)-FSavings(G_nkk-1,h_l2,pi_l2))/(a_grid(G_nkk)-a_grid(G_nkk-1))
                                const=FSavings(G_nkk-1,h_l2,pi_l2)-alpha*a_grid(G_nkk-1)
                                a2=alpha*cash_on_hand+const
                                ExpMarginalUtility=ExpMarginalUtility+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2)
                                if (ExpMarginalUtility<0.0d0) then
                                    print*,marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2),cash_on_hand,cash_on_hand-a2
                                    pause
                                end if
                                non_zero_mu=non_zero_mu+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*marginal_u_fct(cash_on_hand-a2,n_bar(t_l+1),h_l2)
                            end if 
                        else
                            cash_on_hand=max((1.0d0+r)*a+income_grid(y_l,e_l,t_l,min(h_l2,G_h),pi_l2,ts_l2)-m_grid(e_l,t_l,min(h_l2,G_h),xi_l2),0.0d0)
                            ExpMarginalUtility=ExpMarginalUtility+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*marginal_beq_fct(cash_on_hand)
                        end if
                   end do;end do;end do; end do
            end if   
            
END FUNCTION ExpMarginalUtility    
    
FUNCTION marginal_beq_fct (a)
        use nrtype; use preference_p
        implicit none
            REAL(DP):: marginal_beq_fct
            REAL(DP), INTENT(IN) :: a
            marginal_beq_fct=beq_mu**(-RRA)*(a+beq_cur)**(-RRA)
END FUNCTION marginal_beq_fct
    
    
    
FUNCTION value_of_stat_life (a_today,a_prime,FV,t_l,h_l,pi_l,e_l,y_l,df_l)
        use nrtype; use preference_p;use var_first_step
        implicit none
            REAL(DP):: value_of_stat_life,ExpConVal_alive,ExpConVal_dead
            REAL(DP), INTENT(IN) :: a_today,a_prime
            REAL(DP),dimension(G_nkk,G_h+1,G_PI), INTENT(IN) ::FV
            integer, INTENT(IN) :: t_l,h_l,pi_l,e_l,y_l,df_l
            REAL(DP)::cash_on_hand,alpha,sum
            integer:: xi_l2, h_l2,loc_coh,pi_l2,ts_l2
            step=(a_max-a_min)/dble(G_nkk-1)
            value_of_stat_life=0.0d0
            ExpConVal_alive=0.0d0
            ExpConVal_dead=0.0d0
            sum=0.0d0
            if (t_l>=T_R) then
                do xi_l2=1,G_nzz; do h_l2=1,G_h+1
                    cash_on_hand=max((1.0d0+r)*a_prime+PI_grid(pi_l,e_l,y_l)-m_grid(e_l,t_l,h_l,xi_l2),c_floor)
                    !locate cash on hand in capital grid
                    loc_coh=max(min(int(cash_on_hand/step)+1,G_nkk),1)
                    if (loc_coh<G_nkk) then
                        alpha=(cash_on_hand-a_grid(loc_coh))/step
                        if (h_l2<G_h+1) then
                            ExpConVal_alive=ExpConVal_alive+H_sm(h_l,h_l2,t_l,y_l,e_l)/(1.0d0-H_sm(h_l,G_h+1,t_l,y_l,e_l))*pr0_p(xi_l2,1)*((FV(loc_coh+1,h_l2,pi_l)-FV(loc_coh,h_l2,pi_l))*alpha+FV(loc_coh,h_l2,pi_l))
                        else
                            ExpConVal_dead=ExpConVal_dead+pr0_p(xi_l2,1)*((FV(loc_coh+1,h_l2,pi_l)-FV(loc_coh,h_l2,pi_l))*alpha+FV(loc_coh,h_l2,pi_l))
                        end if
                    else
                        alpha=(cash_on_hand-a_grid(G_nkk-1))/step
                        if (h_l2<G_h+1) then
                            ExpConVal_alive=ExpConVal_alive+H_sm(h_l,h_l2,t_l,y_l,e_l)/(1.0d0-H_sm(h_l,G_h+1,t_l,y_l,e_l))*pr0_p(xi_l2,1)*((FV(G_nkk,h_l2,pi_l)-FV(G_nkk-1,h_l2,pi_l))*alpha+FV(G_nkk-1,h_l2,pi_l))
                        else
                            ExpConVal_dead=ExpConVal_dead+pr0_p(xi_l2,1)*((FV(G_nkk,h_l2,pi_l)-FV(G_nkk-1,h_l2,pi_l))*alpha+FV(G_nkk-1,h_l2,pi_l))
                        end if
                    end if  
                end do;end do
        else
            do xi_l2=1,G_nzz; do h_l2=1,G_h+1; do pi_l2=1,G_PI; do ts_l2=1,G_PI
                cash_on_hand=max((1.0d0+r)*a_prime+income_grid(y_l,e_l,t_l,min(h_l2,G_h),pi_l2,ts_l2)-m_grid(e_l,t_l,min(h_l2,G_h),xi_l2),c_floor) 
                    !locate cash on hand in capital grid
                    loc_coh=max(min(int(cash_on_hand/step)+1,G_nkk),1)
                    if (loc_coh<G_nkk) then 
                        alpha=(cash_on_hand-a_grid(loc_coh))/step
                        if (h_l2<G_h+1) then
                            ExpConVal_alive=ExpConVal_alive+H_sm(h_l,h_l2,t_l,y_l,e_l)/(1.0d0-H_sm(h_l,G_h+1,t_l,y_l,e_l))*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*((FV(loc_coh+1,h_l2,pi_l2)-FV(loc_coh,h_l2,pi_l2))*alpha+FV(loc_coh,h_l2,pi_l2))
                        else
                            ExpConVal_dead=ExpConVal_dead+pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*((FV(loc_coh+1,h_l2,pi_l2)-FV(loc_coh,h_l2,pi_l2))*alpha+FV(loc_coh,h_l2,pi_l2))
                        end if
                        sum=sum+H_sm(h_l,h_l2,t_l,y_l,e_l)*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)
                    else
                        alpha=(cash_on_hand-a_grid(G_nkk-1))/step
                        if (h_l2<G_h+1) then
                            ExpConVal_alive=ExpConVal_alive+H_sm(h_l,h_l2,t_l,y_l,e_l)/(1.0d0-H_sm(h_l,G_h+1,t_l,y_l,e_l))*pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*FV(G_nkk,h_l2,pi_l2)
                        else
                            ExpConVal_dead=ExpConVal_dead+pr0_p(xi_l2,1)*Pi_p(pi_l,pi_l2,t_l,h_l,e_l)*Pi_t(ts_l2,e_l)*FV(G_nkk,h_l2,pi_l2)
                        end if
                    end if  
            end do;end do;end do;end do
        end if
            
        value_of_stat_life=betas(df_l)*(ExpConVal_alive-ExpConVal_dead)/(1.0d0/n_bar(t_l)*(((a_today-a_prime)/n_bar(t_l))**(-RRA)))!betas(y_l)*(ExpConVal_alive-ExpConVal_dead)/((1.0d0-betas(y_l))*1.0d0/n_bar(t_l)*(((a_today-a_prime)/n_bar(t_l))**(-RRA)))
        
END FUNCTION value_of_stat_life