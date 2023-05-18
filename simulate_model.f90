subroutine simulate_model(e,y,a_policy,VSL,asset_distribution,av_VSL,av_V_ini,lambda_ref)
    use nrtype;use state_space_dim;use var_first_step;use preference_p
    implicit none
    integer,intent(in)::e,y
    real(DP),dimension(G_nkk,T,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::a_policy
    real(DP),dimension(G_nkk,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::VSL
    real(DP),optional::lambda_ref
    integer,parameter::indv_sim=20000
    real(DP),dimension(indv_sim,T)::panel_assets
    real(DP),dimension(indv_sim,T)::med_exp
    real(DP),dimension(G_DF,indv_sim)::VSL_i,V_i
    real(DP),dimension(T)::unemployed
    real(DP)::cash_on_hand,alpha,lambda_c
    real(DP),dimension(G_df)::av_c
    integer::h,xi,xi2,t_l,ind,i_l,h2,pi_l,loc_coh,df,ts_l,pi_l2
    character::continue_k
    integer,dimension(T)::counter_surv
    integer,dimension(T,G_DF)::counter_surv_df
    real(DP),dimension(T,4),intent(out)::asset_distribution
    real(DP),dimension(G_DF),intent(out)::av_VSL,av_V_ini
    integer(8),dimension(1)::seed=321
    double precision,dimension(9,T)::u
    real(DP)::savings
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
    END INTERFACE
    
    lambda_c=0.0d0
    if(present(lambda_ref))lambda_c=lambda_ref
    
    !Call seed number
    call random_seed(PUT=seed)
    
    counter_surv=0.0d0
    counter_surv_df=0.0d0
    panel_assets=-9.0d0
    med_exp=-9.0d0
    unemployed=0.0d0
    VSL_i=0.0d0
    V_i=0.0d0
    av_c=0.0d0

    
    do i_l=1,indv_sim
        savings=0.0d0
        call RANDOM_NUMBER(u)
        if (u(8,1)<pr_betas(y,e)) then
            df=1
        else
            df=2
        end if
        
        do t_l=1,T-1
        if (t_l==1) then
            !Simulate intial health
            ind=1
            h=-9
            do while (h==-9)
                if (u(9,t_l)<sum(fraction_h_ey(1:ind,e,y)))then
                    h=ind
                else
                    ind=ind+1
                end if
            end do
            !Sample  medical shock from uncond distribution
            ind=1
            xi=-9
            do while (xi==-9)
                if (ind>G_nzz) then
                    print*,'error simulating initial persistent shock '
                    read*,continue_k
                end if
                if (u(1,t_l)<sum(pr0_p(1:ind,1))) then
                    xi=ind
                else
                    ind=ind+1
                end if
            end do
            !Sample persisent income shock from uncond distribution
            ind=1
            pi_l=-9
            do while (pi_l==-9)
                if (ind>G_PI) then
                    print*,'error simulating initial persistent shock ',t_l,h,e
                    read*,continue_k
                end if
                if (u(2,t_l)<sum(Pi_p_0(1:ind,t_l,h,e))) then 
                    pi_l=ind
                else
                    ind=ind+1
                end if
            end do
            !Sample transitory income shock from uncond distribution
            ind=1
            ts_l=-9
            do while (ts_l==-9)
                if (ind>G_PI) then
                    print*,'error simulating initial persistent shock '
                    read*,continue_k
                end if
                if (u(3,t_l)<sum(Pi_t(1:ind,1))) then
                    ts_l=ind
                else
                    ind=ind+1
                end if
            end do
            cash_on_hand=(1+r)*savings+income_grid(y,e,t_l,h,pi_l,ts_l)-m_grid(e,t_l,h,xi)
        else
           !Sample  medical shock from uncond distribution
            ind=1
            xi=-9
            do while (xi==-9)
                if (ind>G_nzz) then
                    print*,'error simulating initial persistent shock '
                    read*,continue_k
                end if
                if (u(4,t_l)<sum(pr0_p(1:ind,1))) then
                    xi=ind
                else
                    ind=ind+1
                end if
            end do
            
            !Sample persisent income shock from cond distribution
            if (t_l<=T_R) then
                ind=1
                pi_l2=-9
                do while (pi_l2==-9)
                    if (ind>G_PI) then
                        print*,'error simulating initial persistent shock '
                        read*,continue_k
                    end if
                    if (u(5,t_l)<sum(Pi_p(pi_l,1:ind,t_l,h,e))) then
                        pi_l2=ind
                        pi_l=pi_l2
                    else
                        ind=ind+1
                    end if
                end do
            
                !Sample transitory income shock from uncond distribution

                ind=1
                ts_l=-9
                do while (ts_l==-9)
                    if (ind>G_PI) then
                        print*,'error simulating initial persistent shock '
                        read*,continue_k
                    end if
                    if (u(6,t_l)<sum(Pi_t(1:ind,1))) then
                        ts_l=ind
                    else
                        ind=ind+1
                    end if
                end do
            end if
            
            !Sample health status
            ind=1
            h2=-9
            do while (h2==-9)
                if (ind>G_h+1) then
                    print*,'error simulating initial persistent shock '
                    read*,continue_k
                end if
                if (u(7,t_l)<sum(H_sm(h,1:ind,t_l,y,e))) then !H_sm(h,:,t_l,y,e)
                    h2=ind
                    h=h2
                else
                    ind=ind+1
                end if
            end do
            
            if (t_l>=T_R) then
                cash_on_hand=(1+r)*savings+PI_grid(pi_l,e,y)-m_grid(e,t_l,min(h,G_h),xi)
            else
                cash_on_hand=(1+r)*savings+income_grid(y,e,t_l,min(h,G_h),pi_l,ts_l)-m_grid(e,t_l,min(h,G_h),xi)
            end if
        end if
        
        !If dead go to next individual
        if (h==G_h+1) then
            V_i(df,i_l)=V_i(df,i_l)+betas(df)**t_l*beq_fct(max(cash_on_hand,0.0d0)*(1.0d0-lambda_c))!
            if (isnan(V_i(df,i_l))) then
                print*,''
            end if
            exit
        else           
            counter_surv(t_l)=counter_surv(t_l)+1
            counter_surv_df(t_l,df)=counter_surv_df(t_l,df)+1
            !Simulate savings decisions
            loc_coh=max(min(int(cash_on_hand/step)+1,G_nkk),1) !
            if (loc_coh<G_nkk .and. cash_on_hand>0.0d0) then
                alpha=1-(cash_on_hand-a_grid(loc_coh))/step
                savings=alpha*a_policy(loc_coh,t_l,h,pi_l,e,y,df)+(1.0d0-alpha)*a_policy(loc_coh+1,t_l,h,pi_l,e,y,df) 
            else
                savings=a_policy(loc_coh,t_l,h,pi_l,e,y,df)
            end if
            panel_assets(counter_surv(t_l),t_l)=savings
            av_c(df)=dble(sum(counter_surv_df(:,df))-1.0d0)/dble(sum(counter_surv_df(:,df)))*av_c(df) + &
                 1.0d0/dble(sum(counter_surv_df(:,df)))*max((cash_on_hand-savings),c_floor)
            V_i(df,i_l)=V_i(df,i_l)+betas(df)**(t_l-1)*(1.0d0-betas(df))*u_fct(max((cash_on_hand-savings),c_floor)*(1.0d0-lambda_c),n_bar(t_l),h,y)
            if (isnan(V_i(df,i_l))) then
                print*,''
            end if
            med_exp(counter_surv(t_l),t_l)=m_grid(e,t_l,h,xi)
            if (pi_l==G_PI) then
                unemployed(t_l)=unemployed(t_l)+1.0d0
            end if
            if (t_l==T_svl) then
                VSL_i(df,counter_surv_df(t_l,df))=VSL(loc_coh,h,pi_l,e,y,df)
                !if (VSL(loc_coh,h,pi_l,e,y,df)<0.0d0) then
                !    print*,'negative v',VSL(loc_coh,h,pi_l,e,y,df),loc_coh,h,pi_l,e,y,df
                !end if
            end if

        end if

    end do;
    end do
    
    do t_l=1,T
            if (counter_surv(t_l)>50) then
                call sort(panel_assets(1:counter_surv(t_l),t_l),counter_surv(t_l))
                asset_distribution(t_l,1)=panel_assets(counter_surv(t_l)/4,t_l)
                asset_distribution(t_l,2)=panel_assets(counter_surv(t_l)/2,t_l)
                asset_distribution(t_l,3)=panel_assets(counter_surv(t_l)/4*3,t_l)
                asset_distribution(t_l,4)=sum(panel_assets(1:counter_surv(t_l),t_l))/dble(counter_surv(t_l))
            else
                asset_distribution(t_l,:)=-9.0d0
            end if

        if (t_l==1) then
            do df=1,G_DF
                if (counter_surv_df(t_l,df)>0) then
                    av_V_ini(df)=sum(V_i(df,:))/dble(counter_surv_df(1,df))!counter_surv_df(1,:)
                else
                    av_V_ini(df)=-9.0d0
                end if
                !print '(A3,I2,A3,I2,A3,I2,A3,F7.4,A7,F7.2)','e=',e,';y=',y,';df=',df,';V=',av_V_ini(df),';av c=',av_c(df)
            end do
        end if 
        if (t_l==T_svl) then
            do df=1,G_DF
                if (counter_surv_df(t_l,df)>0) then
                    av_VSL(df)=sum(VSL_i(df,1:counter_surv_df(t_l,df)))/dble(counter_surv_df(t_l,df))
                else
                    av_VSL(df)=-9.0d0
                end if
               ! print*,'e=',e,';y=',y,';VSL=',av_VSL(df)
            end do
        end if 
        !print*,t_l,sum(med_exp(1:counter_surv(t_l),t_l))/counter_surv(t_l), sqrt(sum((med_exp(1:counter_surv(t_l),t_l)-sum(med_exp(1:counter_surv(t_l),t_l))/counter_surv(t_l))**2.0d0)/counter_surv(t_l))
    end do
    


   
    !panel_assets(:,20)
end subroutine