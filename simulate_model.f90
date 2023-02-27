subroutine simulate_model(a_policy,VSL,V_ini,asset_distribution,av_VSL,av_V_ini)
    use nrtype;use state_space_dim;use var_first_step;use preference_p
    implicit none
    real(DP),dimension(G_nkk,T,G_nzz,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::a_policy
    real(DP),dimension(G_nkk,G_nzz,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::VSL
    real(DP),dimension(G_nkk,G_nzz,G_h,G_PI,G_educ,G_types,G_DF),intent(in)::V_ini
    integer,parameter::indv_sim=25000
    real(DP),dimension(indv_sim,T)::panel_assets,med_exp
    real(DP),dimension(indv_sim)::VSL_i,V_i
    real(DP),dimension(T)::unemployed
    double precision::cash_on_hand,alpha
    integer::h,xi,xi2,t_l,ind,i_l,y,e,h2,pi_l,loc_coh,df,ts_l,pi_l2
    character::continue_k
    integer,dimension(T)::counter_surv
    real(DP),dimension(T,G_educ,G_types,4),intent(out)::asset_distribution
    real(DP),dimension(G_educ,G_types),intent(out)::av_VSL,av_V_ini
    integer(8),dimension(1)::seed=321
    double precision,dimension(9,T)::u
    real(DP)::savings
    
    !Call seed number
    call random_seed(PUT=seed)
    
    do e=1,G_educ

    do y=1,G_types
    counter_surv=0.0d0
    panel_assets=-9.0d0
    med_exp=-9.0d0
    unemployed=0.0d0
    VSL_i=-9.0d0
    V_i=-9.0d0
    
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
            !Sample persisent medical shock from uncond distribution
            
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
            cash_on_hand=income_grid(y,e,t_l,h,pi_l,ts_l)-m_grid(e,t_l,h,xi)
        else
            !Sample persisent medical shock from cond distribution
            ind=1
            xi2=-9
            do while (xi2==-9)
                if (ind>G_nzz) then
                    print*,'error simulating initial persistent shock '
                    read*,continue_k
                end if
                if (u(4,t_l)<sum(Pi_m(xi,1:ind))) then
                    xi2=ind
                    xi=xi2
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
            exit
        else
            counter_surv(t_l)=counter_surv(t_l)+1
            !Simulate savings decisions
            loc_coh=max(min(int(cash_on_hand/step)+1,G_nkk),1) !
            if (loc_coh<G_nkk .and. cash_on_hand>0.0d0) then
                alpha=1-(cash_on_hand-a_grid(loc_coh))/step
                savings=alpha*a_policy(loc_coh,t_l,xi,h,pi_l,e,y,df)+(1.0d0-alpha)*a_policy(loc_coh+1,t_l,xi,h,pi_l,e,y,df) 
            else
                savings=a_policy(loc_coh,t_l,xi,h,pi_l,e,y,df)
            end if
            panel_assets(counter_surv(t_l),t_l)=savings
            
            med_exp(counter_surv(t_l),t_l)=m_grid(e,t_l,h,xi)
            if (pi_l==G_PI) then
                unemployed(t_l)=unemployed(t_l)+1.0d0
            end if
            if (t_l==T_svl) then
                VSL_i(counter_surv(t_l))=VSL(loc_coh,xi,h,pi_l,e,y,df)
            end if
            if (t_l==1) then
                V_i(counter_surv(t_l))=V_ini(loc_coh,xi,h,pi_l,e,y,df)
            end if
        end if

    end do;
    end do
    
    do t_l=1,T
        if (counter_surv(t_l)>50) then
            call sort(panel_assets(1:counter_surv(t_l),t_l),counter_surv(t_l))
            asset_distribution(t_l,e,y,1)=panel_assets(counter_surv(t_l)/4,t_l)
            asset_distribution(t_l,e,y,2)=panel_assets(counter_surv(t_l)/2,t_l)
            asset_distribution(t_l,e,y,3)=panel_assets(counter_surv(t_l)/4*3,t_l)
            asset_distribution(t_l,e,y,4)=sum(panel_assets(1:counter_surv(t_l),t_l))/dble(counter_surv(t_l))
        else
            asset_distribution(t_l,e,y,:)=-9.0d0
        end if
        if (t_l==1) then
            av_V_ini(e,y)=sum(V_i(1:counter_surv(t_l)))/dble(counter_surv(t_l))
            print*,'e=',e,';y=',y,';V=',av_V_ini(e,y)
        end if 
        if (t_l==T_svl) then
            av_VSL(e,y)=sum(VSL_i(1:counter_surv(t_l)))/dble(counter_surv(t_l))
            print*,'e=',e,';y=',y,';VSL=',av_VSL(e,y)
        end if  
    end do
    
    
    end do;end do

   
    !panel_assets(:,20)
end subroutine