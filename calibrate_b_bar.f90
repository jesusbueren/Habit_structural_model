subroutine calibrate_b_bar()
    use nrtype; use preference_p; use global_var; use var_first_step
    implicit none
    real(DP)::obj_function,b_bar_max,b_bar_min,VSL_data
    real(DP),dimension(T,G_h,G_educ,G_types,4)::asset_distribution
    real(DP),dimension(G_educ)::av_VSL
    real(DP),dimension(G_DF,G_educ,G_types)::av_V_ini
    real(DP),dimension(G_df,G_educ,G_types)::joint_pr,cost_ey
    real(DP),dimension(G_df)::k
    real(DP),dimension(8,generations,types,L_educ)::dist_assets_data
    real(DP),dimension(G_educ,G_h,T)::p50_delta
    integer::y_l,t_l,e_l,it,df_l
    character::end_k
    
    open(unit=9,file='parameter.txt')
        read(9,*) betas,c_floor,pr_betas,obj_function
    close (9)
 
    
    print*,'----------------------'
    print*,"Parameter"
    print('(A20,<3>F5.2)'),"beta",betas
    print('(A20,F10.2)'),"beq cur",beq_cur
    print('(A20,F10.2)'),"c floor",c_floor
    print('(A20,F10.2)'),"beq mu",beq_mu
    print('(A20,<9>F6.3)'),"pr low beta",pr_betas

    

    if (G_DF==1) then
        pr_betas=1.0d0
    end if
    
    b_bar_min=0.0d0
    b_bar_max=30.0d0
    it=0
    
    VSL_data=2000.0d0
    
    
    
    do while (abs(av_VSL(1)-VSL_data)>100.0d0 .and. it<10) 
        b_bar=(b_bar_max+b_bar_min)/2.0d0
        !Solve model
        call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini,p50_delta)
        if (av_VSL(1)>VSL_data) then
            b_bar_max=b_bar
        else
            b_bar_min=b_bar
        end if
        it=it+1
        print*,it,b_bar,av_VSL(1)
    end do
    
    !Estimate costs:
    !!!!!!!!!!!!!!!!    
    
    !Specification I. 9 cost parameter, var=1.
    open(unit=9,file='costs_e_y.txt')
    do df_l=1,G_DF;
        do e_l=1,G_educ;do y_l=1,types
            if (df_l==1) then
                joint_pr(df_l,e_l,y_l)=fraction_types(e_l,y_l)*fraction_e(e_l)*pr_betas(y_l,e_l)
            else
                joint_pr(df_l,e_l,y_l)=fraction_types(e_l,y_l)*fraction_e(e_l)*(1.0d0-pr_betas(y_l,e_l))
            end if
        end do;end do
        joint_pr(df_l,:,:)=joint_pr(df_l,:,:)/sum(joint_pr(df_l,:,:))
        print*,'cov',sum((log(joint_pr(df_l,:,:))-sum(log(joint_pr(df_l,:,:)))/dble(G_educ*G_types))*(av_V_ini(df_l,:,:)-sum(av_V_ini(df_l,:,:))/dble(G_educ*G_types)))/dble(G_educ*G_types)
        print*,'var x',sum((av_V_ini(df_l,:,:)-sum(av_V_ini(df_l,:,:))/dble(G_educ*G_types))**2.0d0)/dble(G_educ*G_types)
        k(df_l)=sum((log(joint_pr(df_l,:,:))-sum(log(joint_pr(df_l,:,:)))/dble(G_educ*G_types))*(av_V_ini(df_l,:,:)-sum(av_V_ini(df_l,:,:))/dble(G_educ*G_types)))/sum((av_V_ini(df_l,:,:)-sum(av_V_ini(df_l,:,:))/dble(G_educ*G_types))**2.0d0)
        k(df_l)=1.0d0/k(df_l)
        print*,'k=',k(df_l)
        print*,'cost: e x y'
        
        do e_l=1,L_educ;do y_l=1,types
            cost_ey(df_l,e_l,y_l)=av_V_ini(df_l,e_l,y_l)-av_V_ini(df_l,1,3)-k(df_l)*(log(joint_pr(df_l,e_l,y_l))-log(joint_pr(df_l,1,3)))
            print'(A4,I4,A4,I4,A4,I4,A4,F20.5,A5,F20.5)','df_l',df_l,'e_l',e_l,'y_l',y_l,'%',joint_pr(df_l,e_l,y_l),'cost',cost_ey(df_l,e_l,y_l)
            write(9,'(A4,I4,A4,I4,A4,I4,F20.5,F20.5,F20.5,F20.5)'),'df_l',df_l,'e_l',e_l,'y_l',y_l,cost_ey(df_l,e_l,y_l),av_V_ini(df_l,e_l,y_l),joint_pr(df_l,e_l,y_l),exp(1.0d0/k(df_l)*av_V_ini(df_l,e_l,y_l))/sum(exp(1.0d0/k(df_l)*av_V_ini(df_l,:,:)))
        end do;end do
        
    
    end do
close (9)

    open(unit=9,file='b_bar_costs.txt')
        write(9,*) b_bar,cost_ey,k
    close(9)
    
    open(unit=9,file='v_ini.txt')
        write(9,*) av_V_ini
    close(9)
    
end subroutine  
    