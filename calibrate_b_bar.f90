subroutine calibrate_b_bar()
    use nrtype; use preference_p; use global_var; use var_first_step
    implicit none
    real(DP)::obj_function,av_VSL_do,b_bar_max,b_bar_min,VSL_data,k
    real(DP),dimension(T,G_educ,G_types,4)::asset_distribution
    real(DP),dimension(G_educ,G_types)::cost_ey
    real(DP),dimension(G_educ,G_types)::av_VSL,av_V_ini,joint_pr
    real(DP),dimension(8,generations,types,L_educ)::dist_assets_data
    integer::y_l,t_l,e_l,it
    character::end_k
    
    open(unit=9,file='parameter.txt')
        read(9,*) betas,beq_cur,c_floor,beq_mu,obj_function
    close (9)
    
    print*,'----------------------'
    print*,"Parameter"
    print('(A20,<1>F5.2)'),"beta",betas
    print('(A20,F10.2)'),"beq cur",beq_cur
    print('(A20,F10.2)'),"c floor",c_floor
    print('(A20,F10.2)'),"beq mu",beq_mu
    !print('(A20,<9>F6.3)'),"pr low beta",pr_betas
    pr_betas=1.0d0
    
    b_bar_min=0.0d0
    b_bar_max=7.1d0
    av_VSL_do=-9.0d0
    it=0
    
    VSL_data=6000.0d0
    
    !do while (abs(av_VSL_do-VSL_data)>100.0d0 .and. it<10) 
     !   b_bar=(b_bar_max+b_bar_min)/2.0d0
        !Solve model
        call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini)
        !Av VSL for HS dropouts
        av_VSL_do=sum(av_VSL(1,:)*fraction_types(1,:))
        if (av_VSL_do>VSL_data) then
            b_bar_max=b_bar
        else
            b_bar_min=b_bar
        end if
        it=it+1
        print*,it,b_bar,av_VSL_do
    !end do
    
    !Estimate costs:
    !!!!!!!!!!!!!!!!    
    
    
    !Specification I. 9 cost parameter, var=1.
    do e_l=1,3
        joint_pr(e_l,:)=fraction_types(e_l,:)*fraction_e(e_l)
    end do
    print*,'cov',sum((log(joint_pr)-sum(log(joint_pr))/dble(G_educ*G_types))*(av_V_ini-sum(av_V_ini)/dble(G_educ*G_types)))/dble(G_educ*G_types)
    print*,'var x',sum((av_V_ini-sum(av_V_ini)/dble(G_educ*G_types))**2.0d0)/dble(G_educ*G_types)
    k=sum((log(joint_pr)-sum(log(joint_pr))/dble(G_educ*G_types))*(av_V_ini-sum(av_V_ini)/dble(G_educ*G_types)))/sum((av_V_ini-sum(av_V_ini)/dble(G_educ*G_types))**2.0d0)
    k=1.0d0/k
    print*,'k=',k
    print*,'cost: e x y'
    open(unit=9,file='costs_e_y.txt')
    do e_l=1,L_educ;do y_l=1,types
        cost_ey(e_l,y_l)=av_V_ini(e_l,y_l)-av_V_ini(1,3)-k*(log(joint_pr(e_l,y_l))-log(joint_pr(1,3)))
        print'(A4,I4,A4,I4,A4,F7.3,A5,F7.3)','e_l',e_l,'y_l',y_l,'%',fraction_types(e_l,y_l)*fraction_e(e_l),'cost',cost_ey(e_l,y_l)
         write(9,'(A4,I4,A4,I4,F20.3,F7.3,F7.3,F7.3)'),'e_l',e_l,'y_l',y_l,cost_ey(e_l,y_l),av_V_ini(e_l,y_l),joint_pr(e_l,y_l),exp(1.0d0/k*av_V_ini(e_l,y_l))/sum(exp(1.0d0/k*av_V_ini(:,:)))
    end do;end do
    close (9)


    
    !pause
end subroutine  
    