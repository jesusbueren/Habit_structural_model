program main
    use nrtype; use preference_p; use initial_p
    implicit none
    integer,dimension(1)::seed=254 
    real(DP),dimension(PAR+1)::y
    real(DP),dimension(PAR+1,PAR)::p
    integer::p_l,iter
    real(DP)::ftol=1d-4 
    character::end_k
    interface
        function obj_function(parameters)
            use nrtype; use preference_p
            implicit none
            real(DP),dimension(:),intent(in)::parameters
            real(DP)::obj_function
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
    
    call load_first_step_p()
    
    !Calibrate model by matching wealth profile
    !do p_l=1,PAR+1
    !    p(p_l,:)=(/(betas_ini-beta_min)/(beta_max-beta_min),beq_cur_ini,c_floor_ini,beq_mu_ini,pr_betas_ini/) !
    !    if (p_l>1) then
    !        p(p_l,p_l-1)=p(p_l,p_l-1)*0.95d0
    !    end if
    !    p(p_l,:)=(/log(p(p_l,1:2)/(1.0d0-p(p_l,1:2))),log(p(p_l,3)),log(p(p_l,4)),log(p(p_l,5)),log(p(p_l,6:PAR)/(1.0d0-p(p_l,6:PAR)))/) 
    !    !p(p_l,:)=(/log(p(p_l,1)/(1.0d0-p(p_l,1))),log(p(p_l,2)),log(p(p_l,3)),log(p(p_l,4))/) 
    !    y(p_l)=obj_function(p(p_l,:))
    !end do
    !call amoeba(p,y,ftol,obj_function,iter)
    
    !Calibrate b_bar to match the value statistical life of 4,000,000 for the average dropout
    !call calibrate_b_bar()
    
    
    !call welfare_analysis()
    call counterfactuals()
    
    
    
    print*,'End of program'
    print*,'Press any key to close window'
    read*,end_k
    
    
       
end program
    
function obj_function(parameters)
    use nrtype; use preference_p; use global_var
    implicit none
    real(DP),dimension(:),intent(in)::parameters
    real(DP)::obj_function
    real(DP),dimension(T,G_educ,G_types,4)::asset_distribution
    real(DP),dimension(G_DF,G_educ,G_types)::av_VSL,av_V_ini
    real(DP),dimension(9,generations,G_h,types,L_educ)::dist_assets_data
    real(DP),dimension(G_educ,G_h,T)::p50_delta
    integer::y_l,t_l,e_l,h_l,t_ini,t_final,t_l52
    real(DP),dimension(G_educ,G_h,10)::delta_assets_data_MA
    real(DP),dimension(G_educ,G_h,T)::delta_assets_data
    real(DP)::moment_data,moment_model
    character::end_k
    
    !Moments targetted
    open(unit=9,file=path//'metric_model\Results\wealth_moments_data.txt')
            read(9,*) dist_assets_data 
    close(9)
    dist_assets_data=dist_assets_data/1000.0d0
    open(unit=9,file=path//'data\delta_assets.csv')
            read(9,*) delta_assets_data_MA 
    close(9)
    delta_assets_data=-9.0d0
    do t_l=1,10
        t_l52=(52-first_age_sm)/2+1
        delta_assets_data(:,:,t_l52+(t_l-1)*2)=delta_assets_data_MA(:,:,t_l) !delta_assets_data(3,1,:)
    end do

    betas=1.0d0/(1.0d0 + exp(-parameters(1:2)))*(beta_max-beta_min)+ beta_min !1.0d0/(1.0d0 + exp(-parameters(1:2)))*(beta_max-beta_min)+ beta_min
    beq_cur=exp(parameters(3))
    c_floor=exp(parameters(4))
    beq_mu=exp(parameters(5))
    pr_betas=reshape(1.0d0/(1.0d0 + exp(-parameters(6:PAR))),shape(pr_betas) )
    
    print*,'----------------------'
    print*,"Parameter"
    print('(A20,<3>F5.2)'),"beta",betas
    print('(A20,F10.2)'),"beq cur",beq_cur
    print('(A20,F10.2)'),"c floor",c_floor
    print('(A20,F10.2)'),"beq mu",beq_mu

    RRA_beq=RRA


    print('(A20,<9>F6.3)'),"pr low beta",pr_betas
    
    call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini,p50_delta)
    t_ini=5
    t_final=T_R+4
    
    obj_function=0.0d0
    do e_l=1,G_educ;do y_l=1,G_types     

        obj_function=obj_function & !+sum(((asset_distribution(t_ini:t_final,e_l,y_l,1)-dist_assets_data(5,t_ini:t_final,1,y_l,e_l)))**2.0d0) & 
                                 +sum(((asset_distribution(t_ini:t_final,e_l,y_l,2)-dist_assets_data(6,t_ini:t_final,1,y_l,e_l))/dist_assets_data(6,t_ini:t_final,1,y_l,e_l))**2.0d0) !& 
                                 !+sum(((asset_distribution(t_ini:t_final,e_l,y_l,3)-dist_assets_data(7,t_ini:t_final,1,y_l,e_l)))**2.0d0)    
        print('(I4,I4,<3>F20.2)'),e_l,y_l, & !sum(((asset_distribution(t_ini:t_final,e_l,y_l,1)-dist_assets_data(5,t_ini:t_final,1,y_l,e_l))/dist_assets_data(5,t_ini:t_final,1,y_l,e_l))**2.0d0), &
                                     sum(((asset_distribution(t_ini:t_final,e_l,y_l,2)-dist_assets_data(6,t_ini:t_final,1,y_l,e_l))/dist_assets_data(6,t_ini:t_final,1,y_l,e_l))**2.0d0) !, &
                                     !sum(((asset_distribution(t_ini:t_final,e_l,y_l,3)-dist_assets_data(7,t_ini:t_final,1,y_l,e_l))/dist_assets_data(7,t_ini:t_final,1,y_l,e_l))**2.0d0)

    end do;end do
    
    
    print('(A20,F20.3)'),'obj fct',obj_function
    
    if (obj_function<best_obj_fct) then
        best_obj_fct=obj_function
        open(unit=9,file='asset_distribution_model.txt')
        open(unit=10,file='delta_assets_model_vs_data.txt')
            do e_l=1,G_educ;do y_l=1,G_types;do t_l=1,T
                write(9,'(I4,I4,I4,<4>F20.8)'), t_l,e_l,y_l,asset_distribution(t_l,e_l,y_l,1),asset_distribution(t_l,e_l,y_l,2),asset_distribution(t_l,e_l,y_l,3),asset_distribution(t_l,e_l,y_l,4)
            end do;end do;end do
            do e_l=1,G_educ;do t_l=1,T
                write(10,'(I4,I4,<4>F20.8)'), t_l,e_l,p50_delta(e_l,1,t_l),p50_delta(e_l,2,t_l),delta_assets_data(e_l,1,t_l),delta_assets_data(e_l,2,t_l)
            end do;end do
        close (9)
        close (10)
        open(unit=9,file='parameter.txt')
            write(9,'(<PAR>F10.3,F20.5)'),betas,c_floor,beq_cur,beq_mu,pr_betas,obj_function
        close (9)
    end if
    !pause

end function  
    

    
    