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
    !    p(p_l,:)=(/(betas_ini-beta_min)/(beta_max-beta_min),beq_cur_ini,c_floor_ini,beq_mu_ini/)!(betas_ini-beta_min)/(beta_max-beta_min),,reshape(pr_betas_ini,(/G_types*G_educ,1/))
    !    if (p_l>1) then
    !        p(p_l,p_l-1)=p(p_l,p_l-1)*0.8d0
    !    end if
    !    p(p_l,:)=(/log(p(p_l,1)/(1.0d0-p(p_l,1))),log(p(p_l,2)),log(p(p_l,3)),log(p(p_l,4))/)!,log(p(p_l,6:PAR)/(1.0d0-p(p_l,6:PAR)))
    !    y(p_l)=obj_function(p(p_l,:))
    !end do
    !call amoeba(p,y,ftol,obj_function,iter)
    
    !Calibrate b_bar to match the value statistical life of 4,000,000 for the average dropout
    call calibrate_b_bar()
    
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
    real(DP),dimension(G_educ,G_types)::av_VSL,av_V_ini
    real(DP),dimension(8,generations,types,L_educ)::dist_assets_data
    integer::y_l,t_l,e_l
    character::end_k
    
    open(unit=9,file=path//'metric_model\Results\wealth_moments_data.txt')
            read(9,*) dist_assets_data 
    close(9)
    
    betas=1.0d0/(1.0d0 + exp(-parameters(1))) *(beta_max-beta_min)+ beta_min
    beq_cur=exp(parameters(2))
    c_floor=exp(parameters(3))
    beq_mu=exp(parameters(4))
    pr_betas=1.0d0 !reshape(1.0d0/(1.0d0 + exp(-parameters(6:PAR))),shape(pr_betas) )
    
    print*,'----------------------'
    print*,"Parameter"
    print('(A20,<1>F5.2)'),"beta",betas
    print('(A20,F10.2)'),"beq cur",beq_cur
    print('(A20,F10.2)'),"c floor",c_floor
    print('(A20,F10.2)'),"beq mu",beq_mu
    !print('(A20,<9>F6.3)'),"pr low beta",pr_betas
    
    call solve_and_simulate_model(asset_distribution,av_VSL,av_V_ini)
    
    obj_function=0.0d0
    do e_l=1,G_educ;do y_l=1,G_types        
        obj_function=obj_function+sum((asset_distribution(1:27,e_l,y_l,1)-dist_assets_data(4,1:27,y_l,e_l)/1000.0d0)**2.0d0) & 
                                 +sum((asset_distribution(1:27,e_l,y_l,2)-dist_assets_data(5,1:27,y_l,e_l)/1000.0d0)**2.0d0) & 
                                 +sum((asset_distribution(1:27,e_l,y_l,3)-dist_assets_data(6,1:27,y_l,e_l)/1000.0d0)**2.0d0)    
        print('(I4,I4,<3>F8.2)'),e_l,y_l,sum((asset_distribution(1:27,e_l,y_l,1)-dist_assets_data(4,1:27,y_l,e_l)/1000.0d0)**2.0d0)/1e6, &
                                     sum((asset_distribution(1:27,e_l,y_l,2)-dist_assets_data(5,1:27,y_l,e_l)/1000.0d0)**2.0d0)/1e6, &
                                     sum((asset_distribution(1:27,e_l,y_l,3)-dist_assets_data(6,1:27,y_l,e_l)/1000.0d0)**2.0d0) /1e6
    end do;end do

    
    print('(A20,F20.3)'),'obj fct',obj_function/1e6
    
    if (obj_function<best_obj_fct) then
        best_obj_fct=obj_function
        open(unit=9,file='asset_distribution_model.txt')
            do e_l=1,G_educ;do y_l=1,G_types;do t_l=1,T
                write(9,'(I4,I4,I4,<4>F20.8)'), t_l,e_l,y_l,asset_distribution(t_l,e_l,y_l,1),asset_distribution(t_l,e_l,y_l,2),asset_distribution(t_l,e_l,y_l,3),asset_distribution(t_l,e_l,y_l,4)
            end do;end do;end do
        close (9)
        open(unit=9,file='parameter.txt')
            write(9,'(<PAR>F10.3,F20.2)'),betas,beq_cur,c_floor,beq_mu,obj_function/1e6
        close (9)
    end if
    
    pause
end function  
    

    
    