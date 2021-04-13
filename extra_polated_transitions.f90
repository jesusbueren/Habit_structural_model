subroutine extrapolated_transitions(beta,H)
    use global_var; use nrtype;use state_space_dim
    implicit none
    real(DP),dimension(covariates,clusters,clusters+1),intent(in)::beta
    real(DP),dimension(clusters+1,clusters+1,T,types,L_gender,L_educ),intent(out)::H
    integer::e_l,c_l,c_l2,g_l,ge_l,age,ge_d,it,max_loc,d_l,c_l3,t_l
    integer,dimension(types-1)::e_d
    real(DP),dimension(covariates,1)::x
    real(DP),dimension(clusters,clusters,generations,types,L_gender,L_educ)::H_new
    integer,parameter::sims=1000
    integer,dimension(clusters+1)::counter_h
    real(DP),dimension(clusters+1)::h_star
    double precision,dimension(clusters+1,generations)::p
    integer,parameter::nodes=5
    real(DP),dimension(nodes):: xs=(/0.117581320211778,	1.0745620124369,	3.08593744371755,	6.41472973366203,	11.8071894899717/), &
                                weight=(/1.22172526747065,	0.480277222164629,	0.0677487889109621,	0.00268729149356246,	1.52808657104652E-05/),prod1,prod2
    integer::ind
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    !!$OMP PARALLEL  DEFAULT(PRIVATE) SHARED(H,beta)
    !!$OMP  DO collapse(4)
    do t_l=1,types; do c_l=1,clusters; do g_l=1,T; do ge_l=1,2;do e_l=1,L_educ
        age=first_age_sm+(g_l-1)*2-70
        ge_d=ge_l-1
        x(1:5,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),dble(ge_d),dble(age*ge_d)/)
        x(6:9,1)=0.0d0
        if (t_l==2)then
            x(6,1)=1.0
            x(7,1)=dble(age)
        end if
        x(8:13,1)=0.0d0
        if (e_l==2) then
            x(8,1)=1.0d0
            x(9,1)=dble(age)
        elseif (e_l==3) then
            x(10,1)=1.0d0
            x(11,1)=dble(age)
        end if
        if (t_l==2 .and. e_l==2) then
            x(12,1)=1.0d0
        elseif (t_l==2 .and. e_l==3) then
            x(13,1)=1.0d0
        end if
        !By simulation
        !counter_h=0
        !do it=1,sims
        !    do c_l2=1,clusters+1
        !        h_star(c_l2)=sum(x(:,1)*beta(:,c_l,c_l2))+c4_normal_01( )
        !    end do
        !    max_loc=maxloc(h_star,1)
        !    do c_l2=1,clusters+1
        !        if (c_l2==max_loc)then
        !            counter_h(c_l2)=counter_h(c_l2)+1
        !        end if
        !    end do
        !end do
        !H(c_l,:,g_l,e_l,ge_l)=(dble(counter_h))/dble(sims)
        !By numerical integration
        do c_l2=1,clusters+1
            h_star(c_l2)=sum(x(:,1)*beta(:,c_l,c_l2))
        end do
        do c_l2=1,clusters+1
            prod1=1
            prod2=1
            do c_l3=1,clusters+1
                if (c_l2/=c_l3) then
                    prod1=prod1*0.5d0*(1.0d0+erf((-sqrt(2.0d0*xs)-(h_star(c_l3)-h_star(c_l2)))/sqrt(2.0d0)))
                    prod2=prod2*0.5d0*(1.0d0+erf(( sqrt(2.0d0*xs)-(h_star(c_l3)-h_star(c_l2)))/sqrt(2.0d0)))
                end if
            end do
            H(c_l,c_l2,g_l,t_l,ge_l,e_l)=0.5d0/sqrt(pi)*sum(weight*(prod1+prod2))
        end do
    end do; end do; end do; end do;end do
    !!$OMP END DO
    !!$OMP END PARALLEL

    !No resurection
    H(clusters+1,1:clusters,:,:,:,:)=0.0_dp 
    H(clusters+1,clusters+1,:,:,:,:)=1.0_dp

end subroutine