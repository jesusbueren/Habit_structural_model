	SUBROUTINE powell(p,xi,ftol,iter,fret)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : linmin
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: xi
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(DP), INTENT(IN) :: ftol
	REAL(DP), INTENT(OUT) :: fret
    interface
        function obj_function_costs(parameters)
            use nrtype; use preference_p
            implicit none
            real(DP),dimension(:),intent(in)::parameters
            real(DP)::obj_function_costs
        end function  
    end interface
	INTEGER(I4B), PARAMETER :: ITMAX=400
	REAL(DP), PARAMETER :: TINY=1.0e-25_dp
	INTEGER(I4B) :: i,ibig,n
	REAL(DP) :: del,fp,fptt,t
	REAL(DP), DIMENSION(size(p)) :: pt,ptt,xit
	n=assert_eq(size(p),size(xi,1),size(xi,2),'powell')
	fret=obj_function_costs(p)
	pt(:)=p(:)
	iter=0
	do
		iter=iter+1
		fp=fret
		ibig=0
		del=0.0
		do i=1,n
			xit(:)=xi(:,i)
			fptt=fret
			call linmin(p,xit,fret)
			if (fptt-fret > del) then
				del=fptt-fret
				ibig=i
			end if
        end do
        print*,2.0_sp*(fp-fret) , ftol*(abs(fp)+abs(fret))+TINY
		if (2.0_sp*(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY) RETURN
		if (iter == ITMAX) pause !call &
			!nrerror('powell exceeding maximum iterations')
		ptt(:)=2.0_dp*p(:)-pt(:)
		xit(:)=p(:)-pt(:)
		pt(:)=p(:)
		fptt=obj_function_costs(ptt)
		if (fptt >= fp) cycle
		t=2.0_sp*(fp-2.0_dp*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
		if (t >= 0.0) cycle
		call linmin(p,xit,fret)
		xi(:,ibig)=xi(:,n)
		xi(:,n)=xit(:)
	end do
	END SUBROUTINE powell
