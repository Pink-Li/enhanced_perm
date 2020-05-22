!Subroutines and jacobian for converting the PDE to a MOL based ODE and feeding into the solver
!-------------------------------------------------------------------------------------------------------------------------------
!----------------------------------Pritom Sarma---------------------------------------------------------------------------------
!--------------------------------------IITBBS-----------------------------------------------------------------------------------

module util_eperm

    use global_var_dir
    use omp_lib

    implicit none

    contains

    subroutine derivs_odepack(neq,t,y,dydt)
        implicit none
            INTEGER, INTENT(IN)                     :: neq
            REAL(8), INTENT(IN)                     :: t
            REAL(8), DIMENSION(:), INTENT(IN)       :: y(neq)
            REAL(8), DIMENSION(:), INTENT(OUT)      :: dydt(neq)
            integer :: i,chunksize, nthread
            real(8) :: dx
            dx = xvec(2)-xvec(1) ! Spatial discretization
            dydt(:) = 0.d0
            nthread = 2
            chunksize = (nr-2)/nthread
            !call OMP_SET_NUM_THREADS(nthread)
            !!$OMP PARALLEL DO SHARED(y,chunksize,dx) PRIVATE(i) SCHEDULE(static,chunksize)
            do i = 2,nr-1
                dydt(i) = alpha*gamm*y(i)*((y(i+1) - 2.d0*y(i) + y(i-1)))/(dx**2.d0) !MOL based discretization of the non-linear PDE   
            enddo
            !!$OMP END PARALLEL DO
            ! B.C.'s :
            dydt(1)  = 0.d0
            dydt(nr) = 0.d0
    
    end subroutine derivs_odepack
    

    subroutine derivs_dopri5(neq,t,y,f,rpar,ipar)
        
        implicit none
        
        INTEGER, INTENT(IN)  :: neq
        REAL(8), INTENT(IN)  :: t
        REAL(8), INTENT(IN)  :: rpar
        INTEGER, INTENT(IN)  :: ipar
        REAL(8)              :: y(neq), f(neq)
        INTEGER              :: i
        REAL(8)              :: dx  
        dx = xvec(2)-xvec(1); ! Spatial discretization
        f = 0.d0;  
        do i = 2,nr-1
            f(i)=alpha*gamm*y(i)*((y(i+1) - 2.d0*y(i) + y(i-1)))/(dx**2.d0) !MOL based discretization of the non-linear PDE  
        enddo
        ! B.C.'s :
        f(1)  = 0.d0
        f(nr) = 0.d0
    end subroutine derivs_dopri5    

    subroutine solout(nr,told,t,y,n,con,icomp,nd,rpar,ipar,irtrn)

        implicit none                    
        real(8)                           :: y(n),con(5*nd), t, told
        integer                           :: nd, icomp(nd), irtrn, nr, n
        real(8), dimension(:), intent(in) :: rpar
        integer, dimension(:), intent(in) :: ipar
   
    end subroutine solout
    
    subroutine jac(neq, t, y, ml, mu, pd, nrowpd)
        implicit none

        
        
                real(8), intent(in)                               :: t
                integer, intent(in)                               :: ml, mu, nrowpd, neq
                real(8)                                           :: y(neq), dx
                real(8), dimension(nrowpd,neq)                    :: pd
                
                dx = xvec(2)-xvec(1) ! Spatial discretization
                pd(:,:) = 0.d0
        
                pd(mu+1,:) = -2.d0/dx ! Diagonal
                pd(mu,:)   = 1.d0/dx + 1.d0/(2.d0*xvec(:)) ! Upper diagonal
                pd(mu+2,:) = 1.d0/dx - 1.d0/(2.d0*xvec(:)) ! Lower diagonal
                
                pd(mu,2)       = 2.d0/dx
                pd(mu+1,1)     = -2.d0/dx
                pd(mu+2,neq-1) = 2.d0/dx
                pd(mu+1,neq)   = -2.d0/dx
                pd(:,:) = alpha*pd(:,:)/dx
        
        
    
    end subroutine jac
 
    subroutine derivs_vode(neq,t,y,dydt,rpar, ipar)
        implicit none
            INTEGER, INTENT(IN)                     :: neq
            REAL(8), INTENT(IN)                     :: t
            REAL(8), DIMENSION(:), INTENT(IN)       :: y(neq)
            REAL(8), DIMENSION(:), INTENT(OUT)      :: dydt(neq)
            REAL(8), INTENT(IN)  :: rpar 
            INTEGER, INTENT(IN)  :: ipar
            integer :: i,chunksize, nthread
            real(8) :: dx
            dx = xvec(2)-xvec(1) ! Spatial discretization
            dydt(:) = 0.d0
            nthread = 2
            chunksize = (nr-2)/nthread
            !call OMP_SET_NUM_THREADS(nthread)
            !!$OMP PARALLEL DO SHARED(y,chunksize,dx) PRIVATE(i) SCHEDULE(static,chunksize)
            do i = 2,nr-1
                dydt(i) = alpha*gamm*y(i)*((y(i+1) - 2.d0*y(i) + y(i-1)))/(dx**2.d0) !MOL based discretization of the non-linear PDE   
            enddo
            !!$OMP END PARALLEL DO
            ! B.C.'s :
            dydt(1)  = 0.d0
            dydt(nr) = 0.d0
    
    end subroutine derivs_vode

    subroutine jac_vode(neq, t, y, ml, mu, pd, nrowpd, rpar, ipar)
        implicit none

        
        
                real(8), intent(in)                               :: t
                integer, intent(in)                               :: ml, mu, nrowpd, neq
                real(8)                                           :: y(neq), dx
                real(8), dimension(nrowpd,neq)                    :: pd
                REAL(8), INTENT(IN)  :: rpar
                INTEGER, INTENT(IN)  :: ipar
                
                dx = xvec(2)-xvec(1) ! Spatial discretization
                pd(:,:) = 0.d0
        
                pd(mu+1,:) = -2.d0/dx ! Diagonal
                pd(mu,:)   = 1.d0/dx + 1.d0/(2.d0*xvec(:)) ! Upper diagonal
                pd(mu+2,:) = 1.d0/dx - 1.d0/(2.d0*xvec(:)) ! Lower diagonal
                
                pd(mu,2)       = 2.d0/dx
                pd(mu+1,1)     = -2.d0/dx
                pd(mu+2,neq-1) = 2.d0/dx
                pd(mu+1,neq)   = -2.d0/dx
                pd(:,:) = alpha*pd(:,:)/dx
        
        
    
    end subroutine jac_vode
end module util_eperm    