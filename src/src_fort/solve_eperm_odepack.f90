program solve_eperm_vode
    
    
    use global_var_dir
    use util_eperm
    implicit none
    real(8), dimension(:), allocatable    :: yt, dydt, yt_scale, ytsv, dydtsv, yinit, px, nones, ts, cpu1, cpu2, cpu3, cpu4
    real(8), dimension(:,:), allocatable  :: spaceout 
    real(8)                               :: time, timesv, dt_try, dt_next, rtol, dt_did, atol
    real(8)                               :: start, finish, e_time
    integer                               :: i, j, sizespace(2), it, irecspace, outcount
    real(8), dimension(:), allocatable    :: dt_save
!-----------------------------------------------------------------------------------------------------
! Variables specific to dlsode
!-----------------------------------------------------------------------------------------------------
    integer                               :: itol, itask, iopt, istate
    integer                               :: neq, lrw, liw, jt, ml, mu
    real(8), dimension(:), allocatable    :: rwork
    integer, dimension(:), allocatable    :: iwork
    real(8)                               :: tout, acc
    real(8), dimension(:), allocatable    :: rpar
    integer, dimension(:), allocatable    :: ipar
    call cpu_time(start)    !Start Timer
!-----------------------------------------------------------------------------------------------------
! Open formatted files first
!-----------------------------------------------------------------------------------------------------
    call directories
    open(unit=8,file=inparams,STATUS='OLD')
    outcount = 0
!-----------------------------------------------------------------------------------------------------
! Read in input parameters and initialize variables
!-----------------------------------------------------------------------------------------------------
! Reading data from file 'inparams'
    read(8,*) P0            ! Injection Pressure
    read(8,*) Pend          ! Ambient Pressure
    read(8,*) alpha         ! Diffusivity
    read(8,*) gamm          ! Pressure
    read(8,*) Ld            ! Domain Lenght
    read(8,*) dt            ! Timestepping
    read(8,*) nt            ! Timesteps
    read(8,*) discret_lev   ! Spatial Discreetization Level
    close(8)
! Initializing Variables
    dx   = 1/discret_lev    ! Grid Spacing
    nr   = ceiling(Ld/dx)   ! Number of Grid Points
!-----------------------------------------------------------------------------------------------------
! Allocate all major allocatables and initialize them
!-----------------------------------------------------------------------------------------------------
    allocate(yt(nr), dydt(nr), yt_scale(nr), ytsv(nr), dydtsv(nr), nones(nr),xvec(nr),ts(nt), yinit(nr), px(nr))
    allocate(spaceout(nr,3))
    allocate(cpu1(nt-1), cpu2(nt-1), cpu3(nt-1), cpu4(nt-1))
    xvec(:) = 0.d0
    spaceout = 0.d0
    ts      =  0.d0
    do i = 2,nr
        xvec(i) = xvec(i-1) + dx
    enddo
    do i=1,nt
        ts(i) = (i-1)*dt
    enddo
!Specific to DSLODA    
    neq = nr ! Number of equations
    ml = 1 ! Tridiagonal system, upper half-bandwidth
    mu = 1 ! Tridiagonal system, lower half-bandwidth
    jt  = 25
    lrw = 22 + max(16,10+(2*ml+mu))*neq ! Work variables for dlsoda, max(LRN,LRW), refer to documentation for sizes recommended when jt = 4 or 5, and NEQ is constant
    liw = 20 + neq ! Integer work variable, min value: 20 + neq
    itol=1; itask=1; iopt=1; istate=1
    allocate(dt_save(nt),rwork(lrw))
    allocate(iwork(liw))
    iwork = 0
    rwork = 0.d0
    dt_save = 0.d0
    iwork(1) = ml; iwork(2) = mu; iwork(6) = 5000
!-----------------------------------------------------------------------------------------------------
! Open unformatted files now
!-----------------------------------------------------------------------------------------------------
    inquire (iolength = irecspace) spaceout(:,1)
    open(unit=20,file=spacefile,form="unformatted", access="direct", recl=irecspace, status = 'REPLACE')
!---------Main Job----------------------------------------------------------------------
    
    write(*,*) 'Hey!'
    write(*,*) 'Grid Spacing', nr
    write(*,*) 'Diff. =', alpha, ', Gamma ', gamm
   
    dt_try = dt
    time   = ts(1)
    dydt   = 0.d0
    rtol   = 1.d-6
    atol   = 1.d-2
    acc    = 1.d-6
    
    yinit(1)= exp(P0*gamm)/gamm
    yinit(2:nr) = exp(Pend*gamm)/gamm
    yt = yinit
    
    write(*,*) 'ita_start =', yt(1), ', ita_end ', yt(nr)

do it = 1,nt-1
       
    istate = 1 !Initialize error flag
!-----------------------------------------------------------------------------------
! Use odepack dlsoda
!-----------------------------------------------------------------------------------
    time = ts(it)
    tout = ts(it+1)

    
    call dlsode(derivs_odepack,neq,yt,time,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,jt)


    if (istate .lt. 0) then
        write(*,*) 'Error halt.. istate = ', istate
    if (istate==-5) istate=1
    endif
    dt_did = RWORK(11)
    write(6,60)IWORK(11),IWORK(12),IWORK(13),IWORK(19),RWORK(15)
 60 format(/' No. steps =',I4,'  No. f-s =',I6,'  No. J-s =',I4,' Method last used =',I2,' Last switch was at t =',f12.4)
    outcount = outcount + 1
    yt(1) = yinit(1)
    yt(nr) = yinit(nr)
    px     = log(gamm*yt)/gamm   ! Pressure for output

    write(6,'(i9,e19.6,3e15.6,i7)')it,time,dt_did,px(1),px(nr)/Pend, outcount
 ! Save an array of time steps
    dt_save(it) = dt_did

 ! Assign values to output matrix 1
    spaceout(:,1) = xvec; spaceout(:,2) = px; spaceout(1,3) = tout; spaceout(2,3) = alpha; spaceout(3,3) = P0
    spaceout(4,3) = gamm
    sizespace     = shape(spaceout)
    do i = 1,minval(sizespace)
        write(20,rec=i+(outcount-1)*minval(sizespace)) spaceout(:,i) ! Write out spacefile
    enddo
enddo
call cpu_time(finish)   !Stop Timer
e_time = real(finish-start)
write(*,*) 'Model Run Successful'
write(*,*) 'Total clock time for dvode: ', e_time
write(*,*) 'Time steps', nt 
write(*,*) 'Space nodes', nr  
deallocate(yt,dydt,yinit,ytsv,dydtsv,yt_scale,nones,xvec)
close(20)
end program solve_eperm_vode
