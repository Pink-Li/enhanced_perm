!Global Variables and directory definations for enhanced permeability model
!-------------------------------------------------------------------------------------------------------------------------------
!----------------------------------Pritom Sarma---------------------------------------------------------------------------------
!--------------------------------------IITBBS-----------------------------------------------------------------------------------
module global_var_dir
    implicit none

    integer                                    :: nr, nt !Spatial nodes and Time steps
    real(8)                                    :: P0, Pend, alpha, gamm, Ld, dt, dx, discret_lev !IC, BC, diffusivity, pressure sentivity, length of domain, temporal and spatial discretization
    real(8), dimension(:), allocatable         :: xvec
    ! File names:
    character(len=60), parameter               :: outdir = '/output/', indir = '/input/'
    character(len=500)                         :: cwd, inparams, spacefile, processfile, cpupath

    
    contains

    subroutine directories

        implicit none

        call getcwd(cwd)
        
        inparams    = trim(cwd)//trim(indir)//trim('inparams') ! The input parameters
        processfile = trim(cwd)//trim(outdir)//trim('process_log.out') !Process logs
        spacefile   = trim(cwd)//trim(outdir)//trim('spaceout_nonlinear.out')!Output from model
        
    end subroutine directories

end module global_var_dir  