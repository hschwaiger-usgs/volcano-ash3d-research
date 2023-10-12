      module Variable_Diffusivity

      use precis_param

      use io_units

        ! Smagorinsky (1993) constant for LES horizontal diffusivity
        ! (should be 0.2 - 0.9)
      real(kind=ip),parameter :: Coeff_Stab   =   9.2_ip
      real(kind=ip),parameter :: Coeff_UnStab = -13.0_ip
      real(kind=ip),parameter :: Exp_UnStab   =  -0.5_ip

      real(kind=ip) :: KH_SmagC     != 0.9_ip
      real(kind=ip) :: vonKarman    != 0.4_ip   ! von Karman constant
      real(kind=ip) :: LambdaC      != 30.0_ip  ! Asymptotic length scale (m)
      real(kind=ip) :: RI_CRIT      != 0.25_ip  ! Critical Richardson number

      !Coeff_Stab   =   6.9_ip
      !Coeff_UnStab = -22.0_ip
      !Exp_UnStab   = -0.25

      ! Set the number of output variables for this module
      integer, parameter :: nvar_User2d_static_XY_VarDiff = 0
      integer, parameter :: nvar_User2d_XY_VarDiff        = 0
      integer, parameter :: nvar_User3d_XYGs_VarDiff      = 0
      integer, parameter :: nvar_User3d_XYZ_VarDiff       = 2 ! khorz, kvert
      integer, parameter :: nvar_User4d_XYZGs_VarDiff     = 0

      character(len=30),dimension(nvar_User3d_XYZ_VarDiff) :: temp_3d_name_VarDiff
      character(len=30),dimension(nvar_User3d_XYZ_VarDiff) :: temp_3d_unit_VarDiff
      character(len=30),dimension(nvar_User3d_XYZ_VarDiff) :: temp_3d_lname_VarDiff
      real(kind=op),    dimension(nvar_User3d_XYZ_VarDiff) :: temp_3d_MissVal_VarDiff
      real(kind=op),    dimension(nvar_User3d_XYZ_VarDiff) :: temp_3d_FillVal_VarDiff

      ! These are used to keep track of which index in the global list, this
      ! modules output vars corespond to
      integer :: indx_User2d_static_XY_VarDiff
      integer :: indx_User2d_XY_VarDiff
      integer :: indx_User3d_XYGs_VarDiff
      integer :: indx_User3d_XYZ_VarDiff
      integer :: indx_User4d_XYZGs_VarDiff

      real(kind=ip) :: LES_L2ScaleCoeff

      ! Note: RoughLen_z can be related to Land use
      !    From Stohl et al, ACP, v5n9p2461, 2005 Table 3
      !         Grassland       :: 0.10
      !         Arable land     :: 0.15
      !         Permanent crops :: 0.30
      !         Forest          :: 0.60
      !         Inland water    :: Charnock
      !         Urban areas     :: 0.70
      !         Other           :: 0.10
      !         Ocean           :: Charnock
      !   Note: Charnok relation is z_0 = a(u_star^2/g) with a~ 0.018
      !    From Stohl et al, Tech Note FLEXPART 8.2 :: surfdata.t
      !         landuse   comment                               z0        glcf
      !         --------------------------------------------------------
      !          1 Urban land                                   0.7       13
      !          2 Agricultural land                            0.1       11
      !          3 Range land                                   0.1       10
      !          4 Deciduous forest                             1.         3,4
      !          5 Coniferous forest                            1.         1,2
      !          6 Mixed forest including wetland               0.7        5
      !          7 water, both salt and fresh                   0.001      0
      !          8 barren land mostly desert                    0.01      12
      !          9 nonforested wetland                          0.1
      !         10 mixed agricultural and range land            0.1        6,7,8,9
      !         11 rocky open areas with low grow shrubs        0.05 
      !         12 snow and ice                                 0.001
      !         13 rainforest                                   1.
      ! For resuspension cases, friction velocity is needed at every cell and so
      ! will the z0 (RoughLen_z).  Vr (U10 and V10) need to be calculated of
      ! read for the wind grid and regridded to the comp grid.
      
      ! 3d Variables needed on MetP grid
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dVel_dz_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: du_dx_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: du_dy_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dv_dx_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dv_dy_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dV_dz_MetP_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: SurfRoughLen_Met_sp

        ! and at both meso steps (also MetP)
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Ri_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Ri_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_next_step_MetP_sp

      ! 2d variables needed at meso steps on Met grid
      real(kind=sp),dimension(:,:)    ,allocatable :: PBLH_meso_last_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: PBLH_meso_next_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: L_MonOb_meso_last_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: L_MonOb_meso_next_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: FricVel_meso_last_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: FricVel_meso_next_step_Met_sp

      ! Variables needed on Comp grid (kx,y,z are already allocated)
      real(kind=sp),dimension(:,:)    ,allocatable :: FricVel_meso_last_step_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: FricVel_meso_next_step_sp
      real(kind=ip),dimension(:,:)    ,allocatable :: FricVel_ip
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_last_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_next_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_last_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_next_step_sp

      ! Both Khz adn Kv need U and V values on MetP grid so store local copies
      ! Note: The core Ash3d code reads directly into the computational grid
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vx_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vy_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vx_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vy_meso_next_step_MetP_sp

      contains

!******************************************************************************

      subroutine input_data_VarDiff

      use global_param,  only : &
         nmods,useTemperature,useVarDiffH,useVarDiffV

      use io_data,       only : &
         infile

      use MetReader,     only : &
         MR_Save_Velocities

      implicit none

      character(len=3)  :: answer
      character(len=80)  :: linebuffer
      integer :: ios,ioerr
      character(len=20) :: mod_name
      integer :: substr_pos

      open(unit=10,file=infile,status='old',err=1900)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Searching for OPTMOD=VARDIFF"
      endif;enddo
      nmods = 0
      read(10,'(a80)',iostat=ios)linebuffer
      do while(ios.eq.0)
        read(10,'(a80)',iostat=ios)linebuffer

        substr_pos = index(linebuffer,'OPTMOD')
        if(substr_pos.eq.1)then
          ! found an optional module
          !  Parse for the keyword
          read(linebuffer,1104)mod_name
          if(adjustl(trim(mod_name)).eq.'VARDIFF')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo

      useVarDiffH = .false.
      useVarDiffV = .false.
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Continue reading input file for VarDiff block"
      endif;enddo

      !Check if we're going to use variable diffusivity
      read(10,'(a80)',iostat=ios,err=2010)linebuffer
      read(linebuffer,'(a3)',err=2011) answer
      if (answer.eq.'yes') then
        useVarDiffH = .true.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Using horizontal variable diffusivity"
        endif;enddo
      elseif(answer(1:2).eq.'no') then
        useVarDiffH = .false.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Not using horizontal variable diffusivity"
        endif;enddo
      else
        goto 2011
      endif
      read(10,'(a80)',iostat=ios,err=2010)linebuffer
      read(linebuffer,'(a3)',err=2011) answer
      if (answer.eq.'yes') then
        useVarDiffV = .true.
        useTemperature = .true.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Using vertical variable diffusivity"
        endif;enddo
      elseif(answer(1:2).eq.'no') then
        useVarDiffV = .false.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Not using vertical variable diffusivity"
        endif;enddo
      else
        goto 2011
      endif

      if (useVarDiffH.or.useVarDiffV) then
        ! Check if we're using variable diffusivity, then get the constants
        read(10,'(a80)',iostat=ios,err=2010)linebuffer
        read(linebuffer,*,iostat=ioerr) KH_SmagC
        read(10,'(a80)',iostat=ios,err=2010)linebuffer
        read(linebuffer,*,iostat=ioerr) vonKarman
        read(10,'(a80)',iostat=ios,err=2010)linebuffer
        read(linebuffer,*,iostat=ioerr) LambdaC
        read(10,'(a80)',iostat=ios,err=2010)linebuffer
        read(linebuffer,*,iostat=ioerr) RI_CRIT

        !KH_SmagC  = 0.9
        !vonKarman = 0.4
        !LambdaC   = 30.0
        !RI_CRIT   = 0.25
        !do io=1,2;if(VB(io).le.verbosity_info)then
        !  write(outlog(io),*)KH_SmagC
        !  write(outlog(io),*)vonKarman
        !  write(outlog(io),*)LambdaC
        !  write(outlog(io),*)RI_CRIT
        !endif;enddo

        ! We will want to reuse velocities on the metP grid for this module
        MR_Save_Velocities = .true.

      endif

2010  continue
      close(10)

      return

1900  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot find input file: ',infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

2011  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading whether to use variable diffusivity.'
        write(errlog(io),*) 'Answer must be yes or no.'
        write(errlog(io),*) 'You gave:',linebuffer
        write(errlog(io),*) 'Program stopped'
      endif;enddo
      stop 1

      end subroutine input_data_VarDiff

!******************************************************************************

!******************************************************************************

      subroutine Allocate_VarDiff_Met

      use global_param,  only : &
         PI

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User3d_XYGs,nvar_User2d_XY,&
         nvar_User4d_XYZGs,nvar_User3d_XYZ

      use mesh,          only : &
         nxmax,nymax,nzmax

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet

      implicit none


      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOCATE_VARDIFF_MET ------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      allocate(dVel_dz_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(du_dx_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(du_dy_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(dv_dx_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(dv_dy_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(dV_dz_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(SurfRoughLen_Met_sp(nx_submet,ny_submet))
      allocate(Ri_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Ri_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Khz_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Khz_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Kv_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Kv_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(PBLH_meso_last_step_Met_sp(nx_submet,ny_submet))
      allocate(PBLH_meso_next_step_Met_sp(nx_submet,ny_submet))
      allocate(L_MonOb_meso_last_step_Met_sp(nx_submet,ny_submet))
      allocate(L_MonOb_meso_next_step_Met_sp(nx_submet,ny_submet))
      allocate(FricVel_meso_last_step_Met_sp(nx_submet,ny_submet))
      allocate(FricVel_meso_next_step_Met_sp(nx_submet,ny_submet))

      allocate(FricVel_meso_last_step_sp(nxmax,nymax))
      allocate(FricVel_meso_next_step_sp(nxmax,nymax))
      allocate(FricVel_ip(nxmax,nymax))

      allocate(Khz_meso_last_step_sp(nxmax,nymax,nzmax))
      allocate(Khz_meso_next_step_sp(nxmax,nymax,nzmax))
      allocate(Kv_meso_last_step_sp(nxmax,nymax,nzmax))
      allocate(Kv_meso_next_step_sp(nxmax,nymax,nzmax))

      allocate(vx_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(vy_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(vx_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(vy_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))

        ! Precalculate the LES term
      LES_L2ScaleCoeff = (KH_SmagC*KH_SmagC/real(PI*PI,kind=sp))

      ! Set the start indecies
      indx_User2d_static_XY_VarDiff = nvar_User2d_static_XY
      indx_User2d_XY_VarDiff        = nvar_User2d_XY
      indx_User3d_XYGs_VarDiff      = nvar_User3d_XYGs
      indx_User3d_XYZ_VarDiff       = nvar_User3d_XYZ
      indx_User4d_XYZGs_VarDiff     = nvar_User4d_XYZGs

      temp_3d_name_VarDiff(1) = "Kh"
      temp_3d_lname_VarDiff(1) = "Diffusivity_Horizontal"
      temp_3d_unit_VarDiff(1) = "m2/s"
      temp_3d_MissVal_VarDiff(1) = -9999.0_op
      temp_3d_FillVal_VarDiff(1) = -9999.0_op

      temp_3d_name_VarDiff(2) = "Kv"
      temp_3d_lname_VarDiff(2) = "Diffusivity_Vertical"
      temp_3d_unit_VarDiff(2) = "m2/s"
      temp_3d_MissVal_VarDiff(2) = -9999.0_op
      temp_3d_FillVal_VarDiff(2) = -9999.0_op

      nvar_User2d_static_XY = nvar_User2d_static_XY + nvar_User2d_static_XY_VarDiff
      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_VarDiff
      nvar_User3d_XYGs      = nvar_User3d_XYGs      + nvar_User3d_XYGs_VarDiff
      nvar_User3d_XYZ       = nvar_User3d_XYZ       + nvar_User3d_XYZ_VarDiff
      nvar_User4d_XYZGs     = nvar_User4d_XYZGs     + nvar_User4d_XYZGs_VarDiff

        !call MetReader subroutine to set up rdphi_wind rdlambda_wind

      end subroutine Allocate_VarDiff_Met

!******************************************************************************

      subroutine Prep_output_VarDiff

      use mesh,          only : &
         nxmax,nymax,nzmax

      use Diffusion,     only : &
         kx,kz

      use Output_Vars,   only : &
         var_User3d_XYZ_name,var_User3d_XYZ_unit,var_User3d_XYZ_lname,&
         var_User3d_XYZ_MissVal,var_User3d_XYZ_FillVal,var_User3d_XYZ

      implicit none

      integer :: i,indx

      ! Might have to build in some logic for Kh vs Kz

      do i=1,nvar_User3d_XYZ_VarDiff
        indx = indx_User3d_XYZ_VarDiff+i
        var_User3d_XYZ_name(indx)   = temp_3d_name_VarDiff(i)
        var_User3d_XYZ_unit(indx)   = temp_3d_unit_VarDiff(i)
        var_User3d_XYZ_lname(indx)  = temp_3d_lname_VarDiff(i)
        var_User3d_XYZ_MissVal(indx)= temp_3d_MissVal_VarDiff(i)
        var_User3d_XYZ_FillVal(indx)= temp_3d_FillVal_VarDiff(i)
        if(i.eq.1) var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = kx(1:nxmax,1:nymax,1:nzmax)
        if(i.eq.2) var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = kz(1:nxmax,1:nymax,1:nzmax)
      enddo

      end subroutine Prep_output_VarDiff

!******************************************************************************

      subroutine Deallocate_VarDiff_Met

      implicit none

      if(allocated(dVel_dz_MetP_sp))               deallocate(dVel_dz_MetP_sp)
      if(allocated(du_dx_MetP_sp))                 deallocate(du_dx_MetP_sp)
      if(allocated(du_dy_MetP_sp))                 deallocate(du_dy_MetP_sp)
      if(allocated(dv_dx_MetP_sp))                 deallocate(dv_dx_MetP_sp)
      if(allocated(dv_dy_MetP_sp))                 deallocate(dv_dy_MetP_sp)
      if(allocated(dV_dz_MetP_sp))                 deallocate(dV_dz_MetP_sp)
      if(allocated(SurfRoughLen_Met_sp))           deallocate(SurfRoughLen_Met_sp)
      if(allocated(Ri_meso_last_step_MetP_sp))     deallocate(Ri_meso_last_step_MetP_sp)
      if(allocated(Ri_meso_next_step_MetP_sp))     deallocate(Ri_meso_next_step_MetP_sp)
      if(allocated(Khz_meso_last_step_MetP_sp))    deallocate(Khz_meso_last_step_MetP_sp)
      if(allocated(Khz_meso_next_step_MetP_sp))    deallocate(Khz_meso_next_step_MetP_sp)
      if(allocated(Kv_meso_last_step_MetP_sp))     deallocate(Kv_meso_last_step_MetP_sp)
      if(allocated(Kv_meso_next_step_MetP_sp))     deallocate(Kv_meso_next_step_MetP_sp)
      if(allocated(PBLH_meso_last_step_Met_sp))    deallocate(PBLH_meso_last_step_Met_sp)
      if(allocated(PBLH_meso_next_step_Met_sp))    deallocate(PBLH_meso_next_step_Met_sp)
      if(allocated(L_MonOb_meso_last_step_Met_sp)) deallocate(L_MonOb_meso_last_step_Met_sp)
      if(allocated(L_MonOb_meso_next_step_Met_sp)) deallocate(L_MonOb_meso_next_step_Met_sp)
      if(allocated(FricVel_meso_last_step_Met_sp)) deallocate(FricVel_meso_last_step_Met_sp)
      if(allocated(FricVel_meso_next_step_Met_sp)) deallocate(FricVel_meso_next_step_Met_sp)

      if(allocated(FricVel_meso_last_step_sp))     deallocate(FricVel_meso_last_step_sp)
      if(allocated(FricVel_meso_next_step_sp))     deallocate(FricVel_meso_next_step_sp)
      if(allocated(FricVel_ip))                    deallocate(FricVel_ip)

      if(allocated(Khz_meso_last_step_sp))         deallocate(Khz_meso_last_step_sp)
      if(allocated(Khz_meso_next_step_sp))         deallocate(Khz_meso_next_step_sp)
      if(allocated(Kv_meso_last_step_sp))          deallocate(Kv_meso_last_step_sp)
      if(allocated(Kv_meso_next_step_sp))          deallocate(Kv_meso_next_step_sp)

      if(allocated(vx_meso_last_step_MetP_sp))     deallocate(vx_meso_last_step_MetP_sp)
      if(allocated(vy_meso_last_step_MetP_sp))     deallocate(vy_meso_last_step_MetP_sp)
      if(allocated(vx_meso_next_step_MetP_sp))     deallocate(vx_meso_next_step_MetP_sp)
      if(allocated(vy_meso_next_step_MetP_sp))     deallocate(vy_meso_next_step_MetP_sp)

      end subroutine Deallocate_VarDiff_Met

!******************************************************************************
!******************************************************************************

      subroutine Eddy_diff

      use global_param,  only : &
         HR_2_S

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet

      implicit none

      integer :: i,j,k

      real(kind=ip) :: E11,E12,E21,E22
      real(kind=ip) :: D2_tension,D2_strain
      real(kind=ip) :: LES_TimeScale

      do i=1,nx_submet
        do j=1,ny_submet
          do k=1,np_fullmet

      ! Smagorinsky LES horizontal eddy diffusivity is proportional
      ! to sqrt((E12+E21)^2 + (E11-E22)^2) where E is the velocity gradient
      ! tensor (just in x and y)

        ! spatial derivatives of velocity (in 1/s)
      E11=du_dx_MetP_sp(i,j,k)
      E12=du_dy_MetP_sp(i,j,k)
      E21=dv_dx_MetP_sp(i,j,k)
      E22=dv_dy_MetP_sp(i,j,k)
!
      D2_strain  = (E12+E21)**2.0_ip
      D2_tension = (E11-E22)**2.0_ip

      ! Note: Costa et al (2006) use a different form.  Their "tension" is the
      ! sum of the derivatives (essentially the velocity divergence), whereas
      ! Smagorinsky uses the difference.  They also weigh this term by either
      ! 0.5 or 2.0.
        ! in units of 1/s
      LES_TimeScale = sqrt(D2_tension+D2_strain)
        ! in units of 1/hr
      Khz_meso_next_step_MetP_sp(i,j,k) = real(LES_TimeScale*HR_2_S,kind=sp)
          enddo !k
        enddo !j
      enddo !i
      ! This eddy diffusivity will be scaled to the computational grid
      ! and multiplied by (KH_SmagC*KH_SmagC/PI/PI) * DELTA
      ! where DELTA is either dx*dy or sigma_nz
      ! The units will then be that of diffusivity (km^2/hr).
      return

      end subroutine Eddy_diff

!******************************************************************************

      subroutine Calc_Vert_Diff(last_or_next)

      use global_param,  only : &
         KM_2_M,HR_2_S,KM_2_M

      use MetReader,     only : &
        nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,MR_geoH_metP_next

      implicit none

      integer, intent(in) :: last_or_next

      integer :: i,j,k
      real(kind=ip) :: Ri_col(np_fullmet)
      real(kind=ip) :: z_col(np_fullmet)
      real(kind=ip) :: dv_dz_col(np_fullmet)
      real(kind=ip) :: FricVel
      real(kind=ip) :: Kv_col(np_fullmet)
      real(kind=ip) :: L_MonOb
      real(kind=ip) :: PBLz
      real(kind=ip) :: PBL_coeff
      real(kind=ip) :: PBL_profile_fac
      real(kind=ip) :: Kz_tmp
      real(kind=ip) :: Lc

      do i=1,nx_submet
        do j=1,ny_submet
          ! For the calculations, we need:
          !  PBLz, L_MonOb, FricVel, Ri, dv_dz, and z
          If(last_or_next.eq.0)then
              Ri_col(:) = real(Ri_meso_last_step_MetP_sp(i,j,:),kind=ip)    !
               z_col(:) = real(MR_geoH_metP_last(i,j,:),kind=ip)*KM_2_M        ! m
           dv_dz_col(:) = real(dV_dz_MetP_sp(i,j,:),kind=ip)                ! 1/s
                PBLz    = real(PBLH_meso_last_step_Met_sp(i,j),kind=ip)     ! m
             L_MonOb    = real(L_MonOb_meso_last_step_Met_sp(i,j),kind=ip)  ! m
             FricVel    = real(FricVel_meso_last_step_Met_sp(i,j),kind=ip)  ! m/s
          else
              Ri_col(:) = real(Ri_meso_next_step_MetP_sp(i,j,:),kind=ip)    !
               z_col(:) = real(MR_geoH_metP_next(i,j,:),kind=ip)*KM_2_M        ! m
           dv_dz_col(:) = real(dV_dz_MetP_sp(i,j,:),kind=ip)                ! 1/s
                PBLz    = real(PBLH_meso_next_step_Met_sp(i,j),kind=ip)     ! m
             L_MonOb    = real(L_MonOb_meso_next_step_Met_sp(i,j),kind=ip)  ! m
             FricVel    = real(FricVel_meso_next_step_Met_sp(i,j),kind=ip)  ! m/s
          endif

          do k = np_fullmet,1,-1
            ! Determine which form of Kz based on height relative to
            ! atmospheric boundary layer
            if(z_col(k).le.0.0_sp)then
                ! If point is at a nengative gpm, then assign the kz from the
                ! node above
              Kz_tmp = Kv_col(k+1)
            elseif(z_col(k)*KM_2_M.lt.PBLz)then
              ! Within the PBL, use similarity theory
                ! Parabolic profile factor for Kv between 0 and PBL
              PBL_profile_fac = z_col(k)*(1.0_sp-z_col(k)/PBLz)

              !PBL_coeff = PBL_Similarity_Kansas(z_on_L)
              !PBL_coeff = PBL_Similarity_Kansas_Ri(Ri_col_windp(k))
              !PBL_coeff = PBL_Similarity_Ulke(z_col(k),PBLz,L_MonOb)
              PBL_coeff = Phi_WindShear_NonDim(z_col(k),PBLz,L_MonOb, &
                                    Coeff_Stab,Coeff_UnStab,Exp_UnStab)
    
              ! Kz from similarity theory (Eq. 8.48 of Jacobson)
              Kz_tmp = z_col(k)*vonKarman*FricVel*PBL_profile_fac/PBL_coeff
            else
    
                ! In free atmosphere above the PBL, use Prandtl's mixing
                ! length theory for thermally stratified atmosphere.
                ! This is what Costa et al, 2006 uses (Eq. 8)
                !  First get mixing length scale
                !    There are several ways to parameterize the mixing
                !    length (Randerson, p155, 1984; Monin and Yaglom, v1,
                !    p409. Collins et al, NCAR TN-464, 2004, eq. 4.461)
              Lc = MixLen_CAM3(real(z_col(k),kind=ip))
    
                ! calculate eq 8
                ! The Ri-term seems to zero out anything above the PBL
                ! since Ri is too high
              Kz_tmp = Lc*Lc*abs(dv_dz_col(k))!*Fc(Ri_col_windp(k))
            endif
    
            ! assign to array and convert from m2/s to km2/hr
            Kv_col(k) = Kz_tmp * HR_2_S/KM_2_M/KM_2_M

            !if(i.eq.3.and.j.eq.15)then
            !  write(55,*)z_col(k),Ri_col(k),PBLz,L_MonOb,Kz_tmp
            !endif
          enddo

          If(last_or_next.eq.0)then
            Kv_meso_last_step_MetP_sp(i,j,:) = real(Kv_col(:),kind=sp)
          else
            Kv_meso_next_step_MetP_sp(i,j,:) = real(Kv_col(:),kind=sp)
          endif
        enddo
      enddo

      return

      end subroutine Calc_Vert_Diff

!******************************************************************************

      subroutine Set_VarDiffH_Meso(Load_MesoSteps,Interval_Frac)

      use mesh,          only : &
         nxmax,nymax,nzmax,dx,dy,IsLatLon,sigma_nz_pd

      use Diffusion,     only : &
         kx,ky

      use MetReader,     only : &
         MR_dum3d_compH,MR_vx_metP_last,MR_dum3d_metP,MR_dum3d2_metP,MR_iMetStep_Now,&
         MR_vy_metP_last,MR_vy_metP_next,MR_vx_metP_next,&
           MR_DelMetP_Dx,&
           MR_DelMetP_Dy,&
           MR_Regrid_MetP_to_CompH

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=ip),intent(in) :: Interval_Frac

      logical,save :: first_time = .true.

      if(Load_MesoSteps)then
        if(first_time)then
          ! Need to fill _last_step_sp
          !  First fill next step so that outside this 'first_time' loop, the
          !  'next' can be copied to the 'last'
          ! Load U winds on MetP
          !ivar = 2 ! U winds
          !call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
          !vx_meso_next_step_MetP_sp = MR_dum3d_metP
          vx_meso_next_step_MetP_sp = MR_vx_metP_last
          MR_dum3d_metP             = MR_vx_metP_last
            ! Now differentiate
          call MR_DelMetP_Dx
          du_dx_MetP_sp = MR_dum3d2_metP
          call MR_DelMetP_Dy
          du_dy_MetP_sp = MR_dum3d2_metP
          ! Load V winds on MetP
          !ivar = 3 ! V winds
          !call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
          !vy_meso_next_step_MetP_sp = MR_dum3d_metP
          vy_meso_next_step_MetP_sp = MR_vy_metP_last
          MR_dum3d_metP             = MR_vy_metP_last
            ! Now differentiate
          call MR_DelMetP_Dx
          dv_dx_MetP_sp = MR_dum3d2_metP
          call MR_DelMetP_Dy
          dv_dy_MetP_sp = MR_dum3d2_metP
          call Eddy_diff
           ! Now resample onto computational grid
          MR_dum3d_metP = Khz_meso_next_step_MetP_sp
          call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now)
          Khz_meso_next_step_sp = MR_dum3d_compH
          if(IsLatLon) then
            Khz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) = &
             Khz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) * &
                     real(LES_L2ScaleCoeff * sigma_nz_pd(1:nxmax,1:nymax,1:nzmax),kind=sp)
          else
            Khz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) = &
             Khz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) * &
                     real(LES_L2ScaleCoeff * dx * dy,kind=sp)
          endif
          first_time = .false.
        endif ! first_time
        Khz_meso_last_step_MetP_sp = Khz_meso_next_step_MetP_sp
        Khz_meso_last_step_sp      = Khz_meso_next_step_sp
        vx_meso_last_step_MetP_sp  = vx_meso_next_step_MetP_sp
        vy_meso_last_step_MetP_sp  = vy_meso_next_step_MetP_sp

        ! Need to fill _next_step_sp
        ! Load U winds on MetP
        !ivar = 2 ! U winds
        !call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
        vx_meso_next_step_MetP_sp = MR_vx_metP_next
        MR_dum3d_metP             = MR_vx_metP_next
          ! Now differentiate
        call MR_DelMetP_Dx
        du_dx_MetP_sp = MR_dum3d2_metP
        call MR_DelMetP_Dy
        du_dy_MetP_sp = MR_dum3d2_metP
        ! Load V winds on MetP
        !ivar = 3 ! V winds
        !call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
        vy_meso_next_step_MetP_sp = MR_vy_metP_next
        MR_dum3d_metP             = MR_vy_metP_next
          ! Now differentiate
        call MR_DelMetP_Dx
        dv_dx_MetP_sp = MR_dum3d2_metP
        call MR_DelMetP_Dy
        dv_dy_MetP_sp = MR_dum3d2_metP
        call Eddy_diff
         ! Now resample onto computational grid
        MR_dum3d_metP = Khz_meso_next_step_MetP_sp
        call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now+1)
        Khz_meso_next_step_sp = MR_dum3d_compH
        if(IsLatLon) then
          Khz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) = &
           Khz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) * &
                   real(LES_L2ScaleCoeff * sigma_nz_pd(1:nxmax,1:nymax,1:nzmax),kind=sp)
        else
          Khz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) = &
           Khz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) * &
                   real(LES_L2ScaleCoeff * dx * dy,kind=sp)
        endif

      endif

        kx(1:nxmax,1:nymax,1:nzmax) = real( Khz_meso_last_step_sp(:,:,:),kind=ip) + &
                                      real((Khz_meso_next_step_sp(:,:,:) - &
                                            Khz_meso_last_step_sp(:,:,:)),kind=ip) * &
                                                Interval_Frac
        ky = kx

        ! Set boundary kx and ky
          ! Bottom
        kx(0:nxmax+1,0:nymax+1,0) = kx(0:nxmax+1,0:nymax+1,1)
        ky(0:nxmax+1,0:nymax+1,0) = ky(0:nxmax+1,0:nymax+1,1)
          ! Top
        kx(0:nxmax+1,0:nymax+1,nzmax+1) = kx(0:nxmax+1,0:nymax+1,nzmax)
        ky(0:nxmax+1,0:nymax+1,nzmax+1) = ky(0:nxmax+1,0:nymax+1,nzmax)
          ! Left (West)
        kx(0,0:nymax+1,0:nzmax+1) = kx(1,0:nymax+1,0:nzmax+1)
        ky(0,0:nymax+1,0:nzmax+1) = ky(1,0:nymax+1,0:nzmax+1)
          ! Right (East)
        kx(nxmax+1,0:nymax+1,0:nzmax+1) = kx(nxmax,0:nymax+1,0:nzmax+1)
        ky(nxmax+1,0:nymax+1,0:nzmax+1) = ky(nxmax,0:nymax+1,0:nzmax+1)
          ! -y (South)
        kx(0:nxmax+1,0,0:nzmax+1) = kx(0:nxmax+1,1,0:nzmax+1)
        ky(0:nxmax+1,0,0:nzmax+1) = ky(0:nxmax+1,1,0:nzmax+1)
          ! +y (North)
        kx(0:nxmax+1,nymax+1,0:nzmax+1) = kx(0:nxmax+1,nymax,0:nzmax+1)
        ky(0:nxmax+1,nymax+1,0:nzmax+1) = ky(0:nxmax+1,nymax,0:nzmax+1)

      end subroutine Set_VarDiffH_Meso


!******************************************************************************

      subroutine Set_VarDiffV_Meso(Load_MesoSteps,Interval_Frac)

      use mesh,          only : &
         nxmax,nymax,nzmax

      use Diffusion,     only : &
         kz

      use MetReader,     only : &
         MR_iMetStep_Now,MR_dum3d_MetP,MR_dum3d_compH,&
           MR_Regrid_MetP_to_CompH

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=ip),intent(in) :: Interval_Frac

      logical,save :: first_time = .true.

      ! To set the vertical diffusivity, we need to:
      !  1. Calculate the Richardson Number on MetP grid
      !  2. Calculate friction velocity (if not provided)
      !  3. Calculate boundary layer lengths
      !       Atmos. Boundary Layer Height (if not provided by Met file)
      !       surface layer thickness
      !       Monin-Obukhov Length
      !  4. Calculate Kv(Ri,u*,L,PBLz)
      ! We will calculate these values on the MetP grid, then interpolate Kv on
      ! the compH and then onto the current time.  
      ! Note: these are all non-linear functions so the better approach would be
      ! to evaluate everything on the computational grid at each time, but this
      ! is probably overkill.

      if(Load_MesoSteps)then
        if(first_time)then
          !  Populate values for the 'last' step
          call Calc_Ri(0)
          call Calc_SurfaceRoughnessLength
          call Calc_SurfaceFrictionVelocity(0)
          call Calc_Boundary_Lengths(0)

          call Calc_Vert_Diff(0)
          MR_dum3d_MetP = Kv_meso_last_step_MetP_sp
          call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now)
          Kv_meso_last_step_sp = MR_dum3d_compH

          first_time = .false.
        else
          ! If we've already filled 'next', copy 'next' to 'last'
          Ri_meso_last_step_MetP_sp     = Ri_meso_next_step_MetP_sp
          FricVel_meso_last_step_Met_sp = FricVel_meso_next_step_Met_sp
          PBLH_meso_last_step_Met_sp    = PBLH_meso_next_step_Met_sp
          L_MonOb_meso_last_step_Met_sp = L_MonOb_meso_next_step_Met_sp
          Kv_meso_last_step_sp          = Kv_meso_next_step_sp
        endif ! first_time
          ! Populate Ri for the 'next' step
        call Calc_Ri(1)                           ! sets Ri_meso_next_step_MetP_sp
        call Calc_SurfaceFrictionVelocity(1)      ! sets FricVel_meso_next_step_Met_sp
        call Calc_Boundary_Lengths(1)             ! sets PBLH_meso_next_step_Met_sp
                                                  !  and L_MonOb_meso_next_step_Met_sp
        call Calc_Vert_Diff(1)
        MR_dum3d_MetP = Kv_meso_next_step_MetP_sp
        call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now+1)
        Kv_meso_next_step_sp = MR_dum3d_compH

      endif

      kz(1:nxmax,1:nymax,1:nzmax) = real(Kv_meso_last_step_sp(:,:,:),kind=ip) + &
                                    real((Kv_meso_next_step_sp(:,:,:) - &
                                          Kv_meso_last_step_sp(:,:,:)),kind=ip) * &
                                              Interval_Frac
      kz = kz * 3600.0_ip
      !HFS KLUDGE
      kz = 5.0_ip

      ! Set boundary kz
        ! Bottom
      kz(0:nxmax+1,0:nymax+1,0) = kz(0:nxmax+1,0:nymax+1,1)
        ! Top
      kz(0:nxmax+1,0:nymax+1,nzmax+1) = kz(0:nxmax+1,0:nymax+1,nzmax)
        ! Left (West)
      kz(0,0:nymax+1,0:nzmax+1) = kz(1,0:nymax+1,0:nzmax+1)
        ! Right (East)
      kz(nxmax+1,0:nymax+1,0:nzmax+1) = kz(nxmax,0:nymax+1,0:nzmax+1)
        ! -y (South)
      kz(0:nxmax+1,0,0:nzmax+1) = kz(0:nxmax+1,1,0:nzmax+1)
        ! +y (North)
      kz(0:nxmax+1,nymax+1,0:nzmax+1) = kz(0:nxmax+1,nymax,0:nzmax+1)


      end subroutine Set_VarDiffV_Meso

!******************************************************************************
!******************************************************************************

      subroutine Calc_Ri(last_or_next)

      use global_param,  only : &
         GRAV,KM_2_M,MPS_2_KMPHR,useMoistureVars,EPS_SMALL

      use Atmosphere,    only : &
         AirSH_meso_last_step_MetP_sp,AirSH_meso_next_step_MetP_sp,&
         AirTemp_meso_last_step_MetP_sp,AirTemp_meso_next_step_MetP_sp,&
         R_GAS_DRYAIR,CP_AIR

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,p_fullmet_sp,&
         MR_geoH_metP_last,MR_geoH_MetP_next

      implicit none

      integer, intent(in) :: last_or_next

      real(kind=ip),dimension(:),allocatable :: z
      real(kind=ip),dimension(:),allocatable :: u
      real(kind=ip),dimension(:),allocatable :: v
      real(kind=ip),dimension(:),allocatable :: p
      real(kind=ip),dimension(:),allocatable :: T
      real(kind=ip),dimension(:),allocatable :: Q
      !real(kind=sp) :: dv_dz(np_fullmet)

      integer :: i,j,k
      real(kind=ip) :: refP
      real(kind=ip) :: del_z
      real(kind=ip) :: ptemp1,ptemp2,ptemp,delptemp
      real(kind=ip) :: u1,u2,dv1
      real(kind=ip) :: temp_term,mech_term

      allocate(z(np_fullmet))
      allocate(u(np_fullmet))
      allocate(v(np_fullmet))
      allocate(p(np_fullmet))
      allocate(T(np_fullmet))
      allocate(Q(np_fullmet))


      refP = 1.0e5_ip   ! reference pressure for potential temperature

      p = p_fullmet_sp
      do i=1,nx_submet
        do j=1,ny_submet
          if(last_or_next.eq.0)then
            z = MR_geoH_metP_last(i,j,1:np_fullmet) * KM_2_M
            u = vx_meso_last_step_MetP_sp(i,j,1:np_fullmet)/MPS_2_KMPHR
            v = vy_meso_last_step_MetP_sp(i,j,1:np_fullmet)/MPS_2_KMPHR
            T = AirTemp_meso_last_step_MetP_sp(i,j,1:np_fullmet)
            if(useMoistureVars)then
                ! If moisture is enabled, use virtual potential temperatrue
              Q = AirSH_meso_last_step_MetP_sp(i,j,:)
            else
                ! Otherwise, we will just use potential temperature
              Q = 0.0_ip
            endif
          else
            z = MR_geoH_MetP_next(i,j,1:np_fullmet) * KM_2_M
            u = vx_meso_next_step_MetP_sp(i,j,1:np_fullmet)/MPS_2_KMPHR
            v = vy_meso_next_step_MetP_sp(i,j,1:np_fullmet)/MPS_2_KMPHR
            T = AirTemp_meso_next_step_MetP_sp(i,j,1:np_fullmet)
            if(useMoistureVars)then
              Q = AirSH_meso_next_step_MetP_sp(i,j,1:np_fullmet)
            else
              Q = 0.0_ip
            endif
          endif

          do k=1,np_fullmet
            ! We need vertical derivatives of theta_v and u
            !  Use one-sided differences at the top and bottom, two-sided
            !  elsewhere
            If(k.eq.1)then
              del_z  = (z(k+1)- z(k))
              ptemp1 = T(k)*(refP/p(k))**(R_GAS_DRYAIR/CP_AIR)
              ptemp1 = ptemp1 * (1.0_ip + 0.608_ip*Q(k))
              ptemp2 = T(k+1)*(refP/p(k+1))**(R_GAS_DRYAIR/CP_AIR)
              ptemp2 = ptemp2 * (1.0_ip + 0.608_ip*Q(k+1))
              ptemp  = 0.5_ip*(ptemp1+ptemp2)
              u1     = sqrt(u(k  )*u(k  ) + v(k  )*v(k  ))
              u2     = sqrt(u(k+1)*u(k+1) + v(k+1)*v(k+1))
              dv1    = u2-u1
            elseif(k.eq.np_fullmet)then
              del_z  = (z(k  )- z(k-1))
              ptemp1 = T(k-1)*(refP/p(k-1))**(R_GAS_DRYAIR/CP_AIR)
              ptemp1 = ptemp1 * (1.0_ip + 0.608_ip*Q(k-1))
              ptemp2 = T(k  )*(refP/p(k  ))**(R_GAS_DRYAIR/CP_AIR)
              ptemp2 = ptemp2 * (1.0_ip + 0.608_ip*Q(k))
              ptemp  = 0.5_ip*(ptemp1+ptemp2)
              u1     = sqrt(u(k-1)*u(k-1) + v(k-1)*v(k-1))
              u2     = sqrt(u(k  )*u(k  ) + v(k  )*v(k  ))
              dv1    = u2-u1
            else
              del_z  = (z(k+1)- z(k-1))
              ptemp1 = T(k-1)*(refP/p(k-1))**(R_GAS_DRYAIR/CP_AIR)
              ptemp1 = ptemp1 * (1.0_ip + 0.608_ip*Q(k-1))
              ptemp2 = T(k+1)*(refP/p(k+1))**(R_GAS_DRYAIR/CP_AIR)
              ptemp2 = ptemp2 * (1.0_ip + 0.608_ip*Q(k+1))
              ptemp  = 0.5_ip*(ptemp1+ptemp2)
              u1     = sqrt(u(k-1)*u(k-1) + v(k-1)*v(k-1))
              u2     = sqrt(u(k+1)*u(k+1) + v(k+1)*v(k+1))
              dv1    = u2-u1
            endif
            if(abs(dv1).lt.EPS_SMALL) dv1=EPS_SMALL
  
            delptemp = ptemp2-ptemp1
            temp_term = (1.0_ip/ptemp)*(delptemp/del_z)
            dV_dz_MetP_sp(i,j,k) = real(dv1/del_z,kind=sp)
            ! When comparing the Ri calculation below with that from
            ! MERRA, it seems that the magnitude of dv_dz is at least
            ! 3.0e-3 m/s.  Smaller values cause Ri to become singular.
            ! In fact, they may assume a dv_dz = 3.0e-3 m/s
            mech_term = max(real(dV_dz_MetP_sp(i,j,k),kind=ip),3.0e-3_ip)**2.0_ip / &
                            GRAV
            if(last_or_next.eq.0)then
              Ri_meso_last_step_MetP_sp(i,j,k) = real(temp_term / mech_term,kind=sp)
            else
              Ri_meso_next_step_MetP_sp(i,j,k) = real(temp_term / mech_term,kind=sp)
            endif
          enddo
        enddo
      enddo

      end subroutine Calc_Ri

!******************************************************************************

      subroutine Calc_SurfaceRoughnessLength

      use MetReader,     only : &
         Met_var_IsAvailable,MR_iMetStep_Now,MR_dum2d_Met,&
           MR_Read_2d_Met_Variable

      integer :: ivar

      ! Check if the windfile being used provides surface roughness
      ivar = 17 ! Surface_roughness_surface
      if(Met_var_IsAvailable(ivar))then
        ! Surface roughness is provided, read it from the met file
        call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
        SurfRoughLen_Met_sp  = MR_dum2d_Met

      !elseif(useLandCover)then
      !  ! Set SurfRoughLen_Met_sp from Land use classification

      !  if(LandCover_Format.eq.1)then
      !    ! 1 degree resolution
      !    if(LC_grid)!0       Water
      !    !1       Broadleaf Evergreen Forest
      !    !2       Coniferous Evergreen Forest and Woodland
      !    !3       High Latitude Deciduous Forest and Woodland
      !    !4       Tundra
      !    !5       Mixed Coniferous Forest and Woodland
      !    !6       Wooded Grassland
      !    !7       Grassland
      !    !8       Bare Ground
      !    !9       Shrubs and Bare Ground
      !    !10      Cultivated Crops
      !    !11      Broadleaf Deciduous Forest and Woodland
      !    !12      Data Unavailable
      !  elseif(LandCover_Format.eq.2)then
      !    ! 8 km resolution
      !    !0       Water (and Goode's interrupted space)
      !    !1       Evergreen Needleleaf Forest
      !    !2       Evergreen Broadleaf Foreset
      !    !3       Deciduous Needleleaf Forest
      !    !4       Deciduous Broadleaf Forest
      !    !5       Mixed Forest
      !    !6       Woodland
      !    !7       Wooded Grassland
      !    !8       Closed Shrubland
      !    !9       Open Shrubland
      !    !10      Grassland
      !    !11      Cropland
      !    !12      Bare Ground
      !    !13      Urban and Built-up
      !  elseif(LandCover_Format.eq.3)then
      !    ! 1 km resolution
      !    !0       Water (and Goode's interrupted space)
      !    !1       Evergreen Needleleaf Forest
      !    !2       Evergreen Broadleaf Foreset
      !    !3       Deciduous Needleleaf Forest
      !    !4       Deciduous Broadleaf Forest
      !    !5       Mixed Forest
      !    !6       Woodland
      !    !7       Wooded Grassland
      !    !8       Closed Shrubland
      !    !9       Open Shrubland
      !    !10      Grassland
      !    !11      Cropland
      !    !12      Bare Ground
      !    !13      Urban and Built-up
      !  endif
      else
        ! Set SurfRoughLen_Met_sp by assumption
          SurfRoughLen_Met_sp  = 0.1_sp
      endif

      end subroutine Calc_SurfaceRoughnessLength

!******************************************************************************

      subroutine Calc_SurfaceFrictionVelocity(last_or_next)

      use global_param,  only : &
         KM_2_M,MPS_2_KMPHR

      use MetReader,     only : &
         Met_var_IsAvailable,MR_iMetStep_Now,np_fullmet,MR_iMetStep_Now,&
         MR_geoH_metP_last,MR_geoH_metP_next,MR_dum2d_Met,nx_submet,ny_submet,&
           MR_Read_2d_Met_Variable

      implicit none

      integer, intent(in) :: last_or_next

      real(kind=ip) :: U_mag,denom,z0
      integer :: i,j,k
      integer :: ivar

      ! Check if the windfile being used provides friction velocity
      ivar = 13 ! Friction_velocity_surface
      if(Met_var_IsAvailable(ivar))then
        ! Friction velocity is provided, read it from the met file
        if(last_or_next.eq.0)then
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
        else
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now+1)
        endif
        FricVel_meso_next_step_Met_sp = MR_dum2d_Met
      else
        ! friction velocity is not provided by the met file
        ! Calculate it ourselves
        ! Get surface friction velocity using Panofsky/Dutton p376
        do i=1,nx_submet
          do j=1,ny_submet
            z0 = SurfRoughLen_Met_sp(i,j)
            if(last_or_next.eq.0)then
              do k=1,np_fullmet
                if(MR_geoH_metP_last(i,j,k).gt.z0)then
                  exit
                endif
              enddo
              U_mag = sqrt(vx_meso_last_step_MetP_sp(i,j,k)**2.0_ip + &
                           vy_meso_last_step_MetP_sp(i,j,k)**2.0_ip) / MPS_2_KMPHR
              denom = log(MR_geoH_metP_last(i,j,k)*KM_2_M/z0)
              FricVel_meso_last_step_Met_sp(i,j) = real(U_mag*vonKarman/denom,kind=sp)
            else
              do k=1,np_fullmet
                if(MR_geoH_metP_next(i,j,k).gt.z0)then
                  exit
                endif
              enddo
              U_mag = sqrt(vx_meso_next_step_MetP_sp(i,j,k)**2.0_ip + &
                           vy_meso_next_step_MetP_sp(i,j,k)**2.0_ip) / MPS_2_KMPHR
              denom = log(MR_geoH_metP_next(i,j,k)*KM_2_M/z0)
              FricVel_meso_next_step_Met_sp(i,j) = real(U_mag*vonKarman/denom,kind=sp)
            endif
          enddo
        enddo
      endif

      end subroutine Calc_SurfaceFrictionVelocity

!******************************************************************************

      subroutine Calc_Boundary_Lengths(last_or_next)

      use global_param,  only : &
         EPS_SMALL,KM_2_M

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,Met_var_IsAvailable,&
         MR_geoH_metP_last,MR_geoH_metP_next,MR_dum2d_Met,MR_iMetStep_Now,&
           MR_Read_2d_Met_Variable

      implicit none

      integer, intent(in) :: last_or_next

      integer :: ivar
      integer :: i,j,k,k_L
      real(kind=ip) :: denom,tmp
      real(kind=ip) :: Ri_col(np_fullmet)
      real(kind=ip) :: z_col(np_fullmet)
      real(kind=ip) :: L_MonOb
      real(kind=ip) :: PBLz

      ! Check if the windfile being used provides PBLH
      ivar = 10 ! Planetary Boundary Level Height
      if(Met_var_IsAvailable(ivar))then
        ! PBLH is provided, read it from the met file
        if(last_or_next.eq.0)then
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
            PBLH_meso_last_step_Met_sp = MR_dum2d_Met
        else
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now+1)
            PBLH_meso_next_step_Met_sp = MR_dum2d_Met
        endif

      else
        ! We need to calculate PBLH internally.  There are many more involved
        ! methods of determining the PBL, but using Ri and Ri_crit seems to
        ! work fairly well, although it tends to calculate thiner PBLH than
        ! NARR and MERRA report.  Sometimes the difference is dramatic.
          ! Now loop back through and find where Ri exceeds Ri_crit (~0.3)
          ! Troen and Mahrt, Boundary-Layer Meteorology, v37, p129-148, 1986.
        do i=1,nx_submet
          do j=1,ny_submet
            If(last_or_next.eq.0)then
              Ri_col(:) = Ri_meso_last_step_MetP_sp(i,j,:)
              z_col(:)  = MR_geoH_metP_last(i,j,:)*KM_2_M
            else
              Ri_col(:) = Ri_meso_next_step_MetP_sp(i,j,:)
              z_col(:)  = MR_geoH_metP_next(i,j,:)*KM_2_M
            endif
              ! Initialize boundary layer height to sea level
            PBLz = EPS_SMALL
            do k = 2,np_fullmet-1
              if(Ri_col(k).gt.RI_CRIT.and.Ri_col(k-1).le.RI_CRIT)then
                ! This height is above the PBL; interpolate back to
                ! k-1 to get PBLz
                if(abs(Ri_col(k)-Ri_col(k-1)).lt.EPS_SMALL)cycle
                denom = (Ri_col(k)-Ri_col(k-1))
                if(abs(denom).lt.1.0e-3_ip)then
                  tmp = 1.0_ip
                else
                  tmp = (RI_CRIT-Ri_col(k-1)) / denom
                  tmp = min(tmp,1.0_ip)
                endif
                PBLz = z_col(k-1)+tmp*(z_col(k)-z_col(k-1))
                ! if we have a PBLz, then exit the do loop
                exit
              endif
            enddo
    
            ! Make sure that PBLz is not negative
            PBLz = max(PBLz,EPS_SMALL)
            if(last_or_next.eq.0)then
              PBLH_meso_last_step_Met_sp(i,j) = real(PBLz,kind=sp)
            else
              PBLH_meso_next_step_Met_sp(i,j) = real(PBLz,kind=sp)
            endif
          enddo
        enddo
      endif

      ! The surface layer is typically 10% of the planetary boundary layer.
      ! This is the region where Monin-Obukhov theory applies.
      !SurfLayerThick = 0.1_ip*PBLz

      ! Get Monin-Obukhov length from the
      ! Businger-Dyer-Pandolfo empirical result 
      ! using z and Ri at k=2
        ! Eq 6.7.1 and 6.7 2 of "Atmospheric Turbulence";
        ! Panofsky and Dutton,1984
        ! also Eq 11.24 of "Introduction to Micrometeorology";
        ! Arya, 1988
      do i=1,nx_submet
        do j=1,ny_submet
          If(last_or_next.eq.0)then
            Ri_col(:) = Ri_meso_last_step_MetP_sp(i,j,:)
            z_col(:)      = MR_geoH_metP_last(i,j,:)*KM_2_M
          else
            Ri_col(:) = Ri_meso_next_step_MetP_sp(i,j,:)
            z_col(:)      = MR_geoH_metP_next(i,j,:)*KM_2_M
          endif

          ! Pick the bottom (non-zero) z
          do k_L=1,np_fullmet
            if(z_col(k_L).gt.0.0_ip)then
              exit
            endif
          enddo
          if(abs(Ri_col(k_L)).lt.EPS_SMALL)then
              ! For the neutrally stable case, set L to 0
            L_MonOb = 0.0_ip
          elseif(Ri_col(k_L).lt.0.0_ip)then
              ! Unstable (negative L)
            L_MonOb = z_col(k_L)/Ri_col(k_L)
          elseif(Ri_col(k_L).lt.RI_CRIT)then
              ! Stable (positive L)
            L_MonOb = z_col(k_L)/Ri_col(k_L) &
                       * (1.0_ip-Ri_col(k_L)/RI_CRIT)
          else
            L_MonOb = -1.0e-3_ip
          endif
          if(last_or_next.eq.0)then
            L_MonOb_meso_last_step_Met_sp(i,j) = real(L_MonOb,kind=sp)
          else
            L_MonOb_meso_next_step_Met_sp(i,j) = real(L_MonOb,kind=sp)
          endif
        enddo
      enddo

      end subroutine Calc_Boundary_Lengths
!
!!******************************************************************************

      function Fc(Ri)
      ! Stability function for verticle diffusion
      ! see Costa et al, EPSL v241, p634, 2006; eq 10
      ! Originally from Collins et al, NCAR TN-464, 2004
      ! http://www.cesm.ucar.edu/models/atm-cam/docs/description/description.pdf


      implicit none

      real(kind=ip) :: Fc
      real(kind=ip) :: Ri

      if(Ri.ge.0.0_ip)then
          ! Eq. 4.465
        Fc = 1.0_ip/(1.0_ip+10.0_ip*Ri*(1.0_ip+8.0_ip*Ri))
      else
          ! Eq. 4.464
        Fc = sqrt(1.0_ip-18.0_ip*Ri)
      endif

      return

      end function Fc

!!******************************************************************************

      function Fc_Piedelievre(Ri,rho,z)
      ! Stability function for verticle diffusion
      ! see Piedelievre, Jean Philippe, Lue Musson-Genon, François Bompay, 1990:
      ! MEDIA—An Eulerian Model of Atmospheric Dispersion: First Validation on
      ! the Chernobyl Release. J. Appl. Meteor., 29, 1205–1220.
      ! doi: http://dx.doi.org/10.1175/1520-0450(1990)029<1205:MEMOAD>2.0.CO;2 
      ! This is the model used by MDLP0
      ! 

      implicit none

      real(kind=ip) :: Fc_Piedelievre
      real(kind=ip) :: Ri,rho,z

      if(Ri.ge.0.0_ip)then
        Fc_Piedelievre = 1.0_ip/(1.0_ip+15.0_ip*Ri*sqrt(1.0_ip+5.0_ip*Ri) )
      else
        Fc_Piedelievre = 1.0_ip/(1.0_ip+75.0_ip*(rho*rho*sqrt(Ri)/(z*z*5.19615242270663_ip)))
      endif

      return

      end function Fc_Piedelievre


!!******************************************************************************

      function ABL_Similarity_Kansas(z_on_L)

      ! Similarity function for boundary layer from the 1968 Kansas Field
      ! Program
      ! This function returns 1/phi where phi is given by Eq. 8.29 of Jacobson
      ! ( Eq 11.6 of Arya, Eq. 16.75 of Seinfeld and Pandis
      ! As noted in Costa et al, 2006, this is only really valid for the surface
      ! layer, which is about 0.1*h where h is the thickness of the planetary or
      ! atmospheric boundary layer

      implicit none

      real(kind=ip) :: ABL_Similarity_Kansas
      real(kind=ip) :: z_on_L
      real(kind=ip) :: beta,gammam

      ! if vonKarman = 0.35
      ! beta   = 4.7
      ! gammam = 15.0

      ! if vonKarman = 0.4
      beta   = 6.0_ip
      gammam = 19.3_ip

      if(z_on_L.ge.0.0_ip)then
        ! Stable
        ABL_Similarity_Kansas = 1.0_ip/(1.0_ip+beta*z_on_L)
      else
        ! Unstable
        ABL_Similarity_Kansas = (1.0_ip-gammam*z_on_L)**0.25_ip
      endif

      return

      end function ABL_Similarity_Kansas

!!******************************************************************************

      function ABL_Similarity_Kansas_Ri(Ri)

      ! Similarity function for boundary layer from the 1968 Kansas Field
      ! Program, but expressed in terms of the Richardson number
      ! Eq 11.11 of Arya
      ! As noted in Costa et al, 2006, this is only really valid for the surface
      ! layer, which is about 0.1*h where h is the thickness of the planetary or
      ! atmospheric boundary layer

      implicit none

      real(kind=ip) :: ABL_Similarity_Kansas_Ri
      real(kind=ip)::  Ri

      if(Ri.lt.0.0_ip)then
        ! Unstable
        ABL_Similarity_Kansas_Ri = (1.0_ip-15.0_ip*Ri)**0.25_ip
        ! Note: Kramm et al, Precip Scav and Atmos Surf Exch v.2 p1125-1141,1992
        ! in eq 19, use 16.0*Ri instead of 15.0*Ri
      elseif(Ri.le.RI_CRIT)then
        ! Unstable
        ABL_Similarity_Kansas_Ri = 1.0_ip-5.0_ip*Ri
      else
      !  ! shouldn't be calculating boundary layer diffusivity for Ri this high
        ABL_Similarity_Kansas_Ri = 0.0_ip
      endif

      return

      end function ABL_Similarity_Kansas_Ri

!!******************************************************************************

      function Phi_WindShear_NonDim(z,h,L,CStab,CUnStab,EUnStab)

      use global_param,  only : &
         EPS_SMALL

      ! Similarity function for full atmospheric boundary layer given by Ulke,
      ! Atmospheric Environment, v34 p1029-1042, 2000. Eq. 4a and 4b; 5a and 5b
      ! Note: Costa et al, 2006, cite this parameterization

      implicit none

      real(kind=ip) :: Phi_WindShear_NonDim
      real(kind=ip) :: z,h,L

      real(kind=ip) :: z_on_h,z_on_L,h_on_L
      real(kind=ip) :: CStab
      real(kind=ip) :: CUnStab
      real(kind=ip) :: EUnStab

      if(abs(h).lt.EPS_SMALL.or.abs(L).lt.EPS_SMALL)then
        Phi_WindShear_NonDim = EPS_SMALL
      elseif(z.le.0.0_ip)then
        Phi_WindShear_NonDim = EPS_SMALL
      else
        z_on_h = z/h
        z_on_L = z/L
        h_on_L = h/L

        if(h_on_L.ge.0.0_ip)then
          ! Stable
          Phi_WindShear_NonDim = (1.0_ip+CStab*z_on_L)
        else
          ! Unstable
          Phi_WindShear_NonDim = (1.0_ip+CUnStab*z_on_L)**EUnStab
        endif
      endif

      return

      end function Phi_WindShear_NonDim


!!******************************************************************************

      function MixLen_CAM3(z)

      ! Returns the mixing length for Prandtl's turbulent diffusion
      ! See Collins et al, NCAR TN-464, 2004, eq. 4.461
      ! http://www.cesm.ucar.edu/models/atm-cam/docs/description/description.pdf
      ! This is originally from Blackadar (1962)
      ! Louis (2000) found LambdaC to be 100m which Collins used 30m

      implicit none

      real(kind=ip) :: MixLen_CAM3
      real(kind=ip) :: z

      MixLen_CAM3 = 1.0_ip/(1.0_ip/(z*vonKarman) + 1.0_ip/LambdaC)

      return

      end function MixLen_CAM3

!******************************************************************************

      end module Variable_Diffusivity
