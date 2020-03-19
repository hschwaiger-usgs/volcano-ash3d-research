      module Wet_Deposition

      use precis_param

      use io_units

      ! This module calculates the removal of particles in the atmosphere
      ! through the interaction with hydrometeors.  What is needed on the
      ! computational grid are:
      !   Scavenging Coefficient (Liquid) : (nx,ny,nz,ns) :: describes the decay rate of concen
      !   Scavenging Coefficient (Ice)    : (nx,ny,nz,ns) ::  ''
      !   Height of bottom of cloud       : (nx,ny)       :: used for in-cloud v.s below-cloud scavenging
      !   Height of top of cloud          : (nx,ny)       ::  ''
      !
      ! To get cloud position, we can either estimate it from the specific or
      ! relative humidity (each on MetP), or using cloud top/bottom pressure
      ! variables if the windfile provides them.  The method used is specified
      ! by the variable CloudLoc.
      !
      ! To get scavenging coefficients, we simply need the precipitation rate on
      ! the Met grid.  Most precipitation rates are given as 2d surface value.
      ! The 3d field of downward flux is inferred to be constant from the
      ! surface to the cloud bottom, then decrease linearly to the cloud top.
      ! This rate is assigned to either liquid or ice precip depending on
      ! catagorical values (if available).  NASA/MERRA windfiles give the full
      ! 3d field for both liquid and ice.  We then assume a monodisperse rain
      ! drop distribution with a drop size optimally matched to the collection
      ! efficiency of a Gamma drop distribution, calculate CEffic for each
      ! grainsize for this drop size, then calculate the scavenging coefficient.
      ! This is all done in MetP space.  The scavenging coefficient is interpolated
      ! onto the computational grid and onto the current time step.
      !
      ! At each time step, the scavenging is applied by calculating an
      ! exponential decay of concen using the scavenging coefficient.  Because
      ! this will typically be a very stiff system, sub-stepping is done to
      ! integrate over the whole time step.

      real(kind=ip),parameter :: Rain_dens    = 1.0e3_ip  ! density of water in kg/m3
      real(kind=ip),parameter :: Rain_Mu      = 8.9e-4_ip  ! viscosity of water in kg/m/s

      logical :: USE_WETDEPO     = .false.
      logical :: USE_WASHOUT_LIQ = .false.   ! Below-cloud rain
      logical :: USE_RAINOUT_LIQ = .false.  ! In-cloud drop nucleation
      logical :: USE_WASHOUT_ICE = .false.  ! Below-cloud snow
      logical :: USE_RAINOUT_ICE = .false.  ! In-cloud ice nucleation
      integer :: CloudLoc      ! Flag specifying the method for locating cloud position
                               !  0 = Cloud Bottom is set as top of model
                               !      (i.e. activate all comp nodes for WetDepo)
                               !  1 = via Specific Humidity
                               !  2 = via Relative Humidity
                               !  3 = using cloud bot/top variables from met file

      ! Precipitaion rate
      logical            :: USE_3D_PRECIP = .false.
      integer            :: nc_nPrecip              ! 1 for 2d, >1 for 3d
      integer            :: nc_RH

      ! Set the number of output variables for this module
      integer, parameter :: nvar_User2d_static_XY_WetDepo = 0
      integer, parameter :: nvar_User2d_XY_WetDepo        = 3 ! p0
      integer, parameter :: nvar_User3d_XYGs_WetDepo      = 2 ! Dep_liqWO,Dep_liqRO,Dep_iceWO,Dep_iceRO
      integer, parameter :: nvar_User3d_XYZ_WetDepo       = 1 ! precip
      integer, parameter :: nvar_User4d_XYZGs_WetDepo     = 1 ! ScavCo_liq,ScavCo_ice

      character(len=30),dimension(nvar_User2d_static_XY_WetDepo) :: temp_2ds_name_WetDepo
      character(len=30),dimension(nvar_User2d_static_XY_WetDepo) :: temp_2ds_unit_WetDepo
      character(len=30),dimension(nvar_User2d_static_XY_WetDepo) :: temp_2ds_lname_WetDepo
      real(kind=op),    dimension(nvar_User2d_static_XY_WetDepo) :: temp_2ds_MissVal_WetDepo
      real(kind=op),    dimension(nvar_User2d_static_XY_WetDepo) :: temp_2ds_FillVal_WetDepo

      character(len=30),dimension(nvar_User2d_XY_WetDepo)    :: temp_2d_name_WetDepo
      character(len=30),dimension(nvar_User2d_XY_WetDepo)    :: temp_2d_unit_WetDepo
      character(len=30),dimension(nvar_User2d_XY_WetDepo)    :: temp_2d_lname_WetDepo
      real(kind=op),    dimension(nvar_User2d_XY_WetDepo)    :: temp_2d_MissVal_WetDepo
      real(kind=op),    dimension(nvar_User2d_XY_WetDepo)    :: temp_2d_FillVal_WetDepo

      character(len=30),dimension(nvar_User3d_XYGs_WetDepo)  :: temp_3ds_name_WetDepo
      character(len=30),dimension(nvar_User3d_XYGs_WetDepo)  :: temp_3ds_unit_WetDepo
      character(len=30),dimension(nvar_User3d_XYGs_WetDepo)  :: temp_3ds_lname_WetDepo
      real(kind=op),    dimension(nvar_User3d_XYGs_WetDepo)  :: temp_3ds_MissVal_WetDepo
      real(kind=op),    dimension(nvar_User3d_XYGs_WetDepo)  :: temp_3ds_FillVal_WetDepo

      character(len=30),dimension(nvar_User3d_XYZ_WetDepo)   :: temp_3dz_name_WetDepo
      character(len=30),dimension(nvar_User3d_XYZ_WetDepo)   :: temp_3dz_unit_WetDepo
      character(len=30),dimension(nvar_User3d_XYZ_WetDepo)   :: temp_3dz_lname_WetDepo
      real(kind=op),    dimension(nvar_User3d_XYZ_WetDepo)   :: temp_3dz_MissVal_WetDepo
      real(kind=op),    dimension(nvar_User3d_XYZ_WetDepo)   :: temp_3dz_FillVal_WetDepo

      character(len=30),dimension(nvar_User4d_XYZGs_WetDepo) :: temp_4d_name_WetDepo
      character(len=30),dimension(nvar_User4d_XYZGs_WetDepo) :: temp_4d_unit_WetDepo
      character(len=30),dimension(nvar_User4d_XYZGs_WetDepo) :: temp_4d_lname_WetDepo
      real(kind=op),    dimension(nvar_User4d_XYZGs_WetDepo) :: temp_4d_MissVal_WetDepo
      real(kind=op),    dimension(nvar_User4d_XYZGs_WetDepo) :: temp_4d_FillVal_WetDepo

      ! These are used to keep track of which index in the global list, this
      ! modules output vars corespond to
      integer :: indx_User2d_static_XY_WetDepo
      integer :: indx_User2d_XY_WetDepo
      integer :: indx_User3d_XYGs_WetDepo
      integer :: indx_User3d_XYZ_WetDepo
      integer :: indx_User4d_XYZGs_WetDepo


      ! Define the wet removal variables.
      ! These are needed on the computational grid and exist at the current time
      ! step
      real(kind=ip), dimension(:,:)    ,allocatable :: MetCloudBotHeight
      real(kind=ip), dimension(:,:)    ,allocatable :: MetCloudTopHeight
      real(kind=ip), dimension(:,:,:,:),allocatable :: scav_coeff_Liq_3d
      real(kind=ip), dimension(:,:,:,:),allocatable :: scav_coeff_Ice_3d

      ! Define Deposit variables
      ! These are only allocated if the corresponding removal process is
      ! requested
      real(kind=ip),dimension(:,:,:),allocatable :: Deposit_Liq_Washout
      real(kind=ip),dimension(:,:,:),allocatable :: Deposit_Ice_Washout
      real(kind=ip),dimension(:,:,:),allocatable :: Deposit_Liq_Rainout
      real(kind=ip),dimension(:,:,:),allocatable :: Deposit_Ice_Rainout

      real(kind=ip),dimension(:,:)    ,allocatable :: precipitation_rate_2d
      real(kind=ip),dimension(:,:,:)  ,allocatable :: precipitation_rate_3d

      ! Variables currently just needed on the wind grid      
      real(kind=sp),dimension(:,:,:) ,allocatable :: QL_MetP_sp       ! Cloud liquid water mixing ratio
      real(kind=sp),dimension(:,:,:) ,allocatable :: QI_MetP_sp       ! Cloud ice mixing ration
      real(kind=sp),dimension(:,:)   ,allocatable :: Pres_lct_MetP_sp ! Pressure at lower cloud top
      real(kind=sp),dimension(:,:)   ,allocatable :: Pres_lcb_MetP_sp ! Pressure at lower cloud bottom

      integer(kind=sp),dimension(:,:,:),allocatable :: CatagorPrecip_met
      !integer(kind=sp),dimension(:,:,:),allocatable :: CatPrecip_meso_last_step
      !integer(kind=sp),dimension(:,:,:),allocatable :: CatPrecip_meso_next_step

        ! These MetP variables are needed for determining how to calculate
        ! scavenging coefficients
      real(kind=ip),dimension(:,:)    ,allocatable :: MetCloudTop_MetP
      real(kind=ip),dimension(:,:)    ,allocatable :: MetCloudBot_MetP
        ! These CompH variables are needed for interpolating to the current
        ! time, needed for logical blocks for washout v.s. rainout
      real(kind=ip),dimension(:,:)    ,allocatable :: MetCloudTop_meso_last_step
      real(kind=ip),dimension(:,:)    ,allocatable :: MetCloudBot_meso_last_step
      real(kind=ip),dimension(:,:)    ,allocatable :: MetCloudTop_meso_next_step
      real(kind=ip),dimension(:,:)    ,allocatable :: MetCloudBot_meso_next_step

      real(kind=ip),dimension(:,:,:,:),allocatable :: scav_coeff_MetP
      real(kind=ip),dimension(:,:,:,:),allocatable :: scav_coeff_Liq_meso_last_step
      real(kind=ip),dimension(:,:,:,:),allocatable :: scav_coeff_Liq_meso_next_step
      !real(kind=ip),dimension(:,:,:,:),allocatable :: scav_coeff_Ice_meso_last_step
      !real(kind=ip),dimension(:,:,:,:),allocatable :: scav_coeff_Ice_meso_next_step


        ! These are only needed if precipitation is requested as an output
        ! variable
      real(kind=sp),dimension(:,:,:),allocatable :: prate_meso_last_step
      real(kind=sp),dimension(:,:,:),allocatable :: prate_meso_next_step
      real(kind=sp),dimension(:,:)  ,allocatable :: prate_2d_meso_last_step
      real(kind=sp),dimension(:,:)  ,allocatable :: prate_2d_meso_next_step

      real(kind=sp),dimension(:,:,:)   ,allocatable :: prate_Wat_MetP_sp ! These store the 2d or 3d precipitation
      real(kind=sp),dimension(:,:,:)   ,allocatable :: prate_Ice_MetP_sp ! on the wind grid
      real(kind=sp),dimension(:)       ,allocatable :: prate_col_MetP_sp
      real(kind=sp),dimension(:)       ,allocatable :: prate_col_metH_sp
      real(kind=sp),dimension(:,:,:)   ,allocatable :: prate_metH_sp
      real(kind=sp),dimension(:,:)     ,allocatable :: prate_2d_MetP_sp

      contains

!******************************************************************************

      subroutine input_data_WetDepo

      use global_param,  only : &
         useMoistureVars,nmods

      use io_data,       only : &
         infile

      implicit none

      character(len=3)  :: answer
      character(len=80)  :: linebuffer
      integer :: ios,ioerr
      character(len=20) :: mod_name
      integer :: substr_pos

      open(unit=10,file=infile,status='old',err=1900)

      write(global_info,*)"    Searching for OPTMOD=WETDEPO"
      nmods = 0
      read(10,'(a80)',iostat=ios)linebuffer
      do while(ios.eq.0)
        read(10,'(a80)',iostat=ios)linebuffer

        substr_pos = index(linebuffer,'OPTMOD')
        if(substr_pos.eq.1)then
          ! found an optional module
          !  Parse for the keyword
          read(linebuffer,1104)mod_name
          if(adjustl(trim(mod_name)).eq.'WETDEPO')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo

      USE_WETDEPO     = .false. 
      USE_WASHOUT_LIQ = .false.
      USE_RAINOUT_LIQ = .false.
      USE_WASHOUT_ICE = .false.
      USE_RAINOUT_ICE = .false.
      write(global_info,*)"    Continue reading input file for WetDepo block"
       ! Check if we're going to use wet deposition and which type
          ! Below-cloud rain
        read(10,'(a80)',iostat=ios,err=2010)linebuffer
        read(linebuffer,'(a3)',err=2011) answer
        if (answer.eq.'yes') then
          USE_WASHOUT_LIQ = .true.
          write(global_info,*)"    Using liquid washout (below-cloud rain)"
        elseif(answer(1:2).eq.'no') then
          USE_WASHOUT_LIQ = .false.
          write(global_info,*)"    NOT using liquid washout (below-cloud rain)"
        else
          goto 2011
        endif
          ! In-cloud drop nucleation
        read(10,'(a80)',iostat=ios,err=2010)linebuffer
        read(linebuffer,'(a3)',err=2011) answer
        if (answer.eq.'yes') then
          USE_RAINOUT_LIQ = .true.
          write(global_info,*)"    Using liquid rainout (in-cloud rain)"
        elseif(answer(1:2).eq.'no') then
          USE_RAINOUT_LIQ = .false.
          write(global_info,*)"    NOT using liquid rainout (in-cloud rain)"
        else
          goto 2011
        endif
          ! Below-cloud snow
        read(10,'(a80)',iostat=ios,err=2010)linebuffer
        read(linebuffer,'(a3)',err=2011) answer
        if (answer.eq.'yes') then
          USE_WASHOUT_ICE = .true.
          write(global_info,*)"    Using frozen washout (below-cloud snow)"
        elseif(answer(1:2).eq.'no') then
          USE_WASHOUT_ICE = .false.
          write(global_info,*)"    NOT using frozen washout (below-cloud snow)"
        else
          goto 2011
        endif
          ! In-cloud ice nucleation
        read(10,'(a80)',iostat=ios,err=2010)linebuffer
        read(linebuffer,'(a3)',err=2011) answer
        if (answer.eq.'yes') then
          USE_RAINOUT_ICE = .true.
          write(global_info,*)"    Using frozen rainout (in-cloud snow)"
        elseif(answer(1:2).eq.'no') then
          USE_RAINOUT_ICE = .false.
          write(global_info,*)"    NOT using frozen rainout (in-cloud snow)"
        else
          goto 2011
        endif

        if(USE_WASHOUT_LIQ.or.USE_RAINOUT_LIQ.or.&
           USE_WASHOUT_ICE.or.USE_RAINOUT_ICE)then
          USE_WETDEPO  = .true.
          useMoistureVars = .true.
        endif

        if (USE_WETDEPO) then
          ! We're using some form of wet deposition, then get the constants
          read(10,'(a80)',iostat=ios,err=2010)linebuffer
          read(linebuffer,*,iostat=ioerr) CloudLoc
        endif

2010  continue
      close(10)

      return

1900  write(global_info,*)  'error: cannot find input file: ',infile
      write(global_info,*)  'Program stopped'
      write(global_log,*)  'error: cannot find input file: ',infile
      write(global_log,*)  'Program stopped'
      stop 1

2011  write(global_log,*) 'Error reading whether to use wet depositiony.'
      write(global_log,*) 'Answer must be yes or no.'
      write(global_log,*) 'You gave:',linebuffer
      write(global_log,*) 'Program stopped'
      stop 1

      end subroutine input_data_WetDepo

!******************************************************************************


!******************************************************************************

      subroutine Allocate_WetDepo_global

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User2d_XY,nvar_User3d_XYGs,nvar_User3d_XYZ,&
         nvar_User4d_XYZGs

      use MetReader,     only : &
         MR_iwindformat

      implicit none

      integer :: ivar

      write(global_info,*)"--------------------------------------------------"
      write(global_info,*)"---------- ALLOCATE_WETDEPO_GLOBAL ---------------"
      write(global_info,*)"--------------------------------------------------"

      if(MR_iwindformat.eq.24)then
        USE_3D_PRECIP = .true.
      else
        USE_3D_PRECIP = .false.
      endif

       ! Cloud information is needed for distinguishing in-cloud and
       ! below-cloud processes
      allocate(MetCloudBotHeight(nxmax,nymax))
      MetCloudBotHeight = 0.0_ip
      allocate(MetCloudTopHeight(nxmax,nymax))
      MetCloudTopHeight = 0.0_ip

      ! Scavenging coefficients
      if(USE_WASHOUT_LIQ.or.USE_RAINOUT_LIQ)then
        allocate(scav_coeff_Liq_3d(nxmax,nymax,nzmax,nsmax))
        scav_coeff_Liq_3d = 0.0_ip ! scavenging coefficient
      endif
      if(USE_WASHOUT_ICE.or.USE_RAINOUT_ICE)then
        allocate(scav_coeff_Ice_3d(nxmax,nymax,nzmax,nsmax))
        scav_coeff_Ice_3d = 0.0_ip ! scavenging coefficient
      endif

      ! Deposits
      if(USE_WASHOUT_LIQ)then  ! Below-cloud rain
        allocate(Deposit_Liq_Washout(nxmax,nymax,nsmax))
        Deposit_Liq_Washout = 0.0_ip
      endif
      if(USE_WASHOUT_ICE)then  ! Below-cloud snow
        allocate(Deposit_Ice_Washout(nxmax,nymax,nsmax))
        Deposit_Ice_Washout = 0.0_ip
      endif
      if(USE_RAINOUT_LIQ)then  ! In-cloud drop nucleation
        allocate(Deposit_Liq_Rainout(nxmax,nymax,nsmax))
        Deposit_Liq_Rainout = 0.0_ip
      endif
      if(USE_RAINOUT_ICE)then  ! In-cloud ice nucleation
        allocate(Deposit_Ice_Rainout(nxmax,nymax,nsmax))
        Deposit_Ice_Rainout = 0.0_ip
      endif

      ! Precipitation rates
      if(USE_3D_PRECIP)then
!        ! MERRA uses the full 3d precipitation data
        allocate(precipitation_rate_3d(nxmax,nymax,nzmax))
        precipitation_rate_3d = 0.0_ip
      endif
!      ! This is the default for precipitation
      allocate(precipitation_rate_2d(nxmax,nymax))
      precipitation_rate_2d = 0.0_ip

      ! Set the start indecies
      indx_User2d_static_XY_WetDepo = nvar_User2d_static_XY
      indx_User2d_XY_WetDepo        = nvar_User2d_XY
      indx_User3d_XYGs_WetDepo      = nvar_User3d_XYGs
      indx_User3d_XYZ_WetDepo       = nvar_User3d_XYZ
      indx_User4d_XYZGs_WetDepo     = nvar_User4d_XYZGs

      ivar = 1
      temp_2d_name_WetDepo(ivar) = "p0"
      temp_2d_lname_WetDepo(ivar) = "Precipitation rate @ surf"
      temp_2d_unit_WetDepo(ivar) = "kg/m2/s"
      temp_2d_MissVal_WetDepo(ivar) = -9999.0_op
      temp_2d_FillVal_WetDepo(ivar) = -9999.0_op

      ivar = 2
      temp_2d_name_WetDepo(ivar) = "CloudBot_H"
      temp_2d_lname_WetDepo(ivar) = "Cloud Bottom height"
      temp_2d_unit_WetDepo(ivar) = "km"
      temp_2d_MissVal_WetDepo(ivar) = -9999.0_op
      temp_2d_FillVal_WetDepo(ivar) = -9999.0_op

      ivar = 3
      temp_2d_name_WetDepo(ivar) = "CloudTop_H"
      temp_2d_lname_WetDepo(ivar) = "Cloud Top height"
      temp_2d_unit_WetDepo(ivar) = "km"
      temp_2d_MissVal_WetDepo(ivar) = -9999.0_op
      temp_2d_FillVal_WetDepo(ivar) = -9999.0_op

      if(USE_WASHOUT_LIQ)then  ! Below-cloud rain
        ivar = 1
        temp_3ds_name_WetDepo(ivar) = "Dep_liqWO"
        temp_3ds_lname_WetDepo(ivar) = "Deposit Liquid-Washout"
        temp_3ds_unit_WetDepo(ivar) = "kg/m2"
        temp_3ds_MissVal_WetDepo(ivar) = -9999.0_op
        temp_3ds_FillVal_WetDepo(ivar) = -9999.0_op
      endif
 
      if(USE_RAINOUT_LIQ)then  ! In-cloud drop nucleation
        ivar = ivar + 1
        temp_3ds_name_WetDepo(ivar) = "Dep_liqRO"
        temp_3ds_lname_WetDepo(ivar) = "Deposit Liquid-Rainout"
        temp_3ds_unit_WetDepo(ivar) = "kg/m2"
        temp_3ds_MissVal_WetDepo(ivar) = -9999.0_op
        temp_3ds_FillVal_WetDepo(ivar) = -9999.0_op
      endif
 
      if(USE_WASHOUT_ICE)then  ! Below-cloud snow
        ivar = ivar + 1
        temp_3ds_name_WetDepo(ivar) = "Dep_iceWO"
        temp_3ds_lname_WetDepo(ivar) = "Deposit Frozen-Washout"
        temp_3ds_unit_WetDepo(ivar) = "kg/m2"
        temp_3ds_MissVal_WetDepo(ivar) = -9999.0_op
        temp_3ds_FillVal_WetDepo(ivar) = -9999.0_op
      endif
  
      if(USE_RAINOUT_ICE)then  ! In-cloud ice nucleation
        ivar = ivar + 1
        temp_3ds_name_WetDepo(ivar) = "Dep_iceRO"
        temp_3ds_lname_WetDepo(ivar) = "Deposit Frozen-Rainout"
        temp_3ds_unit_WetDepo(ivar) = "kg/m2"
        temp_3ds_MissVal_WetDepo(ivar) = -9999.0_op
        temp_3ds_FillVal_WetDepo(ivar) = -9999.0_op
      endif

      if(USE_3D_PRECIP)then
        ivar = 1
        temp_3dz_name_WetDepo(ivar) = "precip"
        temp_3dz_lname_WetDepo(ivar) = "Precipitation rate"
        temp_3dz_unit_WetDepo(ivar) = "kg/m2/s"
        temp_3dz_MissVal_WetDepo(ivar) = -9999.0_op
        temp_3dz_FillVal_WetDepo(ivar) = -9999.0_op
      endif

      ivar = 1
      temp_4d_name_WetDepo(ivar) = "ScavCo_liq"
      temp_4d_lname_WetDepo(ivar) = "Scavenging coefficient"
      temp_4d_unit_WetDepo(ivar) = "1/hr"
      temp_4d_MissVal_WetDepo(ivar) = -9999.0_op
      temp_4d_FillVal_WetDepo(ivar) = -9999.0_op

      nvar_User2d_static_XY = nvar_User2d_static_XY + nvar_User2d_static_XY_WetDepo
      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_WetDepo
      nvar_User3d_XYGs      = nvar_User3d_XYGs      + nvar_User3d_XYGs_WetDepo
      nvar_User3d_XYZ       = nvar_User3d_XYZ       + nvar_User3d_XYZ_WetDepo
      nvar_User4d_XYZGs     = nvar_User4d_XYZGs     + nvar_User4d_XYZGs_WetDepo

      end subroutine Allocate_WetDepo_global

!******************************************************************************

      subroutine Allocate_WetDepo_Met

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,np_fullmet_P0,np_fullmet_RH

      implicit none

      write(global_info,*)"--------------------------------------------------"
      write(global_info,*)"---------- ALLOCATE_WETDEPO_MET ------------------"
      write(global_info,*)"--------------------------------------------------"

      nc_RH = np_fullmet_RH

      ! lower cloud position
      allocate(MetCloudTop_MetP(nx_submet,ny_submet))
      allocate(MetCloudTop_meso_last_step(nxmax,nymax))
      allocate(MetCloudTop_meso_next_step(nxmax,nymax))
      allocate(MetCloudBot_MetP(nx_submet,ny_submet))
      allocate(MetCloudBot_meso_last_step(nxmax,nymax))
      allocate(MetCloudBot_meso_next_step(nxmax,nymax))

       !HFS  Fix this!!
      !allocate(MetCloudTop_meso_last_step(nxmax,nymax))
      !allocate(MetCloudTop_meso_next_step(nxmax,nymax))
      !allocate(MetCloudBot_meso_last_step(nxmax,nymax))
      !allocate(MetCloudBot_meso_next_step(nxmax,nymax))


      ! Scavenging coefficients
!      allocate(scav_coeff_col_MetP(np_fullmet,ngs))
!      allocate(scav_coeff_col_metH(nz,ngs))
      allocate(scav_coeff_MetP(nx_submet,ny_submet,np_fullmet,nsmax))
      allocate(scav_coeff_Liq_meso_last_step(nxmax,nymax,nzmax,nsmax))
      allocate(scav_coeff_Liq_meso_next_step(nxmax,nymax,nzmax,nsmax))

      ! Precipitation rates (required)
      allocate(prate_Wat_MetP_sp(nx_submet,ny_submet,np_fullmet_P0))
      prate_Wat_MetP_sp = 0.0_sp
!      allocate(prate_Ice_MetP_sp(nx_submet,ny_submet,np_fullmet_P0))
!      prate_Ice_MetP_sp = 0.0_sp

      ! Categorical values
      allocate(CatagorPrecip_met(nx_submet,ny_submet,4))
      !allocate(CatPrecip_meso_last_step(nxmax,nymax,4))
      !allocate(CatPrecip_meso_next_step(nxmax,nymax,4))

      ! Precipitation rates (only needed if prate is an exported variable)
!      if(USE_ADDITIONAL_VARS)then
!        if(USE_3D_PRECIP.and.log_3d2_Precip_rate)then
!          allocate(prate_col_metH_sp(nzmax))
!          allocate(prate_metH_sp(nx_fullmet_per,ny_fullmet,nzmax))
!          allocate(prate_meso_last_step(nxmax,nymax,nzmax))
!          allocate(prate_meso_next_step(nxmax,nymax,nzmax))
!        endif
!        if(log_2d2_Precip_rate)then
!          allocate(prate_2d_MetP_sp(nx_fullmet_per,ny_fullmet))
!          allocate(prate_2d_meso_last_step(nxmax,nymax))
!          allocate(prate_2d_meso_next_step(xmax,nymax))
!        endif
!      endif

!      if(MR_iwindformat.eq.4.or.MR_iwindformat.eq.24)then
!        allocate(QL_MetP_sp(nx_fullmet_per,ny_fullmet,np_fullmet_RH))
!        QL_MetP_sp = 0.0_sp
!        allocate(QI_wind_sp(nx_fullmet_per,ny_fullmet,np_fullmet_RH))
!        QI_wind_sp = 0.0_sp
!      endif

      if(CloudLoc.eq.3)then
        allocate(Pres_lct_MetP_sp(nx_submet,ny_submet))
        allocate(Pres_lcb_MetP_sp(nx_submet,ny_submet))
      endif

      end subroutine Allocate_WetDepo_Met

!******************************************************************************

      subroutine Prep_output_WetDepo

      use Output_Vars,   only : &
        var_User2d_XY_name,var_User2d_XY_unit,var_User2d_XY_lname,&
        var_User2d_XY_MissVal,var_User2d_XY_FillVal,var_User2d_XY,&
        var_User3d_XYGs_name,var_User3d_XYGs_unit,var_User3d_XYGs_lname,&
        var_User3d_XYGs_MissVal,var_User3d_XYGs_FillVal,var_User3d_XYGs,&
        var_User3d_XYZ_name,var_User3d_XYZ_unit,var_User3d_XYZ_lname,&
        var_User3d_XYZ_MissVal,var_User3d_XYZ_FillVal,var_User3d_XYZ,&
        var_User4d_XYZGs_name,var_User4d_XYZGs_unit,var_User4d_XYZGs_lname,&
        var_User4d_XYZGs_MissVal,var_User4d_XYZGs_FillVal,var_User4d_XYZGs

      use mesh,          only : &
         nxmax,nymax,nzmax

      use Tephra,        only : &
         n_gs_max

      implicit none

      integer :: i,indx
      !integer :: ivar

      !integer, parameter :: nvar_User2d_XY_WetDepo        = 1 ! precip2d
      !integer, parameter :: nvar_User3d_XYGs_WetDepo      = 4 ! Deposit (Liq/Ice)/(Wash/Rainout)
      !integer, parameter :: nvar_User3d_XYZ_WetDepo       = 1 ! precip3d
      !integer, parameter :: nvar_User4d_XYZGs_WetDepo     = 1 ! Scav Coeff

      do i=1,nvar_User2d_XY_WetDepo
        indx = indx_User2d_XY_WetDepo+i
        var_User2d_XY_name(indx)   = temp_2d_name_WetDepo(i)
        var_User2d_XY_unit(indx)   = temp_2d_unit_WetDepo(i)
        var_User2d_XY_lname(indx)  = temp_2d_lname_WetDepo(i)
        var_User2d_XY_MissVal(indx)= temp_2d_MissVal_WetDepo(i)
        var_User2d_XY_FillVal(indx)= temp_2d_FillVal_WetDepo(i)
        if(i.eq.1) var_User2d_XY(1:nxmax,1:nymax,indx) = precipitation_rate_2d(1:nxmax,1:nymax)
        if(i.eq.2) var_User2d_XY(1:nxmax,1:nymax,indx) = MetCloudBotHeight(1:nxmax,1:nymax)
        if(i.eq.3) var_User2d_XY(1:nxmax,1:nymax,indx) = MetCloudTopHeight(1:nxmax,1:nymax)
      enddo
      !write(global_info,*)precipitation_rate_2d

      do i=1,nvar_User3d_XYGs_WetDepo
        indx = indx_User3d_XYGs_WetDepo+i
        var_User3d_XYGs_name(indx)   = temp_3ds_name_WetDepo(i)
        var_User3d_XYGs_unit(indx)   = temp_3ds_unit_WetDepo(i)
        var_User3d_XYGs_lname(indx)  = temp_3ds_lname_WetDepo(i)
        var_User3d_XYGs_MissVal(indx)= temp_3ds_MissVal_WetDepo(i)
        var_User3d_XYGs_FillVal(indx)= temp_3ds_FillVal_WetDepo(i)
        if(i.eq.1.and.USE_WASHOUT_LIQ) var_User3d_XYGs(1:nxmax,1:nymax,1:n_gs_max,indx) = &
                   Deposit_Liq_Washout(1:nxmax,1:nymax,1:n_gs_max)
        if(i.eq.2.and.USE_RAINOUT_LIQ) var_User3d_XYGs(1:nxmax,1:nymax,1:n_gs_max,indx) = &
                   Deposit_Liq_Rainout(1:nxmax,1:nymax,1:n_gs_max)
        if(i.eq.3.and.USE_WASHOUT_ICE) var_User3d_XYGs(1:nxmax,1:nymax,1:n_gs_max,indx) = &
                    Deposit_Ice_Washout(1:nxmax,1:nymax,1:n_gs_max)
        if(i.eq.4.and.USE_RAINOUT_ICE) var_User3d_XYGs(1:nxmax,1:nymax,1:n_gs_max,indx) = &
                    Deposit_Ice_Rainout(1:nxmax,1:nymax,1:n_gs_max)
      enddo

      do i=1,nvar_User3d_XYZ_WetDepo
        indx = indx_User3d_XYZ_WetDepo+i
        var_User3d_XYZ_name(indx)   = temp_3dz_name_WetDepo(i)
        var_User3d_XYZ_unit(indx)   = temp_3dz_unit_WetDepo(i)
        var_User3d_XYZ_lname(indx)  = temp_3dz_lname_WetDepo(i)
        var_User3d_XYZ_MissVal(indx)= temp_3dz_MissVal_WetDepo(i)
        var_User3d_XYZ_FillVal(indx)= temp_3dz_FillVal_WetDepo(i)
        if(i.eq.1)then
          if(USE_3D_PRECIP)then
            var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = precipitation_rate_3d(1:nxmax,1:nymax,1:nzmax)
          else
            var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = 0.0
          endif
        endif
      enddo

      do i=1,nvar_User4d_XYZGs_WetDepo
        indx = indx_User4d_XYZGs_WetDepo+i
        var_User4d_XYZGs_name(indx)   = temp_4d_name_WetDepo(i)
        var_User4d_XYZGs_unit(indx)   = temp_4d_unit_WetDepo(i)
        var_User4d_XYZGs_lname(indx)  = temp_4d_lname_WetDepo(i)
        var_User4d_XYZGs_MissVal(indx)= temp_4d_MissVal_WetDepo(i)
        var_User4d_XYZGs_FillVal(indx)= temp_4d_FillVal_WetDepo(i)
        if(i.eq.1) &
          var_User4d_XYZGs(1:nxmax,1:nymax,1:nzmax,1:n_gs_max,indx) = &
            scav_coeff_Liq_3d(1:nxmax,1:nymax,1:nzmax,1:n_gs_max)
      enddo

      end subroutine Prep_output_WetDepo

!******************************************************************************

      subroutine Deallocate_WetDepo_global

      implicit none

      deallocate(MetCloudBotHeight)
      deallocate(MetCloudTopHeight)

      if(USE_WASHOUT_LIQ.or.USE_RAINOUT_LIQ)then
        deallocate(scav_coeff_Liq_3d)
      endif
      if(USE_WASHOUT_ICE.or.USE_RAINOUT_ICE)then
        deallocate(scav_coeff_Ice_3d)
      endif

      if(USE_WASHOUT_LIQ)then  ! Below-cloud rain
        deallocate(Deposit_Liq_Washout)
      endif
      if(USE_WASHOUT_ICE)then  ! Below-cloud snow
        deallocate(Deposit_Ice_Washout)
      endif
      if(USE_RAINOUT_LIQ)then  ! In-cloud drop nucleation
        deallocate(Deposit_Liq_Rainout)
      endif
      if(USE_RAINOUT_ICE)then  ! In-cloud ice nucleation
        deallocate(Deposit_Ice_Rainout)
      endif

      ! Precipitation rates
      if(USE_3D_PRECIP)then
        deallocate(precipitation_rate_3d)
      endif
      deallocate(precipitation_rate_2d)

      end subroutine Deallocate_WetDepo_global

!******************************************************************************

      subroutine Deallocate_WetDepo_Met

      implicit none

      !deallocate(CatPrecip_meso_last_step)
      !deallocate(CatPrecip_meso_next_step)

      deallocate(CatagorPrecip_met)

      ! lower cloud position
      deallocate(MetCloudTop_MetP)
      deallocate(MetCloudTop_meso_last_step)
      deallocate(MetCloudTop_meso_next_step)
      deallocate(MetCloudBot_MetP)
      deallocate(MetCloudBot_meso_last_step)
      deallocate(MetCloudBot_meso_next_step)

      ! Scavenging coefficients
      deallocate(scav_coeff_MetP)
      deallocate(scav_coeff_Liq_meso_last_step)
      deallocate(scav_coeff_Liq_meso_next_step)


      ! Precipitation rates (required)
      deallocate(prate_Wat_MetP_sp)
!      deallocate(prate_Ice_MetP_sp)

      ! Precipitation rates (only needed if prate is an exported variable)
!      if(USE_ADDITIONAL_VARS)then
!        if(USE_3D_PRECIP.and.log_3d2_Precip_rate)then
!          deallocate(prate_col_metH_sp)
!          deallocate(prate_metH_sp)
!          deallocate(prate_meso_last_step)
!          deallocate(prate_meso_next_step)
!        endif
!        if(log_2d2_Precip_rate)then
!          deallocate(prate_2d_MetP_sp)
!          deallocate(prate_2d_meso_last_step)
!          deallocate(prate_2d_meso_next_step)
!        endif
!      endif

!      if(MR_iwindformat.eq.4.or.MR_iwindformat.eq.24)then
!        deallocate(QL_MetP_sp)
!        deallocate(QI_wind_sp)
!      endif

      end subroutine Deallocate_WetDepo_Met

!******************************************************************************

!******************************************************************************

      subroutine Set_WetDepo_Meso(Load_MesoSteps,Interval_Frac)

      use mesh,          only : &
         nsmax

      use solution,      only : &
         IsAloft

      use MetReader,     only : &
         MR_dum3d_compH,MR_dum3d_MetP, MR_dum2d_Met,MR_iMetStep_Now,MR_dum2d_comp,&
           MR_Read_2d_Met_Variable,&
           MR_Regrid_Met2d_to_Comp2d,&
           MR_Read_2d_Met_Variable,&
           MR_Regrid_MetP_to_CompGrid

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=ip),intent(in) :: Interval_Frac

      !logical,save :: first_time = .true.
      integer :: ivar
      integer :: l

      !write(global_info,*)"Inside Set_WetDepo_Meso"
      if(Load_MesoSteps)then

        if(CloudLoc.eq.3)then
          ivar = 20 !  pressure at lower cloud base
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
          Pres_lcb_MetP_sp = MR_dum2d_Met(:,:)
          call MR_Regrid_Met2d_to_Comp2d
            MetCloudBot_meso_last_step = MR_dum2d_comp(:,:)
          ivar = 21 !  pressure at lower cloud top
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
          Pres_lct_MetP_sp = MR_dum2d_Met(:,:)
          call MR_Regrid_Met2d_to_Comp2d
            MetCloudTop_meso_last_step = MR_dum2d_comp(:,:)
          call Set_Cloud_Level(MR_iMetStep_Now)

          ivar = 20 !  pressure at lower cloud base
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now+1)
          Pres_lcb_MetP_sp = MR_dum2d_Met(:,:)
          call MR_Regrid_Met2d_to_Comp2d
            MetCloudBot_meso_next_step = MR_dum2d_comp(:,:)  
          ivar = 21 !  pressure at lower cloud top
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now+1)
          Pres_lct_MetP_sp = MR_dum2d_Met(:,:)
          call MR_Regrid_Met2d_to_Comp2d
            MetCloudTop_meso_next_step = MR_dum2d_comp(:,:)
          call Set_Cloud_Level(MR_iMetStep_Now+1)
        else
          if(CloudLoc.eq.0)then
            MetCloudBot_MetP = 100.0_ip
            MetCloudTop_MetP = 101.0_ip

            MetCloudBot_meso_last_step = 100.0
            MetCloudBot_meso_next_step = 100.0
            MetCloudTop_meso_last_step = 101.0
            MetCloudTop_meso_next_step = 101.0

          else
            call Set_Cloud_Level(MR_iMetStep_Now)
            MR_dum2d_Met(:,:) = real(MetCloudBot_MetP(:,:),kind=sp)
            call MR_Regrid_Met2d_to_Comp2d
              MetCloudBot_meso_last_step = MR_dum2d_comp(:,:)
            MR_dum2d_Met(:,:) = real(MetCloudTop_MetP(:,:),kind=sp)
            call MR_Regrid_Met2d_to_Comp2d
              MetCloudBot_meso_last_step = MR_dum2d_comp(:,:)

            call Set_Cloud_Level(MR_iMetStep_Now+1)
            MR_dum2d_Met(:,:) = real(MetCloudBot_MetP(:,:),kind=sp)
            call MR_Regrid_Met2d_to_Comp2d
              MetCloudBot_meso_next_step = MR_dum2d_comp(:,:)
            MR_dum2d_Met(:,:) = real(MetCloudTop_MetP(:,:),kind=sp)
            call MR_Regrid_Met2d_to_Comp2d
              MetCloudTop_meso_next_step = MR_dum2d_comp(:,:)
          endif
        endif

        ! To get the scavenging coefficient, first read prate on MetP then
        ! convert to scav_coeff
        ivar = 44 !  Precipitation rate large-scale (liquid)
          ! Read data at MR_iMetStep_Now (expects 2d, but reads 3d for MERRA)
        call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
        if(USE_3D_PRECIP)then
          prate_Wat_MetP_sp(:,:,:) = MR_dum3d_MetP(:,:,:)
            ! This part is only needed for writing p0 to outfile
          MR_dum3d_MetP = prate_Wat_MetP_sp
          call MR_Regrid_MetP_to_CompGrid(MR_iMetStep_Now)
          precipitation_rate_3d(:,:,:) = MR_dum3d_compH
          precipitation_rate_2d(:,:)   = precipitation_rate_3d(:,:,5)
        else
          prate_Wat_MetP_sp(:,:,1) = MR_dum2d_Met(:,:)
           ! remap to comp grid if needed for output
           !  Note: the correct variable is already in MR_dum2d_met
          call MR_Regrid_Met2d_to_Comp2d
          precipitation_rate_2d(:,:)   = MR_dum2d_comp(:,:)
        endif
        call Set_Scav_Coeff(MR_iMetStep_Now,1)
        !write(global_info,*)MR_iMetStep_Now,prate_Wat_MetP_sp(5,5,5),scav_coeff_MetP(5,5,5,:)
        do l=1,nsmax
          if(.not.IsAloft(l)) cycle
          MR_dum3d_MetP(:,:,:) = real(scav_coeff_MetP(:,:,:,l),kind=sp)
          call MR_Regrid_MetP_to_CompGrid(MR_iMetStep_Now)
          scav_coeff_Liq_meso_last_step(:,:,:,l) = MR_dum3d_compH(:,:,:)
        enddo

          ! Read data at MR_iMetStep_Now + 1
        call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now+1)
        if(USE_3D_PRECIP)then
          prate_Wat_MetP_sp = MR_dum3d_MetP
            ! This part is only needed for writing p0 to outfile
          !MR_dum3d_MetP = prate_Wat_MetP_sp
          !call MR_Regrid_MetP_to_CompGrid(MR_iMetStep_Now+1)
          !precipitation_rate_3d(:,:,:) = MR_dum3d_compH
          !precipitation_rate_2d(:,:)   = precipitation_rate_3d(:,:,5)
        else
          prate_Wat_MetP_sp(:,:,1) = MR_dum2d_Met(:,:)
           ! remap to comp grid if needed for output
           !  Note: the correct variable is already in MR_dum2d_met
          !call MR_Regrid_Met2d_to_Comp2d
          !precipitation_rate_2d(:,:)   = MR_dum2d_comp(:,:)
        endif
        call Set_Scav_Coeff(MR_iMetStep_Now+1,2)
        do l=1,nsmax
          MR_dum3d_MetP(:,:,:) = real(scav_coeff_MetP(:,:,:,l),kind=sp)
          call MR_Regrid_MetP_to_CompGrid(MR_iMetStep_Now+1)
          scav_coeff_Liq_meso_next_step(:,:,:,l) = MR_dum3d_compH(:,:,:)
        enddo

      endif

      call Interpolate_WetDepo(Interval_Frac)

      end subroutine Set_WetDepo_Meso

!******************************************************************************

      subroutine Set_Cloud_Level(istep)

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,np_fullmet_RH,MR_iMetStep_Now,&
         MR_geoH_MetP_last,MR_geoH_MetP_next,&
         p_fullmet_RH_sp,MR_dum3d_MetP,&
           MR_Read_3d_MetP_Variable

      implicit none

      integer, intent(in) :: istep

      integer       :: i,j
      !integer       :: nPrecip
      real(kind=sp) :: zH(np_fullmet_RH)
      !real(kind=sp) :: inthresh1,inthresh2
      !real(kind=sp) :: invar(npRH)

      integer :: k
      real(kind=ip) :: var_to_check(np_fullmet_RH)
      real(kind=ip) :: thresh1,thresh2
      real(kind=ip) :: cloud_thick_default
      !real(kind=sp) :: prec_rate_liq_sp
      !real(kind=sp) :: prec_rate_ice_sp
      integer :: ivar

        !  1 = via Specific Humidity
        !  2 = via Relative Humidity
        !  3 = using cloud bot/top variables from met file

      If(CloudLoc.eq.1)then
        !  1 = via Specific Humidity
        allocate(QL_MetP_sp(nx_submet,ny_submet,np_fullmet))
        ivar = 32
        call MR_Read_3d_MetP_Variable(ivar,istep)
        QL_MetP_sp  = MR_dum3d_MetP
      endif

      do i=1,nx_submet
        do j=1,ny_submet

          if(istep.eq.MR_iMetStep_Now)then
            zH(1:np_fullmet_RH) = MR_geoH_MetP_last(i,j,1:np_fullmet_RH)
          elseif(istep.eq.MR_iMetStep_Now+1)then
            zH(1:np_fullmet_RH) = MR_geoH_MetP_next(i,j,1:np_fullmet_RH)
          endif

      MetCloudBot_MetP(i,j) = -1.0_ip
      MetCloudTop_MetP(i,j) = -1.0_ip

      cloud_thick_default = 0.0_ip

      if(CloudLoc.eq.1)then
        ! Here's if we are calculating bases on liquid water mixing ratios
        ! (essentially specific humidity)
        var_to_check = real(QL_MetP_sp(i,j,1:np_fullmet_RH),kind=ip)
        thresh1      = 1.0e-5_ip
        thresh2      = 1.0e-5_ip
      elseIf(CloudLoc.eq.2)then
        ! Here's if we are calculating based on relative humidity
        !var_to_check = real(AirRelH_meso_next_step_MetP_sp(i,j,1:np_fullmet_RH) &
        !                     ,kind=ip)
!        var_to_check = real(invar,kind=ip)

        thresh1      = 8.0e-1_ip
        thresh2      = 6.0e-1_ip
      elseif(CloudLoc.eq.3)then
        ! Use cloud top and bottom pressures from wind file
        var_to_check = real(p_fullmet_RH_sp,kind=ip)
        thresh1      = real(Pres_lcb_MetP_sp(i,j),kind=ip)
        thresh2      = real(Pres_lct_MetP_sp(i,j),kind=ip)
      else
        ! Cloud identification not specified
        ! Should not be here
        stop 1
      endif

      ! Get Cloud bottom
      do k=1,np_fullmet_RH-1
        If(CloudLoc.eq.3)then
          ! using pressure as var_to_check
          if(k.eq.1.and.var_to_check(k).le.thresh1)then
              ! Cloud layer goes all the way to the ground
            MetCloudBot_MetP(i,j) = 0.0_ip
            exit
          elseif(var_to_check(k  ).gt.thresh1.and.&
                 var_to_check(k+1).le.thresh1)then
            MetCloudBot_MetP(i,j) = zH(k) + &
               (zH(k+1)-zH(k)) * &
               (thresh1-var_to_check(k))/(var_to_check(k+1)-var_to_check(k))
            exit
          endif
        else
          if(k.eq.1.and.var_to_check(k).ge.thresh1)then
              ! Cloud layer goes all the way to the ground
            MetCloudBot_MetP(i,j) = 0.0_ip
            exit
          elseif(var_to_check(k  ).lt.thresh1.and.&
                 var_to_check(k+1).ge.thresh1)then
            MetCloudBot_MetP(i,j) = zH(k) + &
               (zH(k+1)-zH(k)) * &
               (thresh1-var_to_check(k))/(var_to_check(k+1)-var_to_check(k))
            exit
          endif
        endif
      enddo
      if(MetCloudBot_MetP(i,j).lt.0.0_ip) MetCloudBot_MetP(i,j) = 0.0_ip

      ! Get Cloud top
      do k=1,np_fullmet_RH-1
          ! cycle for values below cloud bottom
        if(zH(k).le.MetCloudBot_MetP(i,j)) cycle

        if(var_to_check(k  ).ge.thresh2.and.&
           var_to_check(k+1).lt.thresh2)then
          MetCloudTop_MetP(i,j) = zH(k) + &
             (zH(k+1)-zH(k)) * &
             (thresh2-var_to_check(k))/(var_to_check(k+1)-var_to_check(k))
          exit
        endif
      enddo

          If(MetCloudTop_MetP(i,j).lt.0.0_ip)then
            MetCloudTop_MetP(i,j) = MetCloudBot_MetP(i,j) + &
                                            cloud_thick_default
          endif

        enddo
      enddo

      If(CloudLoc.eq.1)then
        deallocate(QL_MetP_sp)
      endif

      end subroutine Set_Cloud_Level

!******************************************************************************

!      subroutine RH_to_SH()
!
!        !http://www.cactus2000.de/js/calchum.pdf
!
!        a0 = 6.107799961
!        a1 = 4.436518521e-1
!        a2 = 1.428945805e-2
!        a3 = 2.650648471e-4
!        a4 = 3.031240396e-6
!        a5 = 2.034080948e-8
!        a6 = 6.136820929e-11
!
!        T = temp in C
!        P = pressure in hPa
!
!          !Lowe, P.R. and J.M. Ficke, 1974: The computation of saturation vapor pressure.
!          !Tech. Paper No. 4-74, Environmental Prediction Research Facility, Naval Postgraduate School,
!          !Monterey, CA, 27 pp
!        vappres = a0+T*(a1+T*(a2+T*(a3+T*(a4+T*(a5+T*a6)))));
!        PH2O = RH*vappres
!        VolMixRat = PH2O/P
!        SH = (VolMixRat*MH2O)/(VolMixRat*MH2O +(1.0-VolMixRat)*MH2O)
!
!      end subroutine RH_to_SH

!******************************************************************************

      subroutine Set_Scav_Coeff(istep,last_flag)

      use global_param,  only : &
         MPS_2_KMPHR,KM_2_M

      use Atmosphere,    only : &
         AirTemp_meso_last_step_MetP_sp,AirTemp_meso_next_step_MetP_sp,&
         AirDens_meso_last_step_MetP_sp,AirDens_meso_next_step_MetP_sp,&
         AirVisc_meso_last_step_MetP_sp,AirVisc_meso_next_step_MetP_sp

      use Tephra,        only : &
         n_gs_max,Tephra_rho_m,Tephra_gsdiam,vf_meso_last_step_MetP_sp,&
         vf_meso_next_step_MetP_sp

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,MR_geoH_MetP_last,MR_geoH_MetP_next

      implicit none

      integer, intent(in) :: istep
      integer, intent(in) :: last_flag

!      real(kind=ip) :: z_col(nkmax)

      integer :: i,j,k,l
      real(kind=ip) :: Coll_Effic
      real(kind=ip) :: rain_diam_opt,rain_vel_opt     ! Optimal values returned
                                                      !  from functions
      real(kind=ip) :: prate_liq
      !real(kind=ip) :: prate_ice
      real(kind=ip) :: prate_surface

      real(kind=ip) :: temp_air,dens_air,eta_air
      real(kind=ip),dimension(n_gs_max) :: fall_vel
      real(kind=ip),dimension(np_fullmet) :: prate_Wat_colP
      real(kind=ip),dimension(np_fullmet) :: geoH_colP
      real(kind=ip) :: Cloud_Top,Cloud_Bot,frac

      !  The scavenging coefficient for liquid is based entirely on
      !  the precipitation rate
      !  Loop over Met grid and calculate scavCo at each pressure value

      do i=1,nx_submet
        do j=1,ny_submet
          !  For this i,j, set up the column of prate values, either from 3d
          !  data (if using NASA/MERRA) or from surface values and cloud
          !  height
          if(USE_3D_PRECIP)then
            prate_Wat_colP(1:np_fullmet) = prate_Wat_MetP_sp(i,j,1:np_fullmet)
          else
            if(last_flag.eq.1)then
              geoH_colP(1:np_fullmet) = MR_geoH_MetP_last(i,j,1:np_fullmet)
              Cloud_Top    = MetCloudTop_MetP(i,j)
              Cloud_Bot    = MetCloudBot_MetP(i,j)
            else
              geoH_colP(1:np_fullmet) = MR_geoH_MetP_next(i,j,1:np_fullmet)
              Cloud_Top    = MetCloudTop_MetP(i,j)
              Cloud_Bot    = MetCloudBot_MetP(i,j)
            endif
            prate_surface = real(prate_Wat_MetP_sp(i,j,1),kind=ip)
            do k=1,np_fullmet  ! loop over pressure values
              if(geoH_colP(k).lt.Cloud_Bot)then
                prate_Wat_colP(k) = prate_surface
              elseif(geoH_colP(k).lt.Cloud_Top)then
                frac = (Cloud_Top-geoH_colP(k))/(Cloud_Top-Cloud_Bot)
                prate_Wat_colP(k) = prate_surface*frac
              else
                prate_Wat_colP(k) = 0.0_ip
              endif
            enddo
          endif ! finished column of p0 values
          do k=1,np_fullmet  ! loop over pressure values
            if(last_flag.eq.1)then
              temp_air = real(AirTemp_meso_last_step_MetP_sp(i,j,k),kind=ip)
              dens_air = real(AirDens_meso_last_step_MetP_sp(i,j,k),kind=ip)
              eta_air  = real(AirVisc_meso_last_step_MetP_sp(i,j,k),kind=ip)
              fall_vel(1:n_gs_max) = &
                         real(     vf_meso_last_step_MetP_sp(i,j,k,1:n_gs_max),kind=ip)
            else
              temp_air = real(AirTemp_meso_next_step_MetP_sp(i,j,k),kind=ip)
              dens_air = real(AirDens_meso_next_step_MetP_sp(i,j,k),kind=ip)
              eta_air  = real(AirVisc_meso_next_step_MetP_sp(i,j,k),kind=ip)
              fall_vel(1:n_gs_max) = &
                         real(     vf_meso_next_step_MetP_sp(i,j,k,1:n_gs_max),kind=ip)
            endif
            do l=1,n_gs_max
        ! Calculate scavenging coefficients for particle
        ! Don't bother with expensive collision efficiency
        ! calculations unless it is actually raining
        !  Recall that prate is in m/s
        !fall_vel = real(vf_meso_last_step_MetP_sp(i,j,k,l),kind=ip)
        if(USE_3D_PRECIP)then
          prate_liq = prate_Wat_colP(k)
          if(prate_liq.gt.1.0e-9_ip)then
            rain_diam_opt = Get_Rain_diam_p0(prate_liq)
            rain_vel_opt  = Get_Rain_vel(rain_diam_opt)
            if(Tephra_gsdiam(l).lt.5.0e-5_ip)then
                ! Particles smaller than 50um have a decreased E
                ! (Greenfield gap)
              Coll_Effic = CEffic(rain_diam_opt,rain_vel_opt, &
                                  dens_air,Tephra_rho_m(l),eta_air,Tephra_gsdiam(l), &
                                  temp_air,&
                                  fall_vel(l))
            else
              Coll_Effic = 1.0_ip
            endif
            ! Scavenging coefficient in 1/hr
            scav_coeff_MetP(i,j,k,l) = Coll_Effic * &
                                  1.5_ip*prate_liq*MPS_2_KMPHR &
                                   /(rain_diam_opt/KM_2_M)
          else
            scav_coeff_MetP(i,j,k,l) = 0.0_ip
          endif
        else
!
!          prate_liq = real(prate_Wat_MetP_sp(i,j,1),kind=ip)
!          if(prate_liq.gt.1.0e-7_ip)then
!            if(z_col(k).le.MetCloudBot_MetP(i,j))then
!                ! Below cloud washout
!              rain_diam_opt = Get_Rain_diam_p0(prate_liq)
!              rain_vel_opt  = Get_Rain_vel(rain_diam_opt)
!              if(Tephra_gsdiam(l).lt.5.0e-5_ip)then
!                ! Particles smaller than 50um have a
!                ! decreased E (Greenfield gap)
!                Coll_Effic = CEffic(rain_diam_opt,rain_vel_opt,&
!                                    dens_col_MetP(k),Tephra_rho_m(l),eta_col_MetP(k),Tephra_gsdiam(l), &
!                                  real(T_col_MetP_sp(k),kind=ip),&
!                                  real(vs_col_MetP_sp(k,l),kind=ip))
!
!              else
!                Coll_Effic = 1.0_ip
!              endif
!                ! Scavenging coefficient in 1/hr
!              scav_coeff_col_MetP(k,l) = Coll_Effic * &
!                                    1.5_ip*prate_liq*MPS_2_KMPHR &
!                                     /(rain_diam_opt/KM_2_M)
!            elseif(z_col(k).le.MetCloudTop_MetP(i,j))then
!                ! in-cloud rainout
!              rain_diam_opt = Get_Rain_diam_p0(prate_liq)
!              rain_vel_opt  = Get_Rain_vel(rain_diam_opt)
!              if(Tephra_gsdiam(l).lt.5.0e-5_ip)then
!                  ! Particles smaller than 50um have a
!                  ! decreased E (Greenfield gap)
!                Coll_Effic = CEffic(rain_diam_opt,rain_vel_opt,&
!                                    dens_col_MetP(k),Tephra_rho_m(l),eta_col_MetP(k),Tephra_gsdiam(l),&
!                                  real(T_col_MetP_sp(k),kind=ip),&
!                                  real(vs_col_MetP_sp(k,l),kind=ip))
!              else
!                Coll_Effic = 1.0_ip
!              endif
!              ! Scavenging coefficient in 1/hr
!              scav_coeff_col_MetP(k,l) = Coll_Effic * &
!                                    1.5_ip*prate_liq*MPS_2_KMPHR &
!                                     /(rain_diam_opt/KM_2_M)
!            else
!              scav_coeff_col_MetP(k,l) = 0.0_ip
!           endif
!          else
!            scav_coeff_col_MetP(k,l) = 0.0_ip
!          endif
        endif
            enddo ! n_gs_max
          enddo ! k
        enddo ! y
      enddo ! x

      end subroutine Set_Scav_Coeff

!******************************************************************************

      subroutine Interpolate_WetDepo(Interval_Frac)

      use mesh,          only : &
         nxmax,nymax,nzmax

      use Tephra,        only : &
         n_gs_max

      implicit none

      real(kind=ip), intent(in) :: Interval_Frac

        !scav_coeff_Liq_3d(1:nxmax,1:nymax,1:nzmax,:) = scav_coeff_meso_last_step + &
        !                                               (scav_coeff_meso_next_step - &
        !                                                scav_coeff_meso_last_step) * &
        !                                               Interval_Frac
      scav_coeff_Liq_3d(1:nxmax,1:nymax,1:nzmax,1:n_gs_max) = &
         scav_coeff_Liq_meso_last_step(1:nxmax,1:nymax,1:nzmax,1:n_gs_max) + &
                 (scav_coeff_Liq_meso_next_step(1:nxmax,1:nymax,1:nzmax,1:n_gs_max) - &
                  scav_coeff_Liq_meso_last_step(1:nxmax,1:nymax,1:nzmax,1:n_gs_max)) * Interval_Frac

      ! This needs to be a function of the ice mass flux
      ! KLUDGE
      if(USE_WASHOUT_ICE)then
        scav_coeff_Ice_3d = scav_coeff_Liq_3d
      endif

      !max_scav = maxval(scav_coeff_Liq_3d)
!      if(USE_3D_PRECIP.and.&
!       USE_ADDITIONAL_VARS.and.log_3d2_Precip_rate)then
!        precipitation_rate_3d(1:nxmax,1:nymax,1:nzmax) = &
!                      real(prate_meso_last_step_sp(:,:,:),kind=ip) + &
!                      real((prate_meso_next_step_sp(:,:,:) - &
!                            prate_meso_last_step_sp(:,:,:)),kind=ip) * &
!                      Interval_Frac
!
!      endif
!
!!      if(USE_ADDITIONAL_VARS.and.log_2d2_Precip_rate)then
!        precipitation_rate_2d(1:nxmax,1:nymax) = &
!                      real(prate_2d_meso_last_step_sp(:,:),kind=ip) + &
!                      real((prate_2d_meso_next_step_sp(:,:) - &
!                            prate_2d_meso_last_step_sp(:,:)),kind=ip) * &
!                      Interval_Frac

        MetCloudBotHeight(1:nxmax,1:nymax) = MetCloudBot_meso_last_step + &
                    (MetCloudBot_meso_next_step - &
                     MetCloudBot_meso_last_step)* &
                    Interval_Frac
        MetCloudTopHeight(1:nxmax,1:nymax) = MetCloudTop_meso_last_step + &
                    (MetCloudTop_meso_last_step - &
                     MetCloudTop_meso_last_step)* &
                    Interval_Frac

!      endif

      end subroutine Interpolate_WetDepo

!******************************************************************************

      subroutine Wet_Depo_Rainout

      use global_param,  only : &
         EPS_TINY,CFL

      use mesh,          only : &
        nxmax,nymax,nzmax,ts0,ts1,z_cc_pd,dz_vec_pd

      use time_data,     only : &
         dt

      use solution,      only : &
         concen_pd

      use Tephra,        only : &
         n_gs_max

      implicit none

      integer :: i,j,k,n
      real(kind=ip) :: scrub
      real(kind=ip) :: max_scav,scav_time,dt_sub
      integer :: n_substeps,it

      ! Note: since this subroutine uses sub-stepping, concentration values are
      ! decayed in place rather than calculating concen(:,:,:,:,ts1) and copying
      ! to concen(:,:,:,:,ts0)

      ! Decay (exponential) concentrations due to wet scavenging
      if(USE_WASHOUT_LIQ.or.USE_RAINOUT_LIQ)then
        do i=1,nxmax
          do j=1,nymax
              ! Get the maximum scavenging coeffient over all k and l for this
              ! column at i,j and calculate the required dt
            max_scav = maxval(scav_coeff_Liq_3d(i,j,:,:))
            if(max_scav.lt.EPS_TINY) cycle
            scav_time = 2.0_ip*max_scav
            dt_sub = CFL/scav_time
            If(dt_sub.ge.dt)then
              ! time step is already more restrictive that scavenging decay.
              ! Use the global time step
              dt_sub = dt
              n_substeps = 1
            else
              ! Scavenging requires a smaller time step that dt.
              ! Get the number of sub-steps to take
              n_substeps = ceiling(dt/dt_sub)
              dt_sub = dt/real(n_substeps,kind=ip)
            endif

            !write(global_info,*)MetCloudBotHeight(i,j),MetCloudTopHeight(i,j)

            do k=1,nzmax
              ! First scrub all the BELOW-CLOUD points
              if(USE_WASHOUT_LIQ.and.z_cc_pd(k).le.MetCloudBotHeight(i,j))then
                do n=1,n_gs_max
                  If(scav_coeff_Liq_3d(i,j,k,n).gt.EPS_TINY.and.&
                     concen_pd(i,j,k,n,ts0).gt.EPS_TINY)then
                    do it=1,n_substeps
                        ! Calculate amount transferred to rain phase
                      scrub =  dt_sub*scav_coeff_Liq_3d(i,j,k,n)*concen_pd(i,j,k,n,ts0)
                        ! decrement aloft concentration by scrubbed amount
                      concen_pd(i,j,k,n,ts0) = concen_pd(i,j,k,n,ts0) - scrub
                        ! increment wet deposit by scrubbed amount
                        !  and convert from kg/km3 to kg/m2
                      Deposit_Liq_Washout(i,j,n) = &
                        Deposit_Liq_Washout(i,j,n) + scrub*dz_vec_pd(k)*1.0e-6_ip
                      !write(global_info,*)"Deposit_Liq_Washout ",Deposit_Liq_Washout(i,j,n)
                    enddo
                  endif
                enddo
              ! Next scrub all the IN-CLOUD points
              elseif(USE_RAINOUT_LIQ.and.z_cc_pd(k).gt.MetCloudBotHeight(i,j).and.&
                                         z_cc_pd(k).le.MetCloudTopHeight(i,j))then
                do n=1,n_gs_max
                  If(scav_coeff_Liq_3d(i,j,k,n).gt.EPS_TINY.and.&
                     concen_pd(i,j,k,n,ts0).gt.EPS_TINY)then
                    do it=1,n_substeps
                        ! Calculate amount transferred to rain phase
                      scrub = dt_sub*scav_coeff_Liq_3d(i,j,k,n)*concen_pd(i,j,k,n,ts0)
                        ! decrement aloft concentration by scrubbed amount
                      concen_pd(i,j,k,n,ts0) = concen_pd(i,j,k,n,ts0) - scrub
                        ! increment wet deposit by scrubbed amount
                        !  and convert from kg/km3 to kg/m2
                      Deposit_Liq_Rainout(i,j,n) = &
                        Deposit_Liq_Rainout(i,j,n) + scrub*dz_vec_pd(k)*1.0e-6_ip
                      write(global_info,*)"Deposit_Liq_Rainout ",Deposit_Liq_Rainout(i,j,n)

                    enddo
                  endif
                enddo
              else
                ! This is a place holder for ABOVE-CLOUD processes

              endif
            enddo
          enddo
        enddo
      endif ! USE_WASHOUT_LIQ.or.USE_RAINOUT_LIQ

      ! Decay (exponential) concentrations due to snow scavenging
      !  This is really just a place-holder using the same scheme as wet
      !  scavenging
      !if(USE_WASHOUT_ICE.or.USE_RAINOUT_ICE)then
      !  do i=1,nxmax
      !    do j=1,nymax
      !        ! Get the maximum scavenging coeffient over all k and l for this
      !        ! column at i,j and calculate the required dt
      !      max_scav = maxval(scav_coeff_Ice_3d(i,j,:,:))
      !      if(max_scav.lt.EPS_TINY) cycle
      !      scav_time = 2.0_ip*max_scav
      !      dt_sub = CFL/scav_time
      !      If(dt_sub.ge.dt)then
      !        ! time step is already more restrictive that scavenging decay.
      !        ! Use the global time step
      !        dt_sub = dt
      !      else
      !        ! Scavenging requires a smaller time step that dt.
      !        ! Get the number of sub-steps to take
      !        n_substeps = ceiling(dt/dt_sub)
      !        dt_sub = dt/real(n_substeps,kind=ip)
      !      endif
      !      do k=1,nzmax
      !        ! First scrub all the BELOW-CLOUD points
      !        if(USE_WASHOUT_ICE.and.z_cc(k).le.MetCloudBotHeight(i,j))then
      !          do n=1,n_gs_max
      !            If(scav_coeff_Ice_3d(i,j,k,n).gt.EPS_TINY)then
      !              do it=1,n_substeps
      !                  ! Calculate amount transferred to rain phase
      !                scrub = dt_sub*scav_coeff_Ice_3d(i,j,k,n)*concen(i,j,k,n,ts0)
      !                  ! decrement aloft concentration by scrubbed amount
      !                concen(i,j,k,n,ts0) = concen(i,j,k,n,ts0) - scrub
      !                  ! increment wet deposit by scrubbed amount
      !                  !  and convert from kg/km3 to kg/m2
      !                Deposit_Ice_Washout(i,j,n) = &
      !                  Deposit_Ice_Washout(i,j,n) + scrub*dz*1.0e-6_ip
      !              enddo
      !            endif
      !          enddo
      !        ! Next scrub all the IN-CLOUD points
      !        elseif(USE_RAINOUT_ICE.and.z_cc(k).gt.MetCloudBotHeight(i,j).and.&
      !                                   z_cc(k).le.MetCloudTopHeight(i,j))then
      !          do n=1,n_gs_max
      !            If(scav_coeff_Ice_3d(i,j,k,n).gt.EPS_TINY)then
      !              do it=1,n_substeps
      !                  ! Calculate amount transferred to rain phase
      !                scrub = dt_sub*scav_coeff_Ice_3d(i,j,k,n)*concen(i,j,k,n,ts0)
      !                  ! decrement aloft concentration by scrubbed amount
      !                concen(i,j,k,n,ts0) = concen(i,j,k,n,ts0) - scrub
      !                  ! increment wet deposit by scrubbed amount
      !                  !  and convert from kg/km3 to kg/m2
      !                Deposit_Ice_Rainout(i,j,n) = &
      !                  Deposit_Ice_Rainout(i,j,n) + scrub*dz*1.0e-6_ip
      !              enddo
      !            endif
      !          enddo
      !        else
      !          ! This is a place holder for ABOVE-CLOUD processes
      !          !   e.g.  ice nucleation
      !        endif
      !      enddo
      !    enddo
      !  enddo
      !endif ! USE_WASHOUT_LIQ.or.USE_RAINOUT_LIQ

      return

      end subroutine Wet_Depo_Rainout

!******************************************************************************

      function CEffic(rain_diam,rain_vel,rho_air,rho_ash,eta,diam,Temp_air,vt_ash)

      use global_param,  only : &
         PI,GRAV

      use Atmosphere,    only : &
         BoltzK

      implicit none

      real(kind=ip) CEffic    ! Collision efficiency of rain scrubbing ash
      real(kind=ip) rain_diam
      real(kind=ip) rain_vel
      real(kind=ip) rho_air   ! density of air in kg/m3
      real(kind=ip) rho_ash   ! density of the particle in km/m3
      real(kind=ip) eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) diam      ! diameter of the particle in m
      real(kind=ip) Temp_air  ! temperature of air in K
      real(kind=ip) vt_ash    ! Fall velocity of ash particle

      real(kind=ip) Cc,D,tau,Re,Sc,St,SpecialK,omega,Sstar
      real(kind=ip) term1,term2,term3

      Cc = 1.0_ip+(2.0_ip*6.5e-8_ip/diam) * &
               (1.257_ip+0.4_ip*exp(-1.1_ip*diam/(2.0_ip*6.5e-8_ip)))  !%Eq 8.36 of Seinfeld

      D = BoltzK * Temp_air * Cc/(6.0_ip * PI * 0.5_ip*diam * eta) !Eq 15.29 of Jacobson

      tau = vt_ash/GRAV               ! time scale for terminal velocity

      !write(global_info,*)'Cc,D,tau',Cc,D,tau
      Re = Rain_diam*Rain_Vel*rho_air/(2.0_ip*eta)   ! Reynolds (Note, this is half the Re from vset)
      Sc = eta/rho_air/D          ! Schmidt
      St = 2.0_ip*tau*(Rain_Vel-vt_ash)/Rain_diam      ! Stokes number
      !write(global_info,*)'Re,Sc,St',Re,Sc,St 
      SpecialK = diam/Rain_diam           ! ratio of diameters
      omega = Rain_Mu/eta;          ! ratio of viscosities
      Sstar = (1.2_ip+(1.0_ip/12.0_ip)*log(1.0_ip+Re))/(1.0_ip+log(1.0_ip+Re))
      !write(global_info,*)'SpecialK,omega,Sstar',SpecialK,omega,Sstar
        ! Brownian diffusion
      term1 = (4.0_ip/(Re*Sc))*(1.0_ip+0.4_ip*sqrt(Re) * &
               (Sc)**(1.0_ip/3.0_ip) + 0.16_ip * sqrt(Re*Sc))
        ! interception
      term2 = 4.0_ip*SpecialK*(1.0_ip/omega+(1.0_ip+2.0_ip*sqrt(Re))*SpecialK)
        !  impaction
      term3 = ((St-Sstar)/(St-Sstar+2.0_ip/3.0_ip))
      if(term3.lt.0.0_ip)then
        term3 = 0.0_ip
      else
        !term3 = real(((St-Sstar)/(St-Sstar+2.0/3.0))**1.5 * &
        !        sqrt(Rain_dens/rho_air))
        term3 = (term3**1.5_ip) * sqrt(Rain_dens/rho_ash)
      endif
      !write(global_info,*)term1,term2,term3,St-Sstar
      CEffic = term1 + term2 + term3
      ! Make sure the efficiency doesn't exceed unity (it could with
      ! electrostatic attraction)
      if (CEffic.gt.1.0_ip) then
        CEffic = 1.0_ip
      endif

      return

      end function CEffic


!******************************************************************************

      function Get_Rain_diam_p0(prec_rate)

      implicit none

      real(kind=ip) Get_Rain_diam_p0    ! Rain diam as a function of precipitation rate (m)
      real(kind=ip) prec_rate           ! Precipitation rate (m/s)

      real(kind=ip),dimension(101) :: LogP0, Dpopt
      real(kind=ip) LogPrecRate_mm_hr,frac
      integer i,idx

        ! Get the precip rate in mm/hr and convert to Log10 for LUT
      LogPrecRate_mm_hr = log10(prec_rate*1000.0_ip*3600.0_ip)

      ! These are the optimal droplet sizes for each precip rate based
      ! on matching scavenging coeffiecients
      LogP0(  1)=-2.00_ip ; Dpopt(  1) =  0.000038383094742_ip 
      LogP0(  2)=-1.96_ip ; Dpopt(  2) =  0.000040872079048_ip 
      LogP0(  3)=-1.92_ip ; Dpopt(  3) =  0.000043510794569_ip 
      LogP0(  4)=-1.88_ip ; Dpopt(  4) =  0.000046324015394_ip 
      LogP0(  5)=-1.84_ip ; Dpopt(  5) =  0.000049305098791_ip 
      LogP0(  6)=-1.80_ip ; Dpopt(  6) =  0.000052472568307_ip 
      LogP0(  7)=-1.76_ip ; Dpopt(  7) =  0.000055837465157_ip 
      LogP0(  8)=-1.72_ip ; Dpopt(  8) =  0.000059415683886_ip 
      LogP0(  9)=-1.68_ip ; Dpopt(  9) =  0.000063211280415_ip 
      LogP0( 10)=-1.64_ip ; Dpopt( 10) =  0.000067241080289_ip 
      LogP0( 11)=-1.60_ip ; Dpopt( 11) =  0.000071518627337_ip 
      LogP0( 12)=-1.56_ip ; Dpopt( 12) =  0.000076062662030_ip 
      LogP0( 13)=-1.52_ip ; Dpopt( 13) =  0.000080881829842_ip 
      LogP0( 14)=-1.48_ip ; Dpopt( 14) =  0.000085991450935_ip 
      LogP0( 15)=-1.44_ip ; Dpopt( 15) =  0.000091410164687_ip 
      LogP0( 16)=-1.40_ip ; Dpopt( 16) =  0.000097155199245_ip 
      LogP0( 17)=-1.36_ip ; Dpopt( 17) =  0.000103244589385_ip 
      LogP0( 18)=-1.32_ip ; Dpopt( 18) =  0.000109702804931_ip 
      LogP0( 19)=-1.28_ip ; Dpopt( 19) =  0.000116544818996_ip 
      LogP0( 20)=-1.24_ip ; Dpopt( 20) =  0.000123788273728_ip 
      LogP0( 21)=-1.20_ip ; Dpopt( 21) =  0.000131448272258_ip 
      LogP0( 22)=-1.16_ip ; Dpopt( 22) =  0.000139556375547_ip 
      LogP0( 23)=-1.12_ip ; Dpopt( 23) =  0.000148147216660_ip 
      LogP0( 24)=-1.08_ip ; Dpopt( 24) =  0.000157225007378_ip 
      LogP0( 25)=-1.04_ip ; Dpopt( 25) =  0.000166824711604_ip 
      LogP0( 26)=-1.00_ip ; Dpopt( 26) =  0.000176970554135_ip 
      LogP0( 27)=-0.96_ip ; Dpopt( 27) =  0.000187689494460_ip 
      LogP0( 28)=-0.92_ip ; Dpopt( 28) =  0.000199009401900_ip 
      LogP0( 29)=-0.88_ip ; Dpopt( 29) =  0.000210959040567_ip 
      LogP0( 30)=-0.84_ip ; Dpopt( 30) =  0.000223593181556_ip 
      LogP0( 31)=-0.80_ip ; Dpopt( 31) =  0.000236896858682_ip 
      LogP0( 32)=-0.76_ip ; Dpopt( 32) =  0.000250922540294_ip 
      LogP0( 33)=-0.72_ip ; Dpopt( 33) =  0.000265702421393_ip 
      LogP0( 34)=-0.68_ip ; Dpopt( 34) =  0.000281269425562_ip 
      LogP0( 35)=-0.64_ip ; Dpopt( 35) =  0.000297657142181_ip 
      LogP0( 36)=-0.60_ip ; Dpopt( 36) =  0.000314884552399_ip 
      LogP0( 37)=-0.56_ip ; Dpopt( 37) =  0.000332999276104_ip 
      LogP0( 38)=-0.52_ip ; Dpopt( 38) =  0.000352053375882_ip 
      LogP0( 39)=-0.48_ip ; Dpopt( 39) =  0.000372061875992_ip 
      LogP0( 40)=-0.44_ip ; Dpopt( 40) =  0.000393075910424_ip 
      LogP0( 41)=-0.40_ip ; Dpopt( 41) =  0.000415068172339_ip 
      LogP0( 42)=-0.36_ip ; Dpopt( 42) =  0.000438180992259_ip 
      LogP0( 43)=-0.32_ip ; Dpopt( 43) =  0.000462326871422_ip 
      LogP0( 44)=-0.28_ip ; Dpopt( 44) =  0.000487521352989_ip 
      LogP0( 45)=-0.24_ip ; Dpopt( 45) =  0.000514065645031_ip 
      LogP0( 46)=-0.20_ip ; Dpopt( 46) =  0.000541578414674_ip 
      LogP0( 47)=-0.16_ip ; Dpopt( 47) =  0.000569820841990_ip 
      LogP0( 48)=-0.12_ip ; Dpopt( 48) =  0.000599820984497_ip 
      LogP0( 49)=-0.08_ip ; Dpopt( 49) =  0.000613880150034_ip 
      LogP0( 50)=-0.04_ip ; Dpopt( 50) =  0.000663957705722_ip 
      LogP0( 51)= 0.0_ip  ; Dpopt( 51) =  0.000699141773661_ip 
      LogP0( 52)= 0.04_ip ; Dpopt( 52) =  0.000734967948251_ip 
      LogP0( 53)= 0.08_ip ; Dpopt( 53) =  0.000772131752146_ip 
      LogP0( 54)= 0.12_ip ; Dpopt( 54) =  0.000810558377982_ip 
      LogP0( 55)= 0.16_ip ; Dpopt( 55) =  0.000851483130866_ip 
      LogP0( 56)= 0.20_ip ; Dpopt( 56) =  0.000892819636830_ip 
      LogP0( 57)= 0.24_ip ; Dpopt( 57) =  0.000935617130347_ip 
      LogP0( 58)= 0.28_ip ; Dpopt( 58) =  0.000979748023649_ip 
      LogP0( 59)= 0.32_ip ; Dpopt( 59) =  0.001025218109007_ip 
      LogP0( 60)= 0.36_ip ; Dpopt( 60) =  0.001072110770174_ip 
      LogP0( 61)= 0.40_ip ; Dpopt( 61) =  0.001120321252018_ip 
      LogP0( 62)= 0.44_ip ; Dpopt( 62) =  0.001169853123239_ip 
      LogP0( 63)= 0.48_ip ; Dpopt( 63) =  0.001220705080154_ip 
      LogP0( 64)= 0.52_ip ; Dpopt( 64) =  0.001272844120924_ip 
      LogP0( 65)= 0.56_ip ; Dpopt( 65) =  0.001326056077780_ip 
      LogP0( 66)= 0.60_ip ; Dpopt( 66) =  0.001380384619962_ip 
      LogP0( 67)= 0.64_ip ; Dpopt( 67) =  0.001435721132826_ip 
      LogP0( 68)= 0.68_ip ; Dpopt( 68) =  0.001491998823966_ip 
      LogP0( 69)= 0.72_ip ; Dpopt( 69) =  0.001549100144659_ip 
      LogP0( 70)= 0.76_ip ; Dpopt( 70) =  0.001606892384379_ip 
      LogP0( 71)= 0.80_ip ; Dpopt( 71) =  0.001665293472714_ip 
      LogP0( 72)= 0.84_ip ; Dpopt( 72) =  0.001724086816954_ip 
      LogP0( 73)= 0.88_ip ; Dpopt( 73) =  0.001783156549395_ip 
      LogP0( 74)= 0.92_ip ; Dpopt( 74) =  0.001842356859510_ip 
      LogP0( 75)= 0.96_ip ; Dpopt( 75) =  0.001901581731783_ip 
      LogP0( 76)= 1.00_ip ; Dpopt( 76) =  0.001960474444800_ip 
      LogP0( 77)= 1.04_ip ; Dpopt( 77) =  0.002018976484259_ip 
      LogP0( 78)= 1.08_ip ; Dpopt( 78) =  0.002076786797170_ip 
      LogP0( 79)= 1.12_ip ; Dpopt( 79) =  0.002133685135737_ip 
      LogP0( 80)= 1.16_ip ; Dpopt( 80) =  0.002189460307726_ip 
      LogP0( 81)= 1.20_ip ; Dpopt( 81) =  0.002243865918166_ip 
      LogP0( 82)= 1.24_ip ; Dpopt( 82) =  0.002296610956834_ip 
      LogP0( 83)= 1.28_ip ; Dpopt( 83) =  0.002347521304726_ip 
      LogP0( 84)= 1.32_ip ; Dpopt( 84) =  0.002396222946548_ip 
      LogP0( 85)= 1.36_ip ; Dpopt( 85) =  0.002442537824868_ip 
      LogP0( 86)= 1.40_ip ; Dpopt( 86) =  0.002486161528615_ip 
      LogP0( 87)= 1.44_ip ; Dpopt( 87) =  0.002526774033698_ip 
      LogP0( 88)= 1.48_ip ; Dpopt( 88) =  0.002564067202557_ip 
      LogP0( 89)= 1.52_ip ; Dpopt( 89) =  0.002596858136719_ip 
      LogP0( 90)= 1.56_ip ; Dpopt( 90) =  0.002622439011437_ip 
      LogP0( 91)= 1.60_ip ; Dpopt( 91) =  0.002648151645577_ip 
      LogP0( 92)= 1.64_ip ; Dpopt( 92) =  0.002674278712952_ip 
      LogP0( 93)= 1.68_ip ; Dpopt( 93) =  0.002692308250996_ip 
      LogP0( 94)= 1.72_ip ; Dpopt( 94) =  0.002704934522981_ip 
      LogP0( 95)= 1.76_ip ; Dpopt( 95) =  0.002712249642764_ip 
      LogP0( 96)= 1.80_ip ; Dpopt( 96) =  0.002714365792766_ip 
      LogP0( 97)= 1.84_ip ; Dpopt( 97) =  0.002711204767355_ip 
      LogP0( 98)= 1.88_ip ; Dpopt( 98) =  0.002702571111501_ip 
      LogP0( 99)= 1.92_ip ; Dpopt( 99) =  0.002688363346259_ip 
      LogP0(100)= 1.96_ip ; Dpopt(100) =  0.002668654428230_ip 
      LogP0(101)= 2.00_ip ; Dpopt(101) =  0.002643247260599_ip 

      if(LogPrecRate_mm_hr.lt.-2.0_ip)then
        Get_Rain_diam_p0 = Dpopt(1)
      elseif(LogPrecRate_mm_hr.gt.2.0_ip)then
        Get_Rain_diam_p0 = Dpopt(101)
      else
        do i=1,100
          if(LogPrecRate_mm_hr.ge.LogP0(i).and.LogPrecRate_mm_hr.lt.LogP0(i+1)) &
            idx = i
        enddo
        frac = (LogPrecRate_mm_hr-LogP0(idx))/(LogP0(idx+1)-LogP0(idx))
        Get_Rain_diam_p0 = Dpopt(idx) + frac*(Dpopt(idx+1)-Dpopt(idx))
      endif

      end function

!******************************************************************************

      function Get_Rain_vel(Rain_diam)

      implicit none

      real(kind=ip) Get_Rain_vel    ! Rain drop velocity (m/s)
      real(kind=ip) Rain_diam       ! Rain diam (m)
      
       ! Willis,1984
      Get_Rain_vel = 4854.0_ip*Rain_diam*exp(-195.0_ip*Rain_diam)

      end function

!******************************************************************************

      subroutine ThicknessCalculator_WetDepo

      use global_param,  only : &
         M_2_MM

      use Tephra,        only : &
         n_gs_max,DepositDensity

      use mesh,          only : &
        nxmax,nymax,kappa_pd,dx,dy,dz_vec_pd,IsLatLon,kappa_pd,dz_vec_pd

      use Output_Vars,   only : &
        AreaCovered,DepositThickness

      use solution,      only : &
         DepositGranularity

      implicit none

      integer :: i,j
      real(kind=ip) :: dvol
      real(kind=ip) :: WD1,WD2,WD_contrib
      

      ! add contribution of WetDepo to deposit thickness in mm, and area covered
      AreaCovered               = 0.0_ip
      if(.not.IsLatLon) dvol = dx*dy*dz_vec_pd(1)

      if(USE_WASHOUT_LIQ.or.USE_RAINOUT_LIQ)then
        do i=1,nxmax
          do j=1,nymax
            if(IsLatLon)dvol = kappa_pd(i,j,1)
            WD1 = 0.0_ip
            if(USE_WASHOUT_LIQ)then
              DepositGranularity(i,j,1:n_gs_max) = DepositGranularity(i,j,1:n_gs_max) + &
                                            Deposit_Liq_Washout(i,j,1:n_gs_max)
              WD1 = sum(Deposit_Liq_Washout(i,j,1:n_gs_max))
            endif
            WD2 = 0.0_ip
            if(USE_RAINOUT_LIQ)then
              DepositGranularity(i,j,1:n_gs_max) = DepositGranularity(i,j,1:n_gs_max) + &
                                            Deposit_Liq_Rainout(i,j,1:n_gs_max)
              WD2 = sum(Deposit_Liq_Rainout(i,j,1:n_gs_max))
            endif
            WD_contrib = M_2_MM*(WD1+WD2)*dz_vec_pd(1)*1.0e-6_ip*DepositDensity
            DepositThickness(i,j) = DepositThickness(i,j) + WD_contrib
          enddo
        enddo
      endif

      end subroutine ThicknessCalculator_WetDepo


!******************************************************************************

      end module Wet_Deposition
