      module Source_Gas

      use precis_param

      use io_units

      ! Set the number of output variables for this module
      !integer, parameter :: nvar_User2d_static_XY_SrcGas = 1 ! DepositMask
         ! SO2_surf (ppm)
         ! SO2_PBL  (DU CMA 0.9 km)
         ! SO2_TRL  (DU CMA 2.5 km)
         ! SO2_TRM  (DU CMA 7.5 km)
         ! SO2_STL  (DU CMA 17  km)
         ! SO4_surf (ug/m3)
      integer, parameter :: nvar_User2d_XY_SrcGas        = 3 ! SO2_surf,SO2_TRM,SO4
      !integer, parameter :: nvar_User3d_XYGs_SrcGas      = 0
      !integer, parameter :: nvar_User3d_XYZ_SrcGas       = 0
      !integer, parameter :: nvar_User4d_XYZGs_SrcGas     = 0

      real(kind=ip),parameter :: SO2_surf_thresh = 1.0e-6_ip
      real(kind=ip),parameter :: SO2_TRM_thresh  = 1.0e-6_ip
      real(kind=ip),parameter :: SO4_surf_thresh = 1.0e-6_ip

      !character(len=30),dimension(nvar_User2d_static_XY_SrcGas) :: temp_2ds_name_SrcGas
      !character(len=30),dimension(nvar_User2d_static_XY_SrcGas) :: temp_2ds_unit_SrcGas
      !character(len=30),dimension(nvar_User2d_static_XY_SrcGas) :: temp_2ds_lname_SrcGas
      !real(kind=op),    dimension(nvar_User2d_static_XY_SrcGas) :: temp_2ds_MissVal_SrcGas
      !real(kind=op),    dimension(nvar_User2d_static_XY_SrcGas) :: temp_2ds_FillVal_SrcGas

      character(len=30),dimension(nvar_User2d_XY_SrcGas) :: temp_2d_name_SrcGas
      character(len=30),dimension(nvar_User2d_XY_SrcGas) :: temp_2d_unit_SrcGas
      character(len=30),dimension(nvar_User2d_XY_SrcGas) :: temp_2d_lname_SrcGas
      real(kind=op),    dimension(nvar_User2d_XY_SrcGas) :: temp_2d_MissVal_SrcGas
      real(kind=op),    dimension(nvar_User2d_XY_SrcGas) :: temp_2d_FillVal_SrcGas

      ! These are used to keep track of which index in the global list, this
      ! modules output vars corespond to
      !integer :: indx_User2d_static_XY_SrcGas
      integer :: indx_User2d_XY_SrcGas
      !integer :: indx_User3d_XYGs_SrcGas
      !integer :: indx_User3d_XYZ_SrcGas
      !integer :: indx_User4d_XYZGs_SrcGas

      !character(len=30),dimension(nvar_User2d_static_XY_SrcGas) :: SourceGas_2dvarname
      !character(len=30),dimension(nvar_User2d_static_XY_SrcGas) :: SourceGas_2dvarunit
      !real(kind=ip)    ,dimension(nvar_User2d_static_XY_SrcGas) :: SourceGas_2dvar_FillValue
      !integer          ,dimension(nvar_User2d_static_XY_SrcGas) :: SourceGas_2dvar_id

      !character(len=30),dimension(nvar_User3d_XYZ_SrcGas) :: SourceGas_3dvarname
      !character(len=30),dimension(nvar_User3d_XYZ_SrcGas) :: SourceGas_3dvarunit
      !real(kind=ip)    ,dimension(nvar_User3d_XYZ_SrcGas) :: SourceGas_3dvar_FillValue
      !integer          ,dimension(nvar_User3d_XYZ_SrcGas) :: SourceGas_3dvar_id


      !real(kind=sp),dimension(:,:),allocatable :: SnD_meso_last_step_sp
      !real(kind=sp),dimension(:,:),allocatable :: SnD_meso_next_step_sp
      !real(kind=sp),dimension(:,:),allocatable :: P0_meso_last_step_sp
      !real(kind=sp),dimension(:,:),allocatable :: P0_meso_next_step_sp

      character (len=130) :: DepPerimInfile
      integer :: DepMaskCount
      !real(kind=ip) :: ustar
      !real(kind=ip) :: u_star_thresh
      !real(kind=ip) :: Fv_coeff
      !integer       :: FvID

      integer :: iconcen_gas_start
      integer :: ngas_max
      real(kind=ip),dimension(:),allocatable        :: EruptGasMassRate
      integer      ,dimension(:),allocatable        :: EruptGasSpeciesID
      integer      ,dimension(:),allocatable,public :: EruptGasSrcStruc
      real(kind=ip),dimension(:),allocatable        :: EruptGasVentLon
      real(kind=ip),dimension(:),allocatable        :: EruptGasVentLat

      integer,dimension(:),allocatable              :: EruptGasVentLon_i
      integer,dimension(:),allocatable              :: EruptGasVentLat_j

      integer,dimension(:,:)  ,allocatable :: DepositMask
      real(kind=op),dimension(:,:)    ,allocatable :: SO2_surf
      real(kind=op),dimension(:,:)    ,allocatable :: SO2_TRM
      real(kind=op),dimension(:,:)    ,allocatable :: SO4_surf

      integer,dimension(:)    ,allocatable :: GS_GasSpeciesID

      logical       :: Gas_SO2_SO4_convert
      integer       :: Gas_SO2_SO4_convert_ID
      real(kind=ip) :: Gas_SO2_SO4_conversion_HLife
      integer,public :: Gas_H2O_index = 0
      integer,public :: Gas_CO2_index = 0
      integer,public :: Gas_SO2_index = 0
      integer,public :: Gas_SO4_index = 0
      integer,public :: Gas_H2S_index = 0
      integer,public :: Gas_HCL_index = 0
      integer,public :: Gas_Cl_index  = 0
      integer,public :: Gas_ClO_index = 0
      integer,public :: Gas_HF_index  = 0

      logical :: USE_GAS = .false.

      contains

!******************************************************************************

      subroutine input_data_Source_Gas

      use global_param,  only : &
         nmods

      use io_data,       only : &
         infile

      use mesh,          only : &
         nxmax,nymax,nsmax,de,dn,lonLL,latLL

      use Source,        only : &
         neruptions,SourceType,e_Volume,e_Duration

      implicit none

      character(len=3)  :: answer
      character(len=80)  :: linebuffer
      character(len=120) :: llinebuffer
      character(len=130) :: lllinebuffer
      character :: testkey
      integer :: ios !,ioerr
      character(len=20) :: mod_name
      integer :: substr_pos
      integer         :: i
      integer         :: dum1_int,dum2_int,dum3_int
      real(kind=ip)   :: dum1_ip,dum2_ip,dum3_ip

      ! First check if the requested source belongs to this module
      if ((SourceType.eq.'gas').or. &
          (SourceType.eq.'Gas').or. &
          (SourceType.eq.'GAS')) then
        SourceType='gas'
        USE_GAS = .true.
      else
        USE_GAS = .false.
        return
      endif

      ! allocate the variables needed for this custom source
      allocate(DepositMask(nxmax,nymax))
      allocate(EruptGasSpeciesID(neruptions))
      allocate(EruptGasSrcStruc(neruptions))
      allocate(EruptGasMassRate(neruptions))
      allocate(EruptGasVentLon(neruptions))
      allocate(EruptGasVentLat(neruptions))
      allocate(EruptGasVentLon_i(neruptions))
      allocate(EruptGasVentLat_j(neruptions))

      ! If we've made it here, the requested source is a gas source, open
      ! the input file again to get needed info
      open(unit=10,file=infile,status='old',err=1900)

      ! For custom sources, we want to read down to the source block
      !   Read first comment block
      read(10,'(a80)')linebuffer
      read(linebuffer,*)testkey
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer
        !write(global_info,*)linebuffer
        read(linebuffer,*)testkey
      enddo
      ! Read through block 1
      do while(testkey.ne.'#'.and.testkey.ne.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer
        !write(global_info,*)linebuffer
        read(linebuffer,*)testkey
      enddo
      ! Read next block header
      !write(global_info,*)linebuffer
      read(10,'(a80)')llinebuffer
      read(llinebuffer,*)testkey
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')llinebuffer
        !write(global_info,*)llinebuffer
        read(llinebuffer,*)testkey
      enddo

      write(global_info,*)"Start reading gas source."
      write(global_info,*)neruptions," eruptions"
      ! Note: SpeciesID =>   1 = H2O
      !                      2 = CO2
      !                      3 = SO2
      !                      4 = SO4
      !                      5 = H2S
      !                      6 = HCl
      !                      7 = Cl
      !                      8 = ClO
      !                      9 = HF
      !       EruptGasSrcStruc => 1 = distributed surface source (like a region)
      !                      2 = point on surface (vent)
      !                      3 = line on surface (fissure)
      !                      4 = vertical line source
      !                      5 = vertical profile source
      !       EruptGasMassRate      = tonnes/day

      do i=1,neruptions
        !read start time, duration, plume height, volume of each pulse
        read(llinebuffer,*,err=1910) dum1_int,dum2_int,dum3_int,dum1_ip, &
                              dum2_ip,dum3_ip,&
                              EruptGasSpeciesID(i),EruptGasSrcStruc(i),EruptGasMassRate(i)
        !write(global_info,*)i,llinebuffer
        if(EruptGasSrcStruc(i).eq.1)then
          if(i.gt.1)then
            write(global_info,*)"ERROR: Only one distributed source region allowed,"
            write(global_info,*)"       which must be listed first."
            stop 1
          endif
          ! Read name of file outlining the contour of the region
          read(10,'(a130)')lllinebuffer
          DepPerimInfile = adjustl(trim(lllinebuffer))
        elseif(EruptGasSrcStruc(i).eq.2)then
          ! Read the lon/lat of the point on surface
          read(llinebuffer,*,err=1910) dum1_int,dum2_int,dum3_int,dum1_ip, &
                                dum2_ip,dum3_ip,&
                                EruptGasSpeciesID(i),EruptGasSrcStruc(i),EruptGasMassRate(i), &
                                EruptGasVentLon(i),EruptGasVentLat(i)

          if(EruptGasVentLon(i).lt.-360.0_ip)then
            write(global_info,*)"ERROR: Vent longitude is less than -360.0"
            stop 1
          elseif(EruptGasVentLon(i).gt.360.0_ip)then
            write(global_info,*)"ERROR: Vent longitude is greater than 360.0"
            stop 1
          elseif(EruptGasVentLon(i).lt.0.0_ip)then
            EruptGasVentLon(i) = EruptGasVentLon(i) + 360.0_ip
          endif
          if(EruptGasVentLat(i).lt.-90.0_ip)then
            write(global_info,*)"ERROR: Vent latitude is less than -90.0"
            stop 1
          elseif(EruptGasVentLat(i).gt.90.0_ip)then
            write(global_info,*)"ERROR: Vent latitude is greater than 90.0"
            stop 1
          endif
          ! Need to do some checking here on input coordinates, mapping to
          ! compuataional grid, etc.
          EruptGasVentLon_i(i) = int((EruptGasVentLon(i)-lonLL)/de) + 1
          EruptGasVentLat_j(i) = int((EruptGasVentLat(i)-latLL)/dn) + 1
          write(global_info,*)EruptGasVentLon(i),lonLL,de,EruptGasVentLon_i(i)
          write(global_info,*)EruptGasVentLat(i),latLL,dn,EruptGasVentLat_j(i)
          !stop 1

        elseif(EruptGasSrcStruc(i).eq.3)then
          ! Read the start lon/lat and end lon/lat of fissure
          write(global_info,*)"Gas Fissure source not yet implementes"
          stop 1
        elseif(EruptGasSrcStruc(i).eq.4)then
          ! Read the lon/lat and height of line
          write(global_info,*)"Gas Line source not yet implementes"
          stop 1
        elseif(EruptGasSrcStruc(i).eq.5)then
          ! Read the profile of Mass fractions
          write(global_info,*)"Gas Profile source not yet implementes"
          stop 1
        endif
          ! Read the next line 
        read(10,'(a130)')llinebuffer

      enddo
      ! Initialize some eruption values
      e_Duration  = 0.0_ip
      e_Volume    = 0.0_ip

      ! Now read to the end of the input file and read the Optional Modudle
      ! block
      write(global_info,*)"    Searching for OPTMOD=SRC_GAS"
      ! Example block
      !OPTMOD=SRC_GAS
      !2                                       # ngas_max  (number of gas species)
      !3 4                                     # list of gas indicies
      !yes                                     # Convert SO2 to SO4
      !1                                       # Convertion scheme (1=fixed decay rate)
      !6.0                                     # Half-life in hours

      nmods = 0
      read(10,'(a80)',iostat=ios)linebuffer
      do while(ios.eq.0)
        read(10,'(a80)',iostat=ios)linebuffer

        !read(linebuffer,*)testkey
        substr_pos = index(linebuffer,'OPTMOD')
        if(substr_pos.eq.1)then
          !write(global_info,*)"Found OPTMOD again"
          ! found an optional module
          !  Parse for the keyword
          read(linebuffer,1104)mod_name
          if(adjustl(trim(mod_name)).eq.'SRC_GAS')then
            !write(global_info,*)"Found SRC_GAS block again"
            read(10,'(a80)',iostat=ios)linebuffer
            read(linebuffer,*)ngas_max
            write(global_info,*)"ngas_max = ",ngas_max
            if(ngas_max.lt.1)then
              write(global_info,*)"No gas species defined."
              ngas_max = 0
            elseif(ngas_max.gt.9)then
              write(global_info,*)"To many gas species."
              stop 1
            else
              write(global_info,*)"Reading list of gas species IDs."
              iconcen_gas_start = nsmax
            endif
            allocate(GS_GasSpeciesID(ngas_max))
            ! Get the list of species to trace
            read(10,'(a80)',iostat=ios)linebuffer
            read(linebuffer,*)GS_GasSpeciesID(1:ngas_max)
              ! HFS do some error-checking here to verify that all the source
              ! terms are accommodated
            do i=1,ngas_max
              if(GS_GasSpeciesID(i).eq.1)then
                Gas_H2O_index = iconcen_gas_start+i
              elseif(GS_GasSpeciesID(i).eq.2)then
                Gas_CO2_index = iconcen_gas_start+i
              elseif(GS_GasSpeciesID(i).eq.3)then
                Gas_SO2_index = iconcen_gas_start+i
              elseif(GS_GasSpeciesID(i).eq.4)then
                Gas_SO4_index = iconcen_gas_start+i
              elseif(GS_GasSpeciesID(i).eq.5)then
                Gas_H2S_index = iconcen_gas_start+i
              elseif(GS_GasSpeciesID(i).eq.6)then
                Gas_HCL_index = iconcen_gas_start+i
              elseif(GS_GasSpeciesID(i).eq.7)then
                Gas_Cl_index  = iconcen_gas_start+i
              elseif(GS_GasSpeciesID(i).eq.8)then
                Gas_ClO_index = iconcen_gas_start+i
              elseif(GS_GasSpeciesID(i).eq.9)then
                Gas_HF_index  = iconcen_gas_start+i
              else
                write(global_info,*)"ERROR: unknown gas code."
                write(global_info,*)linebuffer
                write(global_info,*)GS_GasSpeciesID
                stop 1
              endif
            enddo

            ! Now check if we need to do an SO2 to SO4 conversion
            read(10,'(a80)',iostat=ios)linebuffer
            read(linebuffer,'(a3)') answer
            Gas_SO2_SO4_convert = .false.
            if (answer.eq.'yes') then
              Gas_SO2_SO4_convert = .true.
            endif
            if(Gas_SO2_SO4_convert)then
              ! If yes, then read the conversion scheme
              !  1 = fixed decay rate
              !  2 = something more complicated requiring moisture, radiation
              read(10,'(a80)',iostat=ios)linebuffer
              read(linebuffer,*)Gas_SO2_SO4_convert_ID
              ! Now read something about the convertion process
              if(Gas_SO2_SO4_convert_ID.eq.1)then
                read(10,'(a80)',iostat=ios)linebuffer
                read(linebuffer,*)Gas_SO2_SO4_conversion_HLife
              else
                write(global_info,*)"Only fixed SO2 -> SO4 decay rate implemented."
                stop 1
              endif
            endif

            write(global_info,*)"Gas_H2O_index",Gas_H2O_index
            write(global_info,*)"Gas_CO2_index",Gas_CO2_index
            write(global_info,*)"Gas_SO2_index",Gas_SO2_index
            write(global_info,*)"Gas_SO4_index",Gas_SO4_index
            write(global_info,*)"Gas_H2S_index",Gas_H2S_index
            write(global_info,*)"Gas_HCL_index",Gas_HCL_index
            write(global_info,*)"Gas_Cl_index",Gas_Cl_index
            write(global_info,*)"Gas_ClO_index",Gas_ClO_index
            write(global_info,*)"Gas_HF_index",Gas_HF_index
            write(global_info,*)"Gas_SO2_SO4_convert",Gas_SO2_SO4_convert
            write(global_info,*)"Gas_SO2_SO4_convert_ID and rate",&
              Gas_SO2_SO4_convert_ID,Gas_SO2_SO4_conversion_HLife
            exit
          endif
          !OPTMOD_names(nmods) = adjustl(trim(mod_name))
        endif
1104    format(7x,a20)
      enddo

!2010  continue
      close(10)

      nsmax = nsmax + ngas_max

      return

1900  write(global_info,*)  'error: cannot find input file: ',infile
      write(global_info,*)  'Program stopped'
      write(global_log,*)  'error: cannot find input file: ',infile
      write(global_log,*)  'Program stopped'
      stop 1
1910  write(global_info,*)  'error reading start time, duration, height or',&
                  ' volume of an eruptive pulse.  Program stopped'
      stop 1

      end subroutine input_data_Source_Gas

!******************************************************************************

      subroutine Allocate_Source_Gas

      use io_data,       only : &
         nvar_User2d_XY

      use mesh,          only : &
         nxmax,nymax,ts1

      use Source,        only : &
         SourceNodeFlux_Area

      implicit none

      write(global_info,*)"--------------------------------------------------"
      write(global_info,*)"---------- ALLOCATE_GAS -----------------------"
      write(global_info,*)"--------------------------------------------------"

      ! Set the start indecies
      !indx_User2d_static_XY_SrcGas = nvar_User2d_static_XY
      indx_User2d_XY_SrcGas        = nvar_User2d_XY
      !indx_User3d_XYGs_SrcGas      = nvar_User3d_XYGs
      !indx_User3d_XYZ_SrcGas       = nvar_User3d_XYZ
      !indx_User4d_XYZGs_SrcGas     = nvar_User4d_XYZGs

      !! Surface Roughness
      !!temp_2ds_name_SrcGas(1)    = "SurfRough"
      !!temp_2ds_lname_SrcGas(1)   = "Roughness Length"
      !!temp_2ds_unit_SrcGas(1)    = "m"
      !!temp_2ds_MissVal_SrcGas(1) = -9999.0_ip
      !!temp_2ds_FillVal_SrcGas(1) = -9999.0_ip
      !! Deposite Mask
      !temp_2ds_name_SrcGas(1)    = "DepoMask"
      !temp_2ds_lname_SrcGas(1)   = "Source region"
      !temp_2ds_unit_SrcGas(1)    = "flag"
      !temp_2ds_MissVal_SrcGas(1) = -9999.0_ip
      !temp_2ds_FillVal_SrcGas(1) = -9999.0_ip

      !! SO2_surf
      temp_2d_name_SrcGas(1)    = "SO2_surf"
      temp_2d_lname_SrcGas(1)   = "SO2 concentraion at surface"
      temp_2d_unit_SrcGas(1)    = "ppm"
      temp_2d_MissVal_SrcGas(1) = -9999.0_op
      temp_2d_FillVal_SrcGas(1) = -9999.0_op

      !! SO2_TRM
      temp_2d_name_SrcGas(2)    = "SO2_TRM"
      temp_2d_lname_SrcGas(2)   = "SO2 column Mid-tropo"
      temp_2d_unit_SrcGas(2)    = "DU"
      temp_2d_MissVal_SrcGas(2) = -9999.0_op
      temp_2d_FillVal_SrcGas(2) = -9999.0_op

      !! SO4_surf
      temp_2d_name_SrcGas(3)    = "SO4_surf"
      temp_2d_lname_SrcGas(3)   = "SO4 concentraion at surface"
      temp_2d_unit_SrcGas(3)    = "ug/m3"
      temp_2d_MissVal_SrcGas(3) = -9999.0_op
      temp_2d_FillVal_SrcGas(3) = -9999.0_op

      !nvar_User2d_static_XY = nvar_User2d_static_XY + nvar_User2d_static_XY_SrcGas
      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_SrcGas
      !nvar_User3d_XYGs      = nvar_User3d_XYGs      + nvar_User3d_XYGs_SrcGas
      !nvar_User3d_XYZ       = nvar_User3d_XYZ       + nvar_User3d_XYZ_SrcGas
      !nvar_User4d_XYZGs     = nvar_User4d_XYZGs     + nvar_User4d_XYZGs_SrcGas

      allocate(SO2_surf(nxmax,nymax))
      allocate(SO2_TRM(nxmax,nymax))
      allocate(SO4_surf(nxmax,nymax))

      !allocate(SnD_meso_last_step_sp(nxmax,nymax))
      !allocate(SnD_meso_next_step_sp(nxmax,nymax))
      !allocate(P0_meso_last_step_sp(nxmax,nymax))
      !allocate(P0_meso_next_step_sp(nxmax,nymax))

      allocate(SourceNodeFlux_Area(1:nxmax,1:nymax,1:ngas_max))

      end subroutine Allocate_Source_Gas

!******************************************************************************

      subroutine Prep_output_Source_Gas

      use mesh,          only : &
         nxmax,nymax,nzmax,ts1

      use solution,      only : &
         concen_pd

      use Output_Vars,   only : &
         var_User2d_XY_name,var_User2d_XY_unit,var_User2d_XY_lname,&
         var_User2d_XY_MissVal,var_User2d_XY_FillVal,var_User2d_XY

      implicit none

      integer :: i,j,k,indx
      real(kind=op) :: dum_op

      !do i=1,nvar_User2d_static_XY_SrcGas
      !  indx = indx_User2d_static_XY_SrcGas+i
      !  var_User2d_static_XY_name(indx) = temp_2ds_name_SrcGas(i)
      !  var_User2d_static_XY_unit(indx) = temp_2ds_unit_SrcGas(i)
      !  var_User2d_static_XY_lname(indx)= temp_2ds_lname_SrcGas(i)
      !  var_User2d_static_XY_MissVal(indx)= temp_2ds_MissVal_SrcGas(i)
      !  var_User2d_static_XY_FillVal(indx)= temp_2ds_FillVal_SrcGas(i)
      !  !if(i.eq.1)var_User2d_static_XY(:,:,i) = SurfRough
      !  if(i.eq.1)var_User2d_static_XY(:,:,indx) = DepositMask
      !enddo
      SO2_surf(:,:) = -9999.0_op
      SO2_TRM(:,:)  = -9999.0_op
      SO4_surf(:,:) = -9999.0_op
      do i=1,nxmax
        do j=1,nymax

          dum_op = real(concen_pd(i,j,1,Gas_SO2_index,ts1),kind=op) / &
                                       2.62_op/1000.0_op  ! convert to ppm at 1 atm and 25 C
          if(dum_op.ge.SO2_surf_thresh)then
            SO2_surf(i,j) = dum_op
          endif
          dum_op = 0.0_op
          do k = 1,nzmax
            dum_op = dum_op + real(concen_pd(i,j,k,Gas_SO2_index,ts1),kind=op)
          enddo
          if(dum_op.ge.SO2_TRM_thresh)then
            SO2_TRM(i,j) = dum_op
          endif
          dum_op = real(concen_pd(i,j,1,Gas_SO4_index,ts1),kind=op)
          if(dum_op.ge.SO4_surf_thresh)then
            SO4_surf(i,j) = dum_op
          endif
        enddo
      enddo

      do i=1,nvar_User2d_XY_SrcGas
        indx = indx_User2d_XY_SrcGas+i
        var_User2d_XY_name(indx) = temp_2d_name_SrcGas(i)
        var_User2d_XY_unit(indx) = temp_2d_unit_SrcGas(i)
        var_User2d_XY_lname(indx)= temp_2d_lname_SrcGas(i)
        var_User2d_XY_MissVal(indx)= temp_2d_MissVal_SrcGas(i)
        var_User2d_XY_FillVal(indx)= temp_2d_FillVal_SrcGas(i)
        if(i.eq.1)var_User2d_XY(:,:,indx) = SO2_surf(:,:)
        if(i.eq.2)var_User2d_XY(:,:,indx) = SO2_TRM(:,:)
        if(i.eq.3)var_User2d_XY(:,:,indx) = SO4_surf(:,:)
      enddo

      end subroutine Prep_output_Source_Gas

!******************************************************************************

      subroutine Deallocate_Source_Gas

      use Source,        only : &
         SourceNodeFlux_Area

      implicit none

      if(allocated(DepositMask)) deallocate(DepositMask)

#ifdef USEPOINTERS
      if(associated(SourceNodeFlux_Area))   deallocate(SourceNodeFlux_Area)
#else
      if(allocated(SourceNodeFlux_Area))   deallocate(SourceNodeFlux_Area)
#endif


      end subroutine Deallocate_Source_Gas

!******************************************************************************

!******************************************************************************

      subroutine Read_Deposit_Perimeter_Gas

      use mesh,          only : &
         nxmax,nymax,ivent,jvent,x_cc_pd,y_cc_pd,lon_cc_pd,lat_cc_pd,&
         A3d_iprojflag,A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
         A3d_k0_scale,A3d_radius_earth,de,dn,dx,dy,IsLatLon

      use projection,    only : &
           PJ_proj_for

      implicit none

      integer :: npoints
      real(kind=ip),dimension(:),allocatable :: DepPerm_lon
      real(kind=ip),dimension(:),allocatable :: DepPerm_lat

      !real(kind=ip),dimension(:),allocatable :: DepPerm_x
      !real(kind=ip),dimension(:),allocatable :: DepPerm_y

      real(kind=ip) :: DepPerm_x_min,DepPerm_x_max
      real(kind=ip) :: DepPerm_y_min,DepPerm_y_max
      integer       :: DepPerm_i_min,DepPerm_i_max
      integer       :: DepPerm_j_min,DepPerm_j_max

      real(kind=ip),dimension(2) :: testpoint
      real(kind=ip),dimension(:,:),allocatable :: BCpos ! boundary corner position
      integer      ,dimension(:,:),allocatable :: Belem ! boundary element BC IDs

      real(kind=ip),dimension(nxmax) :: loc_x
      real(kind=ip),dimension(nymax) :: loc_y
      real(kind=ip) :: loc_dx, loc_dy
      real(kind=dp) :: lat_in,lon_in,xout,yout

      integer :: i,j

      write(global_info,*)"Opening ",DepPerimInfile
      open(unit=20,file=DepPerimInfile)
      read(20,*)npoints
      allocate(DepPerm_lon(npoints))
      allocate(DepPerm_lat(npoints))
      allocate(BCpos(2,npoints))
      allocate(Belem(2,npoints))

      ! Read in the points and assign boundary elements with the corner IDs
      do i = 1,npoints
        read(20,*)DepPerm_lon(i),DepPerm_lat(i)
        if(DepPerm_lon(i).lt.0.0_ip)DepPerm_lon(i)=DepPerm_lon(i)+360.0_ip
        Belem(1,i)=i
        Belem(2,i)=i+1
      enddo
        ! The end point of the last boundary element is BC #1
      Belem(2,npoints)=1

      if(IsLatLon)then
        BCpos(1,:) = DepPerm_lon(:)
        BCpos(2,:) = DepPerm_lat(:)
        loc_x(1:nxmax) = lon_cc_pd(1:nxmax)
        loc_y(1:nymax) = lat_cc_pd(1:nymax)
        loc_dx = de
        loc_dy = dn
      else
        ! We need to project each point on the deposit perimeter onto the
        ! computational grid
        do i = 1,npoints
          lon_in = real(DepPerm_lon(i),kind=dp)
          lat_in = real(DepPerm_lat(i),kind=dp)
          call PJ_proj_for(lon_in,lat_in,A3d_iprojflag, &
                      A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,A3d_k0_scale,A3d_radius_earth, &
                      xout,yout)
          BCpos(1,i) = real(xout,kind=ip)
          BCpos(2,i) = real(yout,kind=ip)
        enddo
        loc_x(1:nxmax) = x_cc_pd(1:nxmax)
        loc_y(1:nymax) = y_cc_pd(1:nymax)
        loc_dx = dx
        loc_dy = dy
      endif
        ! Get the min and max extents of the deposit
      DepPerm_x_min = minval(BCpos(1,:))
      DepPerm_x_max = maxval(BCpos(1,:))
      DepPerm_y_min = minval(BCpos(2,:))
      DepPerm_y_max = maxval(BCpos(2,:))
        ! And the indicies on the computational grid braketing these points
      DepPerm_i_min =   floor((DepPerm_x_min-loc_x(1))/loc_dx) - 1
      DepPerm_i_max = ceiling((DepPerm_x_max-loc_x(1))/loc_dx) + 1
      DepPerm_j_min =   floor((DepPerm_y_min-loc_y(1))/loc_dy) - 1
      DepPerm_j_max = ceiling((DepPerm_y_max-loc_y(1))/loc_dy) + 1

      !  Loop through all computational grid points near the deposit and flag
      !  those inside as 1 in DepositMask
      DepositMask = 0
      DepMaskCount = 0
      do i = DepPerm_i_min,DepPerm_i_max
        do j = DepPerm_j_min,DepPerm_j_max
          testpoint(1) = loc_x(i)
          testpoint(2) = loc_y(j)
          if(IsIn(testpoint,npoints,BCpos,npoints,Belem))then
            DepositMask(i,j) = 1
            DepMaskCount = DepMaskCount + 1
          elseif(ivent.eq.i.and.jvent.eq.j)then
            DepositMask(i,j) = 1
          else
            DepositMask(i,j) = 0
          endif
        enddo
      enddo

      write(global_info,*)ivent,jvent, DepositMask(ivent,jvent) , DepMaskCount
      close(20)

      end subroutine Read_Deposit_Perimeter_Gas

!******************************************************************************

      subroutine Set_Gas_Meso(Load_MesoSteps,Interval_Frac)

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=ip),intent(in) :: Interval_Frac

      !integer :: ivar
      !integer :: i,j,k

      !logical,save :: first_time = .true.

!      if(Load_MesoSteps)then
!        if(first_time)then
!          ! Need to fill _last_step_sp
!          !  First fill next step so that outside this 'first_time' loop, the
!          !  'next' can be copied to the 'last'
!
!           ! Now resample Friction Velocity onto computational grid
!          !ivar = 13
!          !call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
!          !FricVel_meso_next_step_sp = MR_dum2d_comp
!          !Note:  This was read in variable_diffusivity on the met grid and now
!          !       only needs to be resampled
!          MR_dum2d_met = FricVel_meso_next_step_Met_sp
!          call MR_Regrid_Met2d_to_Comp2d
!          FricVel_meso_next_step_sp = MR_dum2d_comp
!
!          ivar = 15
!          call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
!          SnD_meso_next_step_sp  = MR_dum2d_comp
!
!          ivar = 44
!          call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
!          P0_meso_next_step_sp  = MR_dum2d_comp*3600.0_sp 
!
!          first_time = .false.
!        endif ! first_time
!
!        FricVel_meso_last_step_sp      = FricVel_meso_next_step_sp
!        SnD_meso_last_step_sp          = SnD_meso_next_step_sp
!        P0_meso_last_step_sp           = P0_meso_next_step_sp
!
!        ! Need to fill _next_step_sp
!         ! Now resample onto computational grid
!        !ivar = 13
!        !call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
!        !FricVel_meso_next_step_sp = MR_dum2d_comp
!        MR_dum2d_met = FricVel_meso_next_step_Met_sp
!        call MR_Regrid_Met2d_to_Comp2d
!        FricVel_meso_next_step_sp = MR_dum2d_comp
!
!        ivar = 15
!        call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
!        SnD_meso_next_step_sp  = MR_dum2d_comp*1.0_sp
!
!        ivar = 44
!        call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
!        P0_meso_next_step_sp  = MR_dum2d_comp*3600.0_sp ! convert from kg/m2/s to mm/hr
!
!      endif

      end subroutine Set_Gas_Meso

!******************************************************************************

      subroutine Set_Gas_Flux

      use Source,        only : &
         SourceNodeFlux_Area

      use mesh,          only : &
        nxmax,nymax

      use Source,        only : &
         neruptions

      implicit none

      integer :: i,j,ie,idx
      real(kind=ip) :: Fv

      do ie=1,neruptions
        ! Get the chem species for this eruptive pulse
        select case (EruptGasSpeciesID(ie))
        case(1)
          idx = Gas_H2O_index
        case(2)
          idx = Gas_CO2_index
        case(3)
          idx = Gas_SO2_index
        case(4)
          idx = Gas_SO4_index
        case(5)
          idx = Gas_H2S_index
        case(6)
          idx = Gas_HCL_index
        case(7)
          idx = Gas_Cl_index 
        case(8)
          idx = Gas_ClO_index
        case(9)
          idx = Gas_HF_index 
        end select

        if(EruptGasSrcStruc(ie).eq.1)then
          ! Distributed region
          do i = 1,nxmax
            do j = 1,nymax
              if(DepositMask(i,j).gt.0)then
                Fv = EruptGasMassRate(1)/real(DepMaskCount,kind=ip)  ! Divide by number o cells
                Fv = Fv * 1000.0_ip /24.0_ip          ! Convert to kg/hr for this cell
              else
                Fv = 0.0_ip
              endif
                ! Mass contribution in kg/hr
              SourceNodeFlux_Area(i,j,idx)= Fv
            enddo
          enddo
        elseif(EruptGasSrcStruc(ie).eq.2)then
          ! Point source
            ! Mass contribution in kg/hr
          SourceNodeFlux_Area(EruptGasVentLon_i(ie),EruptGasVentLat_j(ie),idx) = &
            EruptGasMassRate(ie) * 1000.0_ip/24.0_ip
        elseif(EruptGasSrcStruc(ie).eq.3)then
          ! Fissure source
        elseif(EruptGasSrcStruc(ie).eq.4)then
          ! Verftical line source
        elseif(EruptGasSrcStruc(ie).eq.5)then
          ! Vertical profile source
        endif
      enddo !neruptions

      end subroutine Set_Gas_Flux

!******************************************************************************

      subroutine Set_concen_Gas

      use mesh,          only : &
         nxmax,nymax,ts0,ts1,kappa_pd

      use solution,      only : &
         concen_pd

      use Source,        only : &
         SourceNodeFlux_Area

      !use Topography

      use time_data,     only : &
         dt

      implicit none

      integer :: i,j,k

      do i = 1,nxmax
        do j = 1,nymax
          !k = topo_indx(i,j)
          k = 1
          concen_pd(i,j,k,iconcen_gas_start+1:iconcen_gas_start+ngas_max,ts0) =                     &
                   concen_pd(i,j,k,iconcen_gas_start+1:iconcen_gas_start+ngas_max,ts0)  &
                   + dt*SourceNodeFlux_Area(i,j,1:ngas_max)/kappa_pd(i,j,1)
        enddo
      enddo

      end subroutine Set_concen_Gas



!******************************************************************************

      subroutine Gas_Chem_Convert

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,ts0,ts1

      use solution,      only : &
         concen_pd

      use time_data,     only : &
         dt

      implicit none

      integer :: i,j,k
      real(kind=ip) :: cofac
      real(kind=ip) :: scrub

      if(Gas_SO2_SO4_convert)then
        if(Gas_SO2_SO4_convert_ID.eq.1)then
          ! This is the fixed decay rate
          !  Note: half-life should be around 6 hours so we will be far from
          !  this being the most restrictive on the time step.  So just use the
          !  global dt
          cofac = dt*(0.693147180559945_ip/Gas_SO2_SO4_conversion_HLife)
          !write(global_info,*)"Setting up to scrub: ",cofac,Gas_SO2_index,Gas_SO4_index
          do i = 1,nxmax
            do j = 1,nymax
              do k = 1,nzmax
                  ! Get converted mass
                scrub = cofac*concen_pd(i,j,k,Gas_SO2_index,ts0)
                  ! Decrement SO2
                concen_pd(i,j,k,Gas_SO2_index,ts1) = &
                concen_pd(i,j,k,Gas_SO2_index,ts0) - scrub
                  ! Increment SO4
                concen_pd(i,j,k,Gas_SO4_index,ts1) =                     &
                concen_pd(i,j,k,Gas_SO4_index,ts0) + scrub
              enddo
            enddo
          enddo
        else
          ! No other SO2->SO4 conversion scheme implementes
          write(global_info,*)"ERROR: Gas_SO2_SO4_convert_ID.ne.1", &
            Gas_SO2_SO4_convert_ID
          stop 1
        endif
      endif

      concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine Gas_Chem_Convert

!******************************************************************************


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function IsIn(testpoint,Num_Vert,VertPos,Num_Elem,Elem)
        !  This function returns true if testpoint is
        !  located within the Particle Object.

        ! This should be set up with a linked-list for elements so that
        ! only the elements in specific rows and columns are checked.
        implicit none

        logical                             :: IsIn
        real(kind=ip),dimension(2)          :: testpoint
        integer                             :: Num_Vert
        real(kind=ip),dimension(2,Num_Vert) :: VertPos
        integer                             :: Num_Elem
        integer,dimension(2,Num_Elem)       :: Elem

        real(kind=ip),dimension(2) :: minpoint
        real(kind=ip),dimension(2) :: maxpoint
        integer :: num_uvert_cross
        integer :: num_dvert_cross
        integer :: num_lhorz_cross
        integer :: num_rhorz_cross
        integer :: di
        real(kind=ip) :: tol
        real(kind=ip) :: Tol_Dist
        real(kind=ip) :: B_length
        real(kind=ip) :: P_dot_B
        real(kind=ip) :: bp_top
        real(kind=ip) :: bp_bot
        real(kind=ip) :: bp_left
        real(kind=ip) :: bp_right
        real(kind=ip) :: Dist
        real(kind=ip) :: D2
        real(kind=ip),dimension(2) :: bc1_coord
        real(kind=ip),dimension(2) :: bc2_coord
        real(kind=ip),dimension(2) :: Bvec_1
        real(kind=ip),dimension(2) :: Bn
        real(kind=ip),dimension(2) :: P_shift
        real(kind=ip),dimension(2) :: pPos
        real(kind=ip),dimension(2) :: Fp_Pos
        integer :: beli
        logical :: vert_test
        logical :: horz_test

        tol = 0.499_ip
        Tol_Dist = 0.0_ip !tol*AVG_SPACE

        do di = 1,2
          minpoint(di) = minval(VertPos(di,:))
          maxpoint(di) = maxval(VertPos(di,:))
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!  Two-dimensional particle object
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! The strategy here is to send out rays left and right in
            ! x and up and down in y and count the number of boundary
            ! elements each ray intersects.  If testpoint is within
            ! the PO, then each of of these counts should be odd.
          num_uvert_cross = 0
          num_dvert_cross = 0
          num_lhorz_cross = 0
          num_rhorz_cross = 0
          if(testpoint(1).lt.maxpoint(1).and.  &
             testpoint(2).lt.maxpoint(2))then
            do beli = 1,Num_Elem
              bc1_coord = VertPos(:,Elem(1,beli))
              bc2_coord = VertPos(:,Elem(2,beli))

                ! get the top, bottom, left-most and right-most
                ! coordinate values
              bp_top   = max(bc1_coord(2),bc2_coord(2))
              bp_bot   = min(bc1_coord(2),bc2_coord(2))
              bp_left  = min(bc1_coord(1),bc2_coord(1))
              bp_right = max(bc1_coord(1),bc2_coord(1))

                ! Check to see if testpoint can be projected horz or
                ! vert onto beli
              vert_test = testpoint(1) .ge. bp_left   .and.  &
                          testpoint(1) .lt.  bp_right
              horz_test = testpoint(2) .ge. bp_bot    .and.  &
                          testpoint(2) .lt.  bp_top
              if(vert_test.or.horz_test)then
                  ! Get projection particle
                Bvec_1 = bc2_coord - bc1_coord
                B_length = sqrt(dot_product(Bvec_1,Bvec_1))
                  ! norm points to left of boundary element
                Bn(1) = -Bvec_1(2)/B_length
                Bn(2) =  Bvec_1(1)/B_length
                  ! get projection point
                  ! this isn't the intersection since the projection
                  ! is orthogonal to bel (not along x or y), but it is
                  ! sufficient to determine which side (pos x-ray or
                  ! neg one)
                P_shift = testpoint - bc1_coord
                P_dot_B = dot_product(P_shift,Bn)
                Fp_Pos = -P_dot_B*Bn
                pPos = Fp_Pos + testpoint
                D2 = dot_product((testpoint-pPos),(testpoint-pPos))
                Dist = sqrt(D2)
              endif
              if(vert_test .and. Dist .ge. Tol_Dist)then
                  ! Vertical line crosses this Bel either on left or
                  ! right
                if(testpoint(2).gt.pPos(2))then
                    ! element is below point
                  num_dvert_cross = num_dvert_cross + 1
                else
                    ! element is above point
                  num_uvert_cross = num_uvert_cross + 1
                endif
              endif
              if(horz_test .and. Dist .ge. Tol_Dist)then
                  ! Horizontal line crosses this Bel either on top or
                  ! bottom
                if(testpoint(1).gt.pPos(1))then
                    ! element is right of point
                  num_rhorz_cross = num_rhorz_cross + 1
                else
                    ! element is left of point
                  num_lhorz_cross = num_lhorz_cross + 1
                endif
              endif
            enddo

            if(IsOdd(num_uvert_cross) .and. &
               IsOdd(num_dvert_cross) .and. &
               IsOdd(num_lhorz_cross) .and. &
               IsOdd(num_rhorz_cross))then
              IsIn = .true.
            else
              IsIn = .false.
            endif
          else
              ! testpoint is outside of PO bounding box
            IsIn = .False.
          endif
      end function IsIn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function IsOdd(n)
        !  This function returns true if iPos is
        !  located within the Particle Object.
        implicit none

        logical IsOdd
        integer :: n

        if(mod(n,2).eq.1)then
          IsOdd = .true.
        else
          IsOdd = .false.
        endif
        RETURN
      end function IsOdd


      end module Source_Gas

