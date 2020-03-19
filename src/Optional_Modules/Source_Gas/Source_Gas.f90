      module Source_Gas

      use precis_param

      use io_units

      ! Set the number of output variables for this module
      integer :: nvar_User2d_XY_SrcGas
      !integer, parameter :: nvar_User3d_XYGs_SrcGas      = 0
      !integer, parameter :: nvar_User3d_XYZ_SrcGas       = 0
      !integer, parameter :: nvar_User4d_XYZGs_SrcGas     = 0

      real(kind=ip),parameter :: Surf_Conc_tresh     = 1.0e-6_ip  ! kg/km^3 or ug/m^3
      real(kind=ip),parameter :: Surf_Conc_tresh_ppm = 1.0e-3_ip  ! ppm
      real(kind=ip),parameter :: VertColDens_tresh   = 1.0e-6_ip  ! kg/km^2 or mg/m^2
      real(kind=ip),parameter :: VertColDens_tresh_DU= 2.0e-1_ip  ! DU
                                                                  ! 0.4 for OMI, IASI

      character(len=30),dimension(:),allocatable :: temp_2d_name_SrcGas
      character(len=30),dimension(:),allocatable :: temp_2d_unit_SrcGas
      character(len=30),dimension(:),allocatable :: temp_2d_lname_SrcGas
      real(kind=op),    dimension(:),allocatable :: temp_2d_MissVal_SrcGas
      real(kind=op),    dimension(:),allocatable :: temp_2d_FillVal_SrcGas

      ! These are used to keep track of which index in the global list, this
      ! modules output vars correspond to
      integer :: indx_User2d_XY_SrcGas

      logical            :: setSrcGasSurf_Perim
      character (len=60) :: SrcGasSurf_PerimInfile
      integer            :: SrcGasSurf_MaskCount

      integer :: iconcen_gas_start
      integer :: ngas_max
      real(kind=ip),dimension(:),allocatable        :: EruptGasMassRate
      integer      ,dimension(:),allocatable        :: EruptGasSpeciesID
      integer      ,dimension(:),allocatable,public :: EruptGasSrcStruc
      real(kind=ip),dimension(:),allocatable        :: EruptGasVentLon
      real(kind=ip),dimension(:),allocatable        :: EruptGasVentLat

      integer,dimension(:),allocatable              :: EruptGasVentLon_i
      integer,dimension(:),allocatable              :: EruptGasVentLat_j

      integer      ,dimension(:,:)  ,allocatable :: SrcGasSurf_Mask
      real(kind=op),dimension(:,:,:),allocatable :: SrcGas_SurfConc
      real(kind=op),dimension(:,:,:),allocatable :: SrcGas_VertColDens

      integer,dimension(:)    ,allocatable :: GS_GasSpeciesID
      logical,dimension(:)    ,allocatable :: GS_GasSpecies_OutputSurfConc
      logical,dimension(:)    ,allocatable :: GS_GasSpecies_OutputVertColDens

      ! Note: SpeciesID =>   # Primary gas emmission species as listed in Encyc. of Volc. p805
      integer,parameter :: MAXIMUM_NUM_SPECIES = 60  ! current index goes to 53, but this gives a bit more space
      character(len=30),dimension(MAXIMUM_NUM_SPECIES) :: GS_GasSpecies_name
      character(len=30),dimension(MAXIMUM_NUM_SPECIES) :: GS_GasSpecies_lname
      logical          ,dimension(MAXIMUM_NUM_SPECIES) :: GS_GasSpecies_IsAerosol
      real(kind=ip)    ,dimension(MAXIMUM_NUM_SPECIES) :: GS_GasSpecies_MolWeight ! in g/mol

      integer,public :: Gas_H2O_index = 0 ! water vapor 
      integer,public :: Gas_CO2_index = 0 ! carbon dioxide
      integer,public :: Gas_SO2_index = 0 ! sulfur dioxide
      integer,public :: Gas_H2S_index = 0 ! hydrogen sulfide
      integer,public :: Gas_HCl_index = 0 ! hydrogen cloride
      integer,public :: Gas_HF_index  = 0 ! hydrogen fluride
      integer,public :: Gas_NH3_index = 0 ! ammonia
      integer,public :: Gas_He_index  = 0 ! helium
      integer,public :: Gas_Ar_index  = 0 ! argon
      integer,public :: Gas_H2_index  = 0 ! hydrogen
      integer,public :: Gas_N2_index  = 0 ! nitrogen
      integer,public :: Gas_CH4_index = 0 ! methane
      integer,public :: Gas_CO_index  = 0 ! carbon monoxide
      integer,public :: Gas_Rn_index  = 0 ! radon
      integer,public :: Gas_HBr_index = 0 ! hydrogen bromide
      integer,public :: Gas_BrO_index = 0 ! Bromine oxide
      !                     # Primary metals
      integer,public :: Gas_As_index = 0 ! arsenic
      integer,public :: Gas_Hg_index = 0 ! mercury
      integer,public :: Gas_Pb_index = 0 ! lead
      integer,public :: Gas_Al_index = 0 ! aluminium
      !                     # Secondary constituants, reaction products, aerosols
      integer,public :: Gas_SO4_index   = 0 ! particulate sulfate
      integer,public :: Gas_H2SO4_index = 0 ! sulfuric acid
      integer,public :: Gas_COS_index   = 0 ! carbonyl sulfide

      logical :: Have_H2O    = .false.
      logical :: Have_CO2    = .false.
      logical :: Have_SO2    = .false.
      logical :: Have_H2S    = .false.
      logical :: Have_HCl    = .false.
      logical :: Have_HF     = .false.
      logical :: Have_NH3    = .false.
      logical :: Have_He     = .false.
      logical :: Have_Ar     = .false.
      logical :: Have_H2     = .false.
      logical :: Have_N2     = .false.
      logical :: Have_CH4    = .false.
      logical :: Have_CO     = .false.
      logical :: Have_Rn     = .false.
      logical :: Have_HBr    = .false.
      logical :: Have_BrO    = .false.
      logical :: Have_As     = .false.
      logical :: Have_Hg     = .false.
      logical :: Have_Pb     = .false.
      logical :: Have_Al     = .false.
      logical :: Have_SO4    = .false.
      logical :: Have_H2SO4  = .false.
      logical :: Have_COS    = .false.

      ! Reaction variables
      integer       :: nreactions
      integer,dimension(:),allocatable  :: Reaction_ID
      ! Reaction 1  SO2 -> SO4 via constant decay rate
      !                        decay rate is given in s^-1
      real(kind=ip) :: Gas_SO2_SO4_decay_rate
      real(kind=ip) :: Gas_SO2_SO4_conversion_HLife !half-life in hours

      logical :: USE_GAS = .false.

      contains

!******************************************************************************

      subroutine input_data_Source_Gas

      use global_param,  only : &
         nmods,HR_2_S

      use io_data,       only : &
         infile

      use mesh,          only : &
         nxmax,nymax,nsmax,insmax,de,dn,lonLL,latLL

      use Source,        only : &
         neruptions,SourceType,e_Volume,e_Duration

      implicit none

      character(len=3)   :: answer
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
      character(len=2) :: dumstr1,dumstr2
      integer :: ireac

      ! Set up the species names by index
      !   Note: SpeciesID =>   # Primary gas emmission species as listed in Encyc. of Volc. p805
      GS_GasSpecies_name( 1)="H2O"            ;GS_GasSpecies_lname( 1)="water vapor"
        GS_GasSpecies_IsAerosol( 1)=.false.     ;GS_GasSpecies_MolWeight( 1) = 18.015_ip
      GS_GasSpecies_name( 2)="CO2"            ;GS_GasSpecies_lname( 2)="carbon dioxide"
        GS_GasSpecies_IsAerosol( 2)=.false.     ;GS_GasSpecies_MolWeight( 2) = 44.100_ip
      GS_GasSpecies_name( 3)="SO2"            ;GS_GasSpecies_lname( 3)="sulfur dioxide"
        GS_GasSpecies_IsAerosol( 3)=.false.     ;GS_GasSpecies_MolWeight( 3) = 64.066_ip
      GS_GasSpecies_name( 4)="H2S"            ;GS_GasSpecies_lname( 4)="hydrogen sulfide"
        GS_GasSpecies_IsAerosol( 4)=.false.     ;GS_GasSpecies_MolWeight( 4) = 34.100_ip
      GS_GasSpecies_name( 5)="HCl"            ;GS_GasSpecies_lname( 5)="hydrogen cloride"
        GS_GasSpecies_IsAerosol( 5)=.false.     ;GS_GasSpecies_MolWeight( 5) = 36.460_ip
      GS_GasSpecies_name( 6)="HF"             ;GS_GasSpecies_lname( 6)="hydrogen fluride"
        GS_GasSpecies_IsAerosol( 6)=.false.     ;GS_GasSpecies_MolWeight( 6) = 20.010_ip
      GS_GasSpecies_name( 7)="NH3"            ;GS_GasSpecies_lname( 7)="ammonia"
        GS_GasSpecies_IsAerosol( 7)=.false.     ;GS_GasSpecies_MolWeight( 7) = 17.031_ip
      GS_GasSpecies_name( 8)="He"             ;GS_GasSpecies_lname( 8)="helium"
        GS_GasSpecies_IsAerosol( 8)=.false.     ;GS_GasSpecies_MolWeight( 8) = 4.0026_ip
      GS_GasSpecies_name( 9)="Ar"             ;GS_GasSpecies_lname( 9)="argon"
        GS_GasSpecies_IsAerosol( 9)=.false.     ;GS_GasSpecies_MolWeight( 9) = 39.948_ip
      GS_GasSpecies_name(10)="H2"             ;GS_GasSpecies_lname(10)="hydrogen"
        GS_GasSpecies_IsAerosol(10)=.false.     ;GS_GasSpecies_MolWeight(10) = 2.016_ip
      GS_GasSpecies_name(11)="N2"             ;GS_GasSpecies_lname(11)="nitrogen"
        GS_GasSpecies_IsAerosol(11)=.false.     ;GS_GasSpecies_MolWeight(11) = 14.0067_ip
      GS_GasSpecies_name(12)="CH4"            ;GS_GasSpecies_lname(12)="methane"
        GS_GasSpecies_IsAerosol(12)=.false.     ;GS_GasSpecies_MolWeight(12) = 16.04_ip
      GS_GasSpecies_name(13)="CO"             ;GS_GasSpecies_lname(13)="carbon monoxide"
        GS_GasSpecies_IsAerosol(13)=.false.     ;GS_GasSpecies_MolWeight(13) = 28.01_ip
      GS_GasSpecies_name(14)="Rn"             ;GS_GasSpecies_lname(14)="radon"
        GS_GasSpecies_IsAerosol(14)=.false.     ;GS_GasSpecies_MolWeight(14) = 222.018_ip
      GS_GasSpecies_name(15)="HBr"            ;GS_GasSpecies_lname(15)="hydrogen bromide"
        GS_GasSpecies_IsAerosol(15)=.false.     ;GS_GasSpecies_MolWeight(15) = 80.91_ip
      GS_GasSpecies_name(16)="BrO"            ;GS_GasSpecies_lname(16)="bromine oxide"
        GS_GasSpecies_IsAerosol(16)=.false.     ;GS_GasSpecies_MolWeight(16) = 111.903_ip
      GS_GasSpecies_name(30)="As"             ;GS_GasSpecies_lname(30)="arsenic"
        GS_GasSpecies_IsAerosol(30)=.false.     ;GS_GasSpecies_MolWeight(30) = 74.9216_ip
      GS_GasSpecies_name(31)="Hg"             ;GS_GasSpecies_lname(31)="mercury"
        GS_GasSpecies_IsAerosol(31)=.false.     ;GS_GasSpecies_MolWeight(31) = 200.59_ip
      GS_GasSpecies_name(32)="Pb"             ;GS_GasSpecies_lname(32)="lead"
        GS_GasSpecies_IsAerosol(32)=.false.     ;GS_GasSpecies_MolWeight(32) = 207.2_ip
      GS_GasSpecies_name(33)="Al"             ;GS_GasSpecies_lname(33)="aluminium"
        GS_GasSpecies_IsAerosol(33)=.false.     ;GS_GasSpecies_MolWeight(33) = 26.9815_ip
      GS_GasSpecies_name(50)="SO4"            ;GS_GasSpecies_lname(50)="particulate sulfate"
        GS_GasSpecies_IsAerosol(50)=.true.      ;GS_GasSpecies_MolWeight(50) = 96.060_ip
      GS_GasSpecies_name(51)="H2SO4"          ;GS_GasSpecies_lname(51)="sulfuric acid"
        GS_GasSpecies_IsAerosol(51)=.true.      ;GS_GasSpecies_MolWeight(51) = 98.079_ip
      GS_GasSpecies_name(52)="COS"            ;GS_GasSpecies_lname(52)="carbonyl sulfide"
        GS_GasSpecies_IsAerosol(52)=.true.      ;GS_GasSpecies_MolWeight(52) = 60.07_ip

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
      allocate(SrcGasSurf_Mask(nxmax,nymax))
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
      !       EruptGasSrcStruc => 1 = distributed surface source (like a region)
      !                           2 = point on surface (vent)
      !                           3 = line on surface (fissure)
      !                           4 = vertical line source
      !                           5 = vertical profile source
      !       EruptGasMassRate      = tonnes/day
      setSrcGasSurf_Perim = .false.
      do i=1,neruptions
        !read start time, duration, plume height, volume of each pulse
        read(llinebuffer,*,err=1910) dum1_int,dum2_int,dum3_int, &
                                     dum1_ip, dum2_ip, dum3_ip,  &
                              EruptGasSpeciesID(i),EruptGasSrcStruc(i),EruptGasMassRate(i)

        if(EruptGasSrcStruc(i).eq.1)then
          ! This is a distributed source over a surface region
          ! Only one distributed region is allowed.  If this is the first dist. source, then
          ! read the file name
          if(.not.setSrcGasSurf_Perim)then ! This is initialized to .false. and only set to
                                           ! true if a distributed source has be read
            read(llinebuffer,*) dum1_int,dum2_int,dum3_int, &
                                         dum1_ip, dum2_ip, dum3_ip,  &
                              EruptGasSpeciesID(i),EruptGasSrcStruc(i),EruptGasMassRate(i),&
                              SrcGasSurf_PerimInfile
            setSrcGasSurf_Perim = .true.
          endif
          ! Read name of file outlining the contour of the region
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
          write(global_info,*)"Gas Fissure source not yet implemented"
          stop 1
        elseif(EruptGasSrcStruc(i).eq.4)then
          ! Read the lon/lat and height of line
          write(global_info,*)"Gas Line source not yet implemented"
          stop 1
        elseif(EruptGasSrcStruc(i).eq.5)then
          ! Read the profile of Mass fractions
          write(global_info,*)"Gas Profile source not yet implemented"
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
!      OPTMOD=SRC_GAS
!      2                                       # ngas_max
!      3    T T                                # Species #1, + optional output for surf (ppm) and vert col dens (DU)
!      50   T T                                #   :     #2, ....
!      1                                       # nreact: Number of reactions to calculate
!      1  3.0e-5                                  # Reaction equation ID + parameters : Convert SO2 to SO4

      nvar_User2d_XY_SrcGas = 0
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
            ! First line of SRC_GAS block is just the number of species to track
            read(10,'(a80)',iostat=ios)linebuffer
            read(linebuffer,*)ngas_max
            write(global_info,*)"ngas_max = ",ngas_max
            if(ngas_max.lt.1)then
              write(global_info,*)"No gas species defined."
              ngas_max = 0
            elseif(ngas_max.gt.23)then
              write(global_info,*)"To many gas species."
              stop 1
            else
              write(global_info,*)"Reading list of gas species IDs."
              iconcen_gas_start = insmax
            endif
            allocate(GS_GasSpeciesID(ngas_max))
            allocate(GS_GasSpecies_OutputSurfConc(ngas_max))
            allocate(GS_GasSpecies_OutputVertColDens(ngas_max))
            ! Get the list of species to trace
            do i=1,ngas_max
              read(10,'(a80)',iostat=ios)linebuffer
              read(linebuffer,*)GS_GasSpeciesID(i),dumstr1,dumstr2
              testkey = dumstr1(1:1)
              if(testkey.eq.'T'.or.testkey.eq.'t')then
                GS_GasSpecies_OutputSurfConc(i)=.true.
                nvar_User2d_XY_SrcGas = nvar_User2d_XY_SrcGas + 1
              else
                GS_GasSpecies_OutputSurfConc(i)=.false.
              endif
              testkey = dumstr2(1:1)
              if(testkey.eq.'T'.or.testkey.eq.'t')then
                GS_GasSpecies_OutputVertColDens(i)=.true.
                nvar_User2d_XY_SrcGas = nvar_User2d_XY_SrcGas + 1
              else
                GS_GasSpecies_OutputVertColDens(i)=.false.
              endif
              ! HFS do some error-checking here to verify that all the source
              ! terms are accommodated
              if(GS_GasSpeciesID(i).eq.1)then
                Have_H2O = .true.
                Gas_H2O_index = iconcen_gas_start+i    ! 1 = H2O        water vapor
              elseif(GS_GasSpeciesID(i).eq.2)then
                Have_CO2 = .true.
                Gas_CO2_index = iconcen_gas_start+i    ! 2 = CO2        carbon dioxide
              elseif(GS_GasSpeciesID(i).eq.3)then
                Have_SO2 = .true.
                Gas_SO2_index = iconcen_gas_start+i    ! 3 = SO2        sulfur dioxide
              elseif(GS_GasSpeciesID(i).eq.4)then
                Have_H2S = .true.
                Gas_H2S_index = iconcen_gas_start+i    ! 4 = H2S        hydrogen sulfide
              elseif(GS_GasSpeciesID(i).eq.5)then
                Have_HCl = .true.
                Gas_HCl_index = iconcen_gas_start+i    ! 5 = HCl        hydrogen cloride
              elseif(GS_GasSpeciesID(i).eq.6)then
                Have_HF = .true.
                Gas_HF_index  = iconcen_gas_start+i    ! 6 = HF         hydrogen fluride
              elseif(GS_GasSpeciesID(i).eq.7)then
                Have_NH3 = .true.
                Gas_NH3_index = iconcen_gas_start+i    ! 7 = NH3        ammonia
              elseif(GS_GasSpeciesID(i).eq.8)then
                Have_He = .true.
                Gas_He_index  = iconcen_gas_start+i    ! 8 = He         helium
              elseif(GS_GasSpeciesID(i).eq.9)then
                Have_Ar = .true.
                Gas_Ar_index  = iconcen_gas_start+i    ! 9 = Ar         argon
              elseif(GS_GasSpeciesID(i).eq.10)then
                Have_H2 = .true.
                Gas_H2_index  = iconcen_gas_start+i    ! 10 = H2        hydrogen
              elseif(GS_GasSpeciesID(i).eq.11)then
                Have_N2 = .true.
                Gas_N2_index  = iconcen_gas_start+i    ! 11 = N2        nitrogen
              elseif(GS_GasSpeciesID(i).eq.12)then
                Have_CH4 = .true.
                Gas_CH4_index  = iconcen_gas_start+i    ! 12 = CH4      methane
              elseif(GS_GasSpeciesID(i).eq.13)then
                Have_CO = .true.
                Gas_CO_index  = iconcen_gas_start+i   ! 13 = CO         carbon monoxide
              elseif(GS_GasSpeciesID(i).eq.14)then
                Have_Rn = .true.
                Gas_Rn_index  = iconcen_gas_start+i    ! 14 = Rn        radon
              elseif(GS_GasSpeciesID(i).eq.15)then
                Have_HBr = .true.
                Gas_HBr_index  = iconcen_gas_start+i   ! 15 = HBr       hydrogen bromide
              elseif(GS_GasSpeciesID(i).eq.16)then
                Have_BrO = .true.
                Gas_BrO_index  = iconcen_gas_start+i   ! 16 = BrO       bromine oxide
              elseif(GS_GasSpeciesID(i).eq.30)then
                Have_As = .true.
                Gas_As_index  = iconcen_gas_start+i    ! 30 = As        arsenic
              elseif(GS_GasSpeciesID(i).eq.31)then
                Have_Hg = .true.
                Gas_Hg_index  = iconcen_gas_start+i    ! 31 = Hg        mercury
              elseif(GS_GasSpeciesID(i).eq.32)then
                Have_Pb = .true.
                Gas_Pb_index  = iconcen_gas_start+i    ! 32 = Pb        lead
              elseif(GS_GasSpeciesID(i).eq.33)then
                Have_Al = .true.
                Gas_Al_index  = iconcen_gas_start+i    ! 33 = Al        aluminium
              elseif(GS_GasSpeciesID(i).eq.50)then
                Have_SO4 = .true.
                Gas_SO4_index = iconcen_gas_start+i    ! 50 = SO4       particulate sulfate
              elseif(GS_GasSpeciesID(i).eq.51)then
                Have_H2SO4 = .true.
                Gas_H2SO4_index = iconcen_gas_start+i  ! 51 = H2SO4     sulfuric acid
              elseif(GS_GasSpeciesID(i).eq.52)then
                Have_COS = .true.
                Gas_COS_index = iconcen_gas_start+i    ! 52 = COS       carbonyl sulfide
              else
                write(global_info,*)"ERROR: unknown gas code."
                write(global_info,*)linebuffer
                write(global_info,*)GS_GasSpeciesID
                stop 1
              endif
            enddo

            if(Gas_H2O_index.gt.0)  write(global_info,*)"Gas_H2O_index =",   Gas_H2O_index    ! water vapor 
            if(Gas_CO2_index.gt.0)  write(global_info,*)"Gas_CO2_index =",   Gas_CO2_index    ! carbon dioxide
            if(Gas_SO2_index.gt.0)  write(global_info,*)"Gas_SO2_index =",   Gas_SO2_index    ! sulfur dioxide
            if(Gas_H2S_index.gt.0)  write(global_info,*)"Gas_H2S_index =",   Gas_H2S_index    ! hydrogen sulfide
            if(Gas_HCl_index.gt.0)  write(global_info,*)"Gas_HCl_index =",   Gas_HCl_index    ! hydrogen cloride
            if(Gas_HF_index.gt.0)   write(global_info,*)"Gas_HF_index  =",   Gas_HF_index     ! hydrogen fluride
            if(Gas_NH3_index.gt.0)  write(global_info,*)"Gas_NH3_index =",   Gas_NH3_index    ! ammonia
            if(Gas_He_index.gt.0)   write(global_info,*)"Gas_He_index  =",   Gas_He_index     ! helium
            if(Gas_Ar_index.gt.0)   write(global_info,*)"Gas_Ar_index  =",   Gas_Ar_index     ! argon
            if(Gas_H2_index.gt.0)   write(global_info,*)"Gas_H2_index  =",   Gas_H2_index     ! hydrogen
            if(Gas_N2_index.gt.0)   write(global_info,*)"Gas_N2_index  =",   Gas_N2_index     ! nitrogen
            if(Gas_CH4_index.gt.0)  write(global_info,*)"Gas_CH4_index =",   Gas_CH4_index    ! methane
            if(Gas_CO_index.gt.0)   write(global_info,*)"Gas_CO_index  =",   Gas_CO_index     ! carbon monoxide
            if(Gas_Rn_index.gt.0)   write(global_info,*)"Gas_Rn_index  =",   Gas_Rn_index     ! radon
            if(Gas_HBr_index.gt.0)  write(global_info,*)"Gas_HBr_index =",   Gas_HBr_index    ! hydrogen bromide
            if(Gas_BrO_index.gt.0)  write(global_info,*)"Gas_BrO_index =",   Gas_BrO_index    ! Bromine oxide
            if(Gas_As_index.gt.0)   write(global_info,*)"Gas_As_index  =",   Gas_As_index     ! arsenic
            if(Gas_Hg_index.gt.0)   write(global_info,*)"Gas_Hg_index  =",   Gas_Hg_index     ! mercury
            if(Gas_Pb_index.gt.0)   write(global_info,*)"Gas_Pb_index  =",   Gas_Pb_index     ! lead
            if(Gas_Al_index.gt.0)   write(global_info,*)"Gas_Al_index  =",   Gas_Al_index     ! aluminium
            if(Gas_SO4_index.gt.0)  write(global_info,*)"Gas_SO4_index =",   Gas_SO4_index    ! particulate sulfate
            if(Gas_H2SO4_index.gt.0)write(global_info,*)"Gas_H2SO4_index =", Gas_H2SO4_index  ! sulfuric acid
            if(Gas_COS_index.gt.0)  write(global_info,*)"Gas_COS_index =",   Gas_COS_index    ! carbonyl sulfide


            ! Now read the number of reactions to calculate
            read(10,'(a80)',iostat=ios)linebuffer
            read(linebuffer,*) nreactions
            allocate(Reaction_ID(nreactions))
            if(nreactions.lt.1)then
              write(global_info,*)"No gas reactions identified."
            else
              ! Loop through the reactions, read the reaction ID and parameter list
              do ireac=1,nreactions
                read(10,'(a80)',iostat=ios)linebuffer
                read(linebuffer,*) Reaction_ID(ireac)
                if(Reaction_ID(ireac).eq.1)then
                  ! This is the SO2 -> SO4 reaction via a constant decay rate
                  ! The only parameter needed is the decay rate in seconds
                  read(linebuffer,*) Reaction_ID(ireac),Gas_SO2_SO4_decay_rate
                  Gas_SO2_SO4_conversion_HLife = 0.693147180559945_ip/Gas_SO2_SO4_decay_rate/HR_2_S
                  write(global_info,*)"SO2 to SO4 conversion reaction detected using a constant rate"
                  write(global_info,*)"      rate (s^-1) = ",real(Gas_SO2_SO4_decay_rate,kind=4)
                  write(global_info,*)" half-life (hrs)  = ",real(Gas_SO2_SO4_conversion_HLife,kind=4)
                  ! do some minimal error-checking and make sure we are tracking the needed reactants and products
                  if(Have_SO2.and.Have_SO4)then
                    write(global_info,*)"  All reactants and products are available for this reaction."
                  else
                    if(.not.Have_SO2)then
                      write(global_error,*)"ERROR: Reaction 1 (SO2->SO4) requires SO2 as a tracked species."
                    endif
                    if(.not.Have_SO4)then
                      write(global_error,*)"ERROR: Reaction 1 (SO2->SO4) requires SO4 as a tracked species."
                    endif
                    stop 1
                  endif
                elseif(Reaction_ID(ireac).eq.2)then
                  ! Some other SO2->SO4 conversion, maybe with H2O, sunlight, etc.

                else
                  write(global_info,*)"ERROR: Reaction ID not recognized",Reaction_ID(ireac)
                  stop 1
                endif
              enddo
            endif

          endif ! SRC_GAS
          !OPTMOD_names(nmods) = adjustl(trim(mod_name))
        endif
1104    format(7x,a20)
      enddo

!2010  continue
      close(10)

      insmax = insmax + ngas_max  ! Reset the counter to the max thus far
      nsmax  = insmax             ! and assign the full max value

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
         nxmax,nymax,insmax,ts1

      use solution,      only : &
         SpeciesID,SpeciesSubID

      use Source,        only : &
         SourceNodeFlux_Area

      implicit none

      integer :: gasID,ig,i,indx

      write(global_info,*)"--------------------------------------------------"
      write(global_info,*)"---------- ALLOCATE_GAS -----------------------"
      write(global_info,*)"--------------------------------------------------"

      ! Set the start indecies
      indx_User2d_XY_SrcGas        = nvar_User2d_XY

      allocate(temp_2d_name_SrcGas(nvar_User2d_XY_SrcGas))
      allocate(temp_2d_lname_SrcGas(nvar_User2d_XY_SrcGas))
      allocate(temp_2d_unit_SrcGas(nvar_User2d_XY_SrcGas))
      allocate(temp_2d_MissVal_SrcGas(nvar_User2d_XY_SrcGas))
      allocate(temp_2d_FillVal_SrcGas(nvar_User2d_XY_SrcGas))

      i = 1
      do ig=1,ngas_max
        gasID = GS_GasSpeciesID(ig)
        if(GS_GasSpecies_OutputSurfConc(ig))then
          indx = indx_User2d_XY_SrcGas+i
          temp_2d_name_SrcGas(i)    = adjustl(trim(GS_GasSpecies_name(gasID))) // "_SurfConc"
          if(GS_GasSpecies_IsAerosol(gasID))then
            temp_2d_unit_SrcGas(i)    = "ug/m3"
          else
            temp_2d_unit_SrcGas(i)    = "ppm"
          endif
          temp_2d_lname_SrcGas(i)   = adjustl(trim(GS_GasSpecies_lname(gasID))) // " concentration at surface"
          temp_2d_MissVal_SrcGas(i) = -9999.0_op
          temp_2d_FillVal_SrcGas(i) = -9999.0_op
          i=i+1
        endif
        if(GS_GasSpecies_OutputVertColDens(ig))then
          indx = indx_User2d_XY_SrcGas+i
          temp_2d_name_SrcGas(i)    = adjustl(trim(GS_GasSpecies_name(gasID))) // "_VertColDens"
          if(GS_GasSpecies_IsAerosol(gasID))then
            temp_2d_unit_SrcGas(i)    = "ug/m2"
          else
            temp_2d_unit_SrcGas(i)    = "DU"
          endif
          temp_2d_lname_SrcGas(i)   = adjustl(trim(GS_GasSpecies_lname(gasID))) // " vertical column density"
          temp_2d_MissVal_SrcGas(i) = -9999.0_op
          temp_2d_FillVal_SrcGas(i) = -9999.0_op
          i=i+1
        endif
      enddo

      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_SrcGas

      allocate(SrcGas_SurfConc(ngas_max,nxmax,nymax))
      allocate(SrcGas_VertColDens(ngas_max,nxmax,nymax))

      allocate(SourceNodeFlux_Area(1:nxmax,1:nymax,1:ngas_max))

      ! While here, assign some of the global specied indecies now that the arrays
      ! have been allocated
      SpeciesID(insmax+1:insmax+ngas_max)    = 3 ! chem bins
      SpeciesSubID(insmax+1:insmax+ngas_max) = GS_GasSpeciesID(1:ngas_max)

      end subroutine Allocate_Source_Gas

!******************************************************************************

      subroutine Prep_output_Source_Gas

      use mesh,          only : &
         nxmax,nymax,nzmax,ts1,dz_vec_pd

      use solution,      only : &
         concen_pd

      use Output_Vars,   only : &
         var_User2d_XY_name,var_User2d_XY_unit,var_User2d_XY_lname,&
         var_User2d_XY_MissVal,var_User2d_XY_FillVal,var_User2d_XY

      implicit none

      integer :: i,j,k,ig,indx
      integer :: gasID
      real(kind=op) :: dum_op

      SrcGas_SurfConc(1:ngas_max,1:nxmax,1:nymax)    = -9999.0_op
      SrcGas_VertColDens(1:ngas_max,1:nxmax,1:nymax) = -9999.0_op

      do ig=1,ngas_max
        do i=1,nxmax
          do j=1,nymax

      indx=iconcen_gas_start+ig
      gasID=GS_GasSpeciesID(ig)

      ! Get surface concentration
      dum_op = real(concen_pd(i,j,1,indx,ts1),kind=op) ! in km/km^3
      if(GS_GasSpecies_IsAerosol(gasID))then
        ! If species is an aerosol, then write concentration in kg/km^3
        if(dum_op.gt.Surf_Conc_tresh)    SrcGas_SurfConc(ig,i,j) = dum_op ! in kg/km^3
      else
        ! If a gas (not an aerosol) then convert surface concentration to ppm
        dum_op = dum_op / 2.62_op/1000.0_op  ! convert from kg/km^3 to ppm at 1 atm and 25 C
        if(dum_op.gt.Surf_Conc_tresh_ppm) SrcGas_SurfConc(ig,i,j) = dum_op ! in ppm
      endif

      dum_op = 0.0_op
      do k = 1,nzmax
        dum_op = dum_op + &
                 real(concen_pd(i,j,k,indx,ts1),kind=op) * & ! in kg/km^3
                             dz_vec_pd(k)                             ! convert to kg/km^2
      enddo
      if(GS_GasSpecies_IsAerosol(gasID))then
        ! If species is an aerosol, then write column loading in kg/km^2
        if(dum_op.gt.VertColDens_tresh)    SrcGas_VertColDens(ig,i,j) = dum_op ! in kg/km^2
      else
        ! If a gas (not an aerosol) then convert column loading to Dobson Units (DU)
        dum_op = dum_op / (0.44670_op * GS_GasSpecies_MolWeight(gasID))  ! convert from kg/km^2 to DU
        if(dum_op.gt.VertColDens_tresh_DU) SrcGas_VertColDens(ig,i,j) = dum_op ! in DU
      endif

          enddo
        enddo
      enddo

      i = 1
      do ig=1,ngas_max
        if(GS_GasSpecies_OutputSurfConc(ig))then
          indx = indx_User2d_XY_SrcGas+i
          var_User2d_XY_name(indx) = temp_2d_name_SrcGas(i)
          var_User2d_XY_unit(indx) = temp_2d_unit_SrcGas(i)
          var_User2d_XY_lname(indx)= temp_2d_lname_SrcGas(i)
          var_User2d_XY_MissVal(indx)= temp_2d_MissVal_SrcGas(i)
          var_User2d_XY_FillVal(indx)= temp_2d_FillVal_SrcGas(i)
          var_User2d_XY(:,:,indx) = SrcGas_SurfConc(ig,:,:)
          i=i+1
        endif
        if(GS_GasSpecies_OutputVertColDens(ig))then
          indx = indx_User2d_XY_SrcGas+i
          var_User2d_XY_name(indx) = temp_2d_name_SrcGas(i)
          var_User2d_XY_unit(indx) = temp_2d_unit_SrcGas(i)
          var_User2d_XY_lname(indx)= temp_2d_lname_SrcGas(i)
          var_User2d_XY_MissVal(indx)= temp_2d_MissVal_SrcGas(i)
          var_User2d_XY_FillVal(indx)= temp_2d_FillVal_SrcGas(i)
          var_User2d_XY(:,:,indx) = SrcGas_VertColDens(ig,:,:)
          i=i+1
        endif
      enddo

      !do i=1,nvar_User2d_XY_SrcGas
      !  indx = indx_User2d_XY_SrcGas+i
      !  var_User2d_XY_name(indx) = temp_2d_name_SrcGas(i)
      !  var_User2d_XY_unit(indx) = temp_2d_unit_SrcGas(i)
      !  var_User2d_XY_lname(indx)= temp_2d_lname_SrcGas(i)
      !  var_User2d_XY_MissVal(indx)= temp_2d_MissVal_SrcGas(i)
      !  var_User2d_XY_FillVal(indx)= temp_2d_FillVal_SrcGas(i)
      !  if(i.eq.1)var_User2d_XY(:,:,indx) = SrcGas_SurfConc(1,:,:)
      !  if(i.eq.2)var_User2d_XY(:,:,indx) = SrcGas_VertColDens(1,:,:)
      !  if(i.eq.3)var_User2d_XY(:,:,indx) = SrcGas_SurfConc(2,:,:)
      !  if(i.eq.4)var_User2d_XY(:,:,indx) = SrcGas_VertColDens(2,:,:)
      !enddo

      end subroutine Prep_output_Source_Gas

!******************************************************************************

      subroutine Deallocate_Source_Gas

      use Source,        only : &
         SourceNodeFlux_Area

      implicit none

      if(allocated(SrcGasSurf_Mask)) deallocate(SrcGasSurf_Mask)

#ifdef USEPOINTERS
      if(associated(SourceNodeFlux_Area))   deallocate(SourceNodeFlux_Area)
#else
      if(allocated(SourceNodeFlux_Area))   deallocate(SourceNodeFlux_Area)
#endif


      end subroutine Deallocate_Source_Gas

!******************************************************************************

!******************************************************************************

      subroutine Read_Perimeter_Gas

      use mesh,          only : &
         nxmax,nymax,ivent,jvent,x_cc_pd,y_cc_pd,lon_cc_pd,lat_cc_pd,&
         A3d_iprojflag,A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
         A3d_k0_scale,A3d_radius_earth,de,dn,dx,dy,IsLatLon

      use projection,    only : &
           PJ_proj_for

      implicit none

      integer :: npoints
      real(kind=ip),dimension(:),allocatable :: Perm_lon
      real(kind=ip),dimension(:),allocatable :: Perm_lat

      !real(kind=ip),dimension(:),allocatable :: Perm_x
      !real(kind=ip),dimension(:),allocatable :: Perm_y

      real(kind=ip) :: Perm_x_min,Perm_x_max
      real(kind=ip) :: Perm_y_min,Perm_y_max
      integer       :: Perm_i_min,Perm_i_max
      integer       :: Perm_j_min,Perm_j_max

      real(kind=ip),dimension(2) :: testpoint
      real(kind=ip),dimension(:,:),allocatable :: BCpos ! boundary corner position
      integer      ,dimension(:,:),allocatable :: Belem ! boundary element BC IDs

      real(kind=ip),dimension(nxmax) :: loc_x
      real(kind=ip),dimension(nymax) :: loc_y
      real(kind=ip) :: loc_dx, loc_dy
      real(kind=dp) :: lat_in,lon_in,xout,yout

      integer :: i,j

      write(global_info,*)"Opening ",SrcGasSurf_PerimInfile
      open(unit=20,file=SrcGasSurf_PerimInfile)
      read(20,*)npoints
      allocate(Perm_lon(npoints))
      allocate(Perm_lat(npoints))
      allocate(BCpos(2,npoints))
      allocate(Belem(2,npoints))

      ! Read in the points and assign boundary elements with the corner IDs
      do i = 1,npoints
        read(20,*)Perm_lon(i),Perm_lat(i)
        if(Perm_lon(i).lt.0.0_ip)Perm_lon(i)=Perm_lon(i)+360.0_ip
        Belem(1,i)=i
        Belem(2,i)=i+1
      enddo
        ! The end point of the last boundary element is BC #1
      Belem(2,npoints)=1

      if(IsLatLon)then
        BCpos(1,:) = Perm_lon(:)
        BCpos(2,:) = Perm_lat(:)
        loc_x(1:nxmax) = lon_cc_pd(1:nxmax)
        loc_y(1:nymax) = lat_cc_pd(1:nymax)
        loc_dx = de
        loc_dy = dn
      else
        ! We need to project each point on the deposit perimeter onto the
        ! computational grid
        do i = 1,npoints
          lon_in = real(Perm_lon(i),kind=dp)
          lat_in = real(Perm_lat(i),kind=dp)
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
      Perm_x_min = minval(BCpos(1,:))
      Perm_x_max = maxval(BCpos(1,:))
      Perm_y_min = minval(BCpos(2,:))
      Perm_y_max = maxval(BCpos(2,:))
        ! And the indicies on the computational grid braketing these points
      Perm_i_min =   floor((Perm_x_min-loc_x(1))/loc_dx) - 1
      Perm_i_max = ceiling((Perm_x_max-loc_x(1))/loc_dx) + 1
      Perm_j_min =   floor((Perm_y_min-loc_y(1))/loc_dy) - 1
      Perm_j_max = ceiling((Perm_y_max-loc_y(1))/loc_dy) + 1

      !  Loop through all computational grid points near the deposit and flag
      !  those inside as 1 in SrcGasSurf_Mask
      SrcGasSurf_Mask = 0
      SrcGasSurf_MaskCount = 0
      do i = Perm_i_min,Perm_i_max
        do j = Perm_j_min,Perm_j_max
          testpoint(1) = loc_x(i)
          testpoint(2) = loc_y(j)
          if(IsIn(testpoint,npoints,BCpos,npoints,Belem))then
            SrcGasSurf_Mask(i,j) = 1
            SrcGasSurf_MaskCount = SrcGasSurf_MaskCount + 1
          elseif(ivent.eq.i.and.jvent.eq.j)then
            SrcGasSurf_Mask(i,j) = 1
          else
            SrcGasSurf_Mask(i,j) = 0
          endif
        enddo
      enddo

      write(global_info,*)ivent,jvent, SrcGasSurf_Mask(ivent,jvent) , SrcGasSurf_MaskCount
      close(20)

      end subroutine Read_Perimeter_Gas

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
        case(1)          !                          1 = H2O        water vapor
          idx = Gas_H2O_index
        case(2)          !                          2 = CO2        carbon dioxide
          idx = Gas_CO2_index
        case(3)          !                          3 = SO2        sulfur dioxide
          idx = Gas_SO2_index
        case(4)          !                          4 = H2S        hydrogen sulfide
          idx = Gas_H2S_index
        case(5)          !                          5 = HCl        hydrogen cloride
          idx = Gas_HCl_index
        case(6)          !                          6 = HF         hydrogen fluride
          idx = Gas_HF_index 
        case(7)          !                          7 = NH3        ammonia
          idx = Gas_NH3_index
        case(8)          !                          8 = He         helium
          idx = Gas_He_index 
        case(9)          !                          9 = Ar         argon
          idx = Gas_Ar_index
        case(10)          !                         10 = H2         hydrogen
          idx = Gas_H2_index
        case(11)          !                         11 = N2         nitrogen
          idx = Gas_N2_index
        case(12)          !                         12 = CH4        methane
          idx = Gas_CH4_index
        case(13)          !                         13 = CO         carbon monoxide
          idx = Gas_CO_index
        case(14)          !                         14 = Rn         radon
          idx = Gas_Rn_index
        case(15)          !                         15 = HBr        hydrogen bromide
          idx = Gas_HBr_index
        case(16)          !                         16 = BrO        bromine oxide
          idx = Gas_BrO_index
        case(30)          !                         30 = arsenic
          idx = Gas_As_index
        case(31)          !                         31 = mercury
          idx = Gas_Hg_index
        case(32)          !                         32 = lead
          idx = Gas_Pb_index
        case(33)          !                         33 = aluminium
          idx = Gas_Al_index
        case(50)          !                         50 = SO4        particulate sulfate
          idx = Gas_SO4_index
        case(51)          !                         51 = H2SO4      sulfuric acid
          idx = Gas_H2SO4_index
        case(52)          !                         52 = COS        carbonyl sulfide
          idx = Gas_COS_index
        case default
          write(global_info,*)  'ERROR: Gas species ID not recognized'
          write(global_info,*)  'Program stopped'
          stop 1
        end select

        if(EruptGasSrcStruc(ie).eq.1)then
          ! Distributed region
          do i = 1,nxmax
            do j = 1,nymax
              if(SrcGasSurf_Mask(i,j).gt.0)then
                Fv = EruptGasMassRate(1)/real(SrcGasSurf_MaskCount,kind=ip)  ! Divide by number o cells
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
      integer :: ireac
      real(kind=ip) :: cofac
      real(kind=ip) :: scrub
      real(kind=ip) :: del_molarmass

      ! Loop through all the reactions and integrate.
      ! Note: The SO2 to SO4 conversion is slow and will not be the limiting time-step criterion.
      !       If a faster rate is used, then sub-stepping must be applied for accuracy.
      !       Furthermore, if multiple reactions are to be invoked, the sub-stepping must be in sync
      !       so that simultaneous evolution equations can be properly integrated.

      do ireac = 1,nreactions
        if(Reaction_ID(ireac).eq.1)then
          ! This is the fixed decay rate
          !  Note: half-life should be around 6 hours so we will be far from
          !  this being the most restrictive criterion on the time step.  So
          !  just use the global dt
          cofac = dt*(0.693147180559945_ip/Gas_SO2_SO4_conversion_HLife)
          ! A number of molecules will be scrubed from SO2 and that same number
          ! added to SO4, but the added mass will also include the extra O2
          del_molarmass = GS_GasSpecies_MolWeight(50)/GS_GasSpecies_MolWeight( 3)
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
                concen_pd(i,j,k,Gas_SO4_index,ts0) + scrub*del_molarmass
              enddo
            enddo
          enddo
        elseif(Reaction_ID(ireac).eq.2)then
          ! This is the more sophisticated SO2->SO4 conversion

        endif
      enddo

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

