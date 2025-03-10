      module Source_Satellite

      use precis_param

      use io_units

      ! Set the number of output variables for this module
      integer, parameter :: nvar_User2d_static_XY_SrcSat = 3 ! SatSrc_HeightMean
                                                             ! SatSrc_MeanVar
                                                             ! SatSrc_PixCell
      integer, parameter :: nvar_User2d_XY_SrcSat        = 0
      integer, parameter :: nvar_User3d_XYGs_SrcSat      = 0
      integer, parameter :: nvar_User3d_XYZ_SrcSat       = 0
      integer, parameter :: nvar_User4d_XYZGs_SrcSat     = 0

      character(len=30),dimension(nvar_User2d_static_XY_SrcSat) :: temp_2ds_name_SrcSat
      character(len=30),dimension(nvar_User2d_static_XY_SrcSat) :: temp_2ds_unit_SrcSat
      character(len=30),dimension(nvar_User2d_static_XY_SrcSat) :: temp_2ds_lname_SrcSat
      real(kind=op),    dimension(nvar_User2d_static_XY_SrcSat) :: temp_2ds_MissVal_SrcSat
      real(kind=op),    dimension(nvar_User2d_static_XY_SrcSat) :: temp_2ds_FillVal_SrcSat

      ! These are used to keep track of which index in the global list, this
      ! modules output vars corespond to
      integer :: indx_User2d_static_XY_SrcSat
      integer :: indx_User2d_XY_SrcSat
      integer :: indx_User3d_XYGs_SrcSat
      integer :: indx_User3d_XYZ_SrcSat
      integer :: indx_User4d_XYZGs_SrcSat

      character(len=30),dimension(nvar_User2d_static_XY_SrcSat) :: SourceSat_2dvarname
      character(len=30),dimension(nvar_User2d_static_XY_SrcSat) :: SourceSat_2dvarunit
      real(kind=ip)    ,dimension(nvar_User2d_static_XY_SrcSat) :: SourceSat_2dvar_FillValue
      integer          ,dimension(nvar_User2d_static_XY_SrcSat) :: SourceSat_2dvar_id

      character(len=30),dimension(nvar_User3d_XYZ_SrcSat) :: SourceSat_3dvarname
      character(len=30),dimension(nvar_User3d_XYZ_SrcSat) :: SourceSat_3dvarunit
      real(kind=ip)    ,dimension(nvar_User3d_XYZ_SrcSat) :: SourceSat_3dvar_FillValue
      integer          ,dimension(nvar_User3d_XYZ_SrcSat) :: SourceSat_3dvar_id

      real(kind=ip), dimension(:,:),allocatable :: SatHeightMean
      real(kind=ip), dimension(:,:),allocatable :: SatHeightVar
      real(kind=ip), dimension(:,:),allocatable :: SatAshMassLoad
      !real(kind=ip), dimension(:,:),allocatable :: SatPixCell
      integer, dimension(:,:),allocatable :: nSatPixels_in_CompCell

      ! Parameters needed for satellite source
      real(kind=ip) :: AshCloudThickness
      real(kind=ip) :: AshCloudPixArea

      character (len=130) :: SatInfile
      character (len=130) :: DepPerimInfile


      contains


!******************************************************************************

      subroutine input_data_Source_Satellite

      use global_param,  only : &
         nmods

      use Source,        only : &
        SourceType,e_Duration,e_Volume,MassFluxRate,e_EndTime,e_PlumeHeight

      use io_data,       only : &
        infile

      implicit none

      character(len=80)  :: linebuffer080
      character(len=120) :: linebuffer120
      character(len=130) :: linebuffer130
      character :: testkey 
      integer :: ios!,ioerr
      character(len=20) :: mod_name
      integer :: substr_pos
      integer         :: iyear  ! time data read from files
      integer         :: imonth
      integer         :: iday
      real(kind=ip)   :: hour   ! Start time of eruption in

      ! First check if the requested soure belongs to this module
      if ((SourceType.eq.'satellite').or. &
          (SourceType.eq.'Satellite').or. &
          (SourceType.eq.'SATELLITE')) then
        SourceType='satellite'
      else
        return
      endif

      ! If we've made it here, the requested source is a satellite source, open
      ! the input file again to get needed info
      open(unit=10,file=infile,status='old',err=1900)

      ! For custom sources, we want to read down to the source block
      !   Read first comment block
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
      enddo
      ! Read through block 1
      do while(testkey.ne.'#'.and.testkey.ne.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
      enddo
      ! Read next block header
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
      enddo

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Start reading satellite mass loading."
      endif;enddo

      ! read the satellite-specific source parameters
      ! Note: some of these values may be reset depending on the settings
      !       of the block of the input file OPTMOD=SRC_SAT
      read(linebuffer120,*,err=1910) iyear,imonth,iday,hour, e_Duration(1),&
                                   e_PlumeHeight(1),AshCloudThickness, &
                                   AshCloudPixArea

      read(10,'(a130)')linebuffer130
      SatInfile = adjustl(trim(linebuffer130))

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     Satellite source time = ",iyear,imonth,iday,real(hour,kind=sp)
        write(outlog(io),*)"                  Duration = ",e_Duration(1)
        write(outlog(io),*)"    Maximum height allowed = ",e_PlumeHeight(1)
        write(outlog(io),*)"       Ash cloud thickness = ",AshCloudThickness
        write(outlog(io),*)"                Pixel area = ",AshCloudPixArea
        write(outlog(io),*)"  Satellite data file name = ",SatInfile
      endif;enddo

      ! Initialize some eruption values
      e_Duration(1)  = 0.0_ip
      e_Volume(1)    = 0.0_ip
      MassFluxRate(1)  = 0.0_ip
      e_EndTime(1) = 0.0_ip

      ! Now read to the end of the input file and read the Optional Modudle
      ! block
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Searching for OPTMOD=SRC_SAT"
      endif;enddo
      nmods = 0
      read(10,'(a80)',iostat=ios)linebuffer080
      do while(ios.eq.0)
        read(10,'(a80)',iostat=ios)linebuffer080

        substr_pos = index(linebuffer080,'OPTMOD')
        if(substr_pos.eq.1)then
          ! found an optional module
          !  Parse for the keyword
          read(linebuffer080,1104)mod_name
          if(adjustl(trim(mod_name)).eq.'SRC_SAT')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo

!2010  continue
      close(10)

      return

1900  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot find input file: ',infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1
1910  do io=1,2;if(VB(io).le.verbosity_error)then
        write(outlog(io),*)  'error reading start time, duration, height or',&
                    ' volume of an eruptive pulse.  Program stopped'
      endif;enddo
      stop 1

      end subroutine input_data_Source_Satellite

!******************************************************************************

      subroutine Allocate_Source_Satellite

      use mesh,          only : &
         nxmax,nymax

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User2d_XY,nvar_User3d_XYGs,nvar_User3d_XYZ,&
         nvar_User4d_XYZGs

      implicit none

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOCATE_SrcSat ------------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      ! This only allocates a few array, but it also sets the
      ! variables for the output file

      allocate(nSatPixels_in_CompCell(nxmax,nymax))
      allocate(SatHeightMean(nxmax,nymax))
      allocate(SatHeightVar(nxmax,nymax))
      allocate(SatAshMassLoad(nxmax,nymax))

      ! Set the start indecies
      indx_User2d_static_XY_SrcSat = nvar_User2d_static_XY
      indx_User2d_XY_SrcSat        = nvar_User2d_XY
      indx_User3d_XYGs_SrcSat      = nvar_User3d_XYGs
      indx_User3d_XYZ_SrcSat       = nvar_User3d_XYZ
      indx_User4d_XYZGs_SrcSat     = nvar_User4d_XYZGs

      ! Height_Mean
      temp_2ds_name_SrcSat(1)    = "SatAshHeight_Mean"
      temp_2ds_lname_SrcSat(1)   = "Mean height of retrieved ash"
      temp_2ds_unit_SrcSat(1)    = "km"
      temp_2ds_MissVal_SrcSat(1) = -9999.0_op
      temp_2ds_FillVal_SrcSat(1) = -9999.0_op
      ! Height_Var
      temp_2ds_name_SrcSat(2)    = "SatAshHeight_Var"
      temp_2ds_lname_SrcSat(2)   = "Var of height of retrived ash"
      temp_2ds_unit_SrcSat(2)    = "km^2"
      temp_2ds_MissVal_SrcSat(2) = -9999.0_op
      temp_2ds_FillVal_SrcSat(2) = -9999.0_op
      ! Cell_Count
      temp_2ds_name_SrcSat(3)    = "SatCellCount"
      temp_2ds_lname_SrcSat(3)   = "Sat pixels in comp cell"
      temp_2ds_unit_SrcSat(3)    = "number_of_pixels"
      temp_2ds_MissVal_SrcSat(3) = -9999.0_op
      temp_2ds_FillVal_SrcSat(3) = -9999.0_op

      nvar_User2d_static_XY = nvar_User2d_static_XY + nvar_User2d_static_XY_SrcSat
      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_SrcSat
      nvar_User3d_XYGs      = nvar_User3d_XYGs      + nvar_User3d_XYGs_SrcSat
      nvar_User3d_XYZ       = nvar_User3d_XYZ       + nvar_User3d_XYZ_SrcSat
      nvar_User4d_XYZGs     = nvar_User4d_XYZGs     + nvar_User4d_XYZGs_SrcSat

      end subroutine Allocate_Source_Satellite

!******************************************************************************

      subroutine Prep_output_SrcSat

      use Output_Vars,   only : &
         var_User2d_static_XY_name,var_User2d_static_XY_unit,var_User2d_static_XY_lname,&
         var_User2d_static_XY_MissVal,var_User2d_static_XY_FillVal,var_User2d_static_XY

      implicit none

      integer :: i,indx

      do i=1,nvar_User2d_static_XY_SrcSat
        indx = indx_User2d_static_XY_SrcSat+i
        var_User2d_static_XY_name(indx) = temp_2ds_name_SrcSat(i)
        var_User2d_static_XY_unit(indx) = temp_2ds_unit_SrcSat(i)
        var_User2d_static_XY_lname(indx)= temp_2ds_lname_SrcSat(i)
        var_User2d_static_XY_MissVal(indx)= temp_2ds_MissVal_SrcSat(i)
        var_User2d_static_XY_FillVal(indx)= temp_2ds_FillVal_SrcSat(i)
        if(i.eq.1)var_User2d_static_XY(:,:,indx) = real(SatHeightMean,kind=op)
        if(i.eq.2)var_User2d_static_XY(:,:,indx) = real(SatHeightVar,kind=op)
        if(i.eq.3)var_User2d_static_XY(:,:,indx) = real(nSatPixels_in_CompCell,kind=op)
      enddo

      end subroutine Prep_output_SrcSat

!******************************************************************************

      subroutine Deallocate_Source_Satellite

      implicit none

      if(allocated(nSatPixels_in_CompCell)) deallocate(nSatPixels_in_CompCell)
      if(allocated(SatHeightMean))          deallocate(SatHeightMean)
      if(allocated(SatHeightVar))           deallocate(SatHeightVar)
      if(allocated(SatAshMassLoad))         deallocate(SatAshMassLoad)

      end subroutine Deallocate_Source_Satellite

!******************************************************************************




!******************************************************************************

      subroutine Read_SatMassLoading

      use global_param,    only : &
         EPS_SMALL

      use solution,      only : &
         concen_pd

      use mesh,          only : &
         nxmax,nymax,nzmax,lon_cc_pd,lat_cc_pd,de,dn,ts1,kappa_pd,z_cc_pd,dz_vec_pd

      use netcdf

      implicit none

      integer :: nSTAT
      integer :: ncid

      character(len = nf90_max_name) :: name_dum

      integer :: x_dim_id          = 0 ! x or lon
      integer :: y_dim_id          = 0 ! y or lat
      integer :: s_dim_id          = 0 ! dummy dim (or single)

      integer :: sat_nx
      integer :: sat_ny
      integer :: sat_ns

        ! 2d variables
      integer :: lat_var_id          = 0
      integer :: lon_var_id          = 0
      integer :: AshProb_var_id      = 0
      integer :: AshHeight_var_id    = 0
      integer :: AshMassLoad_var_id  = 0
      integer :: AshEffecRad_var_id  = 0

        ! 1-d variables
      integer :: TotalArea_var_id    = 0
      integer :: TotalMass_var_id    = 0

      real(kind=sp) :: LAT_MissValue
      real(kind=sp) :: LON_MissValue
      real(kind=sp) :: AshProb_MissValue
      real(kind=sp) :: AshHeight_MissValue
      real(kind=sp) :: AshMassLoad_MissValue
      real(kind=sp) :: AshEffecRad_MissValue

        ! Variables that live on the satellite grid
      real(kind=sp), dimension(:,:) ,allocatable :: dum2d_in
      real(kind=sp), dimension(:,:) ,allocatable :: lon_SatIn
      real(kind=sp), dimension(:,:) ,allocatable :: lat_SatIn
      real(kind=sp), dimension(:,:) ,allocatable :: AshProb_SatIn
      real(kind=sp), dimension(:,:) ,allocatable :: AshHeight_SatIn
      real(kind=sp), dimension(:,:) ,allocatable :: AshMassLoad_SatIn
      real(kind=sp), dimension(:,:) ,allocatable :: AshEffecRad_SatIn

        ! Variables that live on the computational grid
      !integer, dimension(:,:),allocatable :: nSatPixels_in_CompCell
      real(kind=sp), dimension(:,:,:) ,allocatable :: AshHeights_in_CompCell
      real(kind=sp), dimension(:,:,:) ,allocatable :: AshMass_in_CompCell

      integer, dimension(:,:) ,allocatable :: x_index
      integer, dimension(:,:) ,allocatable :: y_index
      integer, dimension(:,:) ,allocatable :: z_index

      real(kind=sp) :: TotalArea
      real(kind=sp) :: TotalMass

      integer :: SatData_Type = 1 ! 1 = M.Pav data for Kasatochi

      integer :: Max_Cell_Count
      integer :: dum_int
      real(kind=sp) :: dum_sp
      integer :: i,j,k
      integer :: ic
      integer :: ix,iy,iz,ig,iparts,npix
      real(kind=ip) :: pixlon,pixlat,pixmass,frac
      real(kind=ip) :: mean_height,var_height,denom,numer

      real(kind=ip) :: dz

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"About to open sat file : ",SatInfile
      endif;enddo
      nSTAT=nf90_open(SatInfile,NF90_NOWRITE, ncid)

      if(SatData_Type.eq.1)then
      ! Get dim ids and sizes
      nSTAT = nf90_inq_dimid(ncid,"nx",x_dim_id)
      nSTAT = nf90_Inquire_Dimension(ncid,x_dim_id,name=name_dum,len=sat_nx)

      nSTAT = nf90_inq_dimid(ncid,"ny",y_dim_id)
      nSTAT = nf90_Inquire_Dimension(ncid,y_dim_id,name=name_dum,len=sat_ny)

      nSTAT = nf90_inq_dimid(ncid,"single",s_dim_id)
      nSTAT = nf90_Inquire_Dimension(ncid,s_dim_id,name=name_dum,len=sat_ns)

      ! Get variable ids
      nSTAT = nf90_inq_varid(ncid,"LAT",lat_var_id)
      !nSTAT = nf90_Inquire_Attribute(ncid, lat_var_id,"missing_value",    &
      !                       xtype, length, attnum)
      nSTAT = nf90_get_att(ncid,lat_var_id,"missing_value",LAT_MissValue)


      nSTAT = nf90_inq_varid(ncid,"LON",lon_var_id)
      nSTAT = nf90_get_att(ncid,lon_var_id,"missing_value",LON_MissValue)

      nSTAT = nf90_inq_varid(ncid,"ASH_PROBABILITY",AshProb_var_id)
      nSTAT = nf90_get_att(ncid,AshProb_var_id,"missing_value",AshProb_MissValue)

      nSTAT = nf90_inq_varid(ncid,"ASH_HEIGHT",AshHeight_var_id)
      nSTAT = nf90_get_att(ncid,AshHeight_var_id,"missing_value",AshHeight_MissValue)

      nSTAT = nf90_inq_varid(ncid,"ASH_MASS_LOADING",AshMassLoad_var_id)
      nSTAT = nf90_get_att(ncid,AshMassLoad_var_id,"missing_value",AshMassLoad_MissValue)
      nSTAT = nf90_inq_varid(ncid,"ASH_EFFECTIVE_RADIUS",AshEffecRad_var_id)
      nSTAT = nf90_get_att(ncid,AshEffecRad_var_id,"missing_value",AshEffecRad_MissValue)

      nSTAT = nf90_inq_varid(ncid,"TOTAL_AREA",TotalArea_var_id)
      nSTAT = nf90_inq_varid(ncid,"TOTAL_MASS",TotalMass_var_id)

      !write(outlog(io),*)x_dim_id,y_dim_id,s_dim_id
      !write(outlog(io),*)lat_var_id,lon_var_id
      !write(outlog(io),*)AshProb_var_id,AshHeight_var_id,AshMassLoad_var_id,AshEffecRad_var_id
      !write(outlog(io),*)TotalArea_var_id,TotalMass_var_id

      !write(outlog(io),*)sat_nx,sat_ny,sat_ns
      !write(outlog(io),*)LAT_MissValue
      !write(outlog(io),*)LON_MissValue
      !write(outlog(io),*)AshProb_MissValue
      !write(outlog(io),*)AshHeight_MissValue
      !write(outlog(io),*)AshMassLoad_MissValue
      !write(outlog(io),*)AshEffecRad_MissValue

      allocate(dum2d_in(sat_nx,sat_ny))
      allocate(lon_SatIn(sat_nx,sat_ny))
      allocate(lat_SatIn(sat_nx,sat_ny))
      allocate(AshProb_SatIn(sat_nx,sat_ny))
      allocate(AshHeight_SatIn(sat_nx,sat_ny))
      allocate(AshMassLoad_SatIn(sat_nx,sat_ny))
      allocate(AshEffecRad_SatIn(sat_nx,sat_ny))

      allocate(x_index(sat_nx,sat_ny))
      allocate(y_index(sat_nx,sat_ny))
      allocate(z_index(sat_nx,sat_ny))

      !allocate(nSatPixels_in_CompCell(nxmax,nymax))
      !allocate(SatHeightMean(nxmax,nymax))
      !allocate(SatHeightVar(nxmax,nymax))

      nSatPixels_in_CompCell = 0

      ! Get the 1d variables
      nSTAT = nf90_get_var(ncid,TotalArea_var_id,TotalArea)
      nSTAT = nf90_get_var(ncid,TotalMass_var_id,TotalMass)

      ! Get the 2d variables
      nSTAT = nf90_get_var(ncid,lon_var_id,dum2d_in, &
             start = (/1,1/),count = (/sat_nx,sat_ny/))
         lon_SatIn(:,:) = dum2d_in(:,:)
      nSTAT = nf90_get_var(ncid,lat_var_id,dum2d_in, &
             start = (/1,1/),count = (/sat_nx,sat_ny/))
         lat_SatIn(:,:) = dum2d_in(:,:)
      nSTAT = nf90_get_var(ncid,AshProb_var_id,dum2d_in, &
             start = (/1,1/),count = (/sat_nx,sat_ny/))
         AshProb_SatIn(:,:) = dum2d_in(:,:)
      nSTAT = nf90_get_var(ncid,AshHeight_var_id,dum2d_in, &
             start = (/1,1/),count = (/sat_nx,sat_ny/))
         AshHeight_SatIn(:,:) = dum2d_in(:,:)
      nSTAT = nf90_get_var(ncid,AshMassLoad_var_id,dum2d_in, &
             start = (/1,1/),count = (/sat_nx,sat_ny/))
         AshMassLoad_SatIn(:,:) = dum2d_in(:,:)
      nSTAT = nf90_get_var(ncid,AshEffecRad_var_id,dum2d_in, &
             start = (/1,1/),count = (/sat_nx,sat_ny/))
         AshEffecRad_SatIn(:,:) = dum2d_in(:,:)

      ! Double-check the total mass and area
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Total Area Reported = ",TotalArea," km^2"
        write(outlog(io),*)"Total Mass Reported = ",TotalMass," Tg"
      endif;enddo
      dum_sp = 0.0_sp
      dum_int = 0
      do i=1,sat_nx
        do j=1,sat_ny
          if(abs(AshMassLoad_SatIn(i,j)-AshMassLoad_MissValue).gt.EPS_SMALL)then
            dum_sp = dum_sp + (AshMassLoad_SatIn(i,j)              & ! in g/m2
                               * 1.0e6_sp                          & ! in g/km2
                               * real(AshCloudPixArea,kind=sp)      & ! in g
                               * 1.0e-12_sp)                         ! in Tg
            dum_int = dum_int + 1
          endif
        enddo
      enddo
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Num of non-zero values  = ",dum_int
        write(outlog(io),*)"Mass sum from loading   = ",dum_sp," Tg"
        write(outlog(io),*)"Mass ratio (given/calc) = ",dum_sp/TotalMass
        write(outlog(io),*)"Area ratio (given/calc) = ",(dum_int*AshCloudPixArea)/TotalArea
      endif;enddo

      ! Close file
      nSTAT = nf90_close(ncid)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Successfully read satellite data."
      endif;enddo

      ! Now generate the mapping of each satellite pixel to the computational
      ! grid
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"KLUDGE"
      endif;enddo
      dz = dz_vec_pd(1)
      iparts = floor(AshCloudThickness/dz)
      frac   = 1.0/real(iparts)
      do i=1,sat_nx
        do j=1,sat_ny
          pixlon = lon_SatIn(i,j)
          if(pixlon.lt.0.0_ip) pixlon=pixlon+360.0_ip
          ix = floor((pixlon-lon_cc_pd(1)-0.5_ip*de)/de) + 1
          x_index(i,j) = ix

          pixlat = lat_SatIn(i,j)
          iy = floor((pixlat-lat_cc_pd(1)-0.5_ip*dn)/dn) + 1
          y_index(i,j) = iy

          if(AshProb_SatIn(i,j).gt.0.0_ip)then
            iz = floor((AshHeight_SatIn(i,j))/dz)
            iz = max(iz-iparts,1)
            z_index(i,j) = iz
            nSatPixels_in_CompCell(ix,iy) = nSatPixels_in_CompCell(ix,iy) + 1
          endif
        enddo
      enddo
      Max_Cell_Count = maxval(nSatPixels_in_CompCell)
      allocate(AshHeights_in_CompCell(nxmax,nymax,Max_Cell_Count))
      AshHeights_in_CompCell = 0.0_ip
      allocate(AshMass_in_CompCell(nxmax,nymax,Max_Cell_Count))
      AshMass_in_CompCell    = 0.0_ip

      ! And finally, add the mass to the concen array
      !  This branch is for adding a linear 1km thick profile based on
      !  pixel height
      !if(1.eq.0)then
      !concen_pd(:,:,:,:,:) = 0.0_ip
      !do i=1,sat_nx
      !  do j=1,sat_ny
      !    ix = x_index(i,j)
      !    iy = y_index(i,j)
      !    iz = z_index(i,j)
      !    ig = 1            ! Currently, just put everything in the smallest gs bin
      !    if(AshMassLoad_SatIn(i,j).ne.AshMassLoad_MissValue)then
      !      if(ix.ge.1.and.ix.le.nxmax.and.&
      !         iy.ge.1.and.iy.le.nymax.and.&
      !         iz.ge.1.and.iz.le.nzmax-(iparts-1))then
      !        pixmass = AshMassLoad_SatIn(i,j)                & ! in g/m2
      !                    * 1.0e6_sp                          & ! in g/km2
      !                    * real(AshCloudPixArea,kind=sp)      & ! in g
      !                    * 1.0e-3_sp                           ! in kg
      !        concen_pd(ix,iy,iz:iz+iparts-1,ig,ts1) = &
      !             concen_pd(ix,iy,iz:iz+iparts-1,ig,ts1) + &
      !             pixmass                           / & ! in kg
      !             kappa_pd(ix,iy,iz)                   * & ! in kg/km3
      !             frac                                  ! partitioned
      !      else
      !        do io=1,2;if(VB(io).le.verbosity_info)then
      !          write(outlog(io),*)"Satellite load not in domain: ",i,j,AshMassLoad_SatIn(i,j)
      !          write(outlog(io),*)"                              ",lon_SatIn(i,j),lat_SatIn(i,j)
      !        endif;enddo
      !      endif
      !    endif
      !  enddo
      !enddo
      !else
        ! This branch places a Gaussian distribution based on a weighted
        ! mean height and the variance
      concen_pd(:,:,:,:,:) = 0.0_ip

      nSatPixels_in_CompCell = 0
      do i=1,sat_nx
        do j=1,sat_ny
          if(AshProb_SatIn(i,j).gt.0.0)then
            ix = x_index(i,j)
            iy = y_index(i,j)
            ic = nSatPixels_in_CompCell(ix,iy)+1
            nSatPixels_in_CompCell(ix,iy) = ic
            AshHeights_in_CompCell(ix,iy,ic) = AshHeight_SatIn(i,j)
            AshMass_in_CompCell(ix,iy,ic) = AshMassLoad_SatIn(i,j)                & ! in g/m2
                                              * 1.0e6_sp                          & ! in g/km2
                                              * real(AshCloudPixArea,kind=sp)      & ! in g
                                              * 1.0e-3_sp                           ! in kg
          endif
        enddo
      enddo
        ! Now calculate the variance
      do i=1,nxmax
        do j=1,nymax
          npix = nSatPixels_in_CompCell(i,j)
          if(npix.gt.0)then
            ! First the mean
            numer = 0.0_ip
            denom = 0.0_ip
            do ic = 1,nSatPixels_in_CompCell(i,j)
              numer = numer + AshMass_in_CompCell(i,j,ic) * &
                              AshHeights_in_CompCell(i,j,ic)
              denom = denom + AshMass_in_CompCell(i,j,ic)
            enddo
            mean_height = numer/denom
            SatHeightMean(i,j) = mean_height

            ! Now the variance
            numer = 0.0_ip
            denom = 0.0_ip
            do ic = 1,nSatPixels_in_CompCell(i,j)
              numer = numer + (AshHeights_in_CompCell(i,j,ic) - &
                               mean_height)**2.0_ip
            enddo
            if(npix.eq.1)then
              var_height = (numer/(npix))**0.5_ip
            else
              var_height = (numer/(npix-1))**0.5_ip
            endif
            SatHeightVar(i,j) = var_height

            ig = 1            ! Currently, just put everything in the smallest gs bin
            pixmass = sum(AshMass_in_CompCell(i,j,:))
            do k=1,nzmax
              frac = Gaussian_Frac(z_cc_pd(k)-dz*0.5_ip, &
                                   z_cc_pd(k)+dz*0.5_ip, &
                                   mean_height,var_height)
              concen_pd(i,j,k,ig,ts1) = &
                   concen_pd(i,j,k,ig,ts1) + &
                   pixmass                           / & ! in kg
                   kappa_pd(i,j,k)                   * & ! in kg/km3
                   frac                                  ! partitioned

            enddo
          endif
        enddo
      enddo
      !endif

      endif

      SatAshMassLoad(1:nxmax,1:nymax) = 0.0
      do k = 1,nzmax
        SatAshMassLoad(1:nxmax,1:nymax) = SatAshMassLoad(1:nxmax,1:nymax) &
                                      + concen_pd(1:nxmax,1:nymax,k,ig,ts1)
      enddo
      concen_pd(:,:,:,:,:) = 0.0_ip

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)&
         "Successfully loaded satellite data onto computational grid."
      endif;enddo

      end subroutine Read_SatMassLoading

!******************************************************************************

!******************************************************************************
      function Gaussian_Frac(x1,x2,Gauss_Mean, Gauss_Var)

      implicit none

      real(kind=ip) :: Gaussian_Frac

      real(kind=ip) :: x1
      real(kind=ip) :: x2
      real(kind=ip) :: Gauss_Mean
      real(kind=ip) :: Gauss_Var

      real(kind=ip) :: arg_x1,erf_at_x1
      real(kind=ip) :: arg_x2,erf_at_x2
      

      arg_x1 = (x1-Gauss_Mean)/sqrt(2.0_ip*Gauss_Var)
      erf_at_x1 = 0.5_ip*(1.0_ip + erf(arg_x1))
      arg_x2 = (x2-Gauss_Mean)/sqrt(2.0_ip*Gauss_Var)
      erf_at_x2 = 0.5_ip*(1.0_ip + erf(arg_x2))

      Gaussian_Frac = erf_at_x2-erf_at_x1

      return

      end function Gaussian_Frac

!******************************************************************************

      subroutine Sat_Calc_Error(PNORM,Evol_scale,model_comp,data_comp,locLnorm)

      use mesh,          only : &
         nxmax,nymax

      implicit none

      real(kind=ip) :: PNORM
      real(kind=ip) :: Evol_scale
      real(kind=ip), dimension(nxmax,nymax) :: model_comp
      real(kind=ip), dimension(nxmax,nymax) :: data_comp
      real(kind=ip) :: locLnorm

      !real(kind=ip), dimension(nxmax,nymax) :: abs_error
      integer :: ix,iy
      real(kind=ip) :: locerr

      !abs_error = 0.0
      locLnorm = 0.0
      DO ix = 1,nxmax
        DO iy = 1,nymax
          locerr = abs(Evol_scale*model_comp(ix,iy)-data_comp(ix,iy))
          !abs_error(ix,iy) = locerr
          locLnorm = locLnorm + locerr**PNORM
        ENDDO
      ENDDO
      !locLnorm = (locLnorm)**(1.0/PNORM)

      end subroutine Sat_Calc_Error

      end module Source_Satellite

