!*******************************************************************************
!# Topography
!#  Topofile format must be one of
!#    1 : Gridded lon/lat (netcdf)
!#        ETOPO : https://www.ncei.noaa.gov/products/etopo-global-relief-model
!#          ETOPO 1 (deprecated)
!#           lon double
!#           lat float
!#           z(lat,lon) short (m)
!#            60 arcsec in one file (446 Mb)
!#          ETOPO 2022 (netcdf whole or subset)
!#           lon double
!#           lat double
!#           z(lat,lon) float (m)
!#            15 arcsec in 15deg tiles (~20 Mb each)
!#            30 arcsec in 1 file (1.6 Gb)
!#            60 arcsec in 1 file (457 Mb)
!#        GEBCO : https://www.gebco.net/
!#          GEBCO 08 (deprecated)
!#           x_range double
!#           y_range double
!#           z(xysize) short (m)
!#            30 arcsec in one file (1.8 Gb)
!#          GEBCO 14 (deprecated)
!#           lon double
!#           lat double
!#           elevation(lat,lon) short (m)
!#            30 arcsec in one file (1.8 Gb)
!#          GEBCO 2023 (netcdf whole or subset)
!#           lon double
!#           lat double
!#           elevation(lat,lon) short (m)
!#            15 arcsec in 1 file (7.0 Gb for global)
!#
!#    2 : Gridded Binary
!#        NOAA Globe (1-km/30 arcsec) https://www.ngdc.noaa.gov/mgg/topo/globe.html
!#          16 tiles of binary elevation data
!#           10800 col by 6000 row (-50->50lat) or 4800 or polar
!#           signed 16 bit integer binary
!#
!#        GTOPO30 (1-km/30 arcsec)
!#          https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-30-arc-second-elevation-gtopo30
!#          33 tiles of binary elevation data
!#           DEM as signed 16 bit integer binary
!#           HDR with the layout
!#      
!#        GMTED2010 (30,15,7.5 arcsec) https://www.usgs.gov/coastal-changes-and-impacts/gmted2010
!#
!#        US National Elevation Database (NED) 1/3 arcsec (10m)
!#          https://www.usgs.gov/programs/national-geospatial-program/national-map
!#          1x1 deg tiles of binary elevation data
!#          flt as 4 byte float binary
!#
!#        Shuttle Radar Topography Mission (SRTM) 1 arcsec global
!#          https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-shuttle-radar-topography-mission-srtm-1
!#           hgt as signed 16 bit integer binary
!#
!#    3 : ESRI ASCII 
!#
!*******************************************************************************
!OPTMOD=TOPO
!yes 2                         # use topography?; z-mod (0=none,1=shift,2=sigma)
!2 1.0                         # Topofile format, smoothing radius
!GEBCO_08.nc                   # topofile name

      module Topography

!##############################################################################
!     
!     In principle, netcdf would not be needed for including
!     topography, however, we assume/require that this library is
!     included to read ETOPO and GEBCO
!
!##############################################################################

      use precis_param

      use io_units

      ! Set the number of output variables for this module
      integer, parameter :: nvar_User2d_static_XY_Topo = 1 ! topography
      integer, parameter :: nvar_User2d_XY_Topo        = 0
      integer, parameter :: nvar_User3d_XYGs_Topo      = 0
      integer, parameter :: nvar_User3d_XYZ_Topo       = 0
      integer, parameter :: nvar_User4d_XYZGs_Topo     = 0

      character(len=30),dimension(nvar_User2d_static_XY_Topo) :: temp_2ds_name_Topo
      character(len=30),dimension(nvar_User2d_static_XY_Topo) :: temp_2ds_unit_Topo
      character(len=30),dimension(nvar_User2d_static_XY_Topo) :: temp_2ds_lname_Topo
      real(kind=op),    dimension(nvar_User2d_static_XY_Topo) :: temp_2ds_MissVal_Topo
      real(kind=op),    dimension(nvar_User2d_static_XY_Topo) :: temp_2ds_FillVal_Topo

      ! These are used to keep track of which index in the global list, this
      ! modules output vars correspond to
      integer :: indx_User2d_static_XY_Topo
      integer :: indx_User2d_XY_Topo
      integer :: indx_User3d_XYGs_Topo
      integer :: indx_User3d_XYZ_Topo
      integer :: indx_User4d_XYZGs_Topo

      character (len=50)  :: file_topo
      integer, parameter  :: fid_hdrfile     = 600
      integer, parameter  :: fid_datfile     = 601

      integer             :: topoFormat
      real(kind=ip)       :: rad_smooth
      logical :: useTopo         = .false.
      logical :: useSmoothTopo   = .false.

      integer :: nlat_topo_fullgrid
      integer :: nlon_topo_fullgrid
      real(kind=dp), dimension(:)    ,allocatable :: lat_topo_fullgrid
      real(kind=dp), dimension(:)    ,allocatable :: lon_topo_fullgrid

      integer :: nlat_topo_subgrid
      integer :: nlon_topo_subgrid
      real(kind=dp), dimension(:)    ,allocatable :: lat_topo_subgrid
      real(kind=dp), dimension(:)    ,allocatable :: lon_topo_subgrid
      real(kind=dp) :: dlat_topo
      real(kind=dp) :: dlon_topo
      real(kind=dp) :: cleft,cright
      real(kind=sp), dimension(:,:)  ,allocatable :: topo_subgrid

      ! These are on the computational grid
      real(kind=sp),dimension(:,:)  ,allocatable :: topo_comp ! Used if useTopo=.true.
      integer      ,dimension(:,:)  ,allocatable :: topo_indx ! kindex of topo
      integer,private :: lon_shift_flag

      real(kind=ip) :: minlon_Topo,maxlon_Topo
      real(kind=ip) :: minlat_Topo,maxlat_Topo

      contains

!******************************************************************************

      subroutine input_data_Topo

      use global_param,  only : &
         DEG2RAD,RAD_EARTH,nmods

      use mesh,          only : &
         dx,dy,de,dn,IsLatLon,ZScaling_ID

      use io_data,       only : &
         infile

      use MetReader,       only : &
         MR_useTopo,MR_ZScaling_ID

      implicit none

      character(len=3)  :: answer
      integer           :: dum_int
      character(len=80) :: linebuffer080
      integer           :: ios,ioerr
      character(len=20) :: mod_name
      integer           :: substr_pos
      logical           :: IsThere

      open(unit=10,file=infile,status='old',err=1900)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Searching for OPTMOD=TOPO"
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
          if(adjustl(trim(mod_name)).eq.'TOPO')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo

      useTopo = .false.
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Continue reading input file for topo block"
      endif;enddo

      ! Check if we're going to use topography
      read(10,'(a80)',iostat=ios,err=2010)linebuffer080
      read(linebuffer080,'(a3)') answer
      if (answer.eq.'yes') then
        useTopo = .true.
        read(linebuffer080,*,iostat=ios) answer,dum_int
        if(ios.eq.0)then
          ZScaling_ID=dum_int
        else
          ZScaling_ID=0
        endif
        MR_ZScaling_ID = ZScaling_ID
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Using topography"
          if(ZScaling_ID.eq.0)then
            write(outlog(io),*)"     with altitude coordinates :: s=z"
          elseif(ZScaling_ID.eq.1)then
            write(outlog(io),*)"     with shifted-altitude coordinates :: s=z-Zsurf"
          elseif(ZScaling_ID.eq.2)then
            write(outlog(io),*)"     with sigma-altitude coordinates :: s=(z-Zsurf)/(Ztop-Zsurf)"
          endif
        endif;enddo
      elseif(answer(1:2).eq.'no') then
        useTopo = .false.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Not using topography"
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) 'Error reading whether to use topography.'
          write(errlog(io),*) 'Answer must be yes or no.'
          write(errlog(io),*) 'You gave:',linebuffer080
          write(errlog(io),*) 'Program stopped'
        endif;enddo
        stop 1
      endif

      if (useTopo) then
        ! Check if we're using topography, then get the format code
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        read(linebuffer080,*,iostat=ioerr) topoFormat,rad_smooth
        if(topoFormat.eq.1)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Read topoFormat = 1 (NetCDF gridded)"
            write(outlog(io),*)"    Read smoothing radius = ",&
                               real(rad_smooth,kind=4)
          endif;enddo
#ifndef USENETCDF
        ! If we are here, then we expect to read the netcdf output file.  If netcdf
        ! not linked, give an error and exit
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'Expecting to read a netcdf topography file, but the netcdf'
          write(errlog(io),*)'library is not linked.  Please recompile, linking to'
          write(errlog(io),*)'netcdf or use gridded binary or ASCII topography data.'
        endif;enddo
        stop 1
#endif
        elseif(topoFormat.eq.2)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Read topoFormat = 2 (Gridded binary)"
            write(outlog(io),*)"    Read smoothing radius = ",&
                               real(rad_smooth,kind=4)
          endif;enddo
        elseif(topoFormat.eq.3)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Read topoFormat = 3 (ESRI ASCII)"
            write(outlog(io),*)"    Read smoothing radius = ",&
                               real(rad_smooth,kind=4)
          endif;enddo
        endif
        ! And read the file name
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        read(linebuffer080,*) file_topo
        file_topo = trim(adjustl(file_topo))
        do io=1,2;if(VB(io).le.verbosity_info)then           
          write(outlog(io),*)"    Read file_topo = ",file_topo
        endif;enddo

        inquire( file=file_topo, exist=IsThere )
        if(.not.IsThere)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Cannot find topography file"
            write(errlog(io),*)"     ",file_topo
          endif;enddo
          stop 1
        endif

        ! Turn on flag in MetReader library for topography
        MR_useTopo = .true.

      endif

2010  continue
      close(10)

      return

1900  do io=1,2;if(VB(io).le.verbosity_error)then             
        write(errlog(io),*)  'Error: cannot find input file: ',infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

      end subroutine input_data_Topo

!******************************************************************************

      subroutine Allocate_Topo(nx,ny)

      use io_data,       only : &
         nvar_User3d_XYZ,nvar_User3d_XYGs,nvar_User2d_XY,nvar_User2d_static_XY,&
         nvar_User4d_XYZGs

      implicit none

      integer :: nx,ny

      do io=1,2;if(VB(io).le.verbosity_info)then             
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOCATE_TOPO -------------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

!      allocate(topo_comp(0:nx+1,0:ny+1));       topo_comp = 0.0_sp
!      allocate(topo_indx(0:nx+1,0:ny+1));       topo_indx = 0
      allocate(topo_comp(nx,ny));       topo_comp = 0.0_sp
      allocate(topo_indx(nx,ny));       topo_indx = 0


      ! Set the start indecies
      indx_User2d_static_XY_Topo = nvar_User2d_static_XY
      indx_User2d_XY_Topo        = nvar_User2d_XY
      indx_User3d_XYGs_Topo      = nvar_User3d_XYGs
      indx_User3d_XYZ_Topo       = nvar_User3d_XYZ
      indx_User4d_XYZGs_Topo     = nvar_User4d_XYZGs

      temp_2ds_name_Topo(1)  = "Topography"
      temp_2ds_lname_Topo(1) = "Elevation of surface"
      temp_2ds_unit_Topo(1)  = "km"
      temp_2ds_MissVal_Topo(1) = -9999.0_op
      temp_2ds_FillVal_Topo(1) = -9999.0_op

      nvar_User2d_static_XY = nvar_User2d_static_XY + nvar_User2d_static_XY_Topo
      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_Topo
      nvar_User3d_XYGs      = nvar_User3d_XYGs      + nvar_User3d_XYGs_Topo
      nvar_User3d_XYZ       = nvar_User3d_XYZ       + nvar_User3d_XYZ_Topo
      nvar_User4d_XYZGs     = nvar_User4d_XYZGs     + nvar_User4d_XYZGs_Topo

      end subroutine Allocate_Topo

!******************************************************************************

      subroutine Prep_output_Topo

      use mesh,          only : &
         nxmax,nymax

      use Output_Vars,   only : &
         var_User2d_static_XY_name,var_User2d_static_XY_unit,var_User2d_static_XY_lname,&
         var_User2d_static_XY_MissVal,var_User2d_static_XY_FillVal,&
         var_User2d_static_XY

      implicit none

      integer :: i,indx

      do i=1,nvar_User2d_static_XY_Topo
        indx = indx_User2d_static_XY_Topo+i
        var_User2d_static_XY_name(indx) = temp_2ds_name_Topo(i)
        var_User2d_static_XY_unit(indx) = temp_2ds_unit_Topo(i)
        var_User2d_static_XY_lname(indx)= temp_2ds_lname_Topo(i)
        var_User2d_static_XY_MissVal(indx)= temp_2ds_MissVal_Topo(i)
        var_User2d_static_XY_FillVal(indx)= temp_2ds_FillVal_Topo(i)
        if(i.eq.1) &
          var_User2d_static_XY(1:nxmax,1:nymax,indx) = &
           real(topo_comp(1:nxmax,1:nymax),kind=op)
      enddo

      end subroutine Prep_output_Topo

!******************************************************************************

      subroutine Deallocate_Topo

      implicit none

      deallocate(topo_comp)
      deallocate(topo_indx)

      end subroutine Deallocate_Topo

!##############################################################################
!
!    Get_Topo
!
!##############################################################################

      subroutine Get_Topo

      use mesh,          only : &
         IsLatLon,lon_cc_pd,lat_cc_pd,nxmax,nymax

      use MetReader,       only : &
         nx_submet,ny_submet,x_submet_sp,y_submet_sp

      implicit none

      INTERFACE
        subroutine get_minmax_lonlat(lonmin,lonmax,latmin,latmax)
          real(kind=8) :: lonmin
          real(kind=8) :: lonmax
          real(kind=8) :: latmin
          real(kind=8) :: latmax
        end subroutine 
      END INTERFACE

      ! First we need to get the extents of the computational grid
!      if(IsLatLon)then
        !Just get min and max of lat and lon.
        ! These were already calculated in calc_grid under the names
        ! lonmin,lonmax,latmin,latmax
!        minlon_Topo = min(minval(lon_cc_pd(-1:nxmax+1)),real(minval(x_submet_sp(1:nx_submet)),kind=ip))
!        maxlon_Topo = max(maxval(lon_cc_pd(-1:nxmax+1)),real(maxval(x_submet_sp(1:nx_submet)),kind=ip))
!        minlat_Topo = min(minval(lat_cc_pd(-1:nymax+1)),real(minval(y_submet_sp(1:ny_submet)),kind=ip))
!        maxlat_Topo = max(maxval(lat_cc_pd(-1:nymax+1)),real(maxval(y_submet_sp(1:ny_submet)),kind=ip))

        minlon_Topo = real(minval(x_submet_sp(1:nx_submet)),kind=ip)
        maxlon_Topo = real(maxval(x_submet_sp(1:nx_submet)),kind=ip)
        minlat_Topo = real(minval(y_submet_sp(1:ny_submet)),kind=ip)
        maxlat_Topo = real(maxval(y_submet_sp(1:ny_submet)),kind=ip)

!      else
!        ! This function is in Calc_Mesh
!        write(*,*)"HFS Check this."
!        call get_minmax_lonlat(minlon_Topo,maxlon_Topo,minlat_Topo,maxlat_Topo)
!      endif

      ! ETOPO and GEBCO data are provided on -180->180 grid so map to that domain
      if(minlon_Topo.ge.180_ip)then
        minlon_Topo = minlon_Topo-360.0_ip
        maxlon_Topo = maxlon_Topo-360.0_ip
      endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"min max lon = ",real(minlon_Topo,kind=4),real(maxlon_Topo,kind=4)
        write(outlog(io),*)"min max lat = ",real(minlat_Topo,kind=4),real(maxlat_Topo,kind=4)
      endif;enddo


      if(topoFormat.eq.1)then
#ifdef USENETCDF
        call Load_Topo_Gridded_NC
#endif
      elseif(topoFormat.eq.2)then
        call Load_Topo_Gridded_bin
      elseif(topoFormat.eq.3)then
        !call Load_Topo_ESRI_ASCII
      else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Topography input format not recognized."
          write(errlog(io),*)"    Available options:"
          write(errlog(io),*)"      Topo Format ID = 1 : NetCDF"
          write(errlog(io),*)"      Topo Format ID = 2 : Gridded binary"
          write(errlog(io),*)"      Topo Format ID = 3 : ESRI ASCII"
        endif;enddo
        stop 1
      endif

      call Interp_Topo

      call Smooth_Topo

      end subroutine Get_Topo

!##############################################################################
!
!    Load_Topo_Gridded_NC
!
!##############################################################################

#ifdef USENETCDF

      subroutine Load_Topo_Gridded_NC

      use netcdf
      use Ash3d_Netcdf_IO

      implicit none

      integer :: nSTAT
      integer :: ncid
      integer :: ivar,in_var_id,var_ndims

      integer :: lat_dim_id        = 0
      integer :: lon_dim_id        = 0
      integer :: lat_var_id        = 0
      integer :: lon_var_id        = 0
      integer :: topo_var_id       = 0

      character(len=NF90_MAX_NAME)  :: invar,dimname
      integer :: dim_xtype,dimlen,i_dim
      integer,dimension(:),allocatable :: var_dimIDs
      integer :: var_xtype,var_id,idx
      integer :: xtype, length, attnum

      integer(kind=2), dimension(:)   ,allocatable :: dum1d_short
      integer(kind=2), dimension(:,:) ,allocatable :: dum2d_short

      real(kind=sp)   ,dimension(:)   ,allocatable :: temp1d_sp
      real(kind=dp)   ,dimension(:)   ,allocatable :: temp1d_dp
     
      integer :: nlat_tot,nlon_tot
      integer :: start_lat_idx,start_lon_idx
      integer :: end_lat_idx,end_lon_idx
      integer :: ilat,ilon

      logical :: x_inverted     = .false.
      logical :: y_inverted     = .false.
      logical :: wrapgrid       = .false.

      integer(kind=2) ,dimension(:,:) ,allocatable :: temp2d_intS
      integer(kind=4) ,dimension(:,:) ,allocatable :: temp2d_intL
      real(kind=sp)   ,dimension(:,:) ,allocatable :: temp2d_sp
      real(kind=dp)   ,dimension(:,:) ,allocatable :: temp2d_dp
      integer :: i,ict, ileft(2),iright(2)   !if wrapgrid=.true. ict=2 and left & iright have 2 values, otherwise 1
      integer :: iistart(2),iicount(2)     !if (wrapgrid), iistart(1)=istart, iistart(2)=1

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Inside Load_Topo_Gridded_NC"
      endif;enddo
        ! Check to see if the domain straddles the anti-meridian
      if (minlon_Topo.lt.180.0_ip.and.maxlon_Topo.gt.180.0_ip)then
        lon_shift_flag = 1
      else
        lon_shift_flag = 0
      endif

      nSTAT = nf90_open(adjustl(trim(file_topo)),NF90_NOWRITE,ncid)
      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"nf90_open topofile")

      ! First look up the varaible containing the topographic data
      invar = 'elevation'  ! This is the default for GEBCO
      nSTAT = nf90_inq_varid(ncid,invar,topo_var_id)
      if(nSTAT.ne.NF90_NOERR)then
        call NC_check_status(nSTAT,0,"inq_varid elevation")
        do io=1,nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'  Cannot find variable ',trim(adjustl(invar))
          write(outlog(io),*)'  Testing for known synonyms'
        endif;enddo
        invar = 'z'
        nSTAT = nf90_inq_varid(ncid,invar,topo_var_id)  ! get the var_id for topo
        if(nSTAT.ne.NF90_NOERR)then
          call NC_check_status(nSTAT,0,"inq_varid z")
          do io=1,nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'  Cannot find variable elevation or z'
            write(outlog(io),*)'  Unknown topography format'
          endif;enddo
          stop 1
        endif
      endif

      ! We have found either elevation or z; now determine variable type and dimensions
      nSTAT = nf90_inquire_variable(ncid, topo_var_id, invar, &
                    xtype = var_xtype, &
                    ndims = var_ndims)   ! get the number of dimensions

      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"inq_variable")
      if(.not.allocated(var_dimIDs))allocate(var_dimIDs(var_ndims))
      nSTAT = nf90_inquire_variable(ncid, topo_var_id, invar, &
                dimids = var_dimIDs(:var_ndims))
      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"inq_variable")

      if(var_ndims.eq.1)then
        ! The deprecated GEBCO 08 data stores elevation in one long array
        ! Need to use a special reader for this
        write(*,*)"Looks like GEBCO 08"
        write(*,*)"Calling special reader"
        nSTAT = nf90_close(ncid)
        call Load_Topo_Gridded_NC_GEBCO08
        return
      endif

      ! get the dimension information for lon (or x)
      i_dim = 1
      nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                   name =  dimname, &
                   len = dimlen)
      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"nf90_inquire_dimension X")
      nlon_topo_fullgrid = dimlen
      lon_dim_id    = var_dimIDs(i_dim)

      nSTAT = nf90_inq_varid(ncid,dimname,var_id) ! get the variable associated with this dim
      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"inq_variable X")
      ! Check what temporary array to use
      nSTAT = nf90_inquire_variable(ncid, var_id, dimname, xtype = dim_xtype)
      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"nf90_inquire_variable X")
      allocate(lon_topo_fullgrid(1:nlon_topo_fullgrid))
      if(dim_xtype.eq.NF90_FLOAT)then
        allocate(temp1d_sp(dimlen))
        nSTAT = nf90_get_var(ncid,var_id,temp1d_sp, &
                       start = (/1/),count = (/dimlen/))
        if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"get_var X flt")
        ! copy to local variable
        lon_topo_fullgrid(1:nlon_topo_fullgrid) = real(temp1d_sp(1:nlon_topo_fullgrid),kind=dp)
        deallocate(temp1d_sp)
      elseif(dim_xtype.eq.NF90_DOUBLE)then
        allocate(temp1d_dp(dimlen))
        nSTAT = nf90_get_var(ncid,var_id,temp1d_dp, &
                       start = (/1/),count = (/dimlen/))
        if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"get_var X dbl")
        ! copy to local variable
        lon_topo_fullgrid(1:nlon_topo_fullgrid) = temp1d_dp(1:nlon_topo_fullgrid)
        deallocate(temp1d_dp)
      else
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Cannot recognize variable type for X'
          if(dim_xtype.eq.NF90_BYTE)  write(errlog(io),*)"NF90_BYTE = "  ,NF90_BYTE
          if(dim_xtype.eq.NF90_CHAR)  write(errlog(io),*)"NF90_CHAR = "  ,NF90_CHAR
          if(dim_xtype.eq.NF90_SHORT) write(errlog(io),*)"NF90_SHORT = " ,NF90_SHORT
          if(dim_xtype.eq.NF90_INT)   write(errlog(io),*)"NF90_INT = "   ,NF90_INT
          if(dim_xtype.eq.NF90_FLOAT) write(errlog(io),*)"NF90_FLOAT = " ,NF90_FLOAT
          if(dim_xtype.eq.NF90_DOUBLE)write(errlog(io),*)"NF90_DOUBLE = ",NF90_DOUBLE
          if(dim_xtype.eq.NF90_UBYTE) write(errlog(io),*)"NF90_UBYTE = " ,NF90_UBYTE
          if(dim_xtype.eq.NF90_USHORT)write(errlog(io),*)"NF90_USHORT = ",NF90_USHORT
          if(dim_xtype.eq.NF90_UINT)  write(errlog(io),*)"NF90_UINT = "  ,NF90_UINT
          if(dim_xtype.eq.NF90_INT64) write(errlog(io),*)"NF90_INT64 = " ,NF90_INT64
          if(dim_xtype.eq.NF90_UINT64)write(errlog(io),*)"NF90_UINT64 = ",NF90_UINT64
          if(dim_xtype.eq.NF90_STRING)write(errlog(io),*)"NF90_STRING = ",NF90_STRING
        endif;enddo
        stop 1
      endif
      dlon_topo = lon_topo_fullgrid(2) - lon_topo_fullgrid(1)
      if(dlon_topo.lt.0.0_dp)then
        x_inverted = .true.
        dlon_topo = abs(dlon_topo)
      endif

      ! get the dimension information for lat (or y)
      i_dim = 2
      nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                   name =  dimname, &
                   len = dimlen)
      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"nf90_inquire_dimension Y")
      nlat_topo_fullgrid = dimlen
      lat_dim_id    = var_dimIDs(i_dim)

      nSTAT = nf90_inq_varid(ncid,dimname,var_id) ! get the variable associated with this dim
      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"inq_variable Y")
      ! Check what temporary array to use
      nSTAT = nf90_inquire_variable(ncid, var_id, dimname, xtype = dim_xtype)
      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"nf90_inquire_variable Y")
      allocate(lat_topo_fullgrid(1:nlat_topo_fullgrid))
      if(dim_xtype.eq.NF90_FLOAT)then
        allocate(temp1d_sp(dimlen))
        nSTAT = nf90_get_var(ncid,var_id,temp1d_sp, &
                       start = (/1/),count = (/dimlen/))
        if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"get_var Y flt")
        ! copy to local variable
        lat_topo_fullgrid(1:nlat_topo_fullgrid) = real(temp1d_sp(1:nlat_topo_fullgrid),kind=dp)
        deallocate(temp1d_sp)
      elseif(dim_xtype.eq.NF90_DOUBLE)then
        allocate(temp1d_dp(dimlen))
        nSTAT = nf90_get_var(ncid,var_id,temp1d_dp, &
                       start = (/1/),count = (/dimlen/))
        if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"get_var Y dbl")
        ! copy to local variable
        lat_topo_fullgrid(1:nlat_topo_fullgrid) = temp1d_dp(1:nlat_topo_fullgrid)
        deallocate(temp1d_dp)
      else
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Cannot recognize variable type for Y'
          if(dim_xtype.eq.NF90_BYTE)  write(errlog(io),*)"NF90_BYTE = "  ,NF90_BYTE
          if(dim_xtype.eq.NF90_CHAR)  write(errlog(io),*)"NF90_CHAR = "  ,NF90_CHAR
          if(dim_xtype.eq.NF90_SHORT) write(errlog(io),*)"NF90_SHORT = " ,NF90_SHORT
          if(dim_xtype.eq.NF90_INT)   write(errlog(io),*)"NF90_INT = "   ,NF90_INT
          if(dim_xtype.eq.NF90_FLOAT) write(errlog(io),*)"NF90_FLOAT = " ,NF90_FLOAT
          if(dim_xtype.eq.NF90_DOUBLE)write(errlog(io),*)"NF90_DOUBLE = ",NF90_DOUBLE
          if(dim_xtype.eq.NF90_UBYTE) write(errlog(io),*)"NF90_UBYTE = " ,NF90_UBYTE
          if(dim_xtype.eq.NF90_USHORT)write(errlog(io),*)"NF90_USHORT = ",NF90_USHORT
          if(dim_xtype.eq.NF90_UINT)  write(errlog(io),*)"NF90_UINT = "  ,NF90_UINT
          if(dim_xtype.eq.NF90_INT64) write(errlog(io),*)"NF90_INT64 = " ,NF90_INT64
          if(dim_xtype.eq.NF90_UINT64)write(errlog(io),*)"NF90_UINT64 = ",NF90_UINT64
          if(dim_xtype.eq.NF90_STRING)write(errlog(io),*)"NF90_STRING = ",NF90_STRING
        endif;enddo
        stop 1
      endif
      dlat_topo = lat_topo_fullgrid(2) - lat_topo_fullgrid(1)
      if(dlat_topo.lt.0.0_dp)then
        y_inverted = .true.
        dlat_topo = abs(dlat_topo)
      endif

      ! All of the etopo and gebco files have start longitudes between -180 and 180
      ! Sub grids that straddle the anti-meridian have max longitudes > 180
      ! Now allocate space for the sub-grid of the topo file at the natural resolution
      ! of the topography data
      if(x_inverted.or.y_inverted)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Not yet set up for inverted x or y topo dimensions."
        endif;enddo
        stop 1
      endif

      ! Check that the computational grid is in the correct periodic mapping to the topo grid
      if(minlon_Topo.lt.lon_topo_fullgrid(1)-0.5_dp*dlon_topo)then
        minlon_Topo = minlon_Topo + 360.0_dp
        maxlon_Topo = maxlon_Topo + 360.0_dp
      endif
      start_lon_idx = -1
      cleft = lon_topo_fullgrid(1)-0.5_dp*dlon_topo
      do ilon=1,nlon_topo_fullgrid
        cright = lon_topo_fullgrid(ilon)+0.5_dp*dlon_topo
        if(minlon_Topo.ge.cleft.and. &
           minlon_Topo.lt.cright)then
          start_lon_idx = ilon
        endif
        cleft = cright
      enddo
      nlon_topo_subgrid = floor((maxlon_Topo-minlon_Topo)/dlon_topo)
      if(start_lon_idx.lt.1.or.start_lon_idx.gt.nlon_topo_fullgrid)then
        ! Couldn't find start x
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Couldn't find start x of topo sub-grid"
        endif;enddo
        stop 1
      endif
      if(start_lon_idx+nlon_topo_subgrid.gt.nlon_topo_fullgrid)then
        wrapgrid = .true.
      endif
      ! Load data variables for just the subgrid defined above
      if(wrapgrid)then
        ict        = 2
          ! index on the topo-sub-grid
        ileft(1)   = 1
        iright(1)  = nlon_topo_fullgrid - start_lon_idx
        iistart(1) = start_lon_idx
        iicount(1) = nlon_topo_fullgrid - start_lon_idx
        ileft(2)   = iright(1) + 1
        iright(2)  = nlon_topo_subgrid
        iistart(2) = 1
        iicount(2) = nlon_topo_subgrid - iicount(1)
      else
        ict        = 1
        ileft(1)   = 1
        iright(1)  = nlon_topo_subgrid
        iistart(1) = start_lon_idx
        iicount(1) = nlon_topo_subgrid
      endif

      start_lat_idx = -1
      end_lat_idx   = -1
      cleft = lat_topo_fullgrid(1)-0.5_dp*dlat_topo
      do ilat=1,nlat_topo_fullgrid
        cright = lat_topo_fullgrid(ilat)+0.5_dp*dlat_topo
        if(minlat_Topo.ge.cleft.and. &
           minlat_Topo.lt.cright)then
          start_lat_idx = ilat
        endif
        if(maxlat_Topo.gt.cleft.and. &
           maxlat_Topo.le.cright)then
          end_lat_idx = ilat
        endif
        cleft = cright
      enddo
      if(start_lat_idx.lt.1.or.start_lat_idx.gt.nlat_topo_fullgrid)then
        ! Couldn't find start y
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Couldn't find start y of topo sub-grid"
        endif;enddo
        stop 1
      endif
      if(end_lat_idx.lt.1.or.end_lat_idx.gt.nlat_topo_fullgrid)then
        ! Couldn't find end y
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Couldn't find end y of topo sub-grid"
        endif;enddo
        stop 1
      endif

      nlat_topo_subgrid = end_lat_idx - start_lat_idx + 1
      allocate(lon_topo_subgrid(nlon_topo_subgrid))
      allocate(lat_topo_subgrid(nlat_topo_subgrid))
      allocate(topo_subgrid(nlon_topo_subgrid,nlat_topo_subgrid))
      do ilon=1,nlon_topo_subgrid
        lon_topo_subgrid(ilon)=lon_topo_fullgrid(start_lon_idx)+(ilon-1)*dlat_topo
      enddo
      lat_topo_subgrid(1:nlat_topo_subgrid)=lat_topo_fullgrid(start_lat_idx:end_lat_idx)

      ! Now re-inquire about the topo variables and load the data
      nSTAT = nf90_inquire_variable(ncid, topo_var_id, invar, &
                    xtype = var_xtype, &
                    ndims = var_ndims)   ! get the number of dimensions
      if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"nf90_inquire_variable topo")

      ! Check what temporary array to use; this shoult be short, but to keep things
      ! general, check int, float and double too.
      if(var_xtype.eq.NF90_SHORT)then
        do i=1,ict        !read subgrid at current time step
          allocate(temp2d_intS(iicount(i),nlat_topo_subgrid))
          nSTAT = nf90_get_var(ncid,topo_var_id,temp2d_intS(:,:), &
                               start = (/iistart(i),start_lat_idx/), &
                               count = (/iicount(i),nlat_topo_subgrid/))
          if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"get_var topo short")
          ! copy to local variable
          topo_subgrid(ileft(i):iright(i),1:nlat_topo_subgrid) = &
                  real(temp2d_intS(1:iicount(i),1:nlat_topo_subgrid),kind=sp)
          deallocate(temp2d_intS)
        enddo
      elseif(var_xtype.eq.NF90_INT)then
        do i=1,ict        !read subgrid at current time step
          allocate(temp2d_intL(iicount(i),nlat_topo_subgrid))
          nSTAT = nf90_get_var(ncid,topo_var_id,temp2d_intL(:,:), &
                         start = (/iistart(i),start_lat_idx/), &
                         count = (/iicount(i),nlat_topo_subgrid/))
          if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"get_var topo int")
          ! copy to local variable
          topo_subgrid(ileft(i):iright(i),1:nlat_topo_subgrid) = &
                  real(temp2d_intL(1:iicount(i),1:nlat_topo_subgrid),kind=sp)
          deallocate(temp2d_intL)
        enddo
      elseif(var_xtype.eq.NF90_FLOAT)then
        do i=1,ict        !read subgrid at current time step
          allocate(temp2d_sp(iicount(i),nlat_topo_subgrid))
          nSTAT = nf90_get_var(ncid,topo_var_id,temp2d_sp(:,:), &
                         start = (/iistart(i),start_lat_idx/), &
                         count = (/iicount(i),nlat_topo_subgrid/))
          if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"get_var Y flt")
          ! copy to local variable
          topo_subgrid(ileft(i):iright(i),1:nlat_topo_subgrid) = &
                  real(temp2d_sp(1:iicount(i),1:nlat_topo_subgrid),kind=sp)
          deallocate(temp2d_sp)
        enddo
      elseif(var_xtype.eq.NF90_DOUBLE)then
        do i=1,ict        !read subgrid at current time step
          allocate(temp2d_dp(iicount(i),nlat_topo_subgrid))
          nSTAT = nf90_get_var(ncid,topo_var_id,temp2d_dp(:,:), &
                         start = (/iistart(i),start_lat_idx/), &
                         count = (/iicount(i),nlat_topo_subgrid/))
          if(nSTAT.ne.NF90_NOERR)call NC_check_status(nSTAT,1,"get_var Y dbl")
          ! copy to local variable
          topo_subgrid(ileft(i):iright(i),1:nlat_topo_subgrid) = &
                  real(temp2d_dp(1:iicount(i),1:nlat_topo_subgrid),kind=sp)
          deallocate(temp2d_dp)
        enddo
      else
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Cannot recognize variable type for topo'
          if(dim_xtype.eq.NF90_BYTE)  write(errlog(io),*)"NF90_BYTE = "  ,NF90_BYTE
          if(dim_xtype.eq.NF90_CHAR)  write(errlog(io),*)"NF90_CHAR = "  ,NF90_CHAR
          if(dim_xtype.eq.NF90_SHORT) write(errlog(io),*)"NF90_SHORT = " ,NF90_SHORT
          if(dim_xtype.eq.NF90_INT)   write(errlog(io),*)"NF90_INT = "   ,NF90_INT
          if(dim_xtype.eq.NF90_FLOAT) write(errlog(io),*)"NF90_FLOAT = " ,NF90_FLOAT
          if(dim_xtype.eq.NF90_DOUBLE)write(errlog(io),*)"NF90_DOUBLE = ",NF90_DOUBLE
          if(dim_xtype.eq.NF90_UBYTE) write(errlog(io),*)"NF90_UBYTE = " ,NF90_UBYTE
          if(dim_xtype.eq.NF90_USHORT)write(errlog(io),*)"NF90_USHORT = ",NF90_USHORT
          if(dim_xtype.eq.NF90_UINT)  write(errlog(io),*)"NF90_UINT = "  ,NF90_UINT
          if(dim_xtype.eq.NF90_INT64) write(errlog(io),*)"NF90_INT64 = " ,NF90_INT64
          if(dim_xtype.eq.NF90_UINT64)write(errlog(io),*)"NF90_UINT64 = ",NF90_UINT64
          if(dim_xtype.eq.NF90_STRING)write(errlog(io),*)"NF90_STRING = ",NF90_STRING
        endif;enddo
        stop 1
      endif

      nSTAT = nf90_close(ncid)

      end subroutine Load_Topo_Gridded_NC
#endif

!##############################################################################
!
!    Load_Topo_Gridded_NC_GEBCO08
!
!##############################################################################
#ifdef USENETCDF

      subroutine Load_Topo_Gridded_NC_GEBCO08

      use netcdf

      implicit none

      integer :: nSTAT
      integer :: ncid

      integer :: topo_var_id       = 0

      integer(kind=2), dimension(:)   ,allocatable :: dum1d_short
      integer(kind=2), dimension(:,:) ,allocatable :: dum2d_short

      integer :: nlat_tot,nlon_tot
      integer :: start_lat_idx,start_lon_idx,end_lon_idx
      integer :: ilat,ilon,idx

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Inside load_topo"
      endif;enddo
        ! Check to see if the domain straddles the anti-meridian
      if (minlon_Topo.lt.180.0_ip.and.maxlon_Topo.gt.180.0_ip)then
        lon_shift_flag = 1
      else
        lon_shift_flag = 0
      endif

      nSTAT = nf90_open(adjustl(trim(file_topo)),NF90_NOWRITE,ncid)
      if(nSTAT.ne.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'Could not open -',adjustl(trim(file_topo)),'-'
          write(errlog(io),*)NF90_NOWRITE,ncid
          write(errlog(io),*)'Exiting'
        endif;enddo
        stop 1
      endif
      nSTAT = nf90_inq_varid(ncid,'z',topo_var_id)

      ! GEBCO 08 (30-second global topo/batho)
      nlat_topo_fullgrid = 21600
      nlon_topo_fullgrid = 43200
      dlon_topo = 1.0_dp/120.0_dp
      dlat_topo = 1.0_dp/120.0_dp
      start_lat_idx = int((minlat_Topo+90.0_ip)/dlon_topo)
      start_lat_idx = nlat_topo_fullgrid-start_lat_idx

      ! Define the sub-grid holding the topo data we need
      nlon_topo_subgrid = int((maxlon_Topo-minlon_Topo)/dlon_topo) + 4
      nlat_topo_subgrid = int((maxlat_Topo-minlat_Topo)/dlat_topo) + 4
      allocate(lon_topo_subgrid(nlon_topo_subgrid))
      allocate(lat_topo_subgrid(nlat_topo_subgrid))
      allocate(topo_subgrid(nlon_topo_subgrid,nlat_topo_subgrid))
      ! GEBCO08 has data that start at -180, so shift 0 accordingly
      if(minlon_Topo.lt.0.0_dp)then
        start_lon_idx = int((minlon_Topo+180.0_dp)/dlon_topo)
      elseif(minlon_Topo.ge.0.0_dp.and.minlon_Topo.lt.180.0_dp)then
        start_lon_idx = int((minlon_Topo+180.0_dp)/dlon_topo)
      elseif(minlon_Topo.ge.180.0_dp)then
        start_lon_idx = int((minlon_Topo-180.0_dp)/dlon_topo)
      endif

!      if(lon_shift_flag.eq.0)then
!        start_lon_idx = int((minlon_Topo-180.0_dp)/dlon_topo)
!      else
!        start_lon_idx = int((minlon_Topo+180.0_dp)/dlon_topo)
!      endif

      do ilon=1,nlon_topo_subgrid
        lon_topo_subgrid(ilon) = real(start_lon_idx+ilon-1,kind=ip)*dlon_topo &
                         - 0.5_dp*dlon_topo - 180.0_dp
      enddo
      ! now shift the lon values if they start in the western hemisphere

      if(lon_shift_flag.eq.0)then
        end_lon_idx = nlon_topo_subgrid
      else
        end_lon_idx = int((180.0_dp-minlon_Topo)/dlon_topo)
      endif

      ! GEBCO_08 start at the NW corner of the grid
      ! (89d59'45"N,179d59;45"W) and advances eastward, then
      ! southwards in one long string.
      allocate(dum1d_short(nlon_topo_fullgrid))
      do ilat=1,nlat_topo_subgrid
        ! Get the index point of the start of the line at the right
        ! latitude
        lat_topo_subgrid(ilat)=(real(nlat_topo_fullgrid-(start_lat_idx-ilat+1),kind=ip)*dlat_topo &
                       -0.5_dp*dlat_topo)-90.0_dp
        idx = (start_lat_idx-ilat)*nlon_topo_fullgrid +1
        ! Get the full row of topo values that encircle the globe
!        write(*,*)ilat,idx,nlon_topo_fullgrid
        nSTAT = nf90_get_var(ncid,topo_var_id,dum1d_short, &
                           start = (/idx/),        &
                           count = (/nlon_topo_fullgrid/))
        ! Now find just the subset that we need
        if(lon_shift_flag.eq.0)then
          topo_subgrid(1:nlon_topo_subgrid,ilat)=real(dum1d_short(start_lon_idx:start_lon_idx+nlon_topo_subgrid-1),kind=sp)
        else
          ! Copy the part west of the anti-meridian
          topo_subgrid(1:end_lon_idx,ilat) = real(dum1d_short(start_lon_idx:start_lon_idx+end_lon_idx-1),kind=sp)
          ! And copy the part wrapped over the anti-meridian
!          write(*,*)end_lon_idx,nlon_topo_subgrid
          topo_subgrid(end_lon_idx+1:nlon_topo_subgrid,ilat) = real(dum1d_short(1:nlon_topo_subgrid-end_lon_idx),kind=sp)
        endif
      enddo
      deallocate(dum1d_short)
      nSTAT = nf90_close(ncid)

      end subroutine Load_Topo_Gridded_NC_GEBCO08
#endif

!##############################################################################
!
!    Load_Topo_Gridded_bin
!
!##############################################################################

      subroutine Load_Topo_Gridded_bin

      use global_param,  only : &
         IsLitEnd

      use Ash3d_Netcdf_IO

      implicit none

      character (len=50)  :: file_topo_root
      character (len=50)  :: file_topo_hdr
      logical             :: IsThere1,IsThere2

      integer :: substr_pos,substr_pos2

      integer(kind=2), dimension(:)   ,allocatable :: dum1d_short
      integer(kind=2), dimension(:,:) ,allocatable :: dum2d_short

      real(kind=sp)   ,dimension(:)   ,allocatable :: temp1d_sp
      real(kind=dp)   ,dimension(:)   ,allocatable :: temp1d_dp

      integer :: start_lat_idx,start_lon_idx
      integer :: end_lat_idx,end_lon_idx
      integer :: ilat,ilon

      logical :: y_inverted     = .false.
      logical :: wrapgrid       = .false.

      integer(kind=2) ,dimension(:,:) ,allocatable :: temp2d_intS
      integer(kind=4) ,dimension(:,:) ,allocatable :: temp2d_intL
      real(kind=sp)   ,dimension(:,:) ,allocatable :: temp2d_sp
      real(kind=dp)   ,dimension(:,:) ,allocatable :: temp2d_dp
      integer :: i,ict, ileft(2),iright(2)   !if wrapgrid=.true. ict=2 and left & iright have 2 values, otherwise 1
      integer :: iistart(2),iicount(2)     !if (wrapgrid), iistart(1)=istart, iistart(2)=1

      real(kind=dp) :: startx_topo, starty_topo
      real(kind=dp) :: dy

      character(len=80) :: linebuffer080
      character(len=50) :: linebuffer050
      character         :: testkey
      integer           :: iostatus
      integer           :: ioerr
      character(len=120):: iomessage
      character(len=20) :: key_name
      logical :: IsLitEnd_topo
      logical :: IsInt_topo
      integer :: nlat_tot
      integer :: nlon_tot
      integer :: nbits
      integer :: nodata_int
      real(kind=sp) :: nodata_sp
      real(kind=dp) :: ULx_topo
      real(kind=dp) :: ULy_topo
      real(kind=dp) :: LLx_topo
      real(kind=dp) :: LLy_topo
      logical :: key_found

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Inside Load_Topo_Gridded_bin"
      endif;enddo

      ! Existance of topo file was verified earlier, now look for header; searching for either *HDR or *hdr
      ! First, find the root-extension deliminator '.'.  Two possibilities here:
      !  One or more delimitors -> strip off extension and test for hdr or HDR
      !  No delimitor  -> add hdr or HDR to root and test
      substr_pos = index(file_topo,'.',back=.true.)
      if(substr_pos.gt.0)then
        file_topo_root = file_topo(1:substr_pos-1)
      else
        file_topo_root = file_topo
      endif
      file_topo_hdr = trim(adjustl(file_topo_root)) // '.hdr'
      inquire( file=file_topo_hdr, exist=IsThere1 )
      if(.not.IsThere1)then
        file_topo_hdr = trim(adjustl(file_topo_root)) // '.HDR'
        inquire( file=file_topo_hdr, exist=IsThere2 )
        if(.not.IsThere2)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Cannot find topography hdr file, or HDR file"
            write(errlog(io),*)"     ",file_topo_hdr
          endif;enddo
          stop 1
        endif
      endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Found header file for binary topo data. ",file_topo_hdr
      endif;enddo

      ! Opening header file
      open(unit=fid_hdrfile,file=file_topo_hdr,status='old',action='read',err=9001)
      ! First, look for "BYTEORDER'
      substr_pos = -1
      key_name = 'BYTEORDER'
      key_found = .false.
      read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while(iostatus.eq.0)
        substr_pos = index(linebuffer080,trim(adjustl(key_name)))
        if(substr_pos.gt.0)then
          key_found = .true.
          substr_pos2 = index(linebuffer080(substr_pos:),' ')
          ! Found key, now read it
          read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
          linebuffer050 = trim(adjustl(linebuffer050))
          read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)testkey
          if(testkey.eq.'M'.or.testkey.eq.'B')then
            ! data is in Motorola or Big-Endian format
            IsLitEnd_topo = .false.
          elseif(testkey.eq.'L')then
            ! data is in Little-Endian format
            IsLitEnd_topo = .true.
          else
            ! Cannot determine endian
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Found byteorder, but cannot interpret endian.",linebuffer050
              write(outlog(io),*)"Setting endian-ness to local machine."
            endif;enddo
            IsLitEnd_topo = IsLitEnd
          endif
        endif
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo    

      ! NROWS
      rewind(fid_hdrfile)
      substr_pos = -1
      key_name = 'NROWS'
      key_found = .false.
      read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while(iostatus.eq.0)
        substr_pos = index(linebuffer080,trim(adjustl(key_name)))
        if(substr_pos.gt.0)then
          key_found = .true.
          substr_pos2 = index(linebuffer080(substr_pos:),' ')
          ! Found key, now read it
          read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
          linebuffer050 = trim(adjustl(linebuffer050))
          read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)nlon_topo_fullgrid
          linebuffer050 = "Reading line from topo header file, NROWS"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        endif
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo

      !NCOLS
      rewind(fid_hdrfile)
      substr_pos = -1
      key_name = 'NCOLS'
      key_found = .false.
      read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while(iostatus.eq.0)
        substr_pos = index(linebuffer080,trim(adjustl(key_name)))
        if(substr_pos.gt.0)then
          key_found = .true.
          substr_pos2 = index(linebuffer080(substr_pos:),' ')
          ! Found key, now read it
          read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
          linebuffer050 = trim(adjustl(linebuffer050))
          read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)nlat_topo_fullgrid
          linebuffer050 = "Reading line from topo header file, NCOLS"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        endif
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo

      !NBITS
      rewind(fid_hdrfile)
      substr_pos = -1
      key_name = 'NBITS'
      key_found = .false.
      read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while(iostatus.eq.0)
        substr_pos = index(linebuffer080,trim(adjustl(key_name)))
        if(substr_pos.gt.0)then
          key_found = .true.
          substr_pos2 = index(linebuffer080(substr_pos:),' ')
          ! Found key, now read it
          read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
          linebuffer050 = trim(adjustl(linebuffer050))
          read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)nbits
          linebuffer050 = "Reading line from topo header file, NBITS"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          if(nbits.eq.16)then
            IsInt_topo = .true.
          elseif(nbits.eq.32)then
            IsInt_topo = .false.
          else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: Cannot find binary data type of topo file."
              write(errlog(io),*)"       Expecting NBITS = 16 for signed integer or"
              write(errlog(io),*)"                 NBITS = 32 for floating point"
              write(errlog(io),*)"       Read NBITS = ",nbits
            endif;enddo
            stop 1
          endif
        endif
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo

      !NODATA
      rewind(fid_hdrfile)
      nodata_int = -9999
      nodata_sp  = -9999.0_sp
      substr_pos = -1
      key_name = 'NODATA'
      key_found = .false.
      read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while(iostatus.eq.0)
        substr_pos = index(linebuffer080,trim(adjustl(key_name)))
        if(substr_pos.gt.0)then
          key_found = .true.
          substr_pos2 = index(linebuffer080(substr_pos:),' ')
          ! Found key, now read it
          read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
          linebuffer050 = trim(adjustl(linebuffer050))
          if(IsInt_topo)then
            read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)nodata_int
          else
            read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)nodata_sp
          endif
          linebuffer050 = "Reading line from topo header file, NODATA"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        endif
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo

      !XDIM or DX or DLON
      rewind(fid_hdrfile)
      substr_pos = -1
      key_name = 'XDIM'
      key_found = .false.
      read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while(iostatus.eq.0)
        substr_pos = index(linebuffer080,trim(adjustl(key_name)))
        if(substr_pos.gt.0)then
          key_found = .true.
          substr_pos2 = index(linebuffer080(substr_pos:),' ')
          ! Found key, now read it
          read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
          linebuffer050 = trim(adjustl(linebuffer050))
          read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)dlon_topo
          linebuffer050 = "Reading line from topo header file, XDIM"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        endif
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo

      !YDIM or DY or DLAT
      rewind(fid_hdrfile)
      substr_pos = -1
      key_name = 'YDIM'
      key_found = .false.
      read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while(iostatus.eq.0)
        substr_pos = index(linebuffer080,trim(adjustl(key_name)))
        if(substr_pos.gt.0)then
          key_found = .true.
          substr_pos2 = index(linebuffer080(substr_pos:),' ')
          ! Found key, now read it
          read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
          linebuffer050 = trim(adjustl(linebuffer050))
          read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)dlat_topo
          linebuffer050 = "Reading line from topo header file, YDIM"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        endif
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo

      !ULYMAP or ULY or ULLAT
      y_inverted = .false.  ! initialize
      rewind(fid_hdrfile)
      substr_pos = -1
      key_name = 'ULYMAP'
      key_found = .false.
      read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while(iostatus.eq.0)
        substr_pos = index(linebuffer080,trim(adjustl(key_name)))
        if(substr_pos.gt.0)then
          key_found = .true.
          ! If this key has been found, the upper-left is the start coordinate
          y_inverted = .true.
          substr_pos2 = index(linebuffer080(substr_pos:),' ')
          ! Found key, now read it
          read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
          linebuffer050 = trim(adjustl(linebuffer050))
          read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)starty_topo
          linebuffer050 = "Reading line from topo header file, ULYMAP"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        endif
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo

      if(.not.y_inverted)then
        ! If we didn't find ULYMAP or something analogous above, then look for the lower-left LL
        !LLYMAP or LLY or LLLAT
        rewind(fid_hdrfile)
        substr_pos = -1
        key_name = 'LLYMAP'
        key_found = .false.
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        do while(iostatus.eq.0)
          substr_pos = index(linebuffer080,trim(adjustl(key_name)))
          if(substr_pos.gt.0)then
            key_found = .true.
            ! If this key has been found, the upper-left is the start coordinate
            y_inverted = .true.
            substr_pos2 = index(linebuffer080(substr_pos:),' ')
            ! Found key, now read it
            read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
            linebuffer050 = trim(adjustl(linebuffer050))
            read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)starty_topo
            linebuffer050 = "Reading line from topo header file, LLYMAP"
            if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          endif
          read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        enddo
      endif

      !ULXMAP or ULX or ULLON of LLXMAP or LLX or LLLON
      rewind(fid_hdrfile)
      substr_pos = -1
      key_name = 'ULXMAP'
      key_found = .false.
      read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while(iostatus.eq.0)
        substr_pos = index(linebuffer080,trim(adjustl(key_name)))
        if(substr_pos.gt.0)then
          key_found = .true.
          ! If this key has been found, the upper-left is the start coordinate
          y_inverted = .true.
          substr_pos2 = index(linebuffer080(substr_pos:),' ')
          ! Found key, now read it
          read(linebuffer080(substr_pos2+1:substr_pos2+50),'(a50)',iostat=ioerr,iomsg=iomessage)linebuffer050
          linebuffer050 = trim(adjustl(linebuffer050))
          read(linebuffer050,*,iostat=iostatus,iomsg=iomessage)startx_topo
          linebuffer050 = "Reading line from topo header file, ULXMAP"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        endif
        read(fid_hdrfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo

      close(fid_hdrfile)

      ! Now set up grids 
      allocate(lon_topo_fullgrid(1:nlon_topo_fullgrid))
      allocate(lat_topo_fullgrid(1:nlat_topo_fullgrid))
      do ilon=1,nlon_topo_fullgrid
        lon_topo_fullgrid(ilon) = startx_topo + (ilon-1)*dlon_topo
      enddo
      if(y_inverted)then
        dy = -1.0_dp * dlat_topo
      else
        dy = dlat_topo
      endif
      do ilat=1,nlat_topo_fullgrid
        lat_topo_fullgrid(ilat) = starty_topo + (ilat-1)*dy
      enddo

      ! Double-check that this lon/lat range suffices for the requested computational grid
      nlon_topo_subgrid = floor((maxlon_Topo-minlon_Topo)/dlon_topo)
      nlat_topo_subgrid = int((maxlat_Topo-minlat_Topo)/dlat_topo) 

      if(minlon_Topo.lt.lon_topo_fullgrid(1)-0.5_dp*dlon_topo)then
        minlon_Topo = minlon_Topo + 360.0_dp
        maxlon_Topo = maxlon_Topo + 360.0_dp
      endif
      start_lon_idx = -1
      cleft = lon_topo_fullgrid(1)-0.5_dp*dlon_topo
      do ilon=1,nlon_topo_fullgrid
        cright = lon_topo_fullgrid(ilon)+0.5_dp*dlon_topo
        if(minlon_Topo.ge.cleft.and. &
           minlon_Topo.lt.cright)then
          start_lon_idx = ilon
        endif
        cleft = cright
      enddo
      if(start_lon_idx.lt.1.or.start_lon_idx.gt.nlon_topo_fullgrid)then
        ! Couldn't find start x
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Couldn't find start x of topo sub-grid"
        endif;enddo
        stop 1
      endif
      if(start_lon_idx+nlon_topo_subgrid.gt.nlon_topo_fullgrid)then
        wrapgrid = .true.
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Cannot use grid-wrapping for binary data, at the moment."
        endif;enddo
        stop 1
      else
        end_lon_idx = start_lon_idx+nlon_topo_subgrid
      endif
      !! Load data variables for just the subgrid defined above
      !if(wrapgrid)then
      !  ict        = 2
      !    ! index on the topo-sub-grid
      !  ileft(1)   = 1
      !  iright(1)  = nlon_topo_fullgrid - start_lon_idx
      !  iistart(1) = start_lon_idx
      !  iicount(1) = nlon_topo_fullgrid - start_lon_idx
      !  ileft(2)   = iright(1) + 1
      !  iright(2)  = nlon_topo_subgrid
      !  iistart(2) = 1
      !  iicount(2) = nlon_topo_subgrid - iicount(1)
      !else
      !  ict        = 1
      !  ileft(1)   = 1
      !  iright(1)  = nlon_topo_subgrid
      !  iistart(1) = start_lon_idx
      !  iicount(1) = nlon_topo_subgrid
      !endif

      start_lat_idx = -1
      end_lat_idx   = -1
      cleft = lat_topo_fullgrid(1)-0.5_dp*dy
      do ilat=1,nlat_topo_fullgrid
        cright = lat_topo_fullgrid(ilat)+0.5_dp*dy
        if(y_inverted)then
          if(minlat_Topo.lt.cleft.and.minlat_Topo.ge.cright)then
            start_lat_idx = ilat
          endif
          if(maxlat_Topo.le.cleft.and.maxlat_Topo.gt.cright)then
            end_lat_idx = ilat
          endif
        else
          if(minlat_Topo.ge.cleft.and.minlat_Topo.lt.cright)then
            start_lat_idx = ilat
          endif
          if(maxlat_Topo.gt.cleft.and.maxlat_Topo.le.cright)then
            end_lat_idx = ilat
          endif
        endif
        cleft = cright
      enddo
      if(start_lat_idx.lt.1.or.start_lat_idx.gt.nlat_topo_fullgrid)then
        ! Couldn't find start y
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Couldn't find start y of topo sub-grid"
        endif;enddo
        stop 1
      endif
      if(end_lat_idx.lt.1.or.end_lat_idx.gt.nlat_topo_fullgrid)then
        ! Couldn't find end y
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Couldn't find end y of topo sub-grid"
        endif;enddo
        stop 1
      endif

      write(*,*)start_lat_idx,end_lat_idx
      write(*,*)start_lon_idx,end_lon_idx


      ! Define the sub-grid holding the topo data we need
      write(*,*)"Allocating arrays ",nlon_topo_subgrid,nlat_topo_subgrid
      allocate(lon_topo_subgrid(nlon_topo_subgrid))
      allocate(lat_topo_subgrid(nlat_topo_subgrid))
      allocate(topo_subgrid(nlon_topo_subgrid,nlat_topo_subgrid))

      lon_topo_subgrid(1:nlon_topo_subgrid) = lon_topo_fullgrid(start_lon_idx:end_lon_idx-1)
!      if(y_inverted)then
!        do ilat=1,nlat_topo_subgrid
!          
!        enddo
!      else
        lat_topo_subgrid(1:nlat_topo_subgrid) = lat_topo_fullgrid(start_lat_idx:end_lat_idx-1)
!      endif

       write(*,*)lon_topo_subgrid(1),lon_topo_subgrid(nlon_topo_subgrid)
       write(*,*)lat_topo_subgrid(1),lat_topo_subgrid(nlat_topo_subgrid)


!      ! GEBCO08 has data that start at -180, so shift 0 accordingly
!      if(minlon_Topo.lt.0.0_dp)then
!        start_lon_idx = int((minlon_Topo+180.0_dp)/dlon_topo)
!      elseif(minlon_Topo.ge.0.0_dp.and.minlon_Topo.lt.180.0_dp)then
!        start_lon_idx = int((minlon_Topo+180.0_dp)/dlon_topo)
!      elseif(minlon_Topo.ge.180.0_dp)then
!        start_lon_idx = int((minlon_Topo-180.0_dp)/dlon_topo)
!      endif
!
!!      if(lon_shift_flag.eq.0)then
!!        start_lon_idx = int((minlon_Topo-180.0_dp)/dlon_topo)
!!      else
!!        start_lon_idx = int((minlon_Topo+180.0_dp)/dlon_topo)
!!      endif
!
!      do ilon=1,nlon_topo_subgrid
!        lon_topo_subgrid(ilon) = real(start_lon_idx+ilon-1,kind=ip)*dlon_topo &
!                         - 0.5_dp*dlon_topo - 180.0_dp
!      enddo
!      ! now shift the lon values if they start in the western hemisphere
!
!      if(lon_shift_flag.eq.0)then
!        end_lon_idx = nlon_topo_subgrid
!      else
!        end_lon_idx = int((180.0_dp-minlon_Topo)/dlon_topo)
!      endif
!
!      ! GEBCO_08 start at the NW corner of the grid
!      ! (89d59'45"N,179d59;45"W) and advances eastward, then
!      ! southwards in one long string.
!      allocate(dum1d_short(nlon_topo_fullgrid))
!      do ilat=1,nlat_topo_subgrid
!        ! Get the index point of the start of the line at the right
!        ! latitude
!        lat_topo_subgrid(ilat)=(real(nlat_topo_fullgrid-(start_lat_idx-ilat+1),kind=ip)*dlat_topo &
!                       -0.5_dp*dlat_topo)-90.0_dp
!        idx = (start_lat_idx-ilat)*nlon_topo_fullgrid +1






      stop 4

      ! Error traps (starting with 9000)
      ! For this subroutine, the 100's position refers to block # of control file

9001  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot open header file: ',file_topo_hdr
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

      end subroutine Load_Topo_Gridded_bin

!##############################################################################
!
!    Interp_Topo
!
!##############################################################################

      subroutine Interp_Topo

      use global_param,  only : &
        EPS_SMALL

      use mesh,          only : &
         nxmax,nymax,nzmax,lon_cc_pd,lat_cc_pd,z_cc_pd,xy2ll_ylat,&
         xy2ll_xlon,IsLatLon

      implicit none

      integer :: i,j,k
      real(kind=ip) :: ophi,olam
      real(kind=ip) :: a1,a2,a3,a4
      real(kind=ip) :: xc,yc,xfrac,yfrac
      integer       :: ilon,ilat

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Inside interp_topo"
      endif;enddo

      ! Loop over all the computational grid points
!      do i=0,nxmax+1
!        do j=0,nymax+1
      do i=1,nxmax
        do j=1,nymax

          if(IsLatLon)then
            olam = lon_cc_pd(i)
            ophi = lat_cc_pd(j)
          else
            ! Recover the lat / lon calculated earlier for this point
            olam = xy2ll_xlon(i,j)
            ophi = xy2ll_ylat(i,j)
          endif
!          if(lon_shift_flag.eq.0)then
            if(olam.gt. 180.0_ip.and.&
               lon_topo_subgrid(nlon_topo_subgrid).lt.180.0_ip)olam=olam-360.0_ip
            if(olam.lt.-180.0_ip)olam=olam+360.0_ip
!          endif

          ! Now find the corresponding topo point
!          write(*,*)i,j,olam,lon_topo_subgrid(1),dlon_topo,&
!                    olam-lon_topo_subgrid(1),ophi,lat_topo_subgrid(1)
          ilon = floor((olam-lon_topo_subgrid(1))/dlon_topo) + 1
          ilat = floor((ophi-lat_topo_subgrid(1))/dlat_topo) + 1
          if(olam-lon_topo_subgrid(ilon).lt.0.0_ip)then
            ilon=ilon-1
          elseif(olam-lon_topo_subgrid(ilon).ge.dlon_topo)then
            ilon=ilon+1
          endif
          if(ophi-lat_topo_subgrid(ilat).lt.0.0_ip)then
            ilat=ilat-1
          elseif(ophi-lat_topo_subgrid(ilat).ge.dlat_topo)then
            ilat=ilat+1
          endif

          ! No interp; just lower-left corner
          topo_comp(i,j) = topo_subgrid(ilon,ilat)

          ! Bilinear interpolation
          xfrac=(olam-lon_topo_subgrid(ilon))/(lon_topo_subgrid(ilon+1)-lon_topo_subgrid(ilon))
          yfrac=(ophi-lat_topo_subgrid(ilat))/(lat_topo_subgrid(ilat+1)-lat_topo_subgrid(ilat))
          xc = 1.0_ip-xfrac
          yc = 1.0_ip-yfrac
          if(xc.gt.1.0_ip+EPS_SMALL.or.xc.lt.0.0_ip-EPS_SMALL)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)'lon ',i,olam,lon_topo_subgrid(ilon),xc
            endif;enddo
          endif
          if(yc.gt.1.0_ip+EPS_SMALL.or.yc.lt.0.0_ip-EPS_SMALL)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)'lat ',j,ophi,lat_topo_subgrid(ilat),yc
            endif;enddo
          endif

          a1=xc*yc
          a2=xfrac*yc
          a3=xfrac*yfrac
          a4=yfrac*xc

          topo_comp(i,j) = a1*topo_subgrid(ilon  ,ilat  ) + &
                           a2*topo_subgrid(ilon+1,ilat  ) + &
                           a3*topo_subgrid(ilon+1,ilat+1) + &
                           a4*topo_subgrid(ilon  ,ilat+1)
          topo_comp(i,j) = topo_comp(i,j) / 1000.0_sp ! convert to km

          topo_indx(i,j) = 0
          do k=1,nzmax+1
            if (z_cc_pd(k).le.topo_comp(i,j) .and. &
                z_cc_pd(k+1).gt.topo_comp(i,j))then
              topo_indx(i,j) = k
              exit
            endif
          enddo
        enddo
      enddo
      deallocate(lon_topo_subgrid,lat_topo_subgrid,topo_subgrid)

      if(.not.IsLatLon) then
        deallocate(xy2ll_ylat,xy2ll_xlon)

      endif

      end subroutine Interp_Topo

!##############################################################################

      subroutine Smooth_Topo

      use global_param,  only : &
         RAD_EARTH,DEG2RAD,DEG2RAD,PI

      use mesh,          only : &
         nxmax,nymax,nzmax,IsLatLon,dx,dy,de,dn,lat_cc_pd,lon_cc_pd,&
         x_cc_pd,y_cc_pd,z_cc_pd,xLL,yLL,lonLL,latLL

      use MetReader,       only : &
         MR_minlen,x_submet_sp,y_submet_sp,nx_submet,ny_submet,MR_dx_met,MR_dy_met,&
         MR_Topo_comp,MR_Topo_met

      implicit none

      integer :: i,j,k,it
      integer :: iidx,jidx
      integer :: ncells
      real(kind=ip) :: topo_avg,dist,cell_len
      real(kind=ip) :: deltheta,x1,x2,y1,y2,z1,z2
      real(kind=ip) :: rad_smooth_ang
      real(kind=ip) :: Const1_2d     = 0.6820926132509800_ip
      real(kind=ip) :: fac_1,r,temp1,wg,char_len,norm
      real(kind=ip) :: rad
      !real(kind=sp) :: topo_smooth_comp(0:nxmax+1,0:nymax+1)
      real(kind=sp) :: topo_smooth_comp(nxmax,nymax)
      real(kind=sp) :: topo_smooth_met(nx_submet,ny_submet)
      real(kind=sp) :: xin,yin

      integer :: nx,ny
      real(kind=sp),dimension(:,:),allocatable :: topo_smooth
      integer :: ii,ipad,jj,jpad

      ! We smooth the topographic data for both the computational grid
      ! as well as the met grid since it probably doesn't make sense to use
      ! relatively coarse met data mapped onto more highly resolved topographic
      ! data.
      do it=1,2
        if(it.eq.1)then
          ! First, smooth topography on the computational grid, but using the cell size
          ! of the met grid as the smoothing length (or half the width as the radius).
          ! This will be mapped to MR_Topo_met by selecting the corresponding cells
          ! of the smoothed computational grid.
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Smoothing topography on met grid"
          endif;enddo
          nx = nx_submet
          ny = ny_submet
          rad = 0.25_ip/1000.0_ip*real(MR_minlen,kind=ip)  ! Smooth over the width of the met cells
          cell_len = min(minval(MR_dx_met(:)),minval(MR_dy_met(:)))*DEG2RAD*RAD_EARTH
        else
          ! Second time, smooth topography on the computational grid using the smoothing
          ! radius provided in the input file
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Smoothing topography on comp grid"
          endif;enddo
          nx = nxmax
          ny = nymax
          rad = rad_smooth  ! Smooth over the radius specified in the input file
          if(IsLatLon)then
            cell_len = min(de,dn)*DEG2RAD*RAD_EARTH
          else
            cell_len = min(dx,dy)
          endif
          ! Smoothing radius should smooth over something like the scale of the met grid
          ! Issue a warning if it is much smaller
          if(rad.lt.0.25_ip/1000.0_ip*real(MR_minlen,kind=ip))then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)" WARNING: Width of smoothing kernel (4x radius) is less than"
              write(outlog(io),*)"          the grid size of Met data."
              write(outlog(io),*)"    Smoothing radius (km)         = ",real(rad,kind=4)
              write(outlog(io),*)"    Shortest met grid length (km) = ",real(MR_minlen/1000.0_ip,kind=4)
            endif;enddo
          elseif(rad.gt.cell_len*12.5_ip)then
            ! Conversely, issue a warning if it is too big since this will choke up the smoothing kernel
            ! by having too many topo points
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)" WARNING: Width of smoothing kernel (4x radius) greater than"
              write(outlog(io),*)"          50x the grid size of computational data."
              write(outlog(io),*)"    Smoothing radius (km)          = ",real(rad,kind=4)
              write(outlog(io),*)"    Kernel width (km)              = ",real(4.0_ip*rad,kind=4)
              write(outlog(io),*)"    Shortest comp grid length (km) = ",real(cell_len,kind=4)
            endif;enddo
          endif
        endif
        ! Initialize array for smoothed data with original topo data
        topo_smooth_comp = real(topo_comp,kind=sp)

        if(rad.le.cell_len)then
          useSmoothTopo = .false.
          cycle
        else
          useSmoothTopo = .true.
        endif

        ! On first loop, we always smooth to met grid
        ! Second time, exit if smoothing radius is too small
        if(it.eq.2.and..not.useSmoothTopo)exit

        if(IsLatLon)then
          rad_smooth_ang = rad/RAD_EARTH/DEG2RAD
          ipad = floor(rad_smooth_ang/de)
          jpad = floor(rad_smooth_ang/dn)
        else
          ipad = floor(rad/dx)
          jpad = floor(rad/dy)
        endif

!        do i=0,nxmax+1
!          do j=0,nymax+1
        do i=1,nxmax
          do j=1,nymax

            topo_avg = 0.0_ip
            ncells = 0
            norm = 0.0_ip
            if(IsLatLon)then
              ! Get distance by converting to cartesian on unit sphere
              x1=sin(0.5_ip*PI-lat_cc_pd(j)*DEG2RAD)*cos(lon_cc_pd(i)*DEG2RAD)
              y1=sin(0.5_ip*PI-lat_cc_pd(j)*DEG2RAD)*sin(lon_cc_pd(i)*DEG2RAD)
              z1=cos(0.5_ip*PI-lat_cc_pd(j)*DEG2RAD)
            endif

!            do ii=max(0,i-ipad),min(i+ipad,nxmax+1)
!              do jj=max(0,j-jpad),min(j+jpad,nymax+1)
            do ii=max(1,i-ipad),min(i+ipad,nxmax)
              do jj=max(1,j-jpad),min(j+jpad,nymax)
                if (ii.eq.i.and.jj.eq.j)then
                  dist = 0.0_ip
                else
                  if(IsLatLon)then
                    ! Get distance by converting to cartesian on unit sphere
                    x2=sin(0.5_ip*PI-lat_cc_pd(jj)*DEG2RAD)*cos(lon_cc_pd(ii)*DEG2RAD)
                    y2=sin(0.5_ip*PI-lat_cc_pd(jj)*DEG2RAD)*sin(lon_cc_pd(ii)*DEG2RAD)
                    z2=cos(0.5_ip*PI-lat_cc_pd(jj)*DEG2RAD)
                    ! Position vectors were normalized so angle is just acos
                    ! of inner product
                    deltheta = acos(x1*x2 + y1*y2 + z1*z2)
                    dist = deltheta*RAD_EARTH
                  else
                    dist=sqrt((x_cc_pd(ii)-x_cc_pd(i))**2.0_ip + &
                              (y_cc_pd(jj)-y_cc_pd(j))**2.0_ip)
                  endif
                endif

                if(dist.le.rad)then
                  ncells = ncells + 1
                    ! Here is cubic spline weighting                
                  char_len = 0.5_ip*rad
                  r = dist/char_len
                  fac_1 = Const1_2d/(char_len*char_len)
                  if (r.gt.2.0_ip) then
                    wg = 0.0_ip
                  else if (r.gt.1.0_ip) then
                    temp1 = 2.0_ip-r
                    wg = fac_1*((temp1*temp1*temp1)/6.0_ip)
                  else
                    wg = fac_1*(2.0_ip/3.0_ip - r*r + 0.5_ip*r*r*r)
                  endif
                  topo_avg = topo_avg + wg*topo_comp(ii,jj)
                  norm = norm + wg

                    ! This line can be used to check that the integration
                    ! over the spline is accurate (i.e. should
                    ! reconstructe a topography of 1.0 if the weights are
                    ! properly calculated)
                  !topo_avg = topo_avg + wg * 1.0_ip
                endif
              enddo ! loop over jj
            enddo ! loop over ii
            topo_smooth_comp(i,j) = real(topo_avg/norm,kind=sp)
            ! Assume anything lower than 0.0 is bathymetry; reset to 0.0
            if (topo_smooth_comp(i,j).lt.0.0_sp) topo_smooth_comp(i,j) = 0.0_sp
            ! Reset the topo index
            ! This is only needed if we are letting the grid intersect topography
            topo_indx(i,j) = 0
            do k=1,nzmax+1
              if (z_cc_pd(k).le.topo_smooth_comp(i,j) .and. &
                  z_cc_pd(k+1).gt.topo_smooth_comp(i,j))then
                topo_indx(i,j) = k
                exit
              endif
            enddo
          enddo
        enddo
      enddo
      ! Copy the smoothed data onto the main topo array
!      topo_comp(0:nxmax+1,0:nymax+1) = topo_smooth_comp(0:nxmax+1,0:nymax+1)
      topo_comp(1:nxmax,1:nymax) = topo_smooth_comp(1:nxmax,1:nymax)
      ! And copy to the MetReader array
      MR_Topo_comp(1:nxmax,1:nymax) = topo_smooth_comp(1:nxmax,1:nymax)

      ! Now populate the topo array on the met grid
      do i = 1,nx_submet
        xin = x_submet_sp(i)
        iidx = floor((xin-lonLL)/de) + 1
!        iidx = max(iidx,0)
!        iidx = min(iidx,nxmax+1)
        iidx = max(iidx,1)
        iidx = min(iidx,nxmax)
        do j = 1,ny_submet
          yin = y_submet_sp(j)
          jidx = floor((yin-latLL)/dn) + 1
!          jidx = max(jidx,0)
!          jidx = min(jidx,nymax+1)
          jidx = max(jidx,1)
          jidx = min(jidx,nymax)
          MR_Topo_met(i,j) = topo_comp(iidx,jidx)
          ! Assume anything lower than 0.0 is bathymetry; reset to 0.0
          if (MR_Topo_met(i,j).lt.0.0_sp) MR_Topo_met(i,j) = 0.0_sp
        enddo
      enddo

      return

      end subroutine Smooth_Topo

      end module Topography

