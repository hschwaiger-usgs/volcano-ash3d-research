!*******************************************************************************
!# Topography
!*******************************************************************************
!OPTMOD=TOPO
!yes                           # use topography?
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
      integer             :: topoFormat
      real(kind=ip)       :: rad_smooth
      logical :: useTopo         = .false.
      logical :: useSmoothTopo   = .false.

      real(kind=ip),dimension(:)       ,allocatable :: lat
      real(kind=ip),dimension(:)       ,allocatable :: lon
      integer :: nlat
      integer :: nlon

      real(kind=ip), dimension(:)    ,allocatable :: lat_raw
      real(kind=ip), dimension(:)    ,allocatable :: lon_raw
      real(kind=ip), dimension(:,:)  ,allocatable :: topo_raw

      real(kind=ip) :: dlat,dlon

      real(kind=ip),dimension(:,:)  ,allocatable :: topo_comp ! Used if useTopo=.true.
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
         MR_useTopo

      implicit none

      character(len=3)  :: answer
      integer           :: dum_int
      character(len=80) :: linebuffer080
      integer           :: ios,ioerr
      character(len=20) :: mod_name
      integer           :: substr_pos

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
            write(outlog(io),*)"Read topoFormat = 1 (ETOPO1)"
            write(outlog(io),*)"    Read smoothing radius = ",&
                               real(rad_smooth,kind=4)
          endif;enddo
        elseif(topoFormat.eq.2)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Read topoFormat = 2 (GEBCO08)"
            write(outlog(io),*)"    Read smoothing radius = ",&
                               real(rad_smooth,kind=4)
          endif;enddo
        endif
        ! And read the file name
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        read(linebuffer080,*) file_topo
        do io=1,2;if(VB(io).le.verbosity_info)then           
          write(outlog(io),*)"    Read file_topo = ",file_topo
        endif;enddo

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

      allocate(topo_comp(0:nx+1,0:ny+1));       topo_comp = 0.0_ip
      allocate(topo_indx(0:nx+1,0:ny+1));       topo_indx = 0

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
!    get_topo
!
!##############################################################################

      subroutine get_topo

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
      if(IsLatLon)then
        !Just get min and max of lat and lon.
        ! These were already calculated in calc_grid under the names
        ! latmin,latmax,lonmin,lonmax
        minlat_Topo = min(minval(lat_cc_pd(-1:nymax+1)),real(minval(y_submet_sp(1:ny_submet)),kind=ip))
        maxlat_Topo = max(maxval(lat_cc_pd(-1:nymax+1)),real(maxval(y_submet_sp(1:ny_submet)),kind=ip))
        minlon_Topo = min(minval(lon_cc_pd(-1:nxmax+1)),real(minval(x_submet_sp(1:nx_submet)),kind=ip))
        maxlon_Topo = max(maxval(lon_cc_pd(-1:nxmax+1)),real(maxval(x_submet_sp(1:nx_submet)),kind=ip))
      else
        ! This function is in Calc_Mesh
        write(*,*)"HFS Check this."
        call get_minmax_lonlat(minlon_Topo,maxlon_Topo,minlat_Topo,maxlat_Topo)
      endif

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"min max lon = ",real(minlon_Topo,kind=4),real(maxlon_Topo,kind=4)
        write(outlog(io),*)"min max lat = ",real(minlat_Topo,kind=4),real(maxlat_Topo,kind=4)
      endif;enddo

      call load_topo

      call interp_topo

      call SmoothTopo

      end subroutine get_topo

!##############################################################################
!
!    load_topo
!
!##############################################################################

      subroutine load_topo

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

      if(topoFormat.eq.1)then
        ! ETOPO1 (1-minute global topo/batho)
        nlat_tot = 10800 ! -89.99167          to  89.99167
        nlon_tot = 21600 ! -179.991666666667  to 179.991666666667
        dlon = 1.0_ip/60.0_ip
        dlat = 1.0_ip/60.0_ip
        start_lat_idx = int((minlat_Topo+90.0_ip)/dlon)
      elseif(topoFormat.eq.2)then
        ! GEBCO 08 (30-second global topo/batho)
        nlat_tot = 21600
        nlon_tot = 43200
        dlon = 1.0_ip/120.0_ip
        dlat = 1.0_ip/120.0_ip
        start_lat_idx = int((minlat_Topo+90.0_ip)/dlon)
        start_lat_idx = nlat_tot-start_lat_idx
      endif
      ! Define the sub-grid holding the topo data we need
      nlon = int((maxlon_Topo-minlon_Topo)/dlon) + 4
      nlat = int((maxlat_Topo-minlat_Topo)/dlat) + 4
      allocate(lon_raw(nlon))
      allocate(lat_raw(nlat))
      allocate(topo_raw(nlon,nlat))
      ! ETOPO1 and GEBCO08 have data that start at -180, so shift 0 accordingly
      if(minlon_Topo.lt.0.0_ip)then
        start_lon_idx = int((minlon_Topo+180.0_ip)/dlon)
      elseif(minlon_Topo.ge.0.0_ip.and.minlon_Topo.lt.180.0_ip)then
        start_lon_idx = int((minlon_Topo+180.0_ip)/dlon)
      elseif(minlon_Topo.ge.180.0_ip)then
        start_lon_idx = int((minlon_Topo-180.0_ip)/dlon)
      endif

!      if(lon_shift_flag.eq.0)then
!        start_lon_idx = int((minlon_Topo-180.0_ip)/dlon)
!      else
!        start_lon_idx = int((minlon_Topo+180.0_ip)/dlon)
!      endif

      do ilon=1,nlon
        lon_raw(ilon) = real(start_lon_idx+ilon-1,kind=ip)*dlon &
                         - 0.5_ip*dlon - 180.0_ip
      enddo
      ! now shift the lon values if they start in the western hemisphere

      if(lon_shift_flag.eq.0)then
        end_lon_idx = nlon
      else
        end_lon_idx = int((180.0_ip-minlon_Topo)/dlon)
      endif

      if(topoFormat.eq.1)then
        ! Since this file is about 450Mb, just load the whole thing
        allocate(dum2d_short(nlon_tot,nlat_tot))
        do ilat=1,nlat
          lat_raw(ilat) = real(start_lat_idx+ilat-1,kind=ip)*dlat &
                           - 0.5_ip*dlat - 90.0_ip
        enddo
          ! load the whole array
        nSTAT = nf90_get_var(ncid,topo_var_id,dum2d_short)
        do ilat=1,nlat
          ! Get the index point of the start of the line
          idx = start_lat_idx+ilat - 1
          if(lon_shift_flag.eq.0)then
            topo_raw(:,ilat) = real(dum2d_short(start_lon_idx:start_lon_idx+nlon-1,idx),kind=ip)
          else
            ! Copy the part west of the anti-meridian
            topo_raw(1:end_lon_idx,ilat) = real(dum2d_short(start_lon_idx:start_lon_idx+end_lon_idx-1,idx),kind=ip)
            ! And copy the part wrapped over the anti-meridian
            topo_raw(end_lon_idx+1:nlon,ilat) = real(dum2d_short(1:nlon-end_lon_idx,idx),kind=ip)
          endif
        enddo

        deallocate(dum2d_short)
      elseif(topoFormat.eq.2)then
        ! GEBCO_08 start at the NW corner of the grid
        ! (89d59'45"N,179d59;45"W) and advances eastward, then
        ! southwards in one long string.
        allocate(dum1d_short(nlon_tot))
        do ilat=1,nlat
          ! Get the index point of the start of the line at the right
          ! latitude
          lat_raw(ilat)=(real(nlat_tot-(start_lat_idx-ilat+1),kind=ip)*dlat &
                         -0.5_ip*dlat)-90.0_ip
          idx = (start_lat_idx-ilat)*nlon_tot +1
          ! Get the full row of topo values that encircle the globe
          nSTAT = nf90_get_var(ncid,topo_var_id,dum1d_short, &
                             start = (/idx/),        &
                             count = (/nlon_tot/))
          ! Now find just the subset that we need
          if(lon_shift_flag.eq.0)then
            topo_raw(1:nlon,ilat)=real(dum1d_short(start_lon_idx:start_lon_idx+nlon-1),kind=ip)
          else
            ! Copy the part west of the anti-meridian
            topo_raw(1:end_lon_idx,ilat) = real(dum1d_short(start_lon_idx:start_lon_idx+end_lon_idx-1),kind=ip)
            ! And copy the part wrapped over the anti-meridian
            topo_raw(end_lon_idx+1:nlon,ilat) = real(dum1d_short(1:nlon-end_lon_idx),kind=ip)
          endif
        enddo
        deallocate(dum1d_short)
      else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"topoFormat must equal 1 or 2."
        endif;enddo
        stop 1
      endif
      nSTAT = nf90_close(ncid)

      end subroutine load_topo

!##############################################################################
!
!    interp_topo
!
!##############################################################################

      subroutine interp_topo

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
      do i=0,nxmax+1
        do j=0,nymax+1

          if(IsLatLon)then
            olam = lon_cc_pd(i)
            ophi = lat_cc_pd(j)
          else
            ! Recover the lat / lon calculated earlier for this point
            olam = xy2ll_xlon(i,j)
            ophi = xy2ll_ylat(i,j)
          endif
          if(lon_shift_flag.eq.0)then
            if(olam.gt. 180.0_ip)olam=olam-360.0_ip
            if(olam.lt.-180.0_ip)olam=olam+360.0_ip
          endif

          ! Now find the corresponding topo point
          ilon = floor((olam-lon_raw(1))/dlon) + 1
          ilat = floor((ophi-lat_raw(1))/dlat) + 1
          if(olam-lon_raw(ilon).lt.0.0_ip)then
            ilon=ilon-1
          elseif(olam-lon_raw(ilon).ge.dlon)then
            ilon=ilon+1
          endif
          if(ophi-lat_raw(ilat).lt.0.0_ip)then
            ilat=ilat-1
          elseif(ophi-lat_raw(ilat).ge.dlat)then
            ilat=ilat+1
          endif

          ! No interp; just lower-left corner
          topo_comp(i,j) = topo_raw(ilon,ilat)

          ! Bilinear interpolation
          xfrac=(olam-lon_raw(ilon))/(lon_raw(ilon+1)-lon_raw(ilon))
          yfrac=(ophi-lat_raw(ilat))/(lat_raw(ilat+1)-lat_raw(ilat))
          xc = 1.0_ip-xfrac
          yc = 1.0_ip-yfrac
          if(xc.gt.1.0_ip+EPS_SMALL.or.xc.lt.0.0_ip-EPS_SMALL)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)'lon ',i,olam,lon_raw(ilon),xc
            endif;enddo
          endif
          if(yc.gt.1.0_ip+EPS_SMALL.or.yc.lt.0.0_ip-EPS_SMALL)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)'lat ',j,ophi,lat_raw(ilat),yc
            endif;enddo
          endif

          a1=xc*yc
          a2=xfrac*yc
          a3=xfrac*yfrac
          a4=yfrac*xc

          topo_comp(i,j) = a1*topo_raw(ilon  ,ilat  ) + &
                           a2*topo_raw(ilon+1,ilat  ) + &
                           a3*topo_raw(ilon+1,ilat+1) + &
                           a4*topo_raw(ilon  ,ilat+1)
          topo_comp(i,j) = topo_comp(i,j) / 1000.0_ip ! convert to km

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
      deallocate(lon_raw,lat_raw,topo_raw)

      if(.not.IsLatLon) then
        deallocate(xy2ll_ylat,xy2ll_xlon)

      endif

      end subroutine interp_topo

!##############################################################################

      subroutine SmoothTopo

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
      real(kind=sp) :: topo_smooth_comp(0:nxmax+1,0:nymax+1)
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

        do i=0,nxmax+1
          do j=0,nymax+1
            topo_avg = 0.0_ip
            ncells = 0
            norm = 0.0_ip
            if(IsLatLon)then
              ! Get distance by converting to cartesian on unit sphere
              x1=sin(0.5_ip*PI-lat_cc_pd(j)*DEG2RAD)*cos(lon_cc_pd(i)*DEG2RAD)
              y1=sin(0.5_ip*PI-lat_cc_pd(j)*DEG2RAD)*sin(lon_cc_pd(i)*DEG2RAD)
              z1=cos(0.5_ip*PI-lat_cc_pd(j)*DEG2RAD)
            endif

            do ii=max(0,i-ipad),min(i+ipad,nxmax+1)
              do jj=max(0,j-jpad),min(j+jpad,nymax+1)
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
            ! Reset the topo index
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
      topo_comp(0:nxmax+1,0:nymax+1) = topo_smooth_comp(0:nxmax+1,0:nymax+1)
      ! And copy to the MetReader array
      MR_Topo_comp(1:nxmax,1:nymax) = topo_smooth_comp(1:nxmax,1:nymax)

      ! Now populate the topo array on the met grid
      do i = 1,nx_submet
        xin = x_submet_sp(i)
        iidx = floor((xin-lonLL)/de) + 1
        iidx = max(iidx,0)
        iidx = min(iidx,nxmax+1)
        do j = 1,ny_submet
          yin = y_submet_sp(j)
          jidx = floor((yin-latLL)/dn) + 1
          jidx = max(jidx,0)
          jidx = min(jidx,nymax+1)
          MR_Topo_met(i,j) = topo_comp(iidx,jidx)
          ! Assume anything lower than 0.0 is bathymetry; reset to 0.0
          if (MR_Topo_met(i,j).lt.0.0_sp) MR_Topo_met(i,j) = 0.0_sp
        enddo
      enddo

      return

      end subroutine SmoothTopo

      END module Topography

