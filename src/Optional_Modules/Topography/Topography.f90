      module Topography

!##############################################################################
!     
!     In principle, netcdf would not be needed for including
!     topography, however, we assume/require that this library is
!     included
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
      ! modules output vars corespond to
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

      real(kind=ip),dimension(:,:)  ,allocatable :: topo_grid !Used if useTopo=.true.
      integer      ,dimension(:,:)  ,allocatable :: topo_indx !kindex of topo
      integer,private :: lon_shift_flag

      real(kind=ip) :: minlon_Topo,maxlon_Topo
      real(kind=ip) :: minlat_Topo,maxlat_Topo

      contains

!******************************************************************************

      subroutine input_data_Topo

      use global_param,  only : &
         DEG2RAD,RAD_EARTH,nmods

      use mesh,          only : &
         dx,dy,de,dn,IsLatLon

      use io_data,       only : &
         infile

      implicit none

      character(len=3)  :: answer
      character(len=80)  :: linebuffer
      integer :: ios,ioerr
      character(len=20) :: mod_name
      integer :: substr_pos

      open(unit=10,file=infile,status='old',err=1900)

      write(global_info,*)"    Searching for OPTMOD=TOPO"
      nmods = 0
      read(10,'(a80)',iostat=ios)linebuffer
      do while(ios.eq.0)
        read(10,'(a80)',iostat=ios)linebuffer

        substr_pos = index(linebuffer,'OPTMOD')
        if(substr_pos.eq.1)then
          ! found an optional module
          !  Parse for the keyword
          read(linebuffer,1104)mod_name
          if(adjustl(trim(mod_name)).eq.'TOPO')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo

      useTopo = .false.
      write(global_info,*)&
        "    Continue reading input file for topo block"
       ! Check if we're going to use topography
        read(10,'(a80)',iostat=ios,err=2010)linebuffer

        read(linebuffer,'(a3)',err=2011) answer
        if (answer.eq.'yes') then
          useTopo = .true.
          write(global_info,*)"    Using topography"
        elseif(answer(1:2).eq.'no') then
          useTopo = .false.
          write(global_info,*)"    Not using topography"
        else
          go to 2011
        endif
          
        if (useTopo) then
          ! Check if we're using topography, then get the format code
          read(10,'(a80)',iostat=ios,err=2010)linebuffer
          read(linebuffer,*,iostat=ioerr) topoFormat,rad_smooth
          if(topoFormat.eq.1)&
            write(global_info,*)"Read topoFormat = 1 (ETOPO1)"
          if(topoFormat.eq.2)&
            write(global_info,*)"Read topoFormat = 2 (GEBCO08)"
          write(global_info,*)"    Read smoothing radius = ",&
                               real(rad_smooth,kind=4)
          if(IsLatLon)then
            if(rad_smooth.le.min(de,dn)*DEG2RAD*RAD_EARTH)then
              useSmoothTopo = .false.
            else
              useSmoothTopo = .true.
            endif
          else
            if(rad_smooth.le.min(dx,dy))then
              useSmoothTopo = .false.
            else
              useSmoothTopo = .true.
            endif
          endif
          ! And read the file name
          read(10,'(a80)',iostat=ios,err=2010)linebuffer
          read(linebuffer,*) file_topo
          write(global_info,*)"    Read file_topo = ",file_topo
        endif

2010  continue
      close(10)

      return

1900  write(global_info,*)  'Error: cannot find input file: ',infile
      write(global_info,*)  'Program stopped'
      write(global_log,*)  'error: cannot find input file: ',infile
      write(global_log,*)  'Program stopped'
      stop 1

2011  write(global_log,*) 'Error reading whether to use topography.'
      write(global_log,*) 'Answer must be yes or no.'
      write(global_log,*) 'You gave:',linebuffer
      write(global_log,*) 'Program stopped'
      stop 1


      end subroutine input_data_Topo

!******************************************************************************

      subroutine Allocate_Topo(nx,ny)

      use io_data,       only : &
         nvar_User3d_XYZ,nvar_User3d_XYGs,nvar_User2d_XY,nvar_User2d_static_XY,&
         nvar_User4d_XYZGs

      implicit none

      integer :: nx,ny
      integer :: ngridnode

      write(global_info,*)"--------------------------------------------------"
      write(global_info,*)"---------- ALLOCATE_TOPO -------------------------"
      write(global_info,*)"--------------------------------------------------"
      ngridnode = (nx+3)*(ny+3)

      allocate(topo_grid(0:nx+2,0:ny+2));       topo_grid = 0.0_ip
      allocate(topo_indx(0:nx+2,0:ny+2));       topo_indx = 0

      ! Set the start indecies
      indx_User2d_static_XY_Topo = nvar_User2d_static_XY
      indx_User2d_XY_Topo        = nvar_User2d_XY
      indx_User3d_XYGs_Topo      = nvar_User3d_XYGs
      indx_User3d_XYZ_Topo       = nvar_User3d_XYZ
      indx_User4d_XYZGs_Topo     = nvar_User4d_XYZGs

      temp_2ds_name_Topo(1) = "Topography"
      temp_2ds_lname_Topo(1) = "Elevation of surface"
      temp_2ds_unit_Topo(1) = "km"
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
           real(topo_grid(1:nxmax,1:nymax),kind=op)
      enddo

      end subroutine Prep_output_Topo

!******************************************************************************

      subroutine Deallocate_Topo

      implicit none

      deallocate(topo_grid)
      deallocate(topo_indx)

      end subroutine Deallocate_Topo

!##############################################################################
!
!    get_topo
!
!##############################################################################

      subroutine get_topo

      use global_param,  only : &
         VERB

      use mesh,          only : &
         IsLatLon,lon_cc_pd,lat_cc_pd

      implicit none

      ! First we need to get the extents of the computational grid
      if(IsLatLon)then
        !Just get min and max of lat and lon.
        ! These were already calculated in calc_grid under the names
        ! latmin,latmax,lonmin,lonmax
        minlat_Topo = minval(lat_cc_pd)
        maxlat_Topo = maxval(lat_cc_pd)
        minlon_Topo = minval(lon_cc_pd)
        maxlon_Topo = maxval(lon_cc_pd)
      else
        ! This function is in Calc_Mesh
        call get_minmax_lonlat(minlon_Topo,maxlon_Topo,minlat_Topo,maxlat_Topo)
      endif

      if(VERB.gt.1)write(global_info,*)"min max lon = ",minlon_Topo,maxlon_Topo
      if(VERB.gt.1)write(global_info,*)"min max lat = ",minlat_Topo,maxlat_Topo

      call load_topo
      call interp_topo

      if(useSmoothTopo)then
        write(global_info,*)"Smoothing topography"
        call SmoothTopo
      else
        write(global_info,*)"Topography will not be smoothed."
      endif

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

      write(global_info,*)"Inside load_topo"
        ! Check to see if the domain straddles the anti-meridian
      !if (minlon.lt.180.0_ip.and.maxlon.gt.180.0_ip)then
      !  lon_shift_flag = 1
      !else
        lon_shift_flag = 0
      !endif

      nSTAT = nf90_open(adjustl(trim(file_topo)),NF90_NOWRITE,ncid)
      if(nSTAT.ne.0)write(global_info,*)'ERROR: nf90_open to read header:', &
                           nf90_strerror(nSTAT)
      if(nSTAT.ne.0)then
        write(global_info,*)'Could not open -',adjustl(trim(file_topo)),'-'
        write(global_info,*)NF90_NOWRITE,ncid
        write(global_info,*)'Exiting'
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

      nlon = int((maxlon_Topo-minlon_Topo)/dlon) + 4
      nlat = int((maxlat_Topo-minlat_Topo)/dlat) + 4
      allocate(lon_raw(nlon))
      allocate(lat_raw(nlat))
      allocate(topo_raw(nlon,nlat))
      if(lon_shift_flag.eq.0)then
        start_lon_idx = int((minlon_Topo-180.0_ip)/dlon)
      else
        start_lon_idx = int((minlon_Topo+180.0_ip)/dlon)
      endif
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
        write(global_info,*)"topoFormat must equal 1 or 2."
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

      write(global_info,*)"Inside interp_topo"

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
            if(olam.gt.180.0_ip)olam=olam-360.0_ip
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
          topo_grid(i,j) = topo_raw(ilon,ilat)
          ! Bilinear interpolation

          xfrac=(olam-lon_raw(ilon))/(lon_raw(ilon+1)-lon_raw(ilon))
          yfrac=(ophi-lat_raw(ilat))/(lat_raw(ilat+1)-lat_raw(ilat))
          xc = 1.0_ip-xfrac
          yc = 1.0_ip-yfrac
          if(xc.gt.1.0_ip+EPS_SMALL.or.xc.lt.0.0_ip-EPS_SMALL)then
            write(global_info,*)'lon ',i,olam,lon_raw(ilon),xc
          endif
          if(yc.gt.1.0_ip+EPS_SMALL.or.yc.lt.0.0_ip-EPS_SMALL)then
            write(global_info,*)'lat ',j,ophi,lat_raw(ilat),yc
          endif

          a1=xc*yc
          a2=xfrac*yc
          a3=xfrac*yfrac
          a4=yfrac*xc

          topo_grid(i,j) = a1*topo_raw(ilon  ,ilat  ) + &
                           a2*topo_raw(ilon+1,ilat  ) + &
                           a3*topo_raw(ilon+1,ilat+1) + &
                           a4*topo_raw(ilon  ,ilat+1)
          topo_grid(i,j) = topo_grid(i,j) / 1000.0_ip ! convert to km

          topo_indx(i,j) = 0
          do k=1,nzmax+1
            if (z_cc_pd(k).le.topo_grid(i,j) .and. &
                z_cc_pd(k+1).gt.topo_grid(i,j))then
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
         x_cc_pd,y_cc_pd,z_cc_pd

      implicit none

      integer :: i,j,k
      integer :: ncells
      real(kind=ip) :: topo_avg,dist
      real(kind=ip) :: deltheta,x1,x2,y1,y2,z1,z2
      real(kind=ip) :: rad_smooth_ang
      real(kind=ip) :: Const1_2d     = 0.6820926132509800_ip
      real(kind=ip) :: fac_1,r,temp1,wg,char_len,norm

      real(kind=ip) :: topo_grid_smooth(0:nxmax+2,0:nymax+2)

      integer :: ii,ipad,jj,jpad

      topo_grid_smooth = topo_grid

      if(IsLatLon)then
        rad_smooth_ang = rad_smooth/RAD_EARTH/DEG2RAD
        ipad = floor(rad_smooth_ang/de)
        jpad = floor(rad_smooth_ang/dn)
      else
        ipad = floor(rad_smooth/dx)
        jpad = floor(rad_smooth/dy)
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

              !if(dist.le.rad_smooth)then
              !    ! Here is simple averaging (no weighting)
              !  wg = 1.0_ip
              !  topo_avg = topo_avg + wg * topo_grid(ii,jj)
              !  norm = norm + 1.0_ip
              !endif

              if(dist.le.rad_smooth)then
                ncells = ncells + 1
                  ! Here is cubic spline weighting                
                char_len = 0.5_ip*rad_smooth
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
                topo_avg = topo_avg + wg*topo_grid(ii,jj)
                norm = norm + wg

                  ! This line can be used to check that the integration
                  ! over the spline is accurate (i.e. should
                  ! reconstructe a topography of 1.0 if the weights are
                  ! properly calculated)
                !topo_avg = topo_avg + wg * 1.0_ip
              endif
            enddo ! loop over jj
          enddo ! loop over ii
          topo_grid_smooth(i,j) = topo_avg / norm
          ! Reset the topo index
          topo_indx(i,j) = 0
          do k=1,nzmax+1
            if (z_cc_pd(k).le.topo_grid_smooth(i,j) .and. &
                z_cc_pd(k+1).gt.topo_grid_smooth(i,j))then
              topo_indx(i,j) = k
              exit
            endif
          enddo
        enddo
      enddo

      topo_grid = topo_grid_smooth

      return

      end subroutine SmoothTopo

      END module Topography

