      module land_cover

!##############################################################################
!     
!http://glcf.umd.edu/data/landcover/
! wget ftp://ftp.glcf.umd.edu/glcf/Global_Land_Cover/Global/1deg/gl-latlong-1deg-landcover.bsq.gz
! wget ftp://ftp.glcf.umd.edu/glcf/Global_Land_Cover/Global/8km/gl-goodes-8km-landcover.bsq.gz
! wget ftp://ftp.glcf.umd.edu/glcf/Global_Land_Cover/Global/1km/gl-goodes-1km-landcover.bsq.gz
!
!##############################################################################

      use precis_param

      use io_units

      ! Set the number of output variables for this module
      integer, parameter :: nvar_User2d_static_XY_LC = 1 ! land use
      integer, parameter :: nvar_User2d_XY_LC        = 0
      integer, parameter :: nvar_User3d_XYGs_LC      = 0
      integer, parameter :: nvar_User3d_XYZ_LC       = 0
      integer, parameter :: nvar_User4d_XYZGs_LC     = 0

      character(len=30),dimension(nvar_User2d_static_XY_LC) :: temp_2ds_name_LC
      character(len=30),dimension(nvar_User2d_static_XY_LC) :: temp_2ds_unit_LC
      character(len=30),dimension(nvar_User2d_static_XY_LC) :: temp_2ds_lname_LC
      real(kind=op),    dimension(nvar_User2d_static_XY_LC) :: temp_2ds_MissVal_LC
      real(kind=op),    dimension(nvar_User2d_static_XY_LC) :: temp_2ds_FillVal_LC

      ! These are used to keep track of which index in the global list, this
      ! modules output vars corespond to
      integer :: indx_User2d_static_XY_LC
      integer :: indx_User2d_XY_LC
      integer :: indx_User3d_XYGs_LC
      integer :: indx_User3d_XYZ_LC
      integer :: indx_User4d_XYZGs_LC

      integer :: nlat_LC
      integer :: nlon_LC

      real(kind=ip), dimension(:)    ,allocatable :: lat_raw_LC
      real(kind=ip), dimension(:)    ,allocatable :: lon_raw_LC

      integer(kind=1), dimension(:,:)  ,allocatable :: LC_raw

      real(kind=ip) :: dlat,dlon
      real(kind=ip) :: LC_size_km
      real(kind=ip) :: LC_size_deg
      integer,private :: lon_shift_flag

      real(kind=ip) :: minlon_LC,maxlon_LC
      real(kind=ip) :: minlat_LC,maxlat_LC

      logical :: useLandCover    = .false.
      character (len=50)  :: file_LandCover
      character (len=50)  :: LC_dir
      integer             :: LandCover_Format

      integer(kind=1), dimension(:,:)  ,allocatable :: LC_grid

      contains
!##############################################################################

!##############################################################################

      subroutine input_data_LC

      use global_param,  only : &
         nmods

      use io_data,       only : &
         infile

      implicit none

      character(len=3)  :: answer
      character(len=80)  :: linebuffer
      integer :: ios !,ioerr
      character(len=20) :: mod_name
      integer :: substr_pos

      open(unit=10,file=infile,status='old',err=1900)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Searching for OPTMOD=LC"
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
          if(adjustl(trim(mod_name)).eq.'LC')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo

      useLandCover = .false.
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Continue reading input file for LandCover block"
      endif;enddo
       ! Check if we're going to use land classification data
        read(10,'(a80)',iostat=ios,err=2010)linebuffer

        read(linebuffer,'(a3)',err=2011) answer
        if (answer.eq.'yes') then
          useLandCover = .true.
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using Land Cover data"
          endif;enddo
        elseif(answer(1:2).eq.'no') then
          useLandCover = .false.
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Not using global Land Cover data."
          endif;enddo
        else
          go to 2011
        endif

        if (useLandCover) then
          read(10,'(a80)',iostat=ios,err=2010)linebuffer
          read(linebuffer,*) LandCover_Format
          read(10,'(a80)',iostat=ios,err=2010)linebuffer
          read(linebuffer,*) LC_dir
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    LandCover data located in : ",LC_dir
          endif;enddo
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
        write(errlog(io),*) 'Error reading whether to use land cover.'
        write(errlog(io),*) 'Answer must be yes or no.'
        write(errlog(io),*) 'You gave:',linebuffer
        write(errlog(io),*) 'Program stopped'
      endif;enddo
      stop 1


      end subroutine input_data_LC

!******************************************************************************


!******************************************************************************

      subroutine Allocate_LC

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User2d_XY,&
         nvar_User3d_XYGs,nvar_User3d_XYZ,    &
         nvar_User4d_XYZGs

      use mesh,          only : &
         nxmax,nymax

      implicit none

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOCATE_LC ---------------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      allocate(LC_grid(0:nxmax+2,0:nymax+2)); LC_grid = 0

      ! Set the start indecies
      indx_User2d_static_XY_LC = nvar_User2d_static_XY
      indx_User2d_XY_LC        = nvar_User2d_XY
      indx_User3d_XYGs_LC      = nvar_User3d_XYGs
      indx_User3d_XYZ_LC       = nvar_User3d_XYZ
      indx_User4d_XYZGs_LC     = nvar_User4d_XYZGs

      temp_2ds_name_LC(1)    = "LandCover"
      temp_2ds_lname_LC(1)   = "Land Cover catagory"
      temp_2ds_unit_LC(1)    = "km"
      temp_2ds_MissVal_LC(1) = -9999.0_op
      temp_2ds_FillVal_LC(1) = -9999.0_op
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"units","ID")
!        if(nSTAT.ne.0) &
!          write(global_log,*)'ERROR: put_att land use: ',nf90_strerror(nSTAT)
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"Source",&
!          "Global Land Cover Facility (1-km or 8-km)")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"Site",&
!          "http://glcf.umiacs.umd.edu/data/landcover/")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=0","Water")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=1","Evergreen Needleleaf")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=2","Evergreen Broadleaf")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=3","Deciduous Needleleaf")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=4","Deciduous Broadleaf")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=5","Mixed Forest")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=6","Woodland")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=7","Wooded Grassland")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=8","Closed Shrubland")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=9","Open Shrubland")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=10","Grassland")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=11","Cropland")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=12","Bare Ground")
!        nSTAT = nf90_put_att(ncid,landuse_var_id,"ID=13","Urban")

      nvar_User2d_static_XY = nvar_User2d_static_XY + nvar_User2d_static_XY_LC
      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_LC
      nvar_User3d_XYGs      = nvar_User3d_XYGs      + nvar_User3d_XYGs_LC
      nvar_User3d_XYZ       = nvar_User3d_XYZ       + nvar_User3d_XYZ_LC
      nvar_User4d_XYZGs     = nvar_User4d_XYZGs     + nvar_User4d_XYZGs_LC

      end subroutine Allocate_LC

!******************************************************************************

      subroutine Prep_output_LC

      use mesh,          only : &
         nxmax,nymax

      use Output_Vars, only : &
         var_User2d_static_XY_name,var_User2d_static_XY_unit,var_User2d_static_XY_lname,&
         var_User2d_static_XY_MissVal,var_User2d_static_XY_FillVal,var_User2d_static_XY

      implicit none

      integer :: i,indx

      do i=1,nvar_User2d_static_XY_LC
        indx = indx_User2d_static_XY_LC+i
        var_User2d_static_XY_name(indx) = temp_2ds_name_LC(i)
        var_User2d_static_XY_unit(indx) = temp_2ds_unit_LC(i)
        var_User2d_static_XY_lname(indx)= temp_2ds_lname_LC(i)
        var_User2d_static_XY_MissVal(indx)= temp_2ds_MissVal_LC(i)
        var_User2d_static_XY_FillVal(indx)= temp_2ds_FillVal_LC(i)
        if(i.eq.1)var_User2d_static_XY(1:nxmax,1:nymax,indx) = real(LC_grid(1:nxmax,1:nymax),kind=op)
      enddo

      end subroutine Prep_output_LC

!******************************************************************************

      subroutine Deallocate_LC

      implicit none

      if(allocated(LC_grid)) deallocate(LC_grid)

      end subroutine Deallocate_LC


!##############################################################################
!
!    load_LC
!
!##############################################################################

      subroutine load_LC

      use mesh,          only : &
         IsLatLon,lat_cc_pd,lon_cc_pd

      implicit none

      integer :: nlat_tot,nlon_tot
      integer :: start_lat_idx,start_lon_idx,end_lon_idx
      integer :: ilat,ilon,idx
      logical :: IsThere

      integer(kind=1),allocatable,dimension(:) :: glc_line
      character(len=130),dimension(3) :: LC_files
      integer :: open_status

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Inside load_LC"
      endif;enddo

      write(LC_files(1),116)trim(adjustl(LC_dir)),&
                               "/gl-latlong-1deg-landcover.bsq"
      write(LC_files(2),117)trim(adjustl(LC_dir)),&
                               "/gl-latlong-8km-landcover.bsq"
      write(LC_files(3),117)trim(adjustl(LC_dir)),&
                               "/gl-latlong-1km-landcover.bsq"
 116    format(a50,a30)
 117    format(a50,a29)

      !LC_files(1)="Landcover/gl-latlong-1deg-landcover.bsq"
      !LC_files(2)="Landcover/gl-latlong-8km-landcover.bsq"
      !LC_files(3)="Landcover/gl-latlong-1km-landcover.bsq"

      if(IsLatLon)then
        minlat_LC = minval(lat_cc_pd)
        maxlat_LC = maxval(lat_cc_pd)
        ! The LC files are on the domain -180->180 so map the values to this
        minlon_LC = minval(lon_cc_pd)
        maxlon_LC = maxval(lon_cc_pd)
        if(minlon_LC.gt.180.0_ip)then
          minlon_LC = minlon_LC - 360.0_ip
          maxlon_LC = maxlon_LC - 360.0_ip
        endif

        ! Find out which file we should use based on the de and dn
        !if(de.ge.0.5_ip.or.dn.ge.0.5_ip)then
        !  LandCover_Format = 1
        !elseif(de.ge.0.05_ip.or.dn.ge.0.05_ip)then
        !  LandCover_Format = 2
        !else
        !  LandCover_Format = 3
        !endif
      else
        call get_minmax_lonlat(minlon_LC,maxlon_LC,minlat_LC,maxlat_LC)
        if(minlon_LC.gt.180.0_ip)then
          minlon_LC = minlon_LC - 360.0_ip
          maxlon_LC = maxlon_LC - 360.0_ip
        endif

        ! Find out which file we should use based on the dx and dy
        !if(dx.gt.50.0_ip.or.dy.gt.50.0_ip)then
        !  LandCover_Format = 1
        !elseif(dx.gt.4.0_ip.or.dy.gt.4.0_ip)then
        !  LandCover_Format = 2
        !else
        !  LandCover_Format = 3
        !endif
      endif
      if(LandCover_Format.eq.1)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Using 1-deg landcover data"
        endif;enddo
      elseif(LandCover_Format.eq.2)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Using 8-km landcover data"
        endif;enddo
      elseif(LandCover_Format.eq.3)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Using 1-km landcover data"
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then

        endif;enddo
      endif

        ! Check to see if the domain straddles the anti-meridian
      if (minlon_LC.lt.180.0_ip.and.maxlon_LC.gt.180.0_ip)then
        lon_shift_flag = 1
      else
        lon_shift_flag = 0
      endif

      if(LandCover_Format.eq.1)then
        ! 1 degree resolution
        nlat_tot = 180 !  -90.0 to  90.0
        nlon_tot = 360 ! -180.0 to 180.0
        LC_size_km  = 111.19_ip ! 1 degree with RAD_EARTH=6371.229
        LC_size_deg = 1.0_ip
        !0       Water
        !1       Broadleaf Evergreen Forest
        !2       Coniferous Evergreen Forest and Woodland
        !3       High Latitude Deciduous Forest and Woodland
        !4       Tundra
        !5       Mixed Coniferous Forest and Woodland
        !6       Wooded Grassland
        !7       Grassland
        !8       Bare Ground
        !9       Shrubs and Bare Ground
        !10      Cultivated Crops
        !11      Broadleaf Deciduous Forest and Woodland
        !12      Data Unavailable
      elseif(LandCover_Format.eq.2)then
        ! 8 km resolution
        ! stored top down
        nlat_tot = 2880 !  -90.0  to  90.0
        nlon_tot = 5760 ! -180.0  to 180.0
        LC_size_km  = 8.0_ip
        LC_size_deg = 0.0625_ip
        !0       Water (and Goode's interrupted space)
        !1       Evergreen Needleleaf Forest
        !2       Evergreen Broadleaf Foreset
        !3       Deciduous Needleleaf Forest
        !4       Deciduous Broadleaf Forest
        !5       Mixed Forest
        !6       Woodland
        !7       Wooded Grassland
        !8       Closed Shrubland
        !9       Open Shrubland
        !10      Grassland
        !11      Cropland
        !12      Bare Ground
        !13      Urban and Built-up
      elseif(LandCover_Format.eq.3)then
        ! 1 km resolution
        ! stored top down
        nlat_tot = 21600 ! -90.0   to  90.0
        nlon_tot = 43200 ! -180.0  to 180.0
        LC_size_km  = 1.0_ip
        LC_size_deg = 8.3333e-3_ip
        !0       Water (and Goode's interrupted space)
        !1       Evergreen Needleleaf Forest
        !2       Evergreen Broadleaf Foreset
        !3       Deciduous Needleleaf Forest
        !4       Deciduous Broadleaf Forest
        !5       Mixed Forest
        !6       Woodland
        !7       Wooded Grassland
        !8       Closed Shrubland
        !9       Open Shrubland
        !10      Grassland
        !11      Cropland
        !12      Bare Ground
        !13      Urban and Built-up
      endif

      ! Before we try to open the file, check its existance
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Looking for LC file: ",adjustl(trim(LC_files(LandCover_Format)))
      endif;enddo
      inquire( file=adjustl(trim(LC_files(LandCover_Format))), exist=IsThere )
      if(.not.IsThere)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"LC ERROR: Could not find LandCover file ",&
                      adjustl(trim(LC_files(LandCover_Format)))
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  Found it!"
        endif;enddo
      endif

      open(unit=30,file=trim(adjustl(LC_files(LandCover_Format))),access='direct',&
           recl=nlon_tot,iostat=open_status,status='old')
      if ( open_status /= 0 ) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'Could not open ',trim(adjustl(LC_files(LandCover_Format))),' for reading.'
        endif;enddo
        stop 1
      endif

      allocate(glc_line(nlon_tot))
      dlon = 360.0_ip/real(nlon_tot,kind=ip)
      dlat = 180.0_ip/real(nlat_tot,kind=ip)
      start_lat_idx = int((minlat_LC+90.0_ip)/dlat)
      start_lat_idx = nlat_tot-start_lat_idx
      nlon_LC = int((maxlon_LC-minlon_LC)/dlon) + 3
      nlat_LC = int((maxlat_LC-minlat_LC)/dlat) + 3
      allocate(lon_raw_LC(nlon_LC))
      allocate(lat_raw_LC(nlat_LC))
      allocate(LC_raw(nlon_LC,nlat_LC))
      !if(lon_shift_flag.eq.0)then
        start_lon_idx = int((minlon_LC+180.0_ip)/dlon)
      !else
      !  start_lon_idx = int((minlon_LC+180.0_ip)/dlon)
      !endif
      do ilon=1,nlon_LC
        lon_raw_LC(ilon) = real(start_lon_idx+ilon-1,kind=ip)*dlon &
                         - 0.5_ip*dlon - 180.0_ip
        if(lon_raw_LC(ilon).lt.0.0_ip) lon_raw_LC(ilon) = lon_raw_LC(ilon) + 360.0_ip
      enddo
      if(lon_shift_flag.eq.0)then
        end_lon_idx = nlon_LC
      else
        end_lon_idx = int((180.0_ip-minlon_LC)/dlon)
      endif

      ! data start at the NW corner of the grid
      ! (89d59'45"N,179d59;45"W) and advances eastward, then
      ! southwards in one long string.
      do ilat=1,nlat_LC
        ! Get the index point of the start of the line at the right
        ! latitude
        lat_raw_LC(ilat)=(real(nlat_tot-(start_lat_idx-ilat+1),kind=ip)*dlat &
                       -0.5_ip*dlat)-90.0_ip
        idx = start_lat_idx-ilat + 1
        ! Get the full row of LC values that encircle the globe
        ! Now find just the subset that we need
        read(30,rec=idx)glc_line
        if(lon_shift_flag.eq.0)then
          LC_raw(1:nlon_LC,ilat)=glc_line(start_lon_idx:start_lon_idx+nlon_LC-1)
        else
          ! Copy the part west of the anti-meridian
          LC_raw(1:end_lon_idx,ilat) = glc_line(start_lon_idx:start_lon_idx+end_lon_idx-1)
          ! And copy the part wrapped over the anti-meridian
          LC_raw(end_lon_idx+1:nlon_LC,ilat) = glc_line(1:nlon_LC-end_lon_idx)
        endif
      enddo

      end subroutine load_LC

!##############################################################################
!
!    assign_LC
!
!##############################################################################

      subroutine assign_LC

      use mesh,          only : &
         nxmax,nymax,xy2ll_ylat,xy2ll_xlon,lonLL,latLL,xLL,yLL, &
         A3d_iprojflag,A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
         A3d_k0_scale,A3d_Re,de,dn,dx,dy,IsLatLon,&
         lon_cc_pd,lat_cc_pd

      use projection,    only : &
         PJ_proj_for

      implicit none

      integer :: i,j
      integer       :: ilon,ilat
      integer       :: ix,iy,ilc

      integer :: LC_sum(1:nxmax,1:nymax,14)
      integer :: ic
      integer :: tot_sum,land_sum
      integer :: max_occur,noccur

      real(kind=ip) :: lon_LC,lat_LC
      real(kind=ip) :: xout,yout
      real(kind=ip) :: clon,clat
      integer       :: iidx,jidx

      ! The LC data are on a lon/lat grid with some resolution, but the
      ! computational grid could be lon/lat or x/y.  Moreover the computational
      ! cells might be bigger or smaller (or some mixture) than the LC cells.
      ! The strategy is to loop over all the LC cells that are in the
      ! computational domain and keep track of the number of LC cells that map
      ! to each comp cell.  Then we'll loop over all the comp cells and assign a
      ! LC value based on the most frequent mapped value.  If no LC cell maps to
      ! a particular comp cell, then the comp cell is mapped to the LC grid and
      ! the corresponding value assigned.

      ! First the LC to comp mapping
      LC_sum(1:nxmax,1:nymax,:) = 0
      do ilon=1,nlon_LC
        lon_LC = lon_raw_LC(ilon)
        if(lon_LC.lt.0.0)lon_LC=lon_LC+360.0_ip
        do ilat=1,nlat_LC
          lat_LC = lat_raw_LC(ilat)
          ilc   = LC_raw(ilon,ilat)
          ! Get the index of this LC point on the comp grid
          if(IsLatLon)then
            ! LC and comp grids are both LL: easy-peasy
            ix = int((lon_LC-lonLL)/de) + 1
            iy = int((lat_LC-latLL)/dn) + 1
          else
            call PJ_proj_for(lon_LC,lat_LC, A3d_iprojflag, &
                       A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,A3d_k0_scale,A3d_Re, &
                       xout,yout)
            ix = int((xout-xLL)/dx) + 1
            iy = int((yout-yLL)/dy) + 1
          endif
            ! make sure this ix.iy maps to the comp grid
          if(ix.lt.1.or.ix.gt.nxmax.or.iy.lt.1.or.iy.gt.nymax)then
            cycle
          !else
          endif
            ! Increment counter for this cell and LC type
          LC_sum(ix,iy,ilc+1) = LC_sum(ix,iy,ilc+1) + 1
        enddo
      enddo

      ! Now loop through the comp points and assign the most frequent land class
      do i=1,nxmax
        do j=1,nymax
          ! See how many LC cells mapped to this comp cell
          !  If just one, assign the corresponding value
          !  If more than one, pick the most frequent
          !  If none, then map the comp cell back to LC grid and use that
          tot_sum  = sum(LC_sum(i,j,:))
          land_sum = sum(LC_sum(i,j,2:14))
          if(tot_sum.eq.1)then
            ! only one LC cell, find which one and use its value
            do ic = 1,14
              if(LC_sum(i,j,ic).eq.1)then
                LC_grid(i,j) = ic-1  ! We decrement here because the LC classes
                                     ! are from 0-13
              else
                cycle
              endif
            enddo
          elseif(tot_sum.eq.0)then
            ! No LC cells mapped, Map this comp cell back to LC grid
            ! Get lon/lat of this cell
            if(IsLatLon)then
              clon = lon_cc_pd(i)
              clat = lat_cc_pd(j)
            else
              clon = xy2ll_xlon(i,j)
              clat = xy2ll_ylat(i,j)
            endif
            iidx = (clon-lon_raw_LC(1))/dlon +1
            jidx = (clat-lat_raw_LC(1))/dlat +1
            LC_grid(i,j) = LC_raw(iidx,jidx)
          elseif(tot_sum.gt.1)then
            ! Multiple LC cells mapped here, find which is most frequent
            ! First check land v.s. sea
            if(land_sum.gt.LC_sum(i,j,1))then
               ! this is a land cell, find which type of land cover
              max_occur = maxval(LC_sum(i,j,2:14))
              noccur = 0
              do ic = 2,14
                if(LC_sum(i,j,ic).eq.max_occur)then
                  noccur=noccur+1
                endif
              enddo
              if(noccur.eq.1)then
                ! only one land use class is most frequent, use it
                LC_grid(i,j) = maxloc(LC_sum(i,j,2:14),DIM=1) - 1
              else
                ! more than one land use class is most frequent
                ! for now, just pick the first (most vegitated), but we might
                ! want to do something more sophisticated here
                LC_grid(i,j) = maxloc(LC_sum(i,j,2:14),DIM=1) - 1
              endif
            else
              ! this is water
              LC_grid(i,j) = 0
            endif
          else
            ! Should not be here
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"Sum is less than 0; shouldn't be here."
            endif;enddo
            stop 1
          endif
        enddo
      enddo

      if(allocated(lon_raw_LC)) deallocate(lon_raw_LC)
      if(allocated(lat_raw_LC)) deallocate(lat_raw_LC)
      if(allocated(LC_raw))     deallocate(LC_raw)
      if(allocated(xy2ll_ylat)) deallocate(xy2ll_ylat)
      if(allocated(xy2ll_xlon)) deallocate(xy2ll_xlon)
 
      end subroutine assign_LC

!##############################################################################

      end module land_cover

