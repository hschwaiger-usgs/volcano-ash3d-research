      module ocean_currents

      use precis_param

      use io_units

      ! Set the number of output variables for this module
      integer, parameter :: nvar_User2d_static_XY_OSCAR = 0
      integer, parameter :: nvar_User2d_XY_OSCAR        = 2
      integer, parameter :: nvar_User3d_XYGs_OSCAR      = 0
      integer, parameter :: nvar_User3d_XYZ_OSCAR       = 0
      integer, parameter :: nvar_User4d_XYZGs_OSCAR     = 0

      character(len=30),dimension(nvar_User2d_XY_OSCAR) :: temp_2d_name_OSCAR
      character(len=30),dimension(nvar_User2d_XY_OSCAR) :: temp_2d_unit_OSCAR
      character(len=30),dimension(nvar_User2d_XY_OSCAR) :: temp_2d_lname_OSCAR
      real(kind=op),    dimension(nvar_User2d_XY_OSCAR) :: temp_2d_MissVal_OSCAR
      real(kind=op),    dimension(nvar_User2d_XY_OSCAR) :: temp_2d_FillVal_OSCAR

      ! These are used to keep track of which index in the global list, this
      ! modules output vars corespond to
      integer :: indx_User2d_static_XY_OSCAR
      integer :: indx_User2d_XY_OSCAR
      integer :: indx_User3d_XYGs_OSCAR
      integer :: indx_User3d_XYZ_OSCAR
      integer :: indx_User4d_XYZGs_OSCAR

      character (len=50)  :: oscar_dir
      logical :: useOceanCurrent = .false.

      integer, parameter :: fluc_l = 1
      integer, parameter :: fluc_r = 2

      real(kind=ip)     ,dimension(:), allocatable :: oscar_Hour
      character(len=130),dimension(:), allocatable :: oscar_files
      integer           ,dimension(:), allocatable :: oscar_year
      integer           ,dimension(:), allocatable :: oscar_tslice

      real(kind=sp),allocatable, dimension(:,:,:,:),private :: temp_var_sp
      real(kind=sp),allocatable, dimension(:),private :: time1_var_sp
      real(kind=ip),allocatable, dimension(:),private :: time1hours_var
      real(kind=sp),allocatable, dimension(:),private :: lat1_var_sp
      real(kind=sp),allocatable, dimension(:),private :: lon1_var_sp

      real(kind=sp),allocatable, dimension(:)   :: osc_x_sp,osc_y_sp
      real(kind=sp),allocatable, dimension(:,:) :: osc_u_sp,osc_v_sp
      integer,private :: end_timestep

      integer :: OSCAR_toggle                            ! toggle used for moving pointers between last and next

      real(kind=sp),allocatable, dimension(:,:),target :: Vx_surf_1_sp,Vy_surf_1_sp
      real(kind=sp),allocatable, dimension(:,:),target :: Vx_surf_2_sp,Vy_surf_2_sp
      real(kind=sp),dimension(:,:)  ,pointer           :: Vx_surf_last_step => null()
      real(kind=sp),dimension(:,:)  ,pointer           :: Vx_surf_next_step => null()
      real(kind=sp),dimension(:,:)  ,pointer           :: Vy_surf_last_step => null()
      real(kind=sp),dimension(:,:)  ,pointer           :: Vy_surf_next_step => null()

      integer :: istepNow
      integer :: istepNext

      character(len=130) :: oscar_filename

      integer :: nSTAT
      integer :: ncid_oscar
      integer :: osc_ntmax,osc_nxmax,osc_nymax
      integer :: osc_t_dim_id
      integer :: osc_x_dim_id
      integer :: osc_y_dim_id
      integer :: osc_t_var_id
      integer :: osc_x_var_id
      integer :: osc_y_var_id
      integer :: osc_u_var_id
      integer :: osc_v_var_id

      integer :: osc_nx,osc_ny

      integer :: slon_indx,elon_indx
      integer :: slat_indx,elat_indx

      contains

!******************************************************************************

      subroutine input_data_OSCAR

      use global_param,  only : &
         nmods

      use io_data,       only : &
         infile

      use time_data,     only : &
         BaseYear,useLeap,SimStartHour

      implicit none

      integer            :: HS_YearOfEvent
      character(len=3)   :: answer
      character(len=80)  :: linebuffer
      integer :: ios
      character(len=20)  :: mod_name
      integer :: substr_pos
      integer :: iyear

      iyear = HS_YearOfEvent(SimStartHour,BaseYear,useLeap)

      open(unit=10,file=infile,status='old',err=1900)

      write(global_info,*)"    Searching for OPTMOD=OSCAR"
      nmods = 0
      read(10,'(a80)',iostat=ios)linebuffer
      do while(ios.eq.0)
        read(10,'(a80)',iostat=ios)linebuffer

        substr_pos = index(linebuffer,'OPTMOD')
        if(substr_pos.eq.1)then
          ! found an optional module
          !  Parse for the keyword
          read(linebuffer,1104)mod_name
          if(adjustl(trim(mod_name)).eq.'OSCAR')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo

      useOceanCurrent = .false.
      write(global_info,*)"    Continue reading input file for OSCAR block"
        !write(global_info,*)linebuffer
       ! Check if we're going to use OceanCurrents
        read(10,'(a80)',iostat=ios,err=2010)linebuffer

        read(linebuffer,'(a3)',err=2011) answer
        if (answer.eq.'yes') then
          useOceanCurrent = .true.
          write(global_info,*)"    Using ocean currents"
        elseif(answer(1:2).eq.'no') then
          useOceanCurrent = .false.
          write(global_info,*)"    Not using ocean currents"
        else
          go to 2011
        endif

        if (useOceanCurrent) then
          ! We're going to use OSCAR Sea Surface currents
          !   https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_third-deg
          !If(iyear.lt.1992.or.&
          !   (iyear.eq.1992.and.imonth(1).lt.10))then
          if(iyear.lt.1993)then
            useOceanCurrent = .false.
            write(global_info,*)"Use of OSCAR currents requested, but the eruption occurs"
            write(global_info,*)"before the data are availible."
          else
            useOceanCurrent = .true.
            write(global_info,*)"Planning to use OSCAR current data"
          endif

            ! And read the file name
          read(10,'(a80)',iostat=ios,err=2010)linebuffer
          read(linebuffer,*) oscar_dir
          write(global_info,*)"OSCAR data located in : ",oscar_dir

        endif

2010  continue
      close(10)

      return

1900  write(global_info,*)  'error: cannot find input file: ',infile
      write(global_info,*)  'Program stopped'
      write(global_log,*)  'error: cannot find input file: ',infile
      write(global_log,*)  'Program stopped'
      stop 1


2011  write(global_log,*) 'Error reading whether to use ocean currents.'
      write(global_log,*) 'Answer must be yes or no.'
      write(global_log,*) 'You gave:',linebuffer
      write(global_log,*) 'Program stopped'
      stop 1


      end subroutine input_data_OSCAR

!******************************************************************************

      subroutine Allocate_OSCAR

      use mesh,          only : &
         nxmax,nymax

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User2d_XY,nvar_User3d_XYGs,nvar_User3d_XYZ,&
         nvar_User4d_XYZGs

      implicit none

      write(global_info,*)"--------------------------------------------------"
      write(global_info,*)"---------- ALLOCATE_OSCAR ------------------------"
      write(global_info,*)"--------------------------------------------------"

      allocate(Vx_surf_1_sp(nxmax,nymax)); Vx_surf_1_sp = 0.0_sp
      allocate(Vy_surf_1_sp(nxmax,nymax)); Vy_surf_1_sp = 0.0_sp
      allocate(Vx_surf_2_sp(nxmax,nymax)); Vx_surf_2_sp = 0.0_sp
      allocate(Vy_surf_2_sp(nxmax,nymax)); Vy_surf_2_sp = 0.0_sp


      ! Set the start indecies
      indx_User2d_static_XY_OSCAR = nvar_User2d_static_XY
      indx_User2d_XY_OSCAR        = nvar_User2d_XY
      indx_User3d_XYGs_OSCAR      = nvar_User3d_XYGs
      indx_User3d_XYZ_OSCAR       = nvar_User3d_XYZ
      indx_User4d_XYZGs_OSCAR     = nvar_User4d_XYZGs

      temp_2d_name_OSCAR(1)    = "U-surface-vel"
      temp_2d_lname_OSCAR(1)   = "U ocean current"
      temp_2d_unit_OSCAR(1)    = "m/s"
      temp_2d_MissVal_OSCAR(1) = -9999.0_op
      temp_2d_FillVal_OSCAR(1) = -9999.0_op

      temp_2d_name_OSCAR(2)    = "V-surface-vel"
      temp_2d_lname_OSCAR(2)   = "V ocean current"
      temp_2d_unit_OSCAR(2)    = "m/s"
      temp_2d_MissVal_OSCAR(2) = -9999.0_op
      temp_2d_FillVal_OSCAR(2) = -9999.0_op

      nvar_User2d_static_XY = nvar_User2d_static_XY + nvar_User2d_static_XY_OSCAR
      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_OSCAR
      nvar_User3d_XYGs      = nvar_User3d_XYGs      + nvar_User3d_XYGs_OSCAR
      nvar_User3d_XYZ       = nvar_User3d_XYZ       + nvar_User3d_XYZ_OSCAR
      nvar_User4d_XYZGs     = nvar_User4d_XYZGs     + nvar_User4d_XYZGs_OSCAR

      end subroutine Allocate_OSCAR

!******************************************************************************

      subroutine Prep_output_OSCAR

      use mesh,          only : &
         nxmax,nymax

      use Output_Vars,   only : &
         var_User2d_XY_name,var_User2d_XY_unit,var_User2d_XY_lname,&
         var_User2d_XY_MissVal,var_User2d_XY_FillVal,&
         var_User2d_XY

      use Solution,      only : &
         vx_pd,vy_pd

      implicit none

      integer :: i,indx

      do i=1,nvar_User2d_XY_OSCAR
        indx = indx_User2d_XY_OSCAR+i
        var_User2d_XY_name(indx) = temp_2d_name_OSCAR(i)
        var_User2d_XY_unit(indx) = temp_2d_unit_OSCAR(i)
        var_User2d_XY_lname(indx)= temp_2d_lname_OSCAR(i)
        var_User2d_XY_MissVal(indx)= temp_2d_MissVal_OSCAR(i)
        var_User2d_XY_FillVal(indx)= temp_2d_FillVal_OSCAR(i)
        if(i.eq.1)var_User2d_XY(1:nxmax,1:nymax,indx) = vx_pd(1:nxmax,1:nymax,0)
        if(i.eq.2)var_User2d_XY(1:nxmax,1:nymax,indx) = vy_pd(1:nxmax,1:nymax,0)
      enddo


      end subroutine Prep_output_OSCAR

!******************************************************************************

      subroutine Deallocate_OSCAR

      implicit none

      deallocate(Vx_surf_1_sp)
      deallocate(Vy_surf_1_sp)
      deallocate(Vx_surf_2_sp)
      deallocate(Vy_surf_2_sp)

      end subroutine Deallocate_OSCAR

!******************************************************************************

      subroutine Check_SurfaceVelocity

      use mesh,          only : &
         lon_cc_pd,lat_cc_pd,IsLatLon

      use time_data,     only : &
         BaseYear,useLeap,Simtime_in_hours,SimStartHour

      use netcdf

      implicit none

      integer :: sindex,eindex
      integer :: ti,i,j
      real(kind=dp)   :: HS_hours_since_baseyear
      real(kind=ip)  :: minlon_oscar,maxlon_oscar,minlat_oscar,maxlat_oscar
      integer        :: HS_DayOfYear
      integer        :: HS_YearOfEvent
      integer        :: iyear

      OSCAR_toggle = 1

 110  format(a50,a1)
      write(oscar_filename,110)trim(adjustl(oscar_dir)),'/'
      iyear = HS_YearOfEvent(SimStartHour,BaseYear,useLeap)
      !iday  = DayOfYear(SimStartHour)
      !oscar_tslice(1) = int((iday-1) * 0.2) -1

      ! open first file
      write(oscar_filename,116)trim(adjustl(oscar_dir)),&
                               "/oscar_",iyear,".nc"
 116  format(a50,a7,i4,a3)
      oscar_filename=trim(adjustl(oscar_filename))
      write(global_info,*)"Opening ",oscar_filename
      nSTAT = nf90_open(oscar_filename,NF90_NOWRITE,ncid_oscar)
      if(nSTAT.ne.0) write(global_info,*)'ERROR open oscar file: ',&
                     oscar_filename,nf90_strerror(nSTAT)
      nSTAT = nf90_inq_dimid(ncid_oscar,"time",osc_t_dim_id)
      nSTAT = nf90_Inquire_Dimension(ncid_oscar,osc_t_dim_id,len=osc_ntmax)
      nSTAT = nf90_inq_dimid(ncid_oscar,"lat",osc_y_dim_id)
      nSTAT = nf90_Inquire_Dimension(ncid_oscar,osc_y_dim_id,len=osc_nymax)
      nSTAT = nf90_inq_dimid(ncid_oscar,"lon",osc_x_dim_id)
      nSTAT = nf90_Inquire_Dimension(ncid_oscar,osc_x_dim_id,len=osc_nxmax)
      nSTAT = nf90_inq_varid(ncid_oscar,"time",osc_t_var_id)
      nSTAT = nf90_inq_varid(ncid_oscar,"lat",osc_y_var_id)
      nSTAT = nf90_inq_varid(ncid_oscar,"lon",osc_x_var_id)
      nSTAT = nf90_inq_varid(ncid_oscar,"u",osc_u_var_id)
      nSTAT = nf90_inq_varid(ncid_oscar,"v",osc_v_var_id)

      allocate(time1_var_sp(osc_ntmax))
      nSTAT = nf90_get_var(ncid_oscar,osc_t_var_id,time1_var_sp)
      allocate(lat1_var_sp(osc_nymax))
      nSTAT = nf90_get_var(ncid_oscar,osc_y_var_id,lat1_var_sp)
      allocate(lon1_var_sp(osc_nxmax))
      nSTAT = nf90_get_var(ncid_oscar,osc_x_var_id,lon1_var_sp)

      allocate(time1hours_var(osc_ntmax))
      time1hours_var = time1_var_sp*24.0_sp + &
              HS_hours_since_baseyear(iyear,1,1,0.0_dp,BaseYear,useLeap)

      ! Get the braketing indicies in time
      sindex = 0
      eindex = 0
      do ti = 1,osc_ntmax-1
        if(time1hours_var(ti)  .le.SimStartHour.and.&
           time1hours_var(ti+1).gt.SimStartHour)then
          sindex = ti
        endif
        if(time1hours_var(ti)  .le.SimStartHour+Simtime_in_hours.and.&
           time1hours_var(ti+1).gt.SimStartHour+Simtime_in_hours)then
          eindex = ti+1
        endif
      enddo
      if(eindex.eq.0)then
        write(global_info,*)"Simulation crosses into the next year."
        write(global_info,*)"Need to code up year wrapping for OSCAR data."
        stop 1
      endif
      end_timestep = eindex-sindex+1
      allocate(oscar_Hour(end_timestep))
      oscar_Hour(1:end_timestep) = time1hours_var(sindex:eindex)
      allocate(oscar_files(end_timestep))
      oscar_files(1:end_timestep) = oscar_filename
      allocate(oscar_year(end_timestep))
      oscar_year(1:end_timestep) = iyear
      allocate(oscar_tslice(end_timestep))
      write(global_info,*)"We will need ",end_timestep," oscar files."
      do ti = 1,end_timestep
        write(global_info,*)ti,oscar_Hour(ti),&
          HS_YearOfEvent(oscar_Hour(ti),BaseYear,useLeap),&
          HS_DayOfYear(oscar_Hour(ti),BaseYear,useLeap)
        oscar_tslice(ti) = ti-1+sindex
      enddo

      if(IsLatLon)then
        !Just get min and max of lat and lon.
        ! These were already calculated in calc_grid under the names
        minlon_oscar = minval(lon_cc_pd)
        maxlon_oscar = maxval(lon_cc_pd)
        minlat_oscar = minval(lat_cc_pd)
        maxlat_oscar = maxval(lat_cc_pd)
      else
        call get_minmax_lonlat(minlon_oscar,maxlon_oscar,minlat_oscar,maxlat_oscar)
      endif

      ! Get the bracketing indicies in lon/lat
      slon_indx = 0
      elon_indx = 0
      do i = 1,osc_nxmax-1
        if(lon1_var_sp(i).le.minlon_oscar.and.lon1_var_sp(i+1).gt.minlon_oscar)then
          slon_indx = i
        endif
        if(lon1_var_sp(i).le.maxlon_oscar.and.lon1_var_sp(i+1).gt.maxlon_oscar)then
          elon_indx = i+1
        endif
      enddo

      slat_indx = 0
      elat_indx = 0

      if(maxlat_oscar.gt.60.0_ip)then
        write(global_info,*)"WARNING:  OSCAR data only goes to 60N"
        slat_indx = 1
      endif
      if(minlat_oscar.lt.-60.0_ip)then
        write(global_info,*)"WARNING:  OSCAR data only goes to -60N"
        elat_indx = osc_nymax
      endif

      do j = 1,osc_nymax-1
        ! we're reading from the top down
        if(lat1_var_sp(j).ge.maxlat_oscar.and.lat1_var_sp(j+1).lt.maxlat_oscar)then
          slat_indx = j
        endif
        if(lat1_var_sp(j+1).le.minlat_oscar.and.lat1_var_sp(j).gt.minlat_oscar)then
          elat_indx = j+1
        endif
      enddo

      osc_nx = elon_indx - slon_indx + 1
      osc_ny = elat_indx - slat_indx + 1
      allocate(osc_x_sp(osc_nx))
      osc_x_sp = lon1_var_sp(slon_indx:elon_indx)
      allocate(osc_y_sp(osc_ny))
      do j = 1,osc_ny
        osc_y_sp(j)=lat1_var_sp(elat_indx+1-j)
      enddo

      nSTAT = nf90_close(ncid_oscar)

      istepNow  = 1
      istepNext = 2
      call Read_SurfaceVelocity(istepNow)
      Vx_surf_last_step => Vx_surf_1_sp
      Vy_surf_last_step => Vy_surf_1_sp

      call Read_SurfaceVelocity(istepNext)
      Vx_surf_next_step => Vx_surf_2_sp
      Vy_surf_next_step => Vy_surf_2_sp


      end subroutine Check_SurfaceVelocity

!******************************************************************************

      subroutine Read_SurfaceVelocity(istep)

      use global_param,  only : &
         VERB,EPS_SMALL

      use mesh,          only : &
         nxmax,nymax

      use netcdf

      implicit none

      integer, intent(in) :: istep

      integer :: i,j
      real(kind=sp),dimension(:,:)     ,allocatable :: tmp_regrid2d_sp
      real(kind=sp),allocatable ::  w2_sp(:)
      integer,allocatable ::  iw2(:)
      integer             ::  lw2
      integer             ::  liw2
      integer             ::  ier2
      integer             ::  intpol2(2)
      integer             ::  status


      write(global_info,*)"Reading next Oscar file."
      ! Set pointer toggle
      if(OSCAR_toggle.eq.1)then
        OSCAR_toggle = 2
      else
        OSCAR_toggle = 1
      endif

        ! Now get velocities for last time step
      allocate(temp_var_sp(osc_nx,osc_ny,1,1))
      allocate(osc_u_sp(osc_nx,osc_ny))
      allocate(osc_v_sp(osc_nx,osc_ny))

      nSTAT = nf90_open(oscar_filename,NF90_NOWRITE,ncid_oscar)

      temp_var_sp = 0.0_sp
      nSTAT = nf90_get_var(ncid_oscar,osc_u_var_id,temp_var_sp,&
                 start = (/slon_indx,slat_indx,1,oscar_tslice(istep)/),       &
                 count = (/osc_nx,osc_ny,1,1/))
      if(nSTAT.ne.0) write(global_info,*)'ERROR getting u: ',&
                     oscar_filename,nf90_strerror(nSTAT) !,&
                     !slon_indx,slat_indx,1,oscar_tslice(istep),&
                     !osc_nx,osc_ny,1,1

      do j = 1,osc_ny
        osc_u_sp(:,j) = temp_var_sp(:,osc_ny+1-j,1,1)
      enddo
      do i = 1,osc_nx
        do j = 1,osc_ny
          if(abs(osc_u_sp(i,j)).lt.EPS_SMALL)osc_u_sp(i,j) = 0.0_sp
        enddo
      enddo

      temp_var_sp = 0.0_sp
      nSTAT = nf90_get_var(ncid_oscar,osc_v_var_id,temp_var_sp,&
                 start = (/slon_indx,slat_indx,1,oscar_tslice(istep)/),       &
                 count = (/osc_nx,osc_ny,1,1/))
      if(nSTAT.ne.0) write(global_info,*)'ERROR getting v: ',&
                     oscar_filename,nf90_strerror(nSTAT) !,& 
                     !slon_indx,slat_indx,1,oscar_tslice(istep),&
                     !osc_nx,osc_ny,1,1
      do j = 1,osc_ny
        osc_v_sp(:,j) = temp_var_sp(:,osc_ny+1-j,1,1)
      enddo
      do i = 1,osc_nx
        do j = 1,osc_ny
          if(abs(osc_v_sp(i,j)).lt.EPS_SMALL)osc_v_sp(i,j) = 0.0_sp
        enddo
      enddo


      allocate(tmp_regrid2d_sp(nxmax,nymax))
      lw2 = nxmax+(nymax+2*nxmax)
      allocate(w2_sp(lw2),stat=status)
      liw2 = nxmax+nymax
      allocate(iw2(liw2),stat=status)
      intpol2 = 1
      ier2 = 0

      if(VERB.gt.2)write(global_info,*)"Interpolating osc_u "
      !call rgrd2_sp(osc_nx,   osc_ny,     &
      !              osc_x_sp,  osc_y_sp,    &
      !              osc_u_sp(:,:),         &
      !              nxmax,      nymax,        &
      !              lon_grid_sp, lat_grid_sp,     &
      !              tmp_regrid2d_sp(:,:),       &
      !              intpol2,w2_sp,lw2,  iw2,liw2,ier2)
      do i = 1,nxmax
        do j = 1,nymax
          if(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
        enddo
      enddo
      if(OSCAR_toggle.eq.1)then
        Vx_surf_1_sp = tmp_regrid2d_sp
      else
        Vx_surf_2_sp = tmp_regrid2d_sp
      endif

      if(VERB.gt.2)write(global_info,*)"Interpolating osc_v "
      !call rgrd2_sp(osc_nx,   osc_ny,     &
      !              osc_x_sp,  osc_y_sp,    &
      !              osc_v_sp(:,:),         &
      !              nxmax,      nymax,        &
      !              lon_grid_sp, lat_grid_sp,     &
      !              tmp_regrid2d_sp(:,:),       &
      !              intpol2,w2_sp,lw2,  iw2,liw2,ier2)
      do i = 1,nxmax
        do j = 1,nymax
          if(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
        enddo
      enddo
      if(OSCAR_toggle.eq.1)then
        Vy_surf_1_sp = tmp_regrid2d_sp
      else
        Vy_surf_2_sp = tmp_regrid2d_sp
      endif  

      nSTAT = nf90_close(ncid_oscar)

      deallocate(temp_var_sp)
      deallocate(osc_u_sp)
      deallocate(osc_v_sp)
      deallocate(tmp_regrid2d_sp)
      deallocate(w2_sp,iw2)

      end subroutine Read_SurfaceVelocity

!******************************************************************************

      subroutine set_SurfaceVelocity(TimeNow)

      use mesh,          only : &
         nxmax,nymax

      use time_data,     only : &
         SimStartHour

      use solution,      only : &
         vx_pd,vy_pd

      implicit none

      real(kind=ip), intent(in) :: TimeNow  ! current time, in hours since start of simulation

      real(kind=ip)  :: Interval_Frac
      real(kind=ip)  :: HoursIntoOSCARInterval
      integer        :: istep = 1
      
      if ((TimeNow+SimStartHour).gt.oscar_Hour(istepNext)) then
          istepNext = istepNext + 1
          istepNow  = istepNow  + 1
        call Read_SurfaceVelocity(istepNext)
        if(OSCAR_toggle.eq.1)then
          Vx_surf_last_step => Vx_surf_2_sp
          Vy_surf_last_step => Vy_surf_2_sp
          Vx_surf_next_step => Vx_surf_1_sp
          Vy_surf_next_step => Vy_surf_1_sp
        else
          Vx_surf_last_step => Vx_surf_1_sp
          Vy_surf_last_step => Vy_surf_1_sp
          Vx_surf_next_step => Vx_surf_2_sp
          Vy_surf_next_step => Vy_surf_2_sp
        endif
      endif

      HoursIntoOSCARInterval = TimeNow + SimStartHour -  oscar_Hour(istep)
      Interval_Frac = HoursIntoOSCARInterval / (oscar_Hour(istep+1)-oscar_Hour(istep))

      vx_pd(1:nxmax,1:nymax,0) = 3.6_ip*(real(Vx_surf_last_step(:,:),kind=ip) + &
                              real((Vx_surf_next_step(:,:) - &
                                     Vx_surf_last_step(:,:)),kind=ip) * &
                                     Interval_Frac)
      vy_pd(1:nxmax,1:nymax,0) = 3.6_ip*(real(Vy_surf_last_step(:,:),kind=ip) + &
                              real((Vy_surf_next_step(:,:) - &
                                     Vy_surf_last_step(:,:)),kind=ip) * &
                                     Interval_Frac)

      end subroutine set_SurfaceVelocity

!******************************************************************************

      subroutine advect_deposit(depo)

      use global_param,  only : &
         EPS_THRESH

      use mesh,          only : &
         IsLatLon,nxmax,nymax,nsmax,dz_vec_pd,sigma_nx_pd,sigma_ny_pd,kappa_pd,&
         dx,dy

      use solution,      only : &
         vx_pd,vy_pd,outflow_yz1_pd,outflow_yz2_pd,outflow_xz1_pd,outflow_xz2_pd

      use time_data,     only : &
         dt

      implicit none

      real(kind=ip),dimension(nxmax,nymax,nsmax),intent(inout) :: depo

      real(kind=ip),dimension(0:nxmax+2,0:nymax+2,1:nsmax) :: dum_depo

      integer :: i,j,n
      !integer :: idx_dum
      real(kind=ip) :: sig_p32,sig_p12,sig_m12    ! area of interface
      real(kind=ip) :: usig_p32,usig_p12,usig_m12 ! vel * sig at interface
      real(kind=ip) :: au
      real(kind=ip) :: dqi,ldq
      real(kind=ip) :: rp2,rp1,r0,rm1
      real(kind=ip) :: vol             ! volume of cell
      real(kind=ip) :: dt_vol
      real(kind=ip) :: order1,order2
      real(kind=ip) :: divu_p, divu_m
      real(kind=ip) :: bflux_w,bflux_e
      real(kind=ip) :: bflux_s,bflux_n
      real(kind=ip),dimension(0:nxmax+2,fluc_l:fluc_r) :: xfs
      real(kind=ip),dimension(0:nxmax+2)               :: xfss
      real(kind=ip),dimension(0:nymax+2,fluc_l:fluc_r) :: yfs
      real(kind=ip),dimension(0:nymax+2)               :: yfss

      real(kind=ip) :: dqu,theta
      real(kind=ip) :: dz

      dz = dz_vec_pd(0)

        ! Copy the input deposit array to a workspace
      dum_depo(0:nxmax+2,0:nymax+2,1:nsmax) = 0.0_ip
      dum_depo(1:nxmax  ,1:nymax  ,1:nsmax) = depo(1:nxmax,1:nymax,1:nsmax)
        ! The input array will also be used for the output, advected deposit so
        ! reinitialize it
      depo = 0.0_ip

      ! First advect in x
      if(.not.IsLatLon)then
        ! Cartesian geometry terms do not vary with position, so
        ! pre-calculate them.
        vol     = dx*dy*dz
        dt_vol  = dt/vol
        sig_p32 = dy*dz
        sig_p12 = sig_p32
        sig_m12 = sig_p32
      endif
      do n=1,nsmax
          do j=1,nymax
            xfs  = 0.0e0_ip
            xfss = 0.0e0_ip

            bflux_w = 0.0_ip
            bflux_e = 0.0_ip

            ! First loop over x and build first and second order fluxes
            ! at interface i (the inteface on the right side of cell i)
            do i=1,nxmax
              ! Get the area of each of the interfaces and other
              ! geometry terms
              if(IsLatLon)then
                vol     = kappa_pd(i,j,0)
                dt_vol  = dt/vol
                sig_p32 = sigma_nx_pd(i+1,j,0)
                sig_p12 = sigma_nx_pd(i  ,j,0)
                sig_m12 = sigma_nx_pd(i-1,j,0)
              endif
              ! Get area*velocities at interfaces i+3/2, i+1/2 and i-1/2
              usig_p32 = sig_p32*0.5e0_ip*(vx_pd(i+1,j,0)+vx_pd(i+2,j,0))
              usig_p12 = sig_p12*0.5e0_ip*(vx_pd(i  ,j,0)+vx_pd(i+1,j,0))
              usig_m12 = sig_m12*0.5e0_ip*(vx_pd(i-1,j,0)+vx_pd(i  ,j,0))

              rp2 = max(0.0_ip,dum_depo(i+2,j,n))
              rp1 = max(0.0_ip,dum_depo(i+1,j,n))
              r0  = max(0.0_ip,dum_depo(i  ,j,n))
              rm1 = max(0.0_ip,dum_depo(i-1,j,n))

                ! This cycling is good for production runs, but causes
                ! problems with convergence tests
#ifndef NOCYCLE
              if (rp1.le.EPS_THRESH.and.r0.le.EPS_THRESH) cycle
#endif
                ! delta Q at interface
              dqi  = rp1-r0
                ! Set flux based on upwind velocity for color equation
                !  (equals conservative form if div.v=0)
              xfs(i,fluc_r) = max(0.0_ip,usig_p12)*dqi ! flux OUT OF i,j cell to the cell to the right (pos x)
              xfs(i,fluc_l) = min(0.0_ip,usig_p12)*dqi ! flux INTO i,j cell from the cell to the right (pos x)
                ! Modification for conservative form in
                ! divergent/convergent velocities
              divu_p   = max(0.0_ip,usig_p12) - max(0.0_ip,usig_m12)
              divu_m   = min(0.0_ip,usig_p32) - min(0.0_ip,usig_p12)
              xfs(i,fluc_r) = xfs(i,fluc_r) + rp1 * divu_m
              xfs(i,fluc_l) = xfs(i,fluc_l) + r0  * divu_p

                ! Only calculate dqu and theta if a limiter is being used
                ! Get delta Q at upwind interface relative to interface i
              if(usig_p12.gt.0.0_ip)then
                dqu = r0-rm1  ! upwind is interface i-1
              else
                dqu = rp2-rp1 ! upwind is interface i+1
              endif

                ! Make sure that theta is not singular
              if(abs(dqu).gt.3.0_ip*abs(dqi).or.abs(dqi).lt.EPS_THRESH)then
                theta = sign(3.0_ip,dqu*dqi)
              else
                theta = dqu/dqi
              endif

              ldq = dqi*max(0.0_ip,min(1.0_ip,2.0_ip*theta),min(2.0_ip,theta))

              au = abs(usig_p12)
              xfss(i) = 0.5e0_ip*au*(1.0e0_ip-dt_vol*au)*ldq

                ! For boundary fluxes, we want the flux-difference form rather
                ! than the fluctuation form given above and just use first-order
                ! upwind
              if(i.eq.1)     bflux_w =-min(0.0_ip,usig_m12)*r0
              if(i.eq.nxmax) bflux_e = max(0.0_ip,usig_p12)*r0

            enddo  ! loop over i

            ! Now loop over x again and apply these fluxes to volume i
            do i=1,nxmax
              if(IsLatLon)then
                vol     = kappa_pd(i,j  ,0)
                dt_vol  = dt/vol
              endif

#ifndef NOCYCLE
              rp1 = dum_depo(i+1,j,n)
              r0  = dum_depo(i  ,j,n)

              if(i.gt.0)then
                !rm1 = concen(i-1,j,k,n,ts0)
                rm1 = dum_depo(i-1,j,n)
              else
                rm1 = 0.0_ip
              endif
              if(rm1+r0+rp1 .lt. EPS_THRESH) cycle
#endif
              ! Note: if i=1, xfs(0,:)=0.0 and xfss(0)=0.0
              order1 =  dt_vol * ( -xfs(i-1,fluc_r) -  xfs(i,fluc_l))    ! Eq 19.19 (line 1)
              order2 =  dt_vol * ( xfss(i-1)        - xfss(i)       )    ! Eq 19.19 (line 3)

              !concen(i,j,k,n,ts1) = concen(i,j,k,n,ts0) + order1 + order2
              depo(i,j,n) = dum_depo(i,j,n) + order1 + order2

              if(i.eq.1)then
                ! Flux out the west side
                outflow_yz1_pd(j,0,n) = outflow_yz1_pd(j,0,n) + bflux_w * dt_vol
              elseif(i.eq.nxmax)then
                ! Flux out the east side
                outflow_yz2_pd(j,0,n) = outflow_yz2_pd(j,0,n) + bflux_e * dt_vol
              endif

            enddo
          enddo ! loop over j=1,nymax
      enddo ! loop over idx_dum


        ! Move the x-advected deposit back to the workspace for y-advection
      dum_depo(0:nxmax+2,0:nymax+2,1:nsmax) = 0.0_ip
      dum_depo(1:nxmax  ,1:nymax  ,1:nsmax) = depo(1:nxmax,1:nymax,1:nsmax)
        ! The input array will also be used for the output, advected deposit so
        ! reinitialize it
      depo = 0.0_ip

      ! Now advect in y
      if(.not.IsLatLon)then
        ! Cartesian geometry terms do not vary with position, so
        ! pre-calculate them.
        vol     = dx*dy*dz
        dt_vol  = dt/vol
        sig_p32 = dx*dz
        sig_p12 = sig_p32
        sig_m12 = sig_p32
      endif

      do n=1,nsmax
        do i=1,nxmax
          yfs  = 0.0e0_ip
          yfss = 0.0e0_ip

          bflux_s = 0.0_ip
          bflux_n = 0.0_ip

          ! First loop over y and build first and second order fluxes
          ! at interface j
          do j=1,nymax
            ! Get the area of each of the interfaces and other
            ! geometry terms
            if(IsLatLon)then
              vol     = kappa_pd(i,j,0)
              dt_vol  = dt/vol
              sig_p32 = sigma_ny_pd(i,j+1,0)
              sig_p12 = sigma_ny_pd(i,j  ,0)
              sig_m12 = sigma_ny_pd(i,j-1,0)
            endif
            ! Get area*velocities at interfaces j+3/2, j+1/2 and j-1/2
            usig_p32 = sig_p32*0.5e0_ip*(vy_pd(i,j+1,0)+vy_pd(i,j+2,0))
            usig_p12 = sig_p12*0.5e0_ip*(vy_pd(i,j  ,0)+vy_pd(i,j+1,0))
            usig_m12 = sig_m12*0.5e0_ip*(vy_pd(i,j-1,0)+vy_pd(i,j  ,0))

            rp2 = max(0.0_ip,dum_depo(i,j+2,n))
            rp1 = max(0.0_ip,dum_depo(i,j+1,n))
            r0  = max(0.0_ip,dum_depo(i,j  ,n))
            rm1 = max(0.0_ip,dum_depo(i,j-1,n))

#ifndef NOCYCLE
            if (rp1.le.EPS_THRESH.and.r0.le.EPS_THRESH) cycle
#endif
              ! delta Q at interface
            dqi  = rp1-r0
              ! Set flux based on upwind velocity for color equation
              !  (equals conservative form if div.v=0)
            yfs(j,fluc_r) = max(0.0_ip,usig_p12)*dqi ! flux OUT OF i,j cell to the cell to the back (pos y)
            yfs(j,fluc_l) = min(0.0_ip,usig_p12)*dqi ! flux INTO i,j cell from the cell to the back (pos y)
              ! Modification for conservative form in
              ! divergent/convergent velocities
            divu_p   = max(0.0_ip,usig_p12) - max(0.0_ip,usig_m12)
            divu_m   = min(0.0_ip,usig_p32) - min(0.0_ip,usig_p12)
            yfs(j,fluc_r) = yfs(j,fluc_r) + rp1 * divu_m
            yfs(j,fluc_l) = yfs(j,fluc_l) + r0  * divu_p

              ! Only calculate dqu and theta if a limiter is being used
              ! Get delta Q at upwind interface relative to interface j
            if(usig_p12.gt.0.0_ip)then
              dqu = r0-rm1  ! upwind is interface j-1
            else
              dqu = rp2-rp1 ! upwind is interface j+1
            endif

              ! Make sure that theta is not singular
            if(abs(dqu).gt.3.0_ip*abs(dqi).or.abs(dqi).lt.EPS_THRESH)then
              theta = sign(3.0_ip,dqu*dqi)
            else
              theta = dqu/dqi
            endif

            ldq = dqi*max(0.0_ip,min(1.0_ip,2.0_ip*theta),min(2.0_ip,theta))

            au = abs(usig_p12)
            yfss(j) = 0.5e0_ip*au*(1.0e0_ip-dt_vol*au)*ldq

              ! For boundary fluxes, we want the flux-difference form rather
              ! than the fluctuation form given above and just use first-order
              ! upwind
            if(j.eq.1)     bflux_s =-min(0.0_ip,usig_m12)*r0
            if(j.eq.nymax) bflux_n = max(0.0_ip,usig_p12)*r0

          enddo  ! loop over j

          ! Now loop over y again and apply these fluxes to volume j
          do j=1,nymax
            if(IsLatLon)then
              vol     = kappa_pd(i,j  ,0)
              dt_vol  = dt/vol
            endif

#ifndef NOCYCLE
            rp1 = dum_depo(i,j+1,n)
            r0  = dum_depo(i,j  ,n)
            if(j.gt.0)then
              !rm1 = concen(i,j-1,k,n,ts0)
              rm1 = dum_depo(i,j-1,n)
            else
              rm1 = 0.0_ip
            endif
            if(rm1+r0+rp1 .lt. EPS_THRESH) cycle
#endif
            ! Note: if j=1, yfs(0,:)=0.0 and yfss(0)=0.0
            order1 =  dt_vol * ( -yfs(j-1,fluc_r) -  yfs(j,fluc_l))    ! Eq 19.19 (line 1)
            order2 =  dt_vol * ( yfss(j-1)        - yfss(j)       )    ! Eq 19.19 (line 3)
            !concen(i,j,k,n,ts1) = concen(i,j,k,n,ts0) + order1 + order2
            depo(i,j,n) = dum_depo(i,j,n) + order1 + order2

            if(j.eq.1)then
              ! Flux out the south side
              outflow_xz1_pd(i,0,n) = outflow_xz1_pd(i,0,n) + bflux_s * dt_vol
            elseif(j.eq.nymax)then
              ! Flux out the north side
              outflow_xz2_pd(i,0,n) = outflow_xz2_pd(i,0,n) + bflux_n * dt_vol
            endif

          enddo
        enddo ! loop over i=1,nxmax
      enddo ! loop over idx_dum

      end subroutine advect_deposit

!******************************************************************************

      end module ocean_currents
