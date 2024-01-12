      module Source_Resuspension

      use precis_param

      use io_units

      ! Set the number of output variables for this module
      integer, parameter :: nvar_User2d_static_XY_SrcResusp = 1 ! DepositMask
      integer, parameter :: nvar_User2d_XY_SrcResusp        = 1 ! FricVel
      integer, parameter :: nvar_User3d_XYGs_SrcResusp      = 0
      integer, parameter :: nvar_User3d_XYZ_SrcResusp       = 0
      integer, parameter :: nvar_User4d_XYZGs_SrcResusp     = 0

      character(len=30),dimension(nvar_User2d_static_XY_SrcResusp) :: temp_2ds_name_SrcResusp
      character(len=30),dimension(nvar_User2d_static_XY_SrcResusp) :: temp_2ds_unit_SrcResusp
      character(len=30),dimension(nvar_User2d_static_XY_SrcResusp) :: temp_2ds_lname_SrcResusp
      real(kind=op),    dimension(nvar_User2d_static_XY_SrcResusp) :: temp_2ds_MissVal_SrcResusp
      real(kind=op),    dimension(nvar_User2d_static_XY_SrcResusp) :: temp_2ds_FillVal_SrcResusp

      character(len=30),dimension(nvar_User2d_XY_SrcResusp) :: temp_2d_name_SrcResusp
      character(len=30),dimension(nvar_User2d_XY_SrcResusp) :: temp_2d_unit_SrcResusp
      character(len=30),dimension(nvar_User2d_XY_SrcResusp) :: temp_2d_lname_SrcResusp
      real(kind=op),    dimension(nvar_User2d_XY_SrcResusp) :: temp_2d_MissVal_SrcResusp
      real(kind=op),    dimension(nvar_User2d_XY_SrcResusp) :: temp_2d_FillVal_SrcResusp

      ! These are used to keep track of which index in the global list, this
      ! modules output vars corespond to
      integer :: indx_User2d_static_XY_SrcResusp
      integer :: indx_User2d_XY_SrcResusp
      integer :: indx_User3d_XYGs_SrcResusp
      integer :: indx_User3d_XYZ_SrcResusp
      integer :: indx_User4d_XYZGs_SrcResusp

      !character(len=30),dimension(nvar_User2d_static_XY_SrcResusp) :: SourceResusp_2dvarname
      !character(len=30),dimension(nvar_User2d_static_XY_SrcResusp) :: SourceResusp_2dvarunit
      !real(kind=ip)    ,dimension(nvar_User2d_static_XY_SrcResusp) :: SourceResusp_2dvar_FillValue
      !integer          ,dimension(nvar_User2d_static_XY_SrcResusp) :: SourceResusp_2dvar_id

      !character(len=30),dimension(nvar_User3d_XYZ_SrcResusp) :: SourceResusp_3dvarname
      !character(len=30),dimension(nvar_User3d_XYZ_SrcResusp) :: SourceResusp_3dvarunit
      !real(kind=ip)    ,dimension(nvar_User3d_XYZ_SrcResusp) :: SourceResusp_3dvar_FillValue
      !integer          ,dimension(nvar_User3d_XYZ_SrcResusp) :: SourceResusp_3dvar_id


      real(kind=sp),dimension(:,:),allocatable :: SnD_meso_last_step_sp
      real(kind=sp),dimension(:,:),allocatable :: SnD_meso_next_step_sp
      real(kind=sp),dimension(:,:),allocatable :: P0_meso_last_step_sp
      real(kind=sp),dimension(:,:),allocatable :: P0_meso_next_step_sp

      character (len=130), public :: DepPerimInfile
      integer :: DepMaskCount
      real(kind=ip) :: ustar
      real(kind=ip) :: u_star_thresh
      real(kind=ip) :: SnD_thresh = 0.1_ip ! in meters
      real(kind=ip) :: P0_thresh  = 0.1_ip ! in mm/hour
      real(kind=ip) :: Fv_coeff
      integer       :: FvID
      integer,dimension(:,:)  ,allocatable :: DepositMask

      logical :: useResuspension = .false.

      contains

!******************************************************************************

      subroutine input_data_Source_Resuspension

      use global_param,      only : &
         nmods,CheckConditions

      use io_data,       only : &
         infile

      use Source,        only : &
         SourceType,e_Duration,e_PlumeHeight,e_Volume

      use time_data,     only : &
         Simtime_in_hours

      implicit none

      !character(len=3)  :: answer
      character(len=80)  :: linebuffer080
      character(len=120) :: linebuffer120
      character(len=130) :: linebuffer130
      character :: testkey !,testkey2
      integer :: ios !,ioerr
      character(len=20) :: mod_name
      integer :: substr_pos
      integer,dimension(3) :: idx
      integer         :: iyear  ! time data read from files
      integer         :: imonth
      integer         :: iday
      real(kind=ip)   :: hour   ! Start time of eruption in

      ! First check if the requested source belongs to this module

      idx(1) = index(SourceType,'resus')
      idx(2) = index(SourceType,'Resus')
      idx(3) = index(SourceType,'RESUS')
      if (sum(idx(1:3)).gt.0)then
        SourceType='resuspens'
        useResuspension = .true.
      else
        useResuspension = .false.
        return
      endif

      ! If we've made it here, the requested source is a resuspension source
      ! open the input file again to get needed info
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
        write(outlog(io),*)"Start reading resuspension source."
      endif;enddo

      !read start time, duration, plume height
      !                Duration
      !                |   Height
      !                |   | EmissionScheme (1-west, 2=Lead, 3=Mort)
      !                |   |  |   Ustar
      !                |   |  |    |   Fv_coeff
      !                |   |  |    |    |
      !                v   v  v    v    v
      ! 0 0 0 0.0   35.0 5.0  2  0.5   6.0e-8

      read(linebuffer080,*,err=1910) iyear,imonth,iday,hour, &
                            e_Duration(1),e_PlumeHeight(1),&
                            FvID,u_star_thresh,Fv_coeff
      read(10,'(a130)')linebuffer130
      DepPerimInfile = adjustl(trim(linebuffer130))

      if( e_Duration(1).gt.Simtime_in_hours)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Source duration is longer than the model simulation time."
          write(outlog(io),*)"Resetting to simulation time."
        endif;enddo
        e_Duration(1) = Simtime_in_hours
      endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)iyear,imonth,iday,real(hour,kind=sp)
        write(outlog(io),*)"eDur   = ",real(e_Duration(1),kind=sp)
        write(outlog(io),*)"Height = ",real(e_PlumeHeight(1),kind=sp)
      endif;enddo
      if(FvID.eq.1)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"FvID = 1: Westphal scheme"
        endif;enddo
      elseif(FvID.eq.2)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"FvID = 2: Leadbetter scheme"
        endif;enddo
      elseif(FvID.eq.3)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"FvID = 3: Morticorena scheme"
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Resuspension scheme not recognized."
          write(errlog(io),*)"  Please use 1: Westphal"
          write(errlog(io),*)"             2: Leadbetter"
          write(errlog(io),*)"             3: Morticorena"
        endif;enddo
        stop 1
      endif
      if(u_star_thresh.gt.0.0_ip)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"u_star_thresh = ",&
                              real(u_star_thresh,kind=sp)
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"u_star_thresh will be calculated from local conditions."
          write(outlog(io),*)"Not yet implemented."
        endif;enddo
      endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"FV_scaling coefficient = ",&
                            real(Fv_coeff,kind=sp)
        write(outlog(io),*)DepPerimInfile
      endif;enddo

      ! Initialize some eruption values
      e_Duration(:)  = 0.0_ip
      e_Volume(:)    = 0.0_ip

      ! Now read to the end of the input file and read the Optional Modudle
      ! block
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Searching for OPTMOD=SRC_RESUSP"
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
          if(adjustl(trim(mod_name)).eq.'SRC_RESUSP')then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Found SRC_RESUSP block again"
            endif;enddo
            exit
          endif
          !OPTMOD_names(nmods) = adjustl(trim(mod_name))
        endif
1104    format(7x,a20)
      enddo

      ! Disable some of the stop conditions
      CheckConditions(1) = .false.
      CheckConditions(3) = .false.
      CheckConditions(4) = .false.
      CheckConditions(5) = .false.

!2010  continue
      close(10)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Finished input_data_Source_Resuspension"
      endif;enddo

      return

1900  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot find input file: ',infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

1910  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading start time, duration, height or',&
                    ' volume of an eruptive pulse.  Program stopped'
      endif;enddo
      stop 1

      end subroutine input_data_Source_Resuspension

!******************************************************************************

      subroutine Allocate_Source_Resuspension

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User2d_XY,nvar_User3d_XYGs,&
         nvar_User3d_XYZ,nvar_User4d_XYZGs

      use mesh,          only : &
         nxmax,nymax,nsmax

      use Source,        only : &
         SourceNodeFlux_Area

      implicit none

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOCATE_RESUSP -----------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      ! Set the start indicies
      indx_User2d_static_XY_SrcResusp = nvar_User2d_static_XY
      indx_User2d_XY_SrcResusp        = nvar_User2d_XY
      indx_User3d_XYGs_SrcResusp      = nvar_User3d_XYGs
      indx_User3d_XYZ_SrcResusp       = nvar_User3d_XYZ
      indx_User4d_XYZGs_SrcResusp     = nvar_User4d_XYZGs

      ! Deposite Mask
      temp_2ds_name_SrcResusp(1)    = "DepoMask"
      temp_2ds_lname_SrcResusp(1)   = "Source region"
      temp_2ds_unit_SrcResusp(1)    = "flag"
      temp_2ds_MissVal_SrcResusp(1) = -9999.0_op
      temp_2ds_FillVal_SrcResusp(1) = -9999.0_op

      ! Friction Velocity
      temp_2d_name_SrcResusp(1)    = "FricVel"
      temp_2d_lname_SrcResusp(1)   = "Friction Velocity"
      temp_2d_unit_SrcResusp(1)    = "m/s"
      temp_2d_MissVal_SrcResusp(1) = -9999.0_op
      temp_2d_FillVal_SrcResusp(1) = -9999.0_op

      nvar_User2d_static_XY = nvar_User2d_static_XY + nvar_User2d_static_XY_SrcResusp
      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_SrcResusp
      nvar_User3d_XYGs      = nvar_User3d_XYGs      + nvar_User3d_XYGs_SrcResusp
      nvar_User3d_XYZ       = nvar_User3d_XYZ       + nvar_User3d_XYZ_SrcResusp
      nvar_User4d_XYZGs     = nvar_User4d_XYZGs     + nvar_User4d_XYZGs_SrcResusp

      allocate(DepositMask(nxmax,nymax))

      allocate(SnD_meso_last_step_sp(nxmax,nymax))
      allocate(SnD_meso_next_step_sp(nxmax,nymax))
      allocate(P0_meso_last_step_sp(nxmax,nymax))
      allocate(P0_meso_next_step_sp(nxmax,nymax))

      allocate(SourceNodeFlux_Area(1:nxmax,1:nymax,1:nsmax))

      end subroutine Allocate_Source_Resuspension

!******************************************************************************

      subroutine Prep_output_Source_Resuspension

      use Variable_Diffusivity, only : &
         FricVel_ip

      use Output_Vars,   only : &
         var_User2d_XY_name,var_User2d_XY_unit,var_User2d_XY_lname,&
         var_User2d_XY_MissVal,var_User2d_XY_FillVal,var_User2d_XY,&
         var_User2d_static_XY_name,var_User2d_static_XY_unit,&
         var_User2d_static_XY_lname,var_User2d_static_XY_MissVal,&
         var_User2d_static_XY_FillVal,var_User2d_static_XY

      implicit none

      integer :: i,indx

      do i=1,nvar_User2d_static_XY_SrcResusp
        indx = indx_User2d_static_XY_SrcResusp+i
        var_User2d_static_XY_name(indx) = temp_2ds_name_SrcResusp(i)
        var_User2d_static_XY_unit(indx) = temp_2ds_unit_SrcResusp(i)
        var_User2d_static_XY_lname(indx)= temp_2ds_lname_SrcResusp(i)
        var_User2d_static_XY_MissVal(indx)= temp_2ds_MissVal_SrcResusp(i)
        var_User2d_static_XY_FillVal(indx)= temp_2ds_FillVal_SrcResusp(i)
        if(i.eq.1)var_User2d_static_XY(:,:,indx) = DepositMask
      enddo

      do i=1,nvar_User2d_XY_SrcResusp
        indx = indx_User2d_XY_SrcResusp+i
        var_User2d_XY_name(indx) = temp_2d_name_SrcResusp(i)
        var_User2d_XY_unit(indx) = temp_2d_unit_SrcResusp(i)
        var_User2d_XY_lname(indx)= temp_2d_lname_SrcResusp(i)
        var_User2d_XY_MissVal(indx)= temp_2d_MissVal_SrcResusp(i)
        var_User2d_XY_FillVal(indx)= temp_2d_FillVal_SrcResusp(i)
        if(i.eq.1)var_User2d_XY(:,:,indx) = real(FricVel_ip(:,:),kind=op)
      enddo

      end subroutine Prep_output_Source_Resuspension

!******************************************************************************

      subroutine Deallocate_Source_Resuspension

      use Source,        only : &
         SourceNodeFlux_Area

      implicit none

      if(allocated(DepositMask))           deallocate(DepositMask)
      if(allocated(SnD_meso_last_step_sp)) deallocate(SnD_meso_last_step_sp)
      if(allocated(SnD_meso_next_step_sp)) deallocate(SnD_meso_next_step_sp)
      if(allocated(P0_meso_last_step_sp))  deallocate(P0_meso_last_step_sp)
      if(allocated(P0_meso_next_step_sp))  deallocate(P0_meso_next_step_sp)

#ifdef USEPOINTERS
      if(associated(SourceNodeFlux_Area))   deallocate(SourceNodeFlux_Area)
#else
      if(allocated(SourceNodeFlux_Area))   deallocate(SourceNodeFlux_Area)
#endif
      end subroutine Deallocate_Source_Resuspension

!******************************************************************************

!******************************************************************************

      subroutine Read_Deposit_Perimeter

      use mesh,          only : &
         nxmax,nymax,x_cc_pd,y_cc_pd,lon_cc_pd,lat_cc_pd,dx,dy,de,dn,IsLatLon,&
         ivent,jvent,A3d_iprojflag, &
         A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,A3d_k0_scale,A3d_Re

      use projection,    only : &
           PJ_proj_for

      implicit none

      integer :: npoints
      real(kind=ip),dimension(:),allocatable :: DepPerm_lon
      real(kind=ip),dimension(:),allocatable :: DepPerm_lat

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
      real(kind=dp)  :: lat_in,lon_in,xout,yout
      logical                 :: IsThere

      integer :: i,j

      ! Test for existence of the deposit file
      inquire( file=adjustl(trim(DepPerimInfile)), exist=IsThere )
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)&
          "Trying to read deposit contour file"
        write(outlog(io),*)"Full file name = ",adjustl(trim(DepPerimInfile))
        write(outlog(io),*)"  Exists = ",IsThere
      endif;enddo
      if(.not.IsThere)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Could not find contour file."
          write(errlog(io),*)"       Please copy file to current directory"
        endif;enddo
        stop 1
      endif

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Opening ",DepPerimInfile
      endif;enddo

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
          !if(lon_in.lt.0.0)lon_in=lon_in+360.0_dp
          lat_in = real(DepPerm_lat(i),kind=dp)
          call PJ_proj_for(lon_in,lat_in,A3d_iprojflag, &
                       A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,A3d_k0_scale,A3d_Re, &
                      xout,yout)
          BCpos(1,i) = real(xout,kind=ip)
          BCpos(2,i) = real(yout,kind=ip)
        enddo
        loc_x(1:nxmax) = x_cc_pd(1:nxmax)
        loc_y(1:nymax) = y_cc_pd(1:nymax)
        loc_dx = dx
        loc_dy = dy
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)dx,loc_x(1),loc_x(nxmax)
          write(outlog(io),*)dy,loc_y(1),loc_y(nymax)
        endif;enddo
      endif
        ! Get the min and max extents of the deposit
      DepPerm_x_min = minval(BCpos(1,:))
      DepPerm_x_max = maxval(BCpos(1,:))
      DepPerm_y_min = minval(BCpos(2,:))
      DepPerm_y_max = maxval(BCpos(2,:))
        ! And the indicies on the computational grid bracketing these points
      !write(outlog(io),*)"HFS: double-check this calculation",loc_dx,loc_dy
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

      close(20)

      end subroutine Read_Deposit_Perimeter

!******************************************************************************

      subroutine Set_Resusp_Meso(Load_MesoSteps,Interval_Frac)

      use mesh,          only : &
         nxmax,nymax

      use Variable_Diffusivity, only : &
         FricVel_ip,FricVel_meso_last_step_sp,FricVel_meso_next_step_sp,&
         FricVel_meso_next_step_Met_sp

      use MetReader,     only : &
         MR_dum2d_met,MR_dum2d_comp,MR_iMetStep_Now, &
           MR_Regrid_Met2d_to_Comp2d, &
           MR_Read_2d_Met_Variable_to_CompGrid

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=ip),intent(in) :: Interval_Frac

      integer :: ivar

      logical,save :: first_time = .true.

      if(Load_MesoSteps)then
        if(first_time)then
          ! Need to fill _last_step_sp
          !  First fill next step so that outside this 'first_time' loop, the
          !  'next' can be copied to the 'last'

           ! Now resample Friction Velocity onto computational grid
          ivar = 13
          call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
          FricVel_meso_next_step_sp = MR_dum2d_comp
          !Note:  This was read in variable_diffusivity on the met grid and now
          !       only needs to be resampled
          if(allocated(FricVel_meso_next_step_Met_sp)) then
            write(*,*)"FricVel_meso_next_step_Met_sp is allocated"
          else
            write(*,*)"FricVel_meso_next_step_Met_sp is NOT allocated"
            stop 5
          endif
          MR_dum2d_met = FricVel_meso_next_step_Met_sp
          call MR_Regrid_Met2d_to_Comp2d
          FricVel_meso_next_step_sp = MR_dum2d_comp

          ivar = 15
          call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
          SnD_meso_next_step_sp  = MR_dum2d_comp

          ivar = 44
          call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
          P0_meso_next_step_sp  = MR_dum2d_comp*3600.0_sp 

          first_time = .false.
        endif ! first_time

        FricVel_meso_last_step_sp      = FricVel_meso_next_step_sp
        SnD_meso_last_step_sp          = SnD_meso_next_step_sp
        P0_meso_last_step_sp           = P0_meso_next_step_sp

        ! Need to fill _next_step_sp
         ! Now resample onto computational grid
        !ivar = 13
        !call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
        !FricVel_meso_next_step_sp = MR_dum2d_comp
        MR_dum2d_met = FricVel_meso_next_step_Met_sp
        call MR_Regrid_Met2d_to_Comp2d
        FricVel_meso_next_step_sp = MR_dum2d_comp

        ivar = 15
        call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
        SnD_meso_next_step_sp  = MR_dum2d_comp*1.0_sp

        ivar = 44
        call MR_Read_2d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
        P0_meso_next_step_sp  = MR_dum2d_comp*3600.0_sp ! convert from kg/m2/s to mm/hr

      endif

      FricVel_ip(1:nxmax,1:nymax) = &
              real(FricVel_meso_last_step_sp(1:nxmax,1:nymax),kind=ip) + &
             real((FricVel_meso_next_step_sp(1:nxmax,1:nymax) - &
                   FricVel_meso_last_step_sp(1:nxmax,1:nymax)),kind=ip) * &
                                       Interval_Frac

      end subroutine Set_Resusp_Meso

!******************************************************************************

      subroutine Set_Resusp_Flux

      use global_param,  only : &
         KM_2_M,KM_2_M,HR_2_S,GRAV

      use mesh,          only : &
         nxmax,nymax

      use Source,        only : &
         SourceNodeFlux_Area

      use Variable_Diffusivity, only : &
         FricVel_ip

      implicit none

      integer :: i,j
      real(kind=ip) :: u_star
      real(kind=ip) :: Fv

      SourceNodeFlux_Area = 0.0_ip
      do i = 1,nxmax
        do j = 1,nymax
          if(SnD_meso_next_step_sp(i,j).gt.SnD_thresh) cycle
          if(P0_meso_next_step_sp(i,j).gt.P0_thresh) cycle

          if(FricVel_ip(i,j).gt.u_star_thresh)then
            if(DepositMask(i,j).gt.0)then
              if(FvID.eq.1)then
                ! Westphal
                !  Fv_coeff should be 1.0e-5
                u_star = FricVel_ip(i,j)**4.0_ip
                Fv = Fv_coeff*u_star
              elseif(FvID.eq.2)then
                ! Leadbetter
                u_star = (FricVel_ip(i,j)-u_star_thresh)**3.0_ip
                Fv = Fv_coeff*u_star
              elseif(FvID.eq.3)then
                ! Morticorena
                u_star = (FricVel_ip(i,j)-u_star_thresh)**2.0_ip
                u_star = u_star * FricVel_ip(i,j)
                Fv = Fv_coeff*u_star*2000.0_ip/GRAV
              else
                do io=1,2;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"Unknown resuspension scheme: ",FvID
                endif;enddo
                stop 1
              endif
              ! assign to array and convert from m-2 s-1 to km-2 hr-1
              SourceNodeFlux_Area(i,j,:)= FV *KM_2_M*KM_2_M*HR_2_S
            else
              SourceNodeFlux_Area(i,j,:) = 0.0_ip
            endif
          else
            SourceNodeFlux_Area(i,j,:) = 0.0_ip
          endif
        enddo
      enddo

      end subroutine Set_Resusp_Flux

!******************************************************************************

      subroutine Set_concen_Resusp

      use mesh,          only : &
         nxmax,nymax,nsmax,kappa_pd,ts0

      use solution,      only : &
         concen_pd

      use time_data,     only : &
         dt

      use Source,       only : &
         SourceNodeFlux_Area

      implicit none

      integer :: i,j,k

      do i = 1,nxmax
        do j = 1,nymax
          k = 1
          concen_pd(i,j,k,1:nsmax,ts0) =                     &
                   concen_pd(i,j,k,1:nsmax,ts0)  &
                   + dt*SourceNodeFlux_Area(i,j,1:nsmax)/kappa_pd(i,j,1)
        enddo
      enddo

      end subroutine Set_concen_Resusp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SourceVolInc_Resusp(dt)
!
!  Called from: Ash3d_res.F90
!  Arguments:
!    dt = time step in hours
!
!  This function calculates the total tephra volume inserted in this time step.
!  It is used only for mass-conservation error-checking.  It is similar to
!  SourceVolInc but act on an area source.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function SourceVolInc_Resusp(dt)

      use global_param,  only : &
         KM3_2_M3

      use mesh,          only : &
         kappa_pd,nxmax,nymax


      use Tephra,        only : &
         n_gs_max,MagmaDensity

      use Source,        only : &
         SourceNodeFlux_Area

      real(kind=ip) :: SourceVolInc_Resusp
      real(kind=dp) :: dt

      real(kind=ip) :: tmp
      integer :: i,j,ii,jj,k,isize

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine SourceVolInc_Resusp"
      endif;enddo

      tmp = 0.0_ip

      do isize=1,n_gs_max
        do i = 1,nxmax
          do j = 1,nymax
            k = 1
            tmp = tmp                             + & ! final units is km3
                  real(dt,kind=ip)                * & ! hr
                  SourceNodeFlux_Area(i,j,isize)  / & ! kg/km3 hr
                  kappa_pd(i,j,1)                 / & ! km3
                  MagmaDensity                    / & ! kg/m3
                  KM3_2_M3                            ! m3/km3
          enddo
        enddo
      enddo

      SourceVolInc_Resusp = tmp

      return

      end function SourceVolInc_Resusp

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
            IsIn = .false.
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

      end module Source_Resuspension

