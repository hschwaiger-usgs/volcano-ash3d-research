      program Ash3d

      use precis_param

      use io_units

      use global_param,  only : &
         useCalcFallVel,useDiffusion,useHorzAdvect,useVertAdvect,&
         HR_2_S,useTemperature,DT_MIN,KM3_2_M3,EPS_TINY,EPS_SMALL,&
         nmods,OPTMOD_names,StopConditions,CheckConditions, &
         useVarDiffH,useVarDiffV

      use mesh,          only : &
         ivent,jvent,nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd

      use solution,      only : &
         concen_pd,DepositGranularity,StopValue,aloft_percent_remaining, &
         SourceCumulativeVol,dep_vol,aloft_vol,outflow_vol,tot_vol

      use Output_Vars,   only : &
         DepositAreaCovered,DepositThickness,LoadVal,CloudLoadArea,&
         Calculated_Cloud_Load,Calculated_AshThickness,Calc_vprofile, &
           Allocate_Output_Vars, &
           Allocate_Output_UserVars, &
           Allocate_NTime,   &
           Allocate_Profile, &
           Gen_Output_Vars,&
           FirstAsh

      use io_data,       only : &
         Called_Gen_Output_Vars,isFinal_TS,LoadConcen,log_step,&
         Output_at_logsteps,Output_at_WriteTimes,Output_every_TS,&
         NextWriteTime,iTimeNext,nvprofiles,nWriteTimes,&
         Write_PT_Data,Write_PR_Data

      use time_data,     only : &
         time,dt,Simtime_in_hours,t0,t1,ntmax

      use Source,        only : &
         ibase,itop,SourceNodeFlux,e_EndTime_final,e_Volume,MassFluxRate_now,&
         SourceType,SourceNodeFlux_Area,&
           Allocate_Source_time,&
           MassFluxCalculator,&
           TephraSourceNodes

      use Tephra,        only : &
         n_gs_max,n_gs_aloft,MagmaDensity,&
           Allocate_Tephra,&
           Allocate_Tephra_Met,&
           Prune_GS
   
      use Atmosphere,    only : &
           Allocate_Atmosphere_Met

      use AdvectionHorz, only : &
           AdvectHorz

      use AdvectionVert_DCU, only : &
           advect_z

      use Diffusion,     only : &
           DiffuseHorz,&
           DiffuseVert

      use Airports,      only : &
         Airport_thickness_TS,Airport_thickness,nairports,&
           ReadAirports

#ifdef USENETCDF
      use Ash3d_Netcdf,  only : &
           NC_RestartFile_LoadConcen
#endif

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert 'use' statements here
!
#ifdef TOPO
      use Topography
#endif
#ifdef LC
      use land_cover
#endif
#ifdef OSCAR
      use ocean_currents
#endif
#ifdef WETDEPO
      use Wet_Deposition
#endif
#ifdef VARDIFF
      use variable_diffusivity
#endif
#ifdef SRC_SAT
      use Source_Satellite
#endif
#ifdef SRC_RESUSP
      use Source_Resuspension
#endif
#ifdef SRC_GAS
      use Source_Gas
#endif
!------------------------------------------------------------------------------

      implicit none

      integer               :: iostatus
      integer               :: itime
      integer               :: j,k
      integer               :: ii,jj,iz,isize
      real(kind=ip)         :: avgcon        ! avg concen of cells in umbrella
      real(kind=ip)         :: Interval_Frac
      logical               :: Load_MesoSteps
      logical               :: StopTimeLoop   = .false.
      logical               :: first_time     = .true.
      character(len=130)    :: tmp_str
      real(kind=ip)         :: MassConsErr

      INTERFACE
        subroutine Set_OS_Env
        end subroutine Set_OS_Env
        subroutine Read_Control_File
        end subroutine Read_Control_File
        subroutine input_data_ResetParams
        end subroutine input_data_ResetParams
        subroutine alloc_arrays
        end subroutine alloc_arrays
        subroutine calc_mesh_params
        end subroutine calc_mesh_params
        subroutine MesoInterpolater(TimeNow,Load_MesoSteps,Interval_Frac,first_time)
          integer,parameter  :: dp         = 8 ! Double precision
          real(kind=dp),intent(in)    :: TimeNow
          real(kind=dp),intent(out)   :: Interval_Frac
          logical      ,intent(inout) :: Load_MesoSteps
          logical      ,intent(in)    :: first_time
        end subroutine MesoInterpolater
        subroutine output_results
        end subroutine output_results
        subroutine Set_BC(bc_code)
          integer,intent(in) :: bc_code ! 1 for advection, 2 for diffusion
        end subroutine Set_BC
        subroutine vprofilewriter(itime)
          integer, intent(in) :: itime
        end subroutine vprofilewriter
        subroutine TimeStepTotals(itime)
          integer, intent(in) :: itime
        end subroutine TimeStepTotals
        subroutine dealloc_arrays
        end subroutine dealloc_arrays
      END INTERFACE

!      ! Before we do anything, start a log file
!      open(unit=global_log,file='Ash3d.lst',status='unknown')

      call Set_OS_Env

      aloft_percent_remaining = 1.0_ip
      SourceCumulativeVol     = 0.0_ip

      call cpu_time(t0) !time is a scaler real

        ! input data for ash transport
      call Read_Control_File

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to custom input blocks here
!
!  Loop through all the optional modules compiled and test against the list
!  from the input file (e.g. OPTMOD=TOPO), then call the special input reader
!  for that block
!  Do a sanity check on optional module requested in the input file v.s. those
!  compiled in this executable and for consistency among modules:
!  e.g. SRC_RESUSP will require the VARDIFF and LC be set
!
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Now looping through optional modules found in input file"
      endif;enddo
      do j=1,nmods
        do io=1,2;if(VB(io).le.verbosity_essential)then
          write(outlog(io),*)"Testing for ",OPTMOD_names(j),j
        endif;enddo
        if(OPTMOD_names(j).eq.'RESETPARAMS')then
          do io=1,2;if(VB(io).le.verbosity_essential)then
            write(outlog(io),*)"  Reading input block for RESETPARAMS"
          endif;enddo
          call input_data_ResetParams
        endif
#ifdef TOPO
        if(OPTMOD_names(j).eq.'TOPO')then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Reading input block for TOPO"
          endif;enddo
          call input_data_Topo
        endif
#endif
#ifdef LC
        if(OPTMOD_names(j).eq.'LC')then
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"  Reading input block for LC"
          endif;enddo
          call input_data_LC
        endif
#endif
#ifdef OSCAR
        if(OPTMOD_names(j).eq.'OSCAR')then
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"  Reading input block for OSCAR"
          endif;enddo
          call input_data_OSCAR
        endif
#endif
#ifdef WETDEPO
        if(OPTMOD_names(j).eq.'WETDEPO')then
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"  Reading input block for WETDEPO"
          endif;enddo
          call input_data_WetDepo
        endif
#endif
#ifdef VARDIFF
        if(OPTMOD_names(j).eq.'VARDIFF')then
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"  Reading input block for VARDIFF"
          endif;enddo
          call input_data_VarDiff
        endif
#endif
#ifdef SRC_RESUSP
        if(OPTMOD_names(j).eq.'SRC_RESUSP')then
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"  Reading input block for SRC_RESUSP"
          endif;enddo
          call input_data_Source_Resuspension
        endif
#endif
#ifdef SRC_GAS
        if(OPTMOD_names(j).eq.'SRC_GAS')then
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"  Reading input block for SRC_GAS"
          endif;enddo
          call input_data_Source_Gas
        endif
#endif
#ifdef SRC_SAT
        if(OPTMOD_names(j).eq.'SRC_SAT')then
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"  Reading input block for SRC_SAT"
          endif;enddo
          call input_data_Source_Satellite
        endif
#endif
      enddo
      do io=1,2;if(VB(io).le.verbosity_info)then    
        write(outlog(io),*)"Finished reading all specialized input blocks"
      endif;enddo
!
!------------------------------------------------------------------------------

        ! Read airports/POI and allocate/initilize arrays
        ! We only need to do this if an output variable demands it since this is
        ! a burden every time step
      if(Output_every_TS) &
        call ReadAirports

      call alloc_arrays
        ! Set up grids for solution and Met data
      call calc_mesh_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Initialize concen and any special source terms here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(LoadConcen)then
        ! We are initializing the concentration and time from an output file
        ! Currently, Ash3d assumes the concentration file is compatible with
        ! the computational grid and grainsize distribution
#ifdef USENETCDF
        call NC_RestartFile_LoadConcen
#else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Loading concentration files requires previous netcdf"
          write(errlog(io),*)"       output.  This Ash3d executable was not compiled with"
          write(errlog(io),*)"       netcdf support.  Please recompile Ash3d with"
          write(errlog(io),*)"       USENETCDF=T, or select another source."
        endif;enddo
        stop 1
#endif
      else
        ! Initialize arrays if we haven't already loaded the concentration from
        ! a previous run
        concen_pd = 0.0_ip
        DepositGranularity = 0.0_ip
      endif
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert special source terms here
!
#ifdef SRC_SAT
      if(SourceType.eq.'satellite')then
        ! Satellite initialized runs have a concen array loaded here
        call Allocate_Source_Satellite
        call Read_SatMassLoading
      endif
#endif

#ifdef SRC_RESUSP
      if(SourceType.eq.'resuspens')then
        ! We will also need to allocate some more arrays
        call Allocate_Source_Resuspension
        ! For resuspension cases, the deposit defines the source region
        call Read_Deposit_Perimeter
      endif
#endif

#ifdef SRC_GAS
      if(SourceType.eq.'gas')then
        ! We will also need to allocate some more arrays
        call Allocate_Source_Gas
        if(EruptGasSrcStruc(1).eq.1)call Read_Perimeter_Gas
      endif
#endif

!------------------------------------------------------------------------------

      if(useTemperature)then
        call Allocate_Atmosphere_Met
      endif
      ! Now we can allocate the variables that live on the met grid
      if(useCalcFallVel)then
        call Allocate_Tephra_Met
      endif

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional variable allocation subroutines here
!
#ifdef TOPO
      if(useTopo)                      call Allocate_Topo(nxmax,nymax)
#endif
#ifdef LC
      if(useLandCover)                 call Allocate_LC
#endif
#ifdef OSCAR
      if(useOceanCurrent)              call Allocate_OSCAR
#endif
#ifdef VARDIFF
      if(useVarDiffV)                  call Allocate_Atmosphere_Met
      if(useVarDiffH.or.useVarDiffV)   call Allocate_VarDiff_Met
#endif
#ifdef WETDEPO
      if(USE_WETDEPO)then
        call Allocate_WetDepo_global
        call Allocate_WetDepo_Met
      endif
#endif
!------------------------------------------------------------------------------
      ! Allocate all the output variables
      call Allocate_Output_UserVars(nxmax,nymax,nzmax,nsmax)

      ! Now that we have the Met grids initialized, get the state variables
      ! interpoated on the start time
      time           = 0.0_ip
      Load_MesoSteps = .true.
      Interval_Frac  = 0.0_ip
      first_time     = .true.
      call MesoInterpolater(time , Load_MesoSteps , Interval_Frac, first_time)

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to special MesoInterpolaters subroutines here
!
#ifdef VARDIFF
        if(useVarDiffH)     call Set_VarDiffH_Meso(Load_MesoSteps,Interval_Frac)
        if(useVarDiffV)     call Set_VarDiffV_Meso(Load_MesoSteps,Interval_Frac)
#endif
#ifdef SRC_RESUSP
        if(useResuspension) call Set_Resusp_Meso(Load_MesoSteps,Interval_Frac)
#endif
#ifdef SRC_GAS
        ! No variables need Meso interpolation for SRC_GAS
#endif
#ifdef WETDEPO
        if(USE_WETDEPO)     call Set_WetDepo_Meso(Load_MesoSteps,Interval_Frac)
#endif

#ifdef TOPO
      if(useTopo) call get_topo
#endif
#ifdef LC
      if(useLandCover)then
        call load_LC
        call assign_LC
        !if(useResuspension) call get_roughlen
      endif
#endif
#ifdef OSCAR
      if(useOceanCurrent)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Ocean current branch disabled until rgrd2 is replaced."
        endif;enddo
        stop 1
        call Check_SurfaceVelocity
        call set_SurfaceVelocity(0.0_ip)
      endif
#endif
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to prep user-specified output
!
#ifdef TOPO
      if(useTopo) call Prep_output_Topo
#endif
#ifdef LC
      if(useLandCover) call Prep_output_LC
#endif
#ifdef OSCAR
      if(useOceanCurrent) call Prep_output_OSCAR
#endif
#ifdef VARDIFF
      if(useVarDiffH.or.useVarDiffV) call Prep_output_VarDiff
#endif
#ifdef WETDEPO
      if(USE_WETDEPO) call Prep_output_WetDepo
#endif
#ifdef SRC_SAT
      if(SourceType.eq.'satellite')then
        call Prep_output_SrcSat
      endif
#endif
#ifdef SRC_RESUSP
      if(SourceType.eq.'resuspens')then
        call Prep_output_Source_Resuspension
      endif
#endif
#ifdef SRC_GAS
      if(SourceType.eq.'gas')then
        call Prep_output_Source_Gas
      endif
#endif
!------------------------------------------------------------------------------

        ! Call output_results before time loop to create output files
      call output_results

      ntmax = max(1,3*int(Simtime_in_hours/dt))
      call Allocate_NTime(ntmax)
      if (Write_PR_Data)then
        call Allocate_Profile(nzmax,ntmax,nvprofiles)
      endif

      ! write "Building time array of plume height & eruption rate"
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),7)
      endif;enddo

      call Allocate_Source_time

      !call MassFluxCalculator          !find current mass flux & plume height

      ! Write out starting volume, max time steps, and headers for the table that follows
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),1) tot_vol,ntmax
      endif;enddo

      ! ************************************************************************
      ! ****** begin time simulation *******************************************
      ! ************************************************************************
      itime = 0
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Starting time loop."
      endif;enddo

      do while (StopTimeLoop.eqv..false.)
        ! Note: stop conditions are evaluated at the end of the time loop

        ! Copy previous t=n+1 to t=n
        concen_pd(:,:,:,:,ts0) = concen_pd(:,:,:,:,ts1)
          ! re-initialize slice (1)
        concen_pd(:,:,:,1:nsmax,ts1) = 0.0_ip

        Called_Gen_Output_Vars  = .false.
        Calculated_Cloud_Load   = .false.
        Calculated_AshThickness = .false.

        itime = itime + 1
        if(itime.gt.ntmax)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"WARNING: The number of time steps attempted exceeds 3x that anticipated."
            write(outlog(io),*)"         Check that the winds are stable"
            write(outlog(io),*)"        Simtime_in_hours = ",Simtime_in_hours
            write(outlog(io),*)"                   ntmax = ",ntmax
            write(outlog(io),*)"            current step = ",itime
          endif;enddo
        endif

          ! find the wind field at the current time
        first_time     = .false.
        call MesoInterpolater(time , Load_MesoSteps , Interval_Frac, first_time)

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to special MesoInterpolaters subroutines here
!
#ifdef VARDIFF
        do io=1,2;if(VB(io).le.verbosity_info)then    
          write(outlog(io),*)"Calling Set_VarDiffH_Meso."
        endif;enddo
        if(useVarDiffH)     call Set_VarDiffH_Meso(Load_MesoSteps,Interval_Frac)
        do io=1,2;if(VB(io).le.verbosity_info)then      
          write(outlog(io),*)"Calling Set_VarDiffV_Meso."
        endif;enddo
        if(useVarDiffV)     call Set_VarDiffV_Meso(Load_MesoSteps,Interval_Frac)
#endif
#ifdef SRC_RESUSP
        do io=1,2;if(VB(io).le.verbosity_info)then      
          write(outlog(io),*)"Calling Set_Resusp_Meso."
        endif;enddo
        if(useResuspension) call Set_Resusp_Meso(Load_MesoSteps,Interval_Frac)
#endif
#ifdef SRC_GAS
        do io=1,2;if(VB(io).le.verbosity_info)then      
          write(outlog(io),*)"Calling Set_Gas_Meso."
        endif;enddo
        if(USE_GAS)          call Set_Gas_Meso(Load_MesoSteps,Interval_Frac)
#endif
#ifdef WETDEPO
        do io=1,2;if(VB(io).le.verbosity_info)then      
          write(outlog(io),*)"Calling Set_WetDepo_Meso."
        endif;enddo
        if(USE_WETDEPO)     call Set_WetDepo_Meso(Load_MesoSteps,Interval_Frac)
#endif
!------------------------------------------------------------------------------

        call MassFluxCalculator         ! call subroutine that determines mass flux & plume height

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to specialized MassFluxRate calculations here
!
#ifdef SRC_RESUSP
        if(useResuspension)then
          do io=1,2;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)"Calling Set_Resusp_Flux."
          endif;enddo
          call Set_Resusp_Flux
          MassFluxRate_now = sum(SourceNodeFlux_Area)
        endif
#endif
#ifdef SRC_GAS
        if(USE_GAS)then
          do io=1,2;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)"Calling Set_Gas_Flux."
          endif;enddo
          call Set_Gas_Flux
          MassFluxRate_now = sum(SourceNodeFlux_Area)  ! kg/hr  
        endif
#endif
!------------------------------------------------------------------------------

        ! Add source term
        ! erupt ash into column
        if(MassFluxRate_now.gt.0.0_ip) then
          ! Check if the source type is one of the standard types
          if ((SourceType.eq.'point')  .or. &
              (SourceType.eq.'line')   .or. &
              (SourceType.eq.'profile').or. &
              (SourceType.eq.'suzuki') .or. &
              (SourceType.eq.'umbrella') .or. &
              (SourceType.eq.'umbrella_air'))then
            ! Calculating the flux at the source nodes
            call TephraSourceNodes

            ! Now integrate the ash concentration with the SourceNodeFlux
            if (SourceType.eq.'umbrella'.or. &
               (SourceType.eq.'umbrella_air')) then
              ! Umbrella clouds have a special integration
              !  Below the umbrella cloud, add ash to vent nodes as above
              concen_pd(ivent,jvent,1:ibase-1,1:n_gs_max,ts0) =          & ! 
                       concen_pd(ivent,jvent,1:ibase-1,1:n_gs_max,ts0) + & ! kg/km3
                       dt                                              * & ! hr
                       SourceNodeFlux(1:ibase-1,1:n_gs_max)                ! kg/km3 hr
              do isize=1,n_gs_max
                do k=1,ibase-1
                  SourceCumulativeVol = SourceCumulativeVol + & ! final units is km3
                    dt                              * & ! hr
                    SourceNodeFlux(k,isize)         * & ! kg/km3 hr
                    kappa_pd(ivent,jvent,k)         / & ! km3
                    MagmaDensity                    / & ! kg/m3
                    KM3_2_M3                            ! m3/km3
                enddo
              enddo
              do iz=ibase,itop
                !Within the cloud: first, average the concentration that curently
                !exists in the 9 cells surrounding the vent
                do isize=1,n_gs_max
                  avgcon=sum(concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,iz,isize,ts0))/9.0_ip
                  concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,iz,isize,ts0)=avgcon
                enddo
              enddo
              !Then, add tephra to the 9 nodes surrounding the vent
              ! TephraSourceNodes has a special line to reduce SourceNodeFlux by a factor 9
              ! because it is applied 9 times here.  We need to be careful about mixing mass
              ! and concentration since cell volume differ in lat, but this should be minor
              do ii=ivent-1,ivent+1
                do jj=jvent-1,jvent+1
                  do iz=ibase,itop
                    concen_pd(ii,jj,iz,1:n_gs_max,ts0) =                &
                              concen_pd(ii,jj,iz,1:n_gs_max,ts0)        &
                                 + dt*SourceNodeFlux(iz,1:n_gs_max)
                    do isize=1,n_gs_max
                      SourceCumulativeVol = SourceCumulativeVol + & ! final units is km3
                        dt                              * & ! hr
                        SourceNodeFlux(iz,isize)         * & ! kg/km3 hr
                        kappa_pd(ivent,jvent,iz)         / & ! km3
                        MagmaDensity                    / & ! kg/m3
                        KM3_2_M3                            ! m3/km3
                    enddo
                  enddo
                enddo
              enddo
            else ! (SourceType.eq.'umbrella' or 'umbrella_air')
              ! All other standard source types (point,line,profile, suzuki) are
              ! integrated as follows.
              concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0) =  &
              concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0)    &
                + dt*SourceNodeFlux(1:nzmax+1,1:n_gs_max)
              ! this part is just for book-keeping and error checking
              do isize=1,n_gs_max
                do k=1,nzmax+1
                  SourceCumulativeVol = SourceCumulativeVol + & ! final units is km3
                    dt                              * & ! hr
                    SourceNodeFlux(k,isize)         * & ! kg/km3 hr
                    kappa_pd(ivent,jvent,k)         / & ! km3
                    MagmaDensity                    / & ! kg/m3
                    KM3_2_M3                            ! m3/km3
                enddo
              enddo
            endif
          !else
            ! This is not a standard source.
            !do io=1,2;if(VB(io).le.verbosity_error)then
            !  write(errlog(io),*)"WARNING: source type is non-standard"
            !endif;enddo
            !stop 1
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional sources here
!         These subroutines need to calculate the mass of a species inserted
!         into specific cells and update concen accordingly
!
!          elseif () then
!
          elseif (SourceType.eq.'resuspens') then
#ifdef SRC_RESUSP
            do io=1,2;if(VB(io).le.verbosity_info)then      
              write(outlog(io),*)"Calling Set_concen_Resusp."
            endif;enddo
            call Set_concen_Resusp
#endif
          elseif (SourceType.eq.'gas') then
#ifdef SRC_GAS
            do io=1,2;if(VB(io).le.verbosity_info)then      
              write(outlog(io),*)"Calling Set_concen_Gas."
            endif;enddo
            call Set_concen_Gas
#endif
!------------------------------------------------------------------------------
          endif
        endif !MassFluxRate_now.gt.0.0_ip

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Set Boundary Conditions
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional boundary conditions here

#ifdef OSCAR
        do io=1,2;if(VB(io).le.verbosity_info)then      
          write(outlog(io),*)"Calling set_SurfaceVelocity."
        endif;enddo
        if(useOceanCurrent) call set_SurfaceVelocity(time)
#endif
!------------------------------------------------------------------------------
        call Set_BC(1)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Advection / Diffusion / Deposition
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional advection/diffusion routines here
!
#ifdef OSCAR
        do io=1,2;if(VB(io).le.verbosity_info)then      
          write(outlog(io),*)"Calling advect_deposit."
        endif;enddo
        if(useOceanCurrent) &
          call advect_deposit(concen_pd(1:nxmax,1:nymax,0,1:nsmax,ts0))
#endif
!------------------------------------------------------------------------------

        if(useHorzAdvect) call AdvectHorz(itime)

        if(useVertAdvect) call advect_z

        if(useDiffusion)then
          call Set_BC(2)
          call DiffuseVert
          call DiffuseHorz(itime)
        endif

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional deposition routines here
!
#ifdef SRC_GAS
        if(USE_GAS)then
          do io=1,2;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)"Calling Gas_Chem_Convert."
          endif;enddo
          call Gas_Chem_Convert
        endif
#endif

#ifdef WETDEPO
        if(USE_WETDEPO)then
          do io=1,2;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)"Calling Wet_Depo_Rainout."
          endif;enddo
          call Wet_Depo_Rainout
        endif
#endif
!------------------------------------------------------------------------------

        ! Advance time
        time = time + dt

        ! If there is any time-series output requiring evaluation at each
        ! time-step, then extract output variables from concen here
        if(Output_every_TS)then
          call Gen_Output_Vars

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines (every timestep) here
!
#ifdef WETDEPO
          do io=1,2;if(VB(io).le.verbosity_info)then      
            write(outlog(io),*)"Calling ThicknessCalculator_WetDepo."
          endif;enddo
          if(USE_WETDEPO) call ThicknessCalculator_WetDepo
#endif
!------------------------------------------------------------------------------

            ! SEE WHETHER THE ASH HAS HIT ANY AIRPORTS
          call FirstAsh

            ! Track ash on vertical profiles
          if (Write_PR_Data)then
            call Calc_vprofile(itime)
            call vprofilewriter(itime)     !write out vertical profiles
          endif
        endif

        ! GO TO OUTPUT RESULTS IF WE'RE AT THE NEXT OUTPUT STAGE
        ! Note that dt was set in Adjust_DT so that it is no larger than
        ! DT_MIN, but may be adjusted down so as to land on the next
        ! output time.  time has already been integrated forward so
        ! NextWriteTime-time should be near zero for output steps.
        if(Output_at_WriteTimes.and.(NextWriteTime-time.lt.DT_MIN))then
            ! Generate output variables if we haven't already
          if(.not.Called_Gen_Output_Vars)then
            call Gen_Output_Vars
          endif
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines (every output-step) here
!
#ifdef WETDEPO
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"Calling ThicknessCalculator_WetDepo."
          endif;enddo
          if(USE_WETDEPO) call ThicknessCalculator_WetDepo
#endif
#ifdef VARDIFF
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"Calling Prep_output_VarDiff."
          endif;enddo
          if(useVarDiffH.or.useVarDiffV) call Prep_output_VarDiff
#endif
#ifdef WETDEPO
          do io=1,2;if(VB(io).le.verbosity_info)then    
            write(outlog(io),*)"Calling Prep_output_WetDepo."
          endif;enddo
          if(USE_WETDEPO) call Prep_output_WetDepo
#endif
#ifdef SRC_RESUSP
          if(SourceType.eq.'resuspens')then
            do io=1,2;if(VB(io).le.verbosity_info)then    
              write(outlog(io),*)"Calling Prep_output_Source_Resuspension."
            endif;enddo
            call Prep_output_Source_Resuspension
          endif
#endif
#ifdef SRC_GAS
          if(SourceType.eq.'gas')then
            do io=1,2;if(VB(io).le.verbosity_info)then    
              write(outlog(io),*)"Calling Prep_output_Source_Gas."
            endif;enddo
            call Prep_output_Source_Gas
          endif
#endif
!------------------------------------------------------------------------------
          call output_results
          !if ((WriteAirportFile_ASCII.or.WriteAirportFile_KML).and. &
          if (Write_PT_Data.and. &
              (iTimeNext.lt.nWriteTimes)) then
            do j=iTimeNext,nWriteTimes
              Airport_Thickness_TS(1:nairports,j) = Airport_Thickness(1:nairports)
            enddo
          endif
        endif

        if(Output_at_logsteps)then  
          !WRITE SUMMARY INFORMATION ON MASS CONSERVATION EVERY log_step TIME STEPS
          if(mod(itime,log_step).eq.0) then
            if(.not.Called_Gen_Output_Vars)then
              call Gen_Output_Vars
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines (every log-step) here
!
#ifdef WETDEPO
              do io=1,2;if(VB(io).le.verbosity_info)then    
                write(outlog(io),*)"Calling ThicknessCalculator_WetDepo."
              endif;enddo
              if(USE_WETDEPO) call ThicknessCalculator_WetDepo
#endif
!------------------------------------------------------------------------------
            endif
            call TimeStepTotals(itime)
          endif

            ! Only consider reducing GS bins at logsteps
          if(time.gt.e_EndTime_final)then
            ! Doesn't make sense to flag bins as flushed out while the eruption is on-going
            if(.not.Called_Gen_Output_Vars)then
              call Gen_Output_Vars
#ifdef WETDEPO
              do io=1,2;if(VB(io).le.verbosity_info)then    
                write(outlog(io),*)"Calling ThicknessCalculator_WetDepo."
              endif;enddo
              if(USE_WETDEPO) call ThicknessCalculator_WetDepo
#endif
            endif
            call Prune_GS
          endif
        else
            ! If we are not monitoring deposits through logsteps, then set
            ! tot_vol to 0 and rely on the simulation ending through input
            ! duration values
          tot_vol = 0.0_ip
        endif

        if(tot_vol.gt.EPS_SMALL)then
          aloft_percent_remaining = aloft_vol/tot_vol
        else
          aloft_percent_remaining = 1.0_ip
        endif

        ! Check stop conditions
        !  If any of these is true, then the time loop will stop
           ! Stops if there is less than 1% of ash aloft in the domain
        StopConditions(1) = (aloft_percent_remaining.lt.(1.0_ip-StopValue))
           ! Normal stop condition if simulation exceeds alloted time
        StopConditions(2) = (time.ge.Simtime_in_hours)
           ! Normal stop condition when nothing is left to advect
        StopConditions(3) = (n_gs_aloft.eq.0)
        if(SourceCumulativeVol.gt.EPS_TINY)then
          MassConsErr = abs(SourceCumulativeVol-tot_vol)/SourceCumulativeVol
        endif
           ! Error stop condition if the concen and outflow do not match the source
        StopConditions(4) = (MassConsErr.gt.1.0e-3_ip)
           ! Error stop condition if any volume measure is negative
        StopConditions(5) = (dep_vol.lt.-1.0_ip*EPS_SMALL).or.&
                            (aloft_vol.lt.-1.0_ip*EPS_SMALL).or.&
                            (outflow_vol.lt.-1.0_ip*EPS_SMALL).or.&
                            (SourceCumulativeVol.lt.-1.0_ip*EPS_SMALL)

        if((CheckConditions(1).eqv..true.).and.&
           (StopConditions(1).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(2).eqv..true.).and.&
               (StopConditions(2).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(3).eqv..true.).and.&
               (StopConditions(3).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(4).eqv..true.).and.&
               (StopConditions(4).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(5).eqv..true.).and.&
               (StopConditions(5).eqv..true.))then
          StopTimeLoop = .true.
        else
          StopTimeLoop = .false.
        endif
      enddo  !loop over itime
              !  ((dep_percent_accumulated.le.StopValue).and. &
              !    (time.lt.Simtime_in_hours)        .and. &
              !    (n_gs_aloft.gt.0))

      ! Reset ntmax to the actual number of time steps
      ntmax = itime

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Time integration completed for the following reason:"
      endif;enddo
      if((CheckConditions(1).eqv..true.).and.&
         (StopConditions(1).eqv..true.))then
        ! Normal stop condition set by user tracking the deposit
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Percent accumulated/exited exceeds ",StopValue
        endif;enddo
      endif
      if((CheckConditions(2).eqv..true.).and.&
         (StopConditions(2).eqv..true.))then
        ! Normal stop condition if simulation exceeds alloted time
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"time.ge.Simtime_in_hours"
          write(outlog(io),*)"              Time = ",real(time,kind=4)
          write(outlog(io),*)"  Simtime_in_hours = ",real(Simtime_in_hours,kind=4)
        endif;enddo
      endif
      if((CheckConditions(3).eqv..true.).and.&
         (StopConditions(3).eqv..true.))then
        ! Normal stop condition when nothing is left to advect
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"No ash species remain aloft."
        endif;enddo
      endif
      if((CheckConditions(4).eqv..true.).and.&
         (StopConditions(4).eqv..true.))then
        ! Error stop condition if the concen and outflow do not match the source
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Cummulative source volume does not match aloft + outflow"
          write(errlog(io),*)" tot_vol = ",tot_vol
          write(errlog(io),*)" SourceCumulativeVol = ",SourceCumulativeVol
          write(errlog(io),*)" Abs. Error = ",&
                               abs((tot_vol-SourceCumulativeVol)/SourceCumulativeVol)
          write(errlog(io),*)" e_Volume = ",e_Volume
        endif;enddo
        stop 1
      endif
      if((CheckConditions(5).eqv..true.).and.&
         (StopConditions(5).eqv..true.))then
        ! Error stop condition if any volume measure is negative
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"One of the volume measures is negative."
          write(errlog(io),*)"        dep_vol = ",dep_vol
          write(errlog(io),*)"        aloft_vol = ",aloft_vol
          write(errlog(io),*)"        outflow_vol = ",outflow_vol
          write(errlog(io),*)"        SourceCumulativeVol = ",SourceCumulativeVol
        endif;enddo
        stop 1
      endif

      ! ************************************************************************
      ! ****** end time simulation *********************************************
      ! ************************************************************************

      isFinal_TS = .true.
      Called_Gen_Output_Vars  = .false.
      Calculated_Cloud_Load   = .false.
      Calculated_AshThickness = .false.

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),12)   !put footnotes below output table
        write(outlog(io),*)'time=',real(time,kind=4),',dt=',real(dt,kind=4)
        write(outlog(io),*)"Mass Conservation Error = ",MassConsErr
      endif;enddo

        ! Make sure we have the latest output variables and go to write routines
      call Gen_Output_Vars
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines here
!
#ifdef WETDEPO
      if(USE_WETDEPO) call ThicknessCalculator_WetDepo
#endif
!------------------------------------------------------------------------------

      call output_results
      !WRITE RESULTS TO LOG AND STANDARD OUTPUT
      !TotalTime_sp = etime(elapsed_sp)
      call cpu_time(t1) !time is a scalar real
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),3) t1-t0, time*HR_2_S
        write(outlog(io),4) time*HR_2_S/(t1-t0)
      endif;enddo
      call TimeStepTotals(itime)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),5) dep_vol
        write(outlog(io),6) tot_vol
        write(outlog(io),9) maxval(DepositThickness), DepositAreaCovered
        write(outlog(io),34)       !write out area of cloud at different thresholds
        do k=1,5
          write(outlog(io),35) LoadVal(k), CloudLoadArea(k)
        enddo
        write(outlog(io),33)    !write "normal completion"
      endif;enddo

      ! Format statements
1     format(/,5x,'Starting volume (km3 DRE)    = ',f11.4,       &
             /,5x,'maximum number of time steps = ',i8,          &
             //,21x,'Time',19x,                                  &
                '|--------------------Volume (km3 DRE)-------------------|', &
                3x,'Cloud Area',                                 &
              /,7x,'step',8x,'(hrs)',2x,'yyyymmddhh:mm',         &
                5x,'Source',6x,'Deposit',7x,'Aloft',5x,'Outflow',&
                7x,'Total',10x,'km2')

3     format(/,5x,'Execution time           = ',f15.4,' seconds',&
             /,5x,'Simulation time          = ',f15.4,' seconds')       
4     format(  5x,'Execution time/CPU time  = ',f15.4)       
5     format(  5x,'Ending deposit volume    = ',f15.4,' km3 DRE')       
6     format(  5x,'Ending total volume      = ',f15.4,' km3 DRE')       
7     format(  5x,'Building time array of plume height & eruption rate')
9     format(/,5x,'Maximum deposit thickness (mm)   = ',f10.4, &
             /,5x,'Area covered by >0.01 mm (km2)   = ',f10.1,/)
12    format(4x,'*=files written out')

33    format(/,5x,'Normal completion')
34    format(/,'  Ash load   cloud area',/, &
               '      T/km2         km2')
35    format(2f11.1)

      ! clean up memory
      call dealloc_arrays

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls deallocation routines here
!
#ifdef TOPO
      if(useTopo)                     call Deallocate_Topo
#endif
#ifdef VARDIFF
      if(useVarDiffH.or.useVarDiffV)  call Deallocate_VarDiff_Met
#endif
#ifdef LC
      if(useLandCover)                call Deallocate_LC
#endif
#ifdef WETDEPO
      if(USE_WETDEPO)then
        call Deallocate_WetDepo_global
        call Deallocate_WetDepo_Met
      endif
#endif
#ifdef SRC_SAT
      call Deallocate_Source_Satellite
#endif
#ifdef SRC_RESUSP
      call Deallocate_Source_Resuspension
#endif
#ifdef SRC_GAS
      call Deallocate_Source_Gas
#endif

!------------------------------------------------------------------------------

      close(global_log)       !close log file 

      end program Ash3d
