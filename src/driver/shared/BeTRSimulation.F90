module BeTRSimulation
  !
  ! !DESCRIPTION:
  !  BeTR simulation base class.
  !
  !  BeTR simulation class are API definitions, mapping data
  !  structures from a specific LSM, e.g. CLM, ALM, into BeTR data
  !  structures.
  !
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog, use_cn
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use tracer_varcon  , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use BeTR_decompMod , only : betr_bounds_type
  use decompMod      , only : bounds_type
#if (defined SBETR)
  use PatchType      , only : patch_type
  use ColumnType     , only : column_type
  use LandunitType   , only : landunit_type
#else
  use ColumnType     , only : column_type => column_physical_properties_type
  use VegetationType , only : patch_type  => vegetation_physical_properties_type
  use LandunitType   , only : landunit_type => landunit_physical_properties_type
#endif

  ! !USES:
  use BetrType                 , only : betr_type, create_betr_type
  use betr_ctrl                , only : max_betr_hist_type, betr_offline, max_betr_rest_type, biulog
  use betr_constants           , only : betr_string_length
  use betr_constants           , only : betr_filename_length
  use betr_regression_module   , only : betr_regression_type
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type, create_betr_biogeophys_input
  use BeTR_biogeoStateType     , only : betr_biogeo_state_type, create_betr_biogeo_state
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type, create_betr_biogeoFlux
  use BeTR_TimeMod             , only : betr_time_type
  use betr_varcon              , only : betr_maxpatch_pft
  use BetrStatusType           , only : betr_status_type, create_betr_status_type
  use BetrStatusSimType        , only : betr_status_sim_type, create_betr_status_sim_type
  use betr_columnType          , only : betr_column_type, create_betr_column_type
  use betr_patchType           , only : betr_patch_type, create_betr_patch_type
  use betr_varcon              , only : spval => bspval
  use BeTRHistVarType          , only : betr_hist_var_type
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: betr_simulation_type
     type(betr_type)                     , public, pointer  :: betr(:)  => null()
     type(betr_biogeophys_input_type)    , public, pointer  :: biophys_forc(:) => null()
     type(betr_biogeo_state_type)        , public, pointer  :: biogeo_state(:) => null()
     type(betr_biogeo_flux_type)         , public, pointer  :: biogeo_flux(:) => null()
     type(betr_status_type)              , public, pointer  :: bstatus(:) => null()
     type(betr_status_sim_type)          , public, pointer  :: bsimstatus
     logical                             , public, pointer :: active_col(:) => null()
     type(betr_column_type)              , public, pointer :: betr_col(:) => null()
     type(betr_patch_type)               , public, pointer :: betr_pft(:) => null()
     character(len=betr_filename_length) , private :: base_filename
     character(len=betr_filename_length) , private :: hist_filename

     type(betr_regression_type), private          :: regression
     type(betr_time_type), public                 :: betr_time
     integer, public                              :: num_soilp
     integer, public, allocatable                 :: filter_soilp(:)
     integer, public                              :: num_jtops
     integer, public, allocatable                 :: jtops(:)
     integer, public                              :: num_soilc
     integer, public, allocatable                 :: filter_soilc(:)

     type(betr_hist_var_type), allocatable :: state_hist1d_var(:)
     type(betr_hist_var_type), allocatable :: state_hist2d_var(:)
     type(betr_hist_var_type), allocatable :: flux_hist1d_var(:)
     type(betr_hist_var_type), allocatable :: flux_hist2d_var(:)
     integer :: num_hist_state1d
     integer :: num_hist_state2d
     integer :: num_hist_flux1d
     integer :: num_hist_flux2d

     integer :: num_rest_state1d
     integer :: num_rest_state2d
     real(r8), pointer :: hist_states_2d(:,:,:)
     real(r8), pointer :: hist_states_1d(:,:)
     real(r8), pointer :: rest_states_2d(:,:,:)
     real(r8), pointer :: rest_states_1d(:,:)
     real(r8), pointer :: hist_fluxes_1d(:,:)
     real(r8), pointer :: hist_fluxes_2d(:,:,:)
     logical,  public :: active_soibgc
     ! FIXME(bja, 201603) most of these types should be private!

     ! NOTE(bja, 201603) BeTR types only, no LSM specific types here!

   contains
     procedure, public  :: BeTRInit
     procedure, public  :: BeTRInitAllocateMemory
     procedure, public  :: BeTRInitializeBiophysForc
     procedure, public  :: BeTRInitializeCols
     procedure, public  :: BeTRInitializePFT
     procedure, public  :: BeTRInitializeActiveCols
     procedure, public  :: BeTRInitialize_1
     procedure, public  :: BeTRSetFilter
     procedure, public  :: SetClock
     procedure, public  :: InitOnline                             => BeTRSimulationInit
     procedure, public  :: Init                                   => BeTRSimulationInitOffline
     procedure, public  :: ConsistencyCheck                       => BeTRSimulationConsistencyCheck
     procedure, public  :: PreDiagSoilColWaterFlux                => BeTRSimulationPreDiagSoilColWaterFlux
     procedure, public  :: DiagnoseDtracerFreezeThaw              => BeTRSimulationDiagnoseDtracerFreezeThaw
     procedure, public  :: DiagAdvWaterFlux                       => BeTRSimulationDiagAdvWaterFlux
     procedure, public  :: DiagDrainWaterFlux                     => BeTRSimulationDiagDrainWaterFlux
     procedure, public  :: BeginSnowLayerAdjst                    => BeTRSimulationBeginTracerSnowLayerAdjst
     procedure, public  :: EndSnowLayerAdjst                      => BeTRSimulationEndTracerSnowLayerAdjst
     procedure, public  :: CombineSnowLayers                      => BeTRSimulationCombineSnowLayers
     procedure, public  :: DvideSnowLayers                        => BeTRSimulationDvideSnowLayers
     procedure, public  :: StepWithoutDrainage                    => BeTRSimulationStepWithoutDrainage
     procedure, public  :: StepWithDrainage                       => BeTRSimulationStepWithDrainage
     procedure, public  :: BeginMassBalanceCheck                  => BeTRSimulationBeginMassBalanceCheck
     procedure, public  :: MassBalanceCheck                       => BeTRSimulationMassBalanceCheck
     procedure, public  :: BeTRSetBiophysForcing                  => BeTRSimulationSetBiophysForcing
     procedure, public  :: BeTRSetBiophysForcingWaterstateVars    => BeTRSimulationSetBiophysForcingWaterstateVars
     procedure, public  :: BeTRSetBiophysForcingCarbonFluxVars    => BeTRSimulationSetBiophysForcingCarbonFluxVars
     procedure, public  :: BeTRSetBiophysForcingWaterFluxVars     => BeTRSimulationSetBiophysForcingWaterFluxVars
     procedure, public  :: BeTRSetBiophysForcingTemperatureVars   => BeTRSimulationSetBiophysForcingTemperatureVars
     procedure, public  :: BeTRSetBiophysForcingSoilHydrologyVars => BeTRSimulationSetBiophysForcingSoilHydrologyVars
     procedure, public  :: BeTRSetBiophysForcingAtm2LndVars       => BeTRSimulationSetBiophysForcingAtm2LndVars
     procedure, public  :: BeTRSetBiophysForcingCanopyStateVars   => BeTRSimulationSetBiophysForcingCanopyStateVars
     procedure, public  :: BeTRSetBiophysForcingChemStateVars     => BeTRSimulationSetBiophysForcingChemStateVars
     procedure, public  :: BeTRSetBiophysForcingSoilStateVars     => BeTRSimulationSetBiophysForcingSoilStateVars
     procedure, public  :: RetrieveBiogeoFlux                     => BeTRSimulationRetrieveBiogeoFlux
     procedure, public  :: DiagnoseLnd2atm                        => BeTRSimulationDiagnoseLnd2atm
     procedure, public  :: CreateOfflineHistory                   => hist_htapes_create
     procedure, public  :: WriteOfflineHistory                    => hist_write
     procedure, public  :: SetBaseFilename                        => BeTRSimulationSetBaseFilename

     procedure, public  :: WriteRegressionOutput
     !the following are used to interact with lsm
     procedure, public  :: do_soibgc
     procedure, public  :: BeTRRestart                            => BeTRSimulationRestart
     procedure, public  :: BeTRRestartOffline                     => BeTRSimulationRestartOffline
     procedure, public  :: BeTRRestartOpen                        => BeTRSimulationRestartOpen
     procedure, public  :: BeTRRestartClose                       => BeTRSimulationRestartClose
     procedure, public  :: BeTRCreateHistory                      => BeTRSimulationCreateHistory
     procedure, public  :: HistRetrieval                          => BeTRSimulationHistRetrieval
     procedure, private :: BeTRRetrieveHistoryState               => BeTRSimulationRetrieveHistoryState
     procedure, private :: BeTRRetrieveHistoryFlux                => BeTRSimulationRetrieveHistoryFlux
     procedure, public  :: BeTRSetcps                             => BeTRSimulationSetcps
     procedure, public  :: BeTRSetcpsCol                          => BeTRSimulationSetcpsCol
     procedure, public  :: BeTRSetcpsPFT                          => BeTRSimulationSetcpsPFT
     procedure, public  :: BeTRSetBounds                          => BeTRSimulationSetBounds
     procedure, private :: hist_create_states
     procedure, private :: hist_create_fluxes
     procedure, private :: hist_output_states
     procedure, private :: hist_output_fluxes
     procedure, public  :: RestAlloc                              => BeTRSimulationRestartAlloc
     procedure, private :: HistAlloc                              => BeTRSimulationHistoryAlloc

  end type betr_simulation_type

  public :: BeTRSimulationInit

contains

  subroutine SetClock(this, dtime, nelapstep)
  implicit none
  class(betr_simulation_type)              , intent(inout) :: this
  real(r8), intent(in) :: dtime
  integer, intent(in) :: nelapstep


  call this%betr_time%setClock(dtime, nelapstep)

  end subroutine SetClock
  !-------------------------------------------------------------------------------
  subroutine BeTRSimulationInit(this, bounds, lun, col, pft, waterstate, namelist_buffer, masterproc)
    !
    ! DESCRIPTIONS
    ! Dummy routine for inheritance purposes. don't use.
    !
    !USES
    use WaterstateType , only : waterstate_type
    use betr_constants , only : betr_namelist_buffer_size, betr_filename_length
    implicit none

    class(betr_simulation_type)              , intent(inout) :: this
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(in) :: col
    type(patch_type)                         , intent(in) :: pft
    type(bounds_type)                        , intent(in)    :: bounds
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    type(waterstate_type)                    , intent(inout) :: waterstate
    logical,                        optional , intent(in)    :: masterproc
    character(len=*), parameter :: subname = 'BeTRSimulationInit'

    call endrun(msg="ERROR "//subname//" unimplemented. "//errmsg(mod_filename, __LINE__))

    if (this%num_soilc > 0)                  continue
    if (bounds%begc > 0)                     continue
    if (size(waterstate%h2osoi_liq_col) > 0) continue
    if (len(namelist_buffer) > 0)            continue
  end subroutine BeTRSimulationInit

  !-------------------------------------------------------------------------------
  subroutine BeTRSimulationInitOffline(this, bounds, lun, col, pft, waterstate, &
       namelist_buffer, base_filename)
    !
    ! DESCRIPTIONS
    ! Dummy routine for inheritance purposes. don't use.
    !
    !USES
    use WaterstateType , only : waterstate_type
    use betr_constants , only : betr_namelist_buffer_size, betr_filename_length
    implicit none

    class(betr_simulation_type)              , intent(inout) :: this
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(in) :: col
    type(patch_type)                         , intent(in) :: pft
    type(bounds_type)                        , intent(in)    :: bounds
    type(waterstate_type)                    , intent(inout) :: waterstate
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    character(len=betr_filename_length)      , intent(in)    :: base_filename
    character(len=*), parameter :: subname = 'BeTRSimulationInit'

    call endrun(msg="ERROR "//subname//" unimplemented. "//errmsg(mod_filename, __LINE__))

    if (this%num_soilc > 0)                  continue
    if (bounds%begc > 0)                     continue
    if (size(waterstate%h2osoi_liq_col) > 0) continue
    if (len(base_filename) > 0)              continue
    if (len(namelist_buffer) > 0)            continue
    call this%betr_time%Init(namelist_buffer)
  end subroutine BeTRSimulationInitOffline
!-------------------------------------------------------------------------------
  subroutine BeTRSetFilter(this, maxpft_per_col, boffline)
  !
  !DESCRIPTION
  ! set betr filter, only used for standalone applicaitons
    use betr_ctrl                , only : betr_offline

  implicit none
  !ARGUMENTS
  class(betr_simulation_type), intent(inout) :: this
  integer, intent(in) :: maxpft_per_col
  logical, intent(in) :: boffline
    integer :: p
    !by default, surface litter layer is off
    this%num_jtops = 1
    allocate(this%jtops(this%num_jtops))
    this%jtops(:) = 1

    this%num_soilc = 1
    allocate(this%filter_soilc(this%num_soilc))
    this%filter_soilc(:) = 1

    this%num_soilp = maxpft_per_col
    allocate(this%filter_soilp(this%num_soilp))
    do p = 1, maxpft_per_col
      this%filter_soilp(p) = p
    enddo
    betr_offline=boffline
    betr_maxpatch_pft = maxpft_per_col
  end subroutine BeTRSetFilter
!-------------------------------------------------------------------------------

  subroutine BeTRInit(this, bounds, lun, col, pft, waterstate, namelist_buffer, &
     base_filename, masterproc)
    !
    ! DESCRIPTION
    ! initialize BeTR
    !
    !!USES
    use WaterStateType , only : waterstate_type
    use betr_constants , only : betr_namelist_buffer_size
    use betr_constants , only : betr_filename_length
    use betr_varcon    , only : betr_maxpatch_pft
    use landunit_varcon, only : istsoil, istcrop
    implicit none
    !ARGUMENTS
    class(betr_simulation_type)              , intent(inout) :: this
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(in) :: col
    type(patch_type)                         , intent(in) :: pft
    type(waterstate_type)                    , intent(in) :: waterstate
    character(len=betr_namelist_buffer_size) , intent(in) :: namelist_buffer
    character(len=betr_filename_length)      , optional, intent(in)    :: base_filename
    logical,                      optional   , intent(in) :: masterproc
    !TEMPORARY VARIABLES
    character(len=*), parameter :: subname = 'BeTRInit'
    type(betr_bounds_type) :: betr_bounds
    integer :: c, l
    logical :: asoibgc

    biulog = iulog

    if (present(base_filename)) then
      call this%SetBaseFilename(base_filename)
    else
      call this%SetBaseFilename()
    endif

    if(present(masterproc))then
      call this%betr_time%Init(namelist_buffer, masterproc)
    else
      call this%betr_time%Init(namelist_buffer)
    endif

    !allocate memory
    call this%BeTRInitAllocateMemory(bounds%begc, bounds%endc)

    !grid horizontal bounds
    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRInitializeBiophysForc(betr_bounds, bounds%begc, bounds%endc)
    call this%BeTRInitializeCols(       betr_bounds, bounds%begc, bounds%endc)
    call this%BeTRInitializePFT(        betr_bounds, bounds%begc, bounds%endc)
    call this%BeTRInitializeActiveCols( betr_bounds, bounds%begc, bounds%endc, &
         col%landunit, lun%itype)


    call this%BeTRSetcps(bounds, col, pft)

    call this%BeTRSetBiophysForcingWaterstateVars(&
         bounds%begc, bounds%endc,                   &
         betr_bounds%lbj, betr_bounds%ubj,           &
         waterstate%finundated_col,                  &
         waterstate%frac_h2osfc_col,                 &
         waterstate%h2osoi_liq_col,                  &
         waterstate%h2osoi_ice_col,                  &
         waterstate%h2osoi_liqvol_col,               &
         waterstate%h2osoi_icevol_col,               &
         waterstate%h2osoi_vol_col,                  &
         waterstate%air_vol_col,                     &
         waterstate%rho_vap_col,                               &
         waterstate%rhvap_soi_col,                             &
         waterstate%smp_l_col)

    call this%BeTRInitialize_1(betr_bounds, bounds%begc, bounds%endc, namelist_buffer)
    
    if (betr_offline) then
       call this%CreateOfflineHistory(bounds, &
            this%num_hist_state1d, this%num_hist_state2d,        &
            this%num_hist_flux1d, this%num_hist_flux2d)
    else
       call this%BeTRCreateHistory(bounds)
    endif

    call this%RestAlloc(bounds%begc, bounds%endc)

    if (len(trim(this%base_filename)) > 0 ) then
       call this%regression%Init(trim(this%base_filename), namelist_buffer, this%bsimstatus)
       if (this%bsimstatus%check_status()) call endrun(msg=this%bsimstatus%print_msg())
    endif

  end subroutine BeTRInit

  !---------------------------------------------------------------------------------
  subroutine BeTRInitAllocateMemory(this, begc, endc)
    !
    ! DESCRIPTION
    ! Allocates memory
    !
    class(betr_simulation_type) , intent(inout) :: this
    integer                     , intent(in)    :: begc, endc

    allocate(this%betr(        begc:endc))
    allocate(this%biophys_forc(begc:endc))
    allocate(this%biogeo_flux( begc:endc))
    allocate(this%biogeo_state(begc:endc))
    allocate(this%bstatus(     begc:endc))
    allocate(this%betr_col(    begc:endc))
    allocate(this%betr_pft(    begc:endc))
    allocate(this%active_col(  begc:endc))

    allocate(this%bsimstatus)

    call this%bsimstatus%reset()

  end subroutine BeTRInitAllocateMemory

  !---------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationRestartAlloc(this, begc, endc)
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  integer                     , intent (in)   :: begc, endc

  !identify restart variables
  call this%betr(begc)%get_restartvar_size(this%num_rest_state1d, this%num_rest_state2d)

  allocate(this%rest_states_1d(begc:endc, 1:this%num_rest_state1d))
  this%rest_states_1d(:,:)=spval
  allocate(this%rest_states_2d(begc:endc, 1:betr_nlevtrc_soil, 1:this%num_rest_state2d))
  this%rest_states_2d(:,:,:)=spval

  end subroutine BeTRSimulationRestartAlloc
  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationHistoryAlloc(this, begc, endc)
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  integer                     , intent(in)    :: begc, endc

  !state variables
  allocate(this%state_hist1d_var(this%num_hist_state1d))
  allocate(this%state_hist2d_var(this%num_hist_state2d))

  allocate(this%hist_states_2d(begc:endc, 1:betr_nlevtrc_soil, 1:this%num_hist_state2d))
  allocate(this%hist_states_1d(begc:endc, 1:this%num_hist_state1d))

  !flux variables
  allocate(this%flux_hist1d_var(this%num_hist_flux1d))
  allocate(this%flux_hist2d_var(this%num_hist_flux2d))

  allocate(this%hist_fluxes_2d(begc:endc, 1:betr_nlevtrc_soil, 1:this%num_hist_flux2d))
  allocate(this%hist_fluxes_1d(begc:endc, 1:this%num_hist_flux1d))

  end subroutine BeTRSimulationHistoryAlloc


  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationRestartOpen(this, fname, flag, ncid)
    !
    !! DESCRIPTION
    ! open restart file, note it is only used in the stand alone mode
    ! note
    ! !USES:
    use bncdio_pio, only : file_desc_t, ncd_nowrite, ncd_pio_openfile, ncd_pio_createfile
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    character(len=*), intent(in) :: fname
    type(file_desc_t)          , intent(out) :: ncid ! netcdf id
    character(len=*), intent(in) :: flag

    integer :: c
    print*,'open restart file ',trim(fname), ' for ',trim(flag)
    if(trim(flag)=='read')then
      call ncd_pio_openfile(ncid, trim(fname), ncd_nowrite)
    elseif(trim(flag)=='write')then
      call ncd_pio_createfile(ncid, trim(fname))
    endif
    print*,'creating file succeeded'

  end subroutine BeTRSimulationRestartOpen


  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationRestartClose(this, ncid)
    !
    !! DESCRIPTION
    ! initialize for restart run
    ! !USES:
    use bncdio_pio, only : file_desc_t, ncd_pio_closefile
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(file_desc_t)          , intent(inout) :: ncid ! netcdf id

    call ncd_pio_closefile(ncid)
  end subroutine BeTRSimulationRestartClose
  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithoutDrainage(this, bounds, col, pft)
  !DESCRPTION
  !interface for StepWithoutDrainage
  !
  ! USES
    use SoilStateType     , only : soilstate_type
    use WaterStateType    , only : Waterstate_Type
    use TemperatureType   , only : temperature_type
    use ChemStateType     , only : chemstate_type
    use WaterfluxType     , only : waterflux_type
    use atm2lndType       , only : atm2lnd_type
    use SoilHydrologyType , only : soilhydrology_type
    use CNCarbonFluxType  , only : carbonflux_type
    use CanopyStateType   , only : canopystate_type
    use BeTR_TimeMod      , only : betr_time_type
    use pftvarcon         , only : crop
    implicit none
  !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds ! bounds
    type(column_type)           , intent(in)    :: col ! column type
    type(patch_type)            , intent(in)    :: pft

    ! remove compiler warnings about unused dummy args
    if (this%num_soilc > 0)                           continue
    if (this%betr_time%tstep > 0)                          continue
    if (bounds%begc > 0)                              continue
    if (size(col%z) > 0)                              continue

  end subroutine BeTRSimulationStepWithoutDrainage
  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationDiagnoseLnd2atm(this, bounds,  col, lnd2atm_vars)
   !
   !DESCRIPTION
   !interface for using diagnose land fluxes to atm and river copmonents
   !
   !USES
    use MathfuncMod   , only : safe_div
    use lnd2atmType    , only : lnd2atm_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    type(column_type)           , intent(in)    :: col ! column type
    type(lnd2atm_type)          , intent(inout) :: lnd2atm_vars

    ! remove compiler warnings about unused dummy args
    if (this%num_soilc > 0) continue
    if (bounds%begc > 0)    continue
    if (size(col%z) > 0)    continue

  end subroutine BeTRSimulationDiagnoseLnd2atm
  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithDrainage(this, bounds,  col)
   !
   !DESCRIPTION
   !interface for using StepWithDrainage
   !
   !USES
    use MathfuncMod   , only : safe_div
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    type(column_type)           , intent(in)    :: col ! column type

    ! remove compiler warnings about unused dummy args
    if (this%num_soilc > 0) continue
    if (bounds%begc > 0)    continue
    if (size(col%z) > 0)    continue

  end subroutine BeTRSimulationStepWithDrainage

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationBeginMassBalanceCheck(this, bounds)
    !DESCRIPTION
    !stage tracer mass balance check
    !
    !USES
    use TracerBalanceMod, only : begin_betr_tracer_massbalance
    implicit none
    !ARGUMENTS
    class(betr_simulation_type), intent(inout)   :: this
    type(bounds_type), intent(in) :: bounds

    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj
    integer :: c

    !set lbj and ubj
    call this%BeTRSetBounds(betr_bounds)

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call begin_betr_tracer_massbalance(betr_bounds,             &
         this%betr_col(c), this%num_soilc, this%filter_soilc,     &
         this%betr(c)%tracers, this%betr(c)%tracerstates,         &
         this%betr(c)%tracerfluxes, this%bstatus(c))
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=trim(this%bsimstatus%print_msg()))
  end  subroutine BeTRSimulationBeginMassBalanceCheck
  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationMassBalanceCheck(this, bounds)
   !DESCRIPTION
   ! do tracer mass balance check
   !
   !USES
    use TracerBalanceMod , only : betr_tracer_massbalance_check
    use BeTR_TimeMod     , only : betr_time_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer :: c

    !set lbj and ubj
    call this%BeTRSetBounds(betr_bounds)

   do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call betr_tracer_massbalance_check(this%betr_time, betr_bounds,           &
         this%betr_col(c), this%num_soilc, this%filter_soilc,           &
         this%betr(c)%tracers, this%betr(c)%tracerstates,                &
         this%betr(c)%tracerfluxes, this%bstatus(c))
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
   enddo
   if(this%bsimstatus%check_status()) then
      print*,this%bsimstatus%cindex
      print*,trim(this%bsimstatus%print_msg())
      call endrun(msg=trim(this%bsimstatus%print_msg()))
   endif
  end subroutine BeTRSimulationMassBalanceCheck

!-------------------------------------------------------------------------------

  subroutine hist_htapes_create(this, bounds, &
     num_hist_state1d, num_hist_state2d, num_hist_flux1d, num_hist_flux2d)
  !
  ! DESCRIPTIONS
  ! create history file and define output variables, only for standalone applicaitons
  !
  ! USES
    use netcdf          , only : nf90_float
    use bncdio_pio       , only : file_desc_t
    use bncdio_pio       , only : ncd_pio_createfile
    use bncdio_pio       , only : ncd_pio_closefile
    use bncdio_pio       , only : ncd_enddef
    use bncdio_pio       , only : ncd_defvar
    use bncdio_pio       , only : ncd_putvar
    use bhistFileMod    , only : hist_file_create, hist_def_fld1d, hist_def_fld2d
    use betr_varcon     , only : bspval
    use betr_columnType , only : betr_column_type
    !
    !ARGUMENTS
    implicit none
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds               ! bounds
    integer           ,     intent(in)   :: num_hist_state1d
    integer           ,     intent(in)   :: num_hist_state2d
    integer           ,     intent(in)   :: num_hist_flux1d
    integer           ,     intent(in)   :: num_hist_flux2d


  !TEMPORARY VARIABLES
    integer                     :: jj, kk, c
    integer :: ncol
    type(file_desc_t)           :: ncid
    character(len=*), parameter :: subname = 'hist_htapes_create'

    c = 1
    associate(                                                     &
         ntracers          => this%betr(c)%tracers%ntracers,          &
         ngwmobile_tracers => this%betr(c)%tracers%ngwmobile_tracers, &
         is_volatile       => this%betr(c)%tracers%is_volatile,       &
         is_h2o            => this%betr(c)%tracers%is_h2o,            &
         is_isotope        => this%betr(c)%tracers%is_isotope,        &
         volatileid        => this%betr(c)%tracers%volatileid,        &
         tracernames       => this%betr(c)%tracers%tracernames        &
         )

    ncol = bounds%endc-bounds%begc + 1
    this%hist_filename = trim(this%base_filename) // '.output.nc'
    call ncd_pio_createfile(ncid, this%hist_filename)

    call hist_file_create(ncid, betr_nlevtrc_soil, ncol)

    call ncd_defvar(ncid, "ZSOI", nf90_float, &
        dim1name="ncol",dim2name="levtrc",    &
        long_name="grid center"           ,    &
        units="m", missing_value=bspval, fill_value=bspval)

    call  hist_def_fld2d(ncid, varname="QFLX_ADV", nf90_type=nf90_float,dim1name="ncol",     &
           dim2name="levtrc",long_name="advective flux / velocity", units="m/s")

    write(iulog,*)'hist_create_states'
    call this%hist_create_states(bounds, betr_nlevtrc_soil, num_hist_state1d, num_hist_state2d, ncid=ncid)

    write(iulog,*)'hist_create_fluxes'
    call this%hist_create_fluxes(bounds, betr_nlevtrc_soil, num_hist_flux1d, num_hist_flux2d, ncid=ncid)

    call ncd_enddef(ncid)
    call ncd_putvar(ncid,"ZSOI",1,this%betr_col(c)%z(1:1,1:betr_nlevtrc_soil))
    call ncd_pio_closefile(ncid)

    end associate
  end subroutine hist_htapes_create

  !-------------------------------------------------------------------------------
  subroutine hist_write(this, bounds, record, numf, filter, time_vars, velocity)
    !
    ! DESCRIPTION
    ! output hist file, only for standalone applications
    !
    ! USES
    use bncdio_pio    , only : file_desc_t
    use bncdio_pio    , only : ncd_pio_openfile_for_write
    use bncdio_pio    , only : ncd_putvar
    use bncdio_pio    , only : ncd_pio_closefile
    use BeTR_TimeMod , only : betr_time_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: record
    integer                     , intent(in)    :: numf
    integer                     , intent(in)    :: filter(:)
    type(betr_time_type)        , intent(in)    :: time_vars
    real(r8)                    , intent(in)    :: velocity(:, :)
    !TEMPORARY VARIABLES
    type(file_desc_t)           :: ncid
    integer                     :: jj
    integer                     :: c
    character(len=*), parameter :: subname='hist_write'
      c = 1
    associate(                                                   &
         ntracers          => this%betr(c)%tracers%ntracers,          &
         ngwmobile_tracers => this%betr(c)%tracers%ngwmobile_tracers, &
         is_volatile       => this%betr(c)%tracers%is_volatile,       &
         is_h2o            => this%betr(c)%tracers%is_h2o,            &
         is_isotope        => this%betr(c)%tracers%is_isotope,        &
         volatileid        => this%betr(c)%tracers%volatileid,        &
         tracernames       => this%betr(c)%tracers%tracernames        &
         )
      call ncd_pio_openfile_for_write(ncid, this%hist_filename)

      if (mod(time_vars%time, 86400._r8)==0) then
         write(iulog,*)'day', time_vars%time/86400._r8
      end if
      call ncd_putvar(ncid, "time", record, time_vars%time)

      do c = bounds%begc, bounds%endc
        call ncd_putvar(ncid, 'QFLX_ADV', record, velocity(c:c, 1:betr_nlevtrc_soil))
      enddo

      call this%HistRetrieval(bounds, numf, filter)

      call this%hist_output_states(ncid, record, bounds, numf, filter, betr_nlevtrc_soil, &
            this%num_hist_state1d, this%num_hist_state2d)


      call this%hist_output_fluxes(ncid, record, bounds, numf, filter, betr_nlevtrc_soil, &
           this%num_hist_flux1d, this%num_hist_flux2d)

      call ncd_pio_closefile(ncid)
    end associate
  end subroutine hist_write

  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationHistRetrieval(this, bounds, numf, filter)
  use tracer_varcon  , only : betr_nlevtrc_soil
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds
   integer, intent(in) :: numf
   integer, intent(in) :: filter(:)
  call this%BeTRRetrieveHistoryState(bounds, numf, filter)

  call this%BeTRRetrieveHistoryFlux(bounds, numf, filter)

  end subroutine BeTRSimulationHistRetrieval

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationConsistencyCheck(this, &
     bounds, ubj, num_soilc, filter_soilc, waterstate_vars)
  ! DESCRIPTION
  ! Do consistency check, can be overwritten for varies purpose
  !
  ! USES
    use WaterStateType, only : Waterstate_Type
  implicit none
   !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc ! number of columns in column filter_soilc
    integer                     , intent(in)    :: filter_soilc(:) ! column filter_soilc
    integer                     , intent(in)    :: ubj
    type(Waterstate_Type)       , intent(in)    :: waterstate_vars ! water state variables

    ! remove compiler warnings
    if (this%num_soilc > 0)                         continue
    if (bounds%begc > 0)                            continue
    if (ubj > 0)                                    continue
    if (num_soilc > 0)                              continue
    if (size(filter_soilc) > 0)                     continue
    if (associated(waterstate_vars%h2osoi_liq_col)) continue

  end subroutine BeTRSimulationConsistencyCheck

  !---------------------------------------------------------------------------------

  subroutine WriteRegressionOutput(this, velocity)
  !DESCRIPTION
  ! output for regression test
  !
  ! USES
    use bshr_kind_mod  , only : r8 => shr_kind_r8
    use betr_constants , only : betr_string_length
    implicit none
   !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    real(r8)                    , intent(in)    :: velocity(:, :)
  !TEMPORARY VARIABLES
    integer                           :: jj, tt, begc, endc, c
    character(len=betr_string_length) :: category
    character(len=betr_string_length) :: name
    integer :: loc_id

    ! FIXME(bja, 201603) should we output units as well...?
    begc = 1
    endc = 1
    c= 1
    if (this%regression%write_regression_output) then
       call this%regression%OpenOutput()
       ! NOTE(bja, 201603) currently we are allocating all tracer
       ! state vars all the time.
       do tt = 1, this%betr(c)%tracers%ntracers
          if (tt <= this%betr(c)%tracers%ngwmobile_tracers) then
             category = 'concentration'
             name = trim(this%betr(c)%tracers%tracernames(tt)) // '_total_aqueous_conc'
             call this%regression%WriteData(category, name, &
                  this%betr(c)%tracerstates%tracer_conc_mobile_col(begc, :, tt))
          end if
          if (tt <= this%betr(c)%tracers%nvolatile_tracers) then
             if(.not. this%betr(c)%tracers%is_volatile(tt))cycle
             category = 'pfraction'
             loc_id=this%betr(c)%tracers%volatileid(tt)
             name = trim(this%betr(c)%tracers%tracernames(tt)) // '_gas_partial_fraction'
             call this%regression%WriteData(category, name, &
                  this%betr(c)%tracerstates%tracer_P_gas_frac_col(begc, :, loc_id))
          end if
       end do

       name = 'total_gas_pressure'
       category = 'pressure'
       call this%regression%WriteData(category, name, &
            this%betr(c)%tracerstates%tracer_P_gas_col(begc, :))

       name = 'advective flux'
       category = 'velocity'
       call this%regression%WriteData(category, name, &
            velocity(begc, :))

       call this%regression%CloseOutput()
    end if
  end subroutine WriteRegressionOutput

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcing(this, bounds,  col, pft, lbj, ubj, carbonflux_vars, waterstate_vars, &
       waterflux_vars, temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
       chemstate_vars, soilstate_vars)
    !
    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !USES
    use SoilStateType     , only : soilstate_type
    use WaterStateType    , only : Waterstate_Type
    use TemperatureType   , only : temperature_type
    use ChemStateType     , only : chemstate_type
    use WaterfluxType     , only : waterflux_type
    use atm2lndType       , only : atm2lnd_type
    use SoilHydrologyType , only : soilhydrology_type
    use CNCarbonFluxType  , only : carbonflux_type
    use CanopyStateType   , only : canopystate_type
    use MathfuncMod       , only : isnan => bisnan
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout)        :: this
    type(bounds_type)           , intent(in)           :: bounds
    type(patch_type)            , intent(in)           :: pft
    type(column_type)           , intent(in)           :: col ! column type
    integer                     , intent(in)           :: lbj, ubj
    type(carbonflux_type)       , optional, intent(in) :: carbonflux_vars
    type(Waterstate_Type)       , optional, intent(in) :: Waterstate_vars
    type(waterflux_type)        , optional, intent(in) :: waterflux_vars
    type(temperature_type)      , optional, intent(in) :: temperature_vars
    type(soilhydrology_type)    , optional, intent(in) :: soilhydrology_vars
    type(atm2lnd_type)          , optional, intent(in) :: atm2lnd_vars
    type(canopystate_type)      , optional, intent(in) :: canopystate_vars
    type(chemstate_type)        , optional, intent(in) :: chemstate_vars
    type(soilstate_type)        , optional, intent(in) :: soilstate_vars

    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    if (present(waterstate_vars)) then
       call this%BeTRSetBiophysForcingWaterStateVars( &
            bounds%begc, bounds%endc,                 &
            lbj, ubj,                                 &
            waterstate_vars%finundated_col,           &
            waterstate_vars%frac_h2osfc_col,          &
            waterstate_vars%h2osoi_liq_col,           &
            waterstate_vars%h2osoi_ice_col,           &
            waterstate_vars%h2osoi_liqvol_col,        &
            waterstate_vars%h2osoi_icevol_col,        &
            waterstate_vars%h2osoi_vol_col,           &
            waterstate_vars%air_vol_col,              &
            waterstate_vars%rho_vap_col,              &
            waterstate_vars%rhvap_soi_col,            &
            waterstate_vars%smp_l_col)
    end if

    if (present(carbonflux_vars)) then
       call this%BeTRSetBiophysForcingCarbonFluxVars(  &
            bounds,                                    &
            lbj, ubj, col%pfti, col%npfts, pft%active, &
            carbonflux_vars%annsum_npp_patch,          &
            carbonflux_vars%agnpp_patch,               &
            carbonflux_vars%bgnpp_patch)
    end if

    if (present(waterflux_vars)) then
       call this%BeTRSetBiophysForcingWaterFluxVars(   &
            bounds,                                    &
            lbj, ubj, col%pfti, col%npfts, pft%active, &
            waterflux_vars%qflx_infl_col,              &
            waterflux_vars%qflx_totdrain_col,          &
            waterflux_vars%qflx_gross_evap_soil_col,   &
            waterflux_vars%qflx_gross_infl_soil_col,   &
            waterflux_vars%qflx_surf_col,              &
            waterflux_vars%qflx_dew_grnd_col,          &
            waterflux_vars%qflx_dew_snow_col,          &
            waterflux_vars%qflx_sub_snow_vol_col,      &
            waterflux_vars%qflx_sub_snow_col,          &
            waterflux_vars%qflx_h2osfc2topsoi_col,     &
            waterflux_vars%qflx_snow2topsoi_col,       &
            waterflux_vars%qflx_rootsoi_col,           &
            waterflux_vars%qflx_adv_col,               &
            waterflux_vars%qflx_drain_vr_col,          &
            waterflux_vars%qflx_tran_veg_patch,        &
            waterflux_vars%qflx_rootsoi_frac_patch)

    end if

    if (present(temperature_vars)) then
       call this%BeTRSetBiophysForcingTemperatureVars( &
            bounds,                                    &
            lbj, ubj, col%pfti, col%npfts, pft%active, &
            temperature_vars%t_soi10cm_col,            &
            temperature_vars%t_soisno_col,             &
            temperature_vars%t_veg_patch )
    end if

    if (present(soilhydrology_vars)) then
       call this%BeTRSetBiophysForcingSoilHydrologyVars( &
            bounds, lbj, ubj,                            &
            soilhydrology_vars%qcharge_col,              &
            soilhydrology_vars%fracice_col )
    end if

    if (present(canopystate_vars)) then
       call this%BeTRSetBiophysForcingCanopyStateVars( &
            bounds,                                    &
            lbj, ubj, col%pfti, col%npfts, pft%active, &
            canopystate_vars%altmax_col,               &
            canopystate_vars%altmax_lastyear_col,      &
            canopystate_vars%lbl_rsc_h2o_patch,        &
            canopystate_vars%elai_patch)
    end if

    if (present(atm2lnd_vars)) then
       call this%BeTRSetBiophysForcingAtm2LndVars( &
            bounds, lbj, ubj,                      &
            atm2lnd_vars%forc_pbot_downscaled_col, &
            atm2lnd_vars%forc_t_downscaled_col)
    end if

    if (present(chemstate_vars)) then
       call this%BeTRSetBiophysForcingChemStateVars( &
            bounds, lbj, ubj,                        &
            chemstate_vars%soil_pH)
    end if

    if (present(soilstate_vars)) then
       call this%BeTRSetBiophysForcingSoilStateVars(   &
            bounds,                                    &
            lbj, ubj, col%pfti, col%npfts, pft%active, &
            soilstate_vars%bsw_col,                    &
            soilstate_vars%watsat_col,                 &
            soilstate_vars%eff_porosity_col,           &
            soilstate_vars%cellorg_col,                &
            soilstate_vars%cellclay_col,               &
            soilstate_vars%cellsand_col,               &
            soilstate_vars%bd_col,                     &
            soilstate_vars%watfc_col,                  &
            soilstate_vars%sucsat_col,                 &
            soilstate_vars%rootfr_patch)
    end if

  end subroutine BeTRSimulationSetBiophysForcing

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcingWaterStateVars(this, begc, endc, &
       lbj, ubj, finundated_col, frac_h2osfc_col, h2osoi_liq_col,              &
       h2osoi_ice_col, h2osoi_liqvol_col, h2osoi_icevol_col,                   &
       h2osoi_vol_col, air_vol_col, rho_vap_col, rhvap_soi_col,                &
       smp_l_col)
    !
    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    integer                     , intent(in)    :: begc, endc
    integer                     , intent(in)    :: lbj, ubj
    real(r8), pointer           , intent(in)    :: finundated_col(:)
    real(r8), pointer           , intent(in)    :: frac_h2osfc_col(:)
    real(r8), pointer           , intent(in)    :: h2osoi_liq_col(:,:)
    real(r8), pointer           , intent(in)    :: h2osoi_ice_col(:,:)
    real(r8), pointer           , intent(in)    :: h2osoi_liqvol_col(:,:)
    real(r8), pointer           , intent(in)    :: h2osoi_icevol_col(:,:)
    real(r8), pointer           , intent(in)    :: h2osoi_vol_col(:,:)
    real(r8), pointer           , intent(in)    :: air_vol_col(:,:)
    real(r8), pointer           , intent(in)    :: rho_vap_col(:,:)
    real(r8), pointer           , intent(in)    :: rhvap_soi_col(:,:)
    real(r8), pointer           , intent(in)    :: smp_l_col(:,:)

    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    cc = 1

    do c = begc, endc

       if (.not. this%active_col(c)) cycle

       this%biophys_forc(c)%finundated_col    (cc)         = finundated_col    (c)
       this%biophys_forc(c)%frac_h2osfc_col   (cc)         = frac_h2osfc_col   (c)
       this%biophys_forc(c)%h2osoi_liq_col    (cc,lbj:ubj) = h2osoi_liq_col    (c,lbj:ubj)
       this%biophys_forc(c)%h2osoi_ice_col    (cc,lbj:ubj) = h2osoi_ice_col    (c,lbj:ubj)
       this%biophys_forc(c)%h2osoi_liqvol_col (cc,lbj:ubj) = h2osoi_liqvol_col (c,lbj:ubj)
       this%biophys_forc(c)%h2osoi_icevol_col (cc,lbj:ubj) = h2osoi_icevol_col (c,lbj:ubj)
       this%biophys_forc(c)%h2osoi_vol_col    (cc,lbj:ubj) = h2osoi_vol_col    (c,lbj:ubj)
       this%biophys_forc(c)%air_vol_col       (cc,lbj:ubj) = air_vol_col       (c,lbj:ubj)
       this%biophys_forc(c)%rho_vap           (cc,lbj:ubj) = rho_vap_col       (c,lbj:ubj)
       this%biophys_forc(c)%rhvap_soi         (cc,lbj:ubj) = rhvap_soi_col     (c,lbj:ubj)
       this%biophys_forc(c)%smp_l_col         (cc,lbj:ubj) = smp_l_col         (c,lbj:ubj)

    enddo

  end subroutine BeTRSimulationSetBiophysForcingWaterStateVars

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcingCarbonFluxVars(this, bounds, &
       lbj, ubj, col_pfti, col_npfts, pft_active, &
       annsum_npp_patch, agnpp_patch, bgnpp_patch)
    !
    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: lbj, ubj
    integer  , pointer          , intent(in)    :: col_pfti         (:)
    integer  , pointer          , intent(in)    :: col_npfts        (:)
    logical  , pointer          , intent(in)    :: pft_active       (:)
    real(r8) , pointer          , intent(in)    :: annsum_npp_patch (:)
    real(r8) , pointer          , intent(in)    :: agnpp_patch      (:)
    real(r8) , pointer          , intent(in)    :: bgnpp_patch      (:)
    !	
    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    cc = 1
    do c = bounds%begc, bounds%endc
       if (.not. this%active_col(c)) cycle

       npft_loc = ubound(annsum_npp_patch,1)-lbound(annsum_npp_patch,1)+1

       if (col_pfti(c) /= lbound(annsum_npp_patch,1) .and. npft_loc/=col_npfts(c)) then
          do pi = 1, betr_maxpatch_pft
             this%biophys_forc(c)%annsum_npp_patch(pi) = 0._r8
             this%biophys_forc(c)%agnpp_patch(pi)      = 0._r8
             this%biophys_forc(c)%bgnpp_patch(pi)      = 0._r8
          enddo
       else
          if (use_cn) then
             pp = 0
             do pi = 1, betr_maxpatch_pft
                if (pi <= col_npfts(c)) then
                   p = col_pfti(c) + pi - 1
                   if (pft_active(p)) then
                      pp = pp + 1
                      this%biophys_forc(c)%annsum_npp_patch(pp) = annsum_npp_patch(p)
                      this%biophys_forc(c)%agnpp_patch(pp)      = agnpp_patch(p)
                      this%biophys_forc(c)%bgnpp_patch(pp)      = bgnpp_patch(p)
                   endif
                endif
             enddo
          else
             npft_loc = ubound(annsum_npp_patch,1)-lbound(annsum_npp_patch,1)+1
             if (col_pfti(c) /= lbound(annsum_npp_patch,1) .and. npft_loc/=col_npfts(c)) then
                do pi = 1, betr_maxpatch_pft
                   this%biophys_forc(c)%annsum_npp_patch(pi) = 0._r8
                   this%biophys_forc(c)%agnpp_patch(pi)      = 0._r8
                   this%biophys_forc(c)%bgnpp_patch(pi)      = 0._r8
                enddo
             endif
          endif
       endif
    enddo

  end subroutine BeTRSimulationSetBiophysForcingCarbonFluxVars

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcingWaterFluxVars(this, bounds, &
       lbj, ubj, col_pfti, col_npfts, pft_active,                       &
       qflx_infl_col,                                                   &
       qflx_totdrain_col,                                               &
       qflx_gross_evap_soil_col,                                        &
       qflx_gross_infl_soil_col,                                        &
       qflx_surf_col,                                                   &
       qflx_dew_grnd_col,                                               &
       qflx_dew_snow_col,                                               &
       qflx_sub_snow_vol_col,                                           &
       qflx_sub_snow_col,                                               &
       qflx_h2osfc2topsoi_col,                                          &
       qflx_snow2topsoi_col,                                            &
       qflx_rootsoi_col,                                                &
       qflx_adv_col,                                                    &
       qflx_drain_vr_col,                                               &
       qflx_tran_veg_patch,                                             &
       qflx_rootsoi_frac_patch)

    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: lbj, ubj
    integer  , pointer          , intent(in)    :: col_pfti                 (:)
    integer  , pointer          , intent(in)    :: col_npfts                (:)
    logical  , pointer          , intent(in)    :: pft_active               (:)
    real(r8) , pointer          , intent(in)    :: qflx_infl_col            (:)
    real(r8) , pointer          , intent(in)    :: qflx_totdrain_col        (:)
    real(r8) , pointer          , intent(in)    :: qflx_gross_evap_soil_col (:)
    real(r8) , pointer          , intent(in)    :: qflx_gross_infl_soil_col (:)
    real(r8) , pointer          , intent(in)    :: qflx_surf_col            (:)
    real(r8) , pointer          , intent(in)    :: qflx_dew_grnd_col        (:)
    real(r8) , pointer          , intent(in)    :: qflx_dew_snow_col        (:)
    real(r8) , pointer          , intent(in)    :: qflx_sub_snow_vol_col    (:)
    real(r8) , pointer          , intent(in)    :: qflx_sub_snow_col        (:)
    real(r8) , pointer          , intent(in)    :: qflx_h2osfc2topsoi_col   (:)
    real(r8) , pointer          , intent(in)    :: qflx_snow2topsoi_col     (:)
    real(r8) , pointer          , intent(in)    :: qflx_rootsoi_col         (:,:)
    real(r8) , pointer          , intent(in)    :: qflx_adv_col             (:,:)
    real(r8) , pointer          , intent(in)    :: qflx_drain_vr_col        (:,:)
    real(r8) , pointer          , intent(in)    :: qflx_tran_veg_patch      (:)
    real(r8) , pointer          , intent(in)    :: qflx_rootsoi_frac_patch  (:,:)

    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    cc = 1
    do c = bounds%begc, bounds%endc

       if(.not. this%active_col(c))cycle

       this%biogeo_flux(c)%qflx_infl_col(cc)             = qflx_infl_col(c)
       this%biogeo_flux(c)%qflx_totdrain_col(cc)         = qflx_totdrain_col(c)
       this%biogeo_flux(c)%qflx_gross_evap_soil_col(cc)  = qflx_gross_evap_soil_col(c)
       this%biogeo_flux(c)%qflx_gross_infl_soil_col(cc)  = qflx_gross_infl_soil_col(c)
       this%biophys_forc(c)%qflx_surf_col(cc)            = qflx_surf_col(c)
       this%biophys_forc(c)%qflx_dew_grnd_col(cc)        = qflx_dew_grnd_col(c)
       this%biophys_forc(c)%qflx_dew_snow_col(cc)        = qflx_dew_snow_col(c)
       this%biophys_forc(c)%qflx_sub_snow_vol_col(cc)    = qflx_sub_snow_vol_col(c)
       this%biophys_forc(c)%qflx_sub_snow_col(cc)        = qflx_sub_snow_col(c)
       this%biophys_forc(c)%qflx_h2osfc2topsoi_col(cc)   = qflx_h2osfc2topsoi_col(c)
       this%biophys_forc(c)%qflx_snow2topsoi_col(cc)     = qflx_snow2topsoi_col(c)
       this%biophys_forc(c)%qflx_rootsoi_col(cc,lbj:ubj) = qflx_rootsoi_col(c,lbj:ubj)*1.e-3_r8

       this%biogeo_flux(c)%qflx_adv_col(cc,lbj-1:ubj)    = qflx_adv_col(c,lbj-1:ubj)
       this%biogeo_flux(c)%qflx_drain_vr_col(cc,lbj:ubj) = qflx_drain_vr_col(c,lbj:ubj)

       pp = 0
       do pi = 1, betr_maxpatch_pft
          if (pi <= col_npfts(c)) then
             p = col_pfti(c) + pi - 1
             if (pft_active(p)) then
                pp = pp + 1
                this%biophys_forc(c)%qflx_tran_veg_patch(pp)             = qflx_tran_veg_patch(p)
                this%biophys_forc(c)%qflx_rootsoi_frac_patch(pp,lbj:ubj) = qflx_rootsoi_frac_patch(p,lbj:ubj)
             endif
          endif
       enddo

    enddo

  end subroutine BeTRSimulationSetBiophysForcingWaterFluxVars

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcingTemperatureVars(this, bounds, &
       lbj, ubj, col_pfti, col_npfts, pft_active,                         &
       t_soi10cm_col, t_soisno_col, t_veg_patch )
    !
    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !USES
    !
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: lbj, ubj
    integer  , pointer          , intent(in)    :: col_pfti      (:)
    integer  , pointer          , intent(in)    :: col_npfts     (:)
    logical  , pointer          , intent(in)    :: pft_active    (:)
    real(r8) , pointer          , intent(in)    :: t_soi10cm_col (:)
    real(r8) , pointer          , intent(in)    :: t_soisno_col  (:,:)
    real(r8) , pointer          , intent(in)    :: t_veg_patch   (:)
    !
    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    cc = 1
    do c = bounds%begc, bounds%endc
       if(.not. this%active_col(c))cycle

       this%biophys_forc(c)%t_soi_10cm(cc)           = t_soi10cm_col(c)
       this%biophys_forc(c)%t_soisno_col(cc,lbj:ubj) = t_soisno_col(c,lbj:ubj)
       pp = 0
       do pi = 1, betr_maxpatch_pft
          if (pi <= col_npfts(c)) then
             p = col_pfti(c) + pi - 1
             if (pft_active(p)) then
                pp = pp + 1
                this%biophys_forc(c)%t_veg_patch(pp) = t_veg_patch(p)
             endif
          endif
       enddo

    enddo
  end subroutine BeTRSimulationSetBiophysForcingTemperatureVars

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcingSoilHydrologyVars(this, bounds, &
       lbj, ubj, qcharge_col, fracice_col)
    !
    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: lbj, ubj
    real(r8) , pointer          , intent(in)    :: qcharge_col (:)
    real(r8) , pointer          , intent(in)    :: fracice_col (:,:)
    !
    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    cc = 1
    do c = bounds%begc, bounds%endc
       if (.not. this%active_col(c)) cycle

       this%biophys_forc(c)%qflx_bot_col(cc)        = qcharge_col(c)
       this%biophys_forc(c)%fracice_col(cc,lbj:ubj) = fracice_col(c,lbj:ubj)
    enddo

  end subroutine BeTRSimulationSetBiophysForcingSoilHydrologyVars

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcingCanopyStateVars(this, bounds, &
       lbj, ubj, col_pfti, col_npfts, pft_active,          &
       altmax_col, altmax_lastyear_col, lbl_rsc_h2o_patch, &
       elai_patch)
    !
    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !
    implicit none
    !
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: lbj, ubj
    integer  , pointer          , intent(in)    :: col_pfti            (:)
    integer  , pointer          , intent(in)    :: col_npfts           (:)
    logical  , pointer          , intent(in)    :: pft_active          (:)
    real(r8) , pointer          , intent(in)    :: altmax_col          (:)
    real(r8) , pointer          , intent(in)    :: altmax_lastyear_col (:)
    real(r8) , pointer          , intent(in)    :: lbl_rsc_h2o_patch   (:)
    real(r8) , pointer          , intent(in)    :: elai_patch          (:)

    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    cc = 1
    do c = bounds%begc, bounds%endc

       if (.not. this%active_col(c)) cycle

       this%biophys_forc(c)%altmax_col(cc)          = altmax_col(c)
       this%biophys_forc(c)%altmax_lastyear_col(cc) = altmax_lastyear_col(c)

       pp = 0
       do pi = 1, betr_maxpatch_pft
          if (pi <= col_npfts(c)) then
             p = col_pfti(c) + pi - 1
             if (pft_active(p)) then
                pp = pp + 1
                this%biophys_forc(c)%lbl_rsc_h2o_patch(pp) = lbl_rsc_h2o_patch(p)
                this%biophys_forc(c)%elai_patch(pp)        = elai_patch(p)
             endif
          endif
       enddo

    enddo
  end subroutine BeTRSimulationSetBiophysForcingCanopyStateVars

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcingAtm2LndVars(this, bounds, lbj, ubj, &
       forc_pbot_downscaled_col, forc_t_downscaled_col )
    !
    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: lbj, ubj
    real(r8) , pointer          , intent(in)    :: forc_pbot_downscaled_col (:)
    real(r8) , pointer          , intent(in)    :: forc_t_downscaled_col    (:)
    !
    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    cc = 1
    do c = bounds%begc, bounds%endc

       if (.not. this%active_col(c)) cycle

       this%biophys_forc(c)%forc_pbot_downscaled_col(cc) = forc_pbot_downscaled_col(c)
       this%biophys_forc(c)%forc_t_downscaled_col(cc)    = forc_t_downscaled_col(c)

    enddo

  end subroutine BeTRSimulationSetBiophysForcingAtm2LndVars

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcingChemStateVars(this, bounds, lbj, ubj, soil_pH)
    !
    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !
    implicit none
    !ARGUMENTS
    !
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: lbj, ubj
    real(r8) , pointer          , intent(in)    :: soil_pH (:,:)

    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    cc = 1
    do c = bounds%begc, bounds%endc

       if (.not. this%active_col(c)) cycle

       this%biophys_forc(c)%soil_pH(cc,lbj:ubj) = soil_pH(c,lbj:ubj)

    enddo

  end subroutine BeTRSimulationSetBiophysForcingChemStateVars

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcingSoilStateVars(this, bounds, &
       lbj, ubj, col_pfti, col_npfts, pft_active,                       &
       bsw_col,                                                         &
       watsat_col,                                                      &
       eff_porosity_col,                                                &
       cellorg_col,                                                     &
       cellclay_col,                                                    &
       cellsand_col,                                                    &
       bd_col,                                                          &
       watfc_col,                                                       &
       sucsat_col,                                                      &
       rootfr_patch)
    !
    !DESCRIPTION
    !pass in biogeophysical variables for running betr
    !
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: lbj, ubj
    integer  , pointer          , intent(in)    :: col_pfti         (:)
    integer  , pointer          , intent(in)    :: col_npfts        (:)
    logical  , pointer          , intent(in)    :: pft_active       (:)
    real(r8) , pointer          , intent(in)    :: bsw_col          (:,:)
    real(r8) , pointer          , intent(in)    :: watsat_col       (:,:)
    real(r8) , pointer          , intent(in)    :: eff_porosity_col (:,:)
    real(r8) , pointer          , intent(in)    :: cellorg_col      (:,:)
    real(r8) , pointer          , intent(in)    :: cellclay_col     (:,:)
    real(r8) , pointer          , intent(in)    :: cellsand_col     (:,:)
    real(r8) , pointer          , intent(in)    :: bd_col           (:,:)
    real(r8) , pointer          , intent(in)    :: watfc_col        (:,:)
    real(r8) , pointer          , intent(in)    :: sucsat_col       (:,:)
    real(r8) , pointer          , intent(in)    :: rootfr_patch     (:,:)

    !TEMPORARY VARIABLES
    integer :: p, pi, cc, c, l, pp
    integer :: npft_loc

    cc = 1
    do c = bounds%begc, bounds%endc

       if (.not. this%active_col(c)) cycle

       this%biophys_forc(c)%bsw_col(cc,lbj:ubj)          = bsw_col(c,lbj:ubj)
       this%biophys_forc(c)%watsat_col(cc,lbj:ubj)       = watsat_col(c,lbj:ubj)
       this%biophys_forc(c)%eff_porosity_col(cc,lbj:ubj) = eff_porosity_col(c,lbj:ubj)
       this%biophys_forc(c)%cellorg_col(cc,lbj:ubj)      = cellorg_col(c,lbj:ubj)
       this%biophys_forc(c)%cellclay_col(cc,lbj:ubj)     = cellclay_col(c,lbj:ubj)
       this%biophys_forc(c)%cellsand_col(cc,lbj:ubj)     = cellsand_col(c,lbj:ubj)
       this%biophys_forc(c)%bd_col(cc,lbj:ubj)           = bd_col(c,lbj:ubj)
       this%biophys_forc(c)%watfc_col(cc,lbj:ubj)        = watfc_col(c,lbj:ubj)
       this%biophys_forc(c)%sucsat_col(cc,lbj:ubj)       = sucsat_col(c,lbj:ubj)

       pp = 0
       do pi = 1, betr_maxpatch_pft
          if (pi <= col_npfts(c)) then
             p = col_pfti(c) + pi - 1
             if (pft_active(p)) then
                pp = pp + 1
                this%biophys_forc(c)%rootfr_patch(pp,lbj:ubj) = rootfr_patch(p,lbj:ubj)
             endif
          endif
       enddo

    enddo

  end subroutine BeTRSimulationSetBiophysForcingSoilStateVars

  !------------------------------------------------------------------------
  subroutine BeTRSimulationRetrieveBiogeoFlux(this, bounds, lbj,ubj, carbonflux_vars,  &
    waterflux_vars)
  ! DESCRIPTIONS
  ! update and return fluxes, this eventually will be expanded to
  ! include other fluxes
  ! USES
    use WaterfluxType    , only : waterflux_type
    use CNCarbonFluxType , only : carbonflux_type
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout)           :: this
  type(bounds_type)           , intent(in)              :: bounds
  integer                     , intent(in)              :: lbj, ubj
  type(carbonflux_type)       , optional, intent(inout) :: carbonflux_vars
  type(waterflux_type)        , optional, intent(inout) :: waterflux_vars

  integer :: begp, begc, endp, endc
  integer :: p, c, cc
  cc = 1
  if(present(carbonflux_vars))then
    !do nothing
  endif
  if(present(waterflux_vars))then
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      waterflux_vars%qflx_infl_col(c)            = this%biogeo_flux(c)%qflx_infl_col(cc)
      waterflux_vars%qflx_adv_col(c,lbj-1:ubj)   = this%biogeo_flux(c)%qflx_adv_col(cc,lbj-1:ubj)
      waterflux_vars%qflx_totdrain_col(c)        = this%biogeo_flux(c)%qflx_totdrain_col(cc)
      waterflux_vars%qflx_gross_evap_soil_col(c) = this%biogeo_flux(c)%qflx_gross_evap_soil_col(cc)
      waterflux_vars%qflx_gross_infl_soil_col(c) = this%biogeo_flux(c)%qflx_gross_infl_soil_col(cc)
      waterflux_vars%qflx_drain_vr_col(c,1:ubj)  = this%biogeo_flux(c)%qflx_drain_vr_col(cc,1:ubj)
    enddo
  endif

  end subroutine BeTRSimulationRetrieveBiogeoFlux

  !------------------------------------------------------------------------
  subroutine BeTRSimulationPreDiagSoilColWaterFlux(this, num_nolakec, filter_nolakec)
  !DESCRIPTION
  !prepare for water flux diagnosis. it is called before diagnosing the advective fluxes and applying
  !freeze-thaw tracer partition.
  !
  !USES

  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
   integer                    , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
   integer                    , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points

  !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c


   call this%BeTRSetBounds(betr_bounds)

   do fc= 1, num_nolakec
     c = filter_nolakec(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%pre_diagnose_soilcol_water_flux(betr_bounds, this%num_soilc, &
       this%filter_soilc, this%biophys_forc(c))
   enddo
  end subroutine BeTRSimulationPreDiagSoilColWaterFlux

  !------------------------------------------------------------------------
  subroutine BeTRSimulationDiagnoseDtracerFreezeThaw(this, bounds, num_nolakec, filter_nolakec)
  !
  ! DESCRIPTION
  ! aqueous tracer partition based on freeze-thaw
  !
  ! USES
  use WaterStateType        , only : waterstate_type
  implicit none
  !
  ! Arguments
  class(betr_simulation_type), intent(inout)   :: this
  type(bounds_type)     , intent(in) :: bounds
  integer               , intent(in) :: num_nolakec                        ! number of column non-lake points in column filter
  integer               , intent(in) :: filter_nolakec(:)                  ! column filter for non-lake points

  !temporary variables
  type(betr_bounds_type)     :: betr_bounds
  integer :: fc, c

  call this%BeTRSetBounds(betr_bounds)

  !call this%BeTRSetcpsCol(bounds, col)

  do fc = 1, num_nolakec
    c = filter_nolakec(fc)
    if(.not. this%active_col(c))cycle
    call this%betr(c)%diagnose_dtracer_freeze_thaw(betr_bounds, this%num_soilc, this%filter_soilc,  &
      this%biophys_forc(c))
  enddo
  end subroutine BeTRSimulationDiagnoseDtracerFreezeThaw

  !------------------------------------------------------------------------
  subroutine BeTRSimulationDiagAdvWaterFlux(this, num_hydrologyc, &
    filter_hydrologyc)

  !DESCRIPTION
  ! diagnose water fluxes for tracer advection
  !
  ! USES
  !
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_hydrologyc                        ! number of column non-lake points in column filter
   integer                     , intent(in)    :: filter_hydrologyc(:)                  ! column filter for non-lake points

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c, j

   call this%BeTRSetBounds(betr_bounds)
   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%diagnose_advect_water_flux(this%betr_time,              &
       betr_bounds, this%num_soilc, this%filter_soilc,                         &
       this%biophys_forc(c), this%biogeo_flux(c))
   enddo

  end subroutine BeTRSimulationDiagAdvWaterFlux
  !------------------------------------------------------------------------
  subroutine BeTRSimulationDiagDrainWaterFlux(this, num_hydrologyc, filter_hydrologyc)
  !DESCRIPTION
  ! diagnose water fluxes due to subsurface drainage
  !
  ! USES
  !
    use WaterfluxType     , only : waterflux_type
    use WaterStateType    , only : Waterstate_Type
    use SoilHydrologyType , only : soilhydrology_type
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_hydrologyc                        ! number of column non-lake points in column filter
   integer                     , intent(in)    :: filter_hydrologyc(:)                  ! column filter for non-lake points

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%diagnose_drainage_water_flux(this%betr_time, &
       betr_bounds, this%num_soilc, this%filter_soilc,      &
       this%biophys_forc(c), this%biogeo_flux(c))
  enddo
  end subroutine BeTRSimulationDiagDrainWaterFlux
  !------------------------------------------------------------------------
  subroutine BeTRSimulationBeginTracerSnowLayerAdjst(this,  num_snowc, filter_snowc)
  !DESCRIPTION
  !prepare for tracer adjustment in snow layers
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_snowc
   integer                     , intent(in)    :: filter_snowc(:)

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%Enter_tracer_LayerAdjustment(betr_bounds, this%betr_col(c), &
       this%num_soilc, this%filter_soilc)
   enddo
  end subroutine BeTRSimulationBeginTracerSnowLayerAdjst
  !------------------------------------------------------------------------
  subroutine BeTRSimulationEndTracerSnowLayerAdjst(this, num_snowc, filter_snowc)
  !DESCRIPTION
  !wrap up tracer adjustment in snow layers
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_snowc
   integer                     , intent(in)    :: filter_snowc(:)

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%Exit_tracer_LayerAdjustment(betr_bounds, this%betr_col(c), &
       this%num_soilc, this%filter_soilc)
   enddo

  end subroutine BeTRSimulationEndTracerSnowLayerAdjst
  !------------------------------------------------------------------------
  subroutine BeTRSimulationDvideSnowLayers(this, bounds, num_snowc, filter_snowc, divide_matrix)
  !DESCRIPTIONS
  !redistribute tracer in snow layers due to division
  !
  !USES
  use clm_varpar, only : nlevsno
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds               ! bounds
   integer                     , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                     , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)                    , intent(in)    :: divide_matrix(bounds%begc:bounds%endc , 1:nlevsno , 1:nlevsno )

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%tracer_DivideSnowLayers(betr_bounds, this%betr_col(c),this%num_soilc, &
       this%filter_soilc, divide_matrix(c:c,:,:), this%bstatus(c))
     if(this%bstatus(c)%check_status())then
       call this%bsimstatus%setcol(c)
       call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
       exit
     endif
   enddo
  if(this%bsimstatus%check_status()) &
    call endrun(msg=trim(this%bsimstatus%print_msg()))
  end subroutine BeTRSimulationDvideSnowLayers

  !------------------------------------------------------------------------
  subroutine BeTRSimulationCombineSnowLayers(this, bounds, num_snowc, filter_snowc, combine_matrix)
  !DESCRIPTIONS
  !redistribute tracer in snow layers due to division
  !
  !USES
  use clm_varpar, only : nlevsno
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds               ! bounds
   integer                     , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                     , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)                    , intent(in)    :: combine_matrix(bounds%begc:bounds%endc,-nlevsno+1:1 ,-nlevsno+1:1 )

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%tracer_CombineSnowLayers(betr_bounds, this%betr_col(c),this%num_soilc,&
       this%filter_soilc, combine_matrix(c:c,:,:),this%bstatus(c))
     if(this%bstatus(c)%check_status())then
       call this%bsimstatus%setcol(c)
       call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
       exit
     endif
   enddo
   if(this%bsimstatus%check_status()) &
    call endrun(msg=trim(this%bsimstatus%print_msg()))
  end subroutine BeTRSimulationCombineSnowLayers

  !------------------------------------------------------------------------
  subroutine hist_create_fluxes(this, bounds, betr_nlevtrc_soil, num_flux1d, num_flux2d, ncid)
  !
  !DESCRIPTION
  !create history file for betr fluxes
  !
  use histFileMod         , only: hist_addfld1d, hist_addfld2d
  use bhistFileMod        , only : hist_def_fld2d , hist_def_fld1d
  use bncdio_pio           , only : file_desc_t, ncd_float
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: betr_nlevtrc_soil
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d
  type(file_desc_t) ,   optional,  intent(inout)   :: ncid
  !local variables
  integer :: jj, begc, endc

  character(len=*), parameter :: subname = 'hist_create_fluxes'

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  if(betr_offline .and. (.not. present(ncid)))then
    call endrun(msg="ncid not defined in "//subname//errmsg(mod_filename, __LINE__))
  endif

  begc=bounds%begc; endc=bounds%endc

  do jj = 1, num_flux2d

    if(betr_offline)then

      call hist_def_fld2d (ncid, varname=this%flux_hist2d_var(jj)%varname, &
            nf90_type=ncd_float, dim1name = "ncol",&
            dim2name="levtrc", long_name=this%flux_hist2d_var(jj)%long_name, &
            units=this%flux_hist2d_var(jj)%units)
    else
      this%hist_fluxes_2d(begc:endc,1:betr_nlevtrc_soil, jj) = spval
      data2dptr => this%hist_fluxes_2d(begc:endc,1:betr_nlevtrc_soil, jj)
      call hist_addfld2d (fname=this%flux_hist2d_var(jj)%varname, &
           units=this%flux_hist1d_var(jj)%units, type2d="levtrc",  &
           avgflag=this%flux_hist2d_var(jj)%avg_flag, &
           long_name=this%flux_hist2d_var(jj)%long_name,  ptr_col=data2dptr, &
           default=this%flux_hist2d_var(jj)%use_default)
    endif
  enddo

  do jj = 1, num_flux1d

    if(betr_offline)then

      call hist_def_fld1d (ncid, varname=this%flux_hist1d_var(jj)%varname, &
        nf90_type=ncd_float, &
        dim1name="ncol", long_name=this%flux_hist1d_var(jj)%long_name,&
        units=this%flux_hist1d_var(jj)%units)
    else
      this%hist_fluxes_1d(begc:endc,jj) = spval
      data1dptr => this%hist_fluxes_1d(begc:endc, jj)
      call hist_addfld1d (fname=this%flux_hist1d_var(jj)%varname, &
        units=this%flux_hist1d_var(jj)%units,  &
        avgflag=this%flux_hist1d_var(jj)%avg_flag, &
        long_name=this%flux_hist1d_var(jj)%long_name, &
        ptr_col=data1dptr, default=this%flux_hist1d_var(jj)%use_default)
    endif

  enddo

  end subroutine hist_create_fluxes

  !------------------------------------------------------------------------
  subroutine hist_create_states(this, bounds, betr_nlevtrc_soil, num_state1d, num_state2d, ncid)
  !
  !create history file for betr states variables
  use histFileMod   , only: hist_addfld1d, hist_addfld2d
  use bhistFileMod        , only : hist_def_fld2d , hist_def_fld1d
  use bncdio_pio , only : file_desc_t, ncd_float
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer, intent(in) :: betr_nlevtrc_soil
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  type(file_desc_t) ,   optional,  intent(inout)   :: ncid
  !local variables
  integer :: begc, endc
  integer :: jj

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  character(len=*), parameter :: subname = 'hist_create_states'


  if(betr_offline .and. (.not. present(ncid)))then
    call endrun(msg="ncid not defined in "//subname//errmsg(mod_filename, __LINE__))
  endif

  begc = bounds%begc; endc = bounds%endc

  do jj = 1, num_state2d
    !read namelist

    if(betr_offline)then

      call hist_def_fld2d (ncid=ncid, varname=this%state_hist2d_var(jj)%varname, &
          nf90_type=ncd_float, dim1name = "ncol",&
          dim2name="levtrc", long_name=this%state_hist2d_var(jj)%long_name,&
          units=this%state_hist2d_var(jj)%units)
    else
      this%hist_states_2d(begc:endc,1:betr_nlevtrc_soil, jj) = spval
      data2dptr => this%hist_states_2d(begc:endc,1:betr_nlevtrc_soil, jj)
      call hist_addfld2d (fname=this%state_hist2d_var(jj)%varname, &
           units=this%state_hist2d_var(jj)%units, type2d="levtrc",  &
           avgflag=this%state_hist2d_var(jj)%avg_flag, &
           long_name=this%state_hist2d_var(jj)%long_name, &
           ptr_col=data2dptr, default=this%state_hist2d_var(jj)%use_default)
    endif
  enddo

  do jj = 1, num_state1d

    if(betr_offline)then

      call hist_def_fld1d (ncid, varname=this%state_hist1d_var(jj)%varname,  nf90_type=ncd_float, &
        dim1name="ncol", long_name=this%state_hist1d_var(jj)%long_name, units=this%state_hist1d_var(jj)%units)
    else
      this%hist_states_1d(begc:endc,jj) = spval
      data1dptr => this%hist_states_1d(begc:endc,jj)
      call hist_addfld1d (fname=this%state_hist1d_var(jj)%varname, &
          units=this%state_hist1d_var(jj)%units,      &
          avgflag=this%state_hist1d_var(jj)%avg_flag, &
          long_name=this%state_hist1d_var(jj)%long_name, &
          ptr_col=data1dptr, default=this%state_hist1d_var(jj)%use_default)
    endif
  enddo

  end subroutine hist_create_states
  !------------------------------------------------------------------------

  subroutine hist_output_fluxes(this,  ncid, record, bounds, numf, filter, &
     betr_nlevtrc_soil, num_flux1d, num_flux2d)
  !
  !DESCRIPTION
  !create history file for betr fluxes
  !
  use histFileMod         , only: hist_addfld1d, hist_addfld2d
  use bncdio_pio           , only :   file_desc_t, ncd_putvar
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: betr_nlevtrc_soil
  integer, intent(in) :: record
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d
  type(file_desc_t) ,     intent(inout)   :: ncid
  !local variables
  integer :: jj, begc, endc

  character(len=*), parameter :: subname = 'hist_output_fluxes'

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays


  begc=bounds%begc; endc=bounds%endc

  do jj = 1, num_flux2d

    data2dptr => this%hist_fluxes_2d(begc:endc,1:betr_nlevtrc_soil, jj)
    call ncd_putvar(ncid, this%flux_hist2d_var(jj)%varname, record, data2dptr)

  enddo

  do jj = 1, num_flux1d

    data1dptr => this%hist_fluxes_1d(begc:endc, jj)
    call ncd_putvar(ncid,this%flux_hist1d_var(jj)%varname, record, data1dptr)

  enddo

  end subroutine hist_output_fluxes
  !------------------------------------------------------------------------
  subroutine hist_output_states(this,  ncid,  record, bounds, numf, filter, &
     betr_nlevtrc_soil, num_state1d, num_state2d)
  !
  !create history file for betr states variables
  use histFileMod   , only: hist_addfld1d, hist_addfld2d
  use bncdio_pio , only : file_desc_t, ncd_putvar
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: record
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)
  integer, intent(in) :: betr_nlevtrc_soil
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  type(file_desc_t) ,     intent(inout)   :: ncid
  !local variables
  integer :: begc, endc
  integer :: jj

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  character(len=*), parameter :: subname = 'hist_output_states'


  begc = bounds%begc; endc = bounds%endc


  do jj = 1, num_state2d

    data2dptr => this%hist_states_2d(begc:endc,1:betr_nlevtrc_soil, jj)

    call ncd_putvar(ncid,this%state_hist2d_var(jj)%varname, record, data2dptr)

  enddo

  do jj = 1, num_state1d

    data1dptr => this%hist_states_1d(begc:endc,jj)

    call ncd_putvar(ncid,this%state_hist1d_var(jj)%varname, record, data1dptr)

  enddo

  end subroutine hist_output_states
  !------------------------------------------------------------------------
  subroutine BeTRSimulationCreateHistory(this, bounds)
  !
  !links the variable to output
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds               ! bounds

  call this%hist_create_states(bounds, betr_nlevtrc_soil, this%num_hist_state1d, this%num_hist_state2d)

  call this%hist_create_fluxes(bounds, betr_nlevtrc_soil, this%num_hist_flux1d, this%num_hist_flux2d)

  end subroutine BeTRSimulationCreateHistory
  !------------------------------------------------------------------------
  subroutine BeTRSimulationRetrieveHistoryState(this, bounds, numf, filter)
  use tracer_varcon  , only : betr_nlevtrc_soil
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !TEMPORARY VARIABLES
  integer :: fc, c
  type(betr_bounds_type)     :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)

  do fc = 1, numf
    c = filter(fc)
    if(.not. this%active_col(c))cycle
    call this%betr(c)%HistRetrieveState(betr_bounds, 1, betr_nlevtrc_soil, &
       this%num_hist_state1d, this%num_hist_state2d,&
       this%hist_states_1d(c:c,:), this%hist_states_2d(c:c,1:betr_nlevtrc_soil,:))
  enddo

  end subroutine BeTRSimulationRetrieveHistoryState
  !------------------------------------------------------------------------
  subroutine BeTRSimulationRetrieveHistoryFlux(this, bounds, numf, filter)
  use tracer_varcon  , only :  betr_nlevtrc_soil
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !TEMPORARY VARIABLES
  integer :: fc, c
  type(betr_bounds_type)     :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)
  this%hist_fluxes_1d(:,:)=spval
  this%hist_fluxes_2d(:,:,:)=spval
  do fc = 1, numf
    c = filter(fc)
    if(.not. this%active_col(c))cycle
    call this%betr(c)%HistRetrieveFlux(betr_bounds, 1, betr_nlevtrc_soil, &
       this%num_hist_flux1d,this%num_hist_flux2d, &
       this%hist_fluxes_1d(c:c,:),this%hist_fluxes_2d(c:c,1:betr_nlevtrc_soil,:))
  enddo

  end subroutine BeTRSimulationRetrieveHistoryFlux
  !------------------------------------------------------------------------
  subroutine BeTRSimulationRestartOffline(this, bounds, ncid, numf, filter, flag)
  !DESCRIPTION
  !create or read restart file
  use restUtilMod    , only : restartvar
  use bncdio_pio      , only : file_desc_t,ncd_double
  use bncdio_pio      , only : ncd_defvar, ncd_defdim
  use bncdio_pio      , only : ncd_enddef, ncd_putvar
  use bncdio_pio      , only : ncd_getvar
  implicit none
  ! !ARGUMENTS:
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)    , intent(in)    :: bounds
  type(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
  character(len=*)     , intent(in)    :: flag ! 'read' or 'write'
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !local variables

  integer :: c, jj, fc
  character(len=255),allocatable :: rest_varname_1d(:)
  character(len=255),allocatable :: rest_varname_2d(:)
  logical :: readvar      ! determine if variable is on initial file
  real(r8), pointer :: ptr1d(:)
  real(r8), pointer :: ptr2d(:,:)
  type(betr_bounds_type)     :: betr_bounds
  integer :: recordDimID


  allocate(rest_varname_1d(this%num_rest_state1d)); rest_varname_1d=''
  allocate(rest_varname_2d(this%num_rest_state2d)); rest_varname_2d=''

  c = bounds%begc
  call this%betr(c)%get_restartvar_info(this%num_rest_state1d, &
    this%num_rest_state2d,rest_varname_1d, rest_varname_2d)

  !assign initial conditions
  call this%BeTRSetBounds(betr_bounds)

  ! print*,'offline restart', flag
  if(flag=='define')then
    ! print*,'define restart file'
    ! define the dimensions
    !the temporal dimension is infinite

    !number of vertical layers
    call ncd_defdim(ncid, 'levtrc', betr_nlevsoi, recordDimID)

    !number of columns
    call ncd_defdim(ncid, 'column', this%num_soilc, recordDimID)

    !define the time dimension
    call ncd_defvar(ncid, 'time',ncd_double, long_name='', &
         units = '',  missing_value=spval, fill_value=spval)

    do jj = 1, this%num_rest_state1d
        !x print*,jj,trim(rest_varname_1d(jj))
      call ncd_defvar(ncid, trim(rest_varname_1d(jj)),ncd_double,dim1name='column',  &
          long_name='', units = '',  missing_value=spval, fill_value=spval)
    enddo

    do jj =1, this%num_rest_state2d
      !x print*,jj,trim(rest_varname_2d(jj))
      call ncd_defvar(ncid, trim(rest_varname_2d(jj)),ncd_double,dim1name='column',  &
        dim2name='levtrc', long_name='', units = '',  missing_value=spval, fill_value=spval)
    enddo
    call ncd_enddef(ncid)

  elseif(flag=='write')then

    do fc = 1, numf
      c = filter(fc)
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
        this%num_rest_state1d, this%num_rest_state2d, &
        this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo

    ! print*,'write restart file'
    do jj = 1, this%num_rest_state1d
       ptr1d => this%rest_states_1d(:, jj)
       call ncd_putvar(ncid, trim(rest_varname_1d(jj)), 1, ptr1d)
    enddo

    do jj = 1, this%num_rest_state2d
      ptr2d => this%rest_states_2d(:, :, jj)
      call ncd_putvar(ncid, trim(rest_varname_2d(jj)), 1, ptr2d)
    enddo
  elseif(flag=='read')then
      ! print*,'read restart file'
    do jj = 1, this%num_rest_state1d
       ptr1d => this%rest_states_1d(:, jj)
       call ncd_getvar(ncid, trim(rest_varname_1d(jj)), ptr1d)
    enddo

    do jj = 1, this%num_rest_state2d
      ptr2d => this%rest_states_2d(:, :, jj)
      call ncd_getvar(ncid, trim(rest_varname_2d(jj)), ptr2d)
    enddo

    ! print*,'assign values to state variables',flag
    do fc = 1, numf
      c = filter(fc)
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
      this%num_rest_state1d,this%num_rest_state2d, &
      this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo
  endif

  deallocate(rest_varname_1d)
  deallocate(rest_varname_2d)

  end subroutine BeTRSimulationRestartOffline

  !------------------------------------------------------------------------
  subroutine BeTRSimulationRestart(this, bounds, ncid, flag)
  !DESCRIPTION
  !create or read restart file
  use restUtilMod    , only : restartvar
  use ncdio_pio      , only : file_desc_t,ncd_double
  implicit none
  ! !ARGUMENTS:
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)    , intent(in)    :: bounds
  type(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
  character(len=*)     , intent(in)    :: flag ! 'read' or 'write'

  !local variables
  integer :: nrest_1d, nrest_2d
  integer :: c, jj, fc
  character(len=255), allocatable :: rest_varname_1d(:)
  character(len=255), allocatable :: rest_varname_2d(:)
  logical :: readvar      ! determine if variable is on initial file
  real(r8), pointer :: ptr1d(:)
  real(r8), pointer :: ptr2d(:,:)
  type(betr_bounds_type)     :: betr_bounds
  integer :: recordDimID

  c = bounds%begc

  allocate(rest_varname_1d(this%num_rest_state1d)); rest_varname_1d=''
  allocate(rest_varname_2d(this%num_rest_state2d)); rest_varname_2d=''

  c = bounds%begc
  call this%betr(c)%get_restartvar_info(this%num_rest_state1d, &
    this%num_rest_state2d,rest_varname_1d, rest_varname_2d)

  call this%BeTRSetBounds(betr_bounds)
  if(trim(flag)=='write')then
    do c = bounds%begc, bounds%endc
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
        this%num_rest_state1d,this%num_rest_state2d, &
        this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo
  endif

  do jj = 1, this%num_rest_state1d
    ptr1d => this%rest_states_1d(:, jj)
    call restartvar(ncid=ncid, flag=flag, varname=trim(rest_varname_1d(jj)), &
       xtype=ncd_double,  dim1name='column', long_name='',  units='', &
       interpinic_flag='interp' , readvar=readvar, data=ptr1d)
  enddo
  do jj = 1, this%num_rest_state2d
    ptr2d => this%rest_states_2d(:, :, jj)
    call restartvar(ncid=ncid, flag=flag, varname=trim(rest_varname_2d(jj)), xtype=ncd_double,  &
      dim1name='column',dim2name='levtrc', switchdim=.true., &
      long_name='',  units='', interpinic_flag='interp',readvar=readvar, data=ptr2d)
  enddo

  if(trim(flag)=='read')then
    !assign initial conditions
    do c = bounds%begc, bounds%endc
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
        this%num_rest_state1d,this%num_rest_state2d, &
        this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo
  endif

  deallocate(rest_varname_1d)
  deallocate(rest_varname_2d)
  end subroutine BeTRSimulationRestart
  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetcps(this, bounds, col, pft)
  !
  !DESCRIPTION
  ! set up columns
  !USES
  use decompMod             , only : bounds_type
  use pftvarcon             , only : noveg, crop
  use tracer_varcon         , only : betr_nlevsoi
  !ARGUMENTS
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type), intent(in) :: bounds
  type(column_type), intent(in) :: col
  type(patch_type), optional, intent(in) :: pft

  call this%BeTRSetcpsCol(bounds%begc, bounds%endc, &
       col%snl, col%zi, col%dz, col%z, col%pfti, col%pftf, col%npfts)

  if (present(pft)) then
     call this%BeTRSetcpsPFT(bounds%begc, bounds%endc, &
          col%pfti, col%npfts, pft%active, pft%itype, pft%wtcol)
  endif

  end subroutine BeTRSimulationSetcps

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetcpsCol(this, begc, endc, snl, zi, dz, z, pfti, pftf, npfts)
  !
  !DESCRIPTION
  ! Sets up BeTR columns
  !USES
  use tracer_varcon         , only : betr_nlevsoi
  !ARGUMENTS
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  integer                     , intent(in)    :: begc, endc
  integer, pointer            , intent(in)    :: snl(:)
  real(r8), pointer           , intent(in)    :: zi(:,:)
  real(r8), pointer           , intent(in)    :: dz(:,:)
  real(r8), pointer           , intent(in)    :: z(:,:)
  integer, pointer            , intent(in)    :: pfti(:)
  integer, pointer            , intent(in)    :: pftf(:)
  integer, pointer            , intent(in)    :: npfts(:)
  !
  ! !LOCAL VARIABLES:
  integer :: c

  do c = begc, endc

    this%betr_col(c)%snl(1)               = snl(c)
    this%betr_col(c)%zi(1,0:betr_nlevsoi) = zi(c,0:betr_nlevsoi)
    this%betr_col(c)%dz(1,1:betr_nlevsoi) = dz(c,1:betr_nlevsoi)
    this%betr_col(c)%z(1,1:betr_nlevsoi)  = z(c,1:betr_nlevsoi)
    this%betr_col(c)%pfti(1)              = pfti(c)
    this%betr_col(c)%pftf(1)              = pftf(c)
    this%betr_col(c)%npfts(1)             = npfts(c)

  enddo

  end subroutine BeTRSimulationSetcpsCol

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetcpsPFT(this, begc, endc, pfti, npfts, active, itype, wtcol)
  !
  !DESCRIPTION
  ! Set up BeTR PFT
  !
  !USES
  use betr_varcon , only : betr_maxpatch_pft
  use pftvarcon   , only : noveg, crop
  !
  !ARGUMENTS
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  integer                     , intent(in)    :: begc, endc
  integer, pointer            , intent(in)    :: pfti(:)
  integer, pointer            , intent(in)    :: npfts(:)
  logical, pointer            , intent(in)    :: active(:)
  integer, pointer            , intent(in)    :: itype(:)
  real(r8)                    , intent(in)    :: wtcol(:)
  !
  ! !LOCAL VARIABLES:
  integer :: c, p, pi, pp

  do c = begc, endc

      this%betr_pft(c)%column(:) = 1
      this%betr_pft(c)%npfts     = 0
      pp = 0
      do pi = 1, betr_maxpatch_pft
        if (pi <= npfts(c)) then
          p = pfti(c) + pi - 1
          if (active(p) .and. (itype(p) .ne. noveg)) then
            pp = pp + 1
            this%betr_pft(c)%wtcol(pp) = wtcol(p)
            this%betr_pft(c)%itype(pp) = itype(p)
            this%betr_pft(c)%crop(pp)  = crop(pi)         !the crop looks weird here, jyt
          endif
        endif
      enddo
      this%betr_pft(c)%npfts = pp

  enddo

  end subroutine BeTRSimulationSetcpsPFT

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBounds(this, betr_bounds)
  !
  !DESCRIPTION
  !set betr_bounds
  !
  use betr_varcon    , only : betr_maxpatch_pft
  use tracer_varcon  , only : betr_nlevsoi
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  type(betr_bounds_type), intent(out)  :: betr_bounds

  betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
  betr_bounds%begp = 1; betr_bounds%endp = betr_maxpatch_pft
  betr_bounds%begc = 1; betr_bounds%endc = 1
  betr_bounds%begl = 1; betr_bounds%endl = 1
  betr_bounds%begg = 1; betr_bounds%endg = 1

  end subroutine BeTRSimulationSetBounds

  !------------------------------------------------------------------------
  function do_soibgc(this)result(yesno)

  implicit none
  class(betr_simulation_type) , intent(inout) :: this

  logical :: yesno
  yesno = this%active_soibgc
  return
  end function do_soibgc

  !------------------------------------------------------------------------
  subroutine BeTRInitializeBiophysForc(this, betr_bounds, begc, endc)
    !
    !DESCRIPTION
    ! set up columns
    !
    implicit none
    !ARGUMENTS
    !
    class(betr_simulation_type) , intent(inout) :: this
    type(betr_bounds_type)                      :: betr_bounds
    integer                     , intent(in)    :: begc, endc
    !
    ! !LOCAL VARIABLES:
    integer :: c, l

    do c = begc, endc
       call this%biophys_forc(c)%Init(betr_bounds)
    enddo

  end subroutine BeTRInitializeBiophysForc

  !------------------------------------------------------------------------
  subroutine BeTRInitializeCols(this, betr_bounds, begc, endc)
    !
    !DESCRIPTION
    ! set up columns
    !
    implicit none
    !ARGUMENTS
    !
    class(betr_simulation_type) , intent(inout) :: this
    type(betr_bounds_type)                      :: betr_bounds
    integer                     , intent(in)    :: begc, endc
    !
    ! !LOCAL VARIABLES:
    integer :: c

    do c = begc, endc
       call this%betr_col(c)%Init(betr_bounds)
    enddo

  end subroutine BeTRInitializeCols

  !------------------------------------------------------------------------
  subroutine BeTRInitializePFT(this, betr_bounds, begc, endc)
    !
    !DESCRIPTION
    ! set up columns
    !
    implicit none
    !ARGUMENTS
    !
    class(betr_simulation_type) , intent(inout) :: this
    type(betr_bounds_type)                      :: betr_bounds
    integer                     , intent(in)    :: begc, endc
    !
    ! !LOCAL VARIABLES:
    integer :: c

    do c = begc, endc
       call this%betr_pft(c)%Init(betr_bounds)
    enddo

  end subroutine BeTRInitializePFT

  !------------------------------------------------------------------------
  subroutine BeTRInitializeActiveCols(this, betr_bounds, begc, endc, landunit, lun_itype)
    !
    !DESCRIPTION
    ! set up columns
    !
    !USES
    use landunit_varcon, only : istsoil, istcrop
    implicit none
    !ARGUMENTS
    !
    class(betr_simulation_type) , intent(inout) :: this
    type(betr_bounds_type)                      :: betr_bounds
    integer                     , intent(in)    :: begc, endc
    integer                     , intent(in)    :: landunit(:)
    integer                     , intent(in)    :: lun_itype(:)
    !
    ! !LOCAL VARIABLES:
    integer :: c, l

    do c = begc, endc

       l = landunit(c)

       if (lun_itype(l) == istsoil .or. lun_itype(l) == istcrop) then
          this%active_col(c) = .true.
       else
          this%active_col(c) = .false.
       endif
    enddo

  end subroutine BeTRInitializeActiveCols

  !------------------------------------------------------------------------
  subroutine BeTRInitialize_1(this, betr_bounds, begc, endc, namelist_buffer)
    !
    !DESCRIPTION
    ! set up columns
    !
    !USES
    use betr_constants , only : betr_namelist_buffer_size
    use tracer_varcon  , only : betr_nlevsoi
    !
    implicit none
    !
    !ARGUMENTS
    class(betr_simulation_type)              , intent(inout) :: this
    type(betr_bounds_type)                                   :: betr_bounds
    integer                                  , intent(in)    :: begc, endc
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    !
    ! !LOCAL VARIABLES:
    integer :: c, l
    logical :: asoibgc

    do c = begc, endc

       call this%betr(c)%Init(namelist_buffer, betr_bounds, &
            this%betr_col(c), this%biophys_forc(c), &
            asoibgc, this%bstatus(c))

       if (c==begc) this%active_soibgc=asoibgc

       if (this%bstatus(c)%check_status())then
          call this%bsimstatus%setcol(c)
          call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
          exit
       endif

    enddo

    if (this%bsimstatus%check_status()) call endrun(msg=this%bsimstatus%print_msg())

    do c = begc, endc

       call this%biogeo_state(c)%Init(betr_bounds, this%active_soibgc)

       call this%biogeo_flux(c)%Init(betr_bounds, this%active_soibgc)
    enddo

    !identify variables that are used for history output
    c = begc
    call this%betr(c)%get_hist_size(this%num_hist_state1d, this%num_hist_state2d, &
         this%num_hist_flux1d, this%num_hist_flux2d)

    call this%HistAlloc(begc, endc)

    call this%betr(c)%get_hist_info(this%num_hist_state1d, &
         this%num_hist_state2d,                            &
         this%num_hist_flux1d, this%num_hist_flux2d,       &
         this%state_hist1d_var(1:this%num_hist_state1d),   &
         this%state_hist2d_var(1:this%num_hist_state2d),   &
         this%flux_hist1d_var(1:this%num_hist_flux1d),     &
         this%flux_hist2d_var(1:this%num_hist_flux2d))

    if (this%bstatus(c)%check_status()) then
       call this%bsimstatus%setcol(c)
       call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())

       if (this%bsimstatus%check_status()) call endrun(msg=this%bsimstatus%print_msg())
    endif

  end subroutine BeTRInitialize_1

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBaseFilename(this,base_filename)
    !
    !DESCRIPTION
    ! Sets base filename
    !ARGUMENTS
    implicit none
    class(betr_simulation_type)         , intent(inout)        :: this
    character(len=betr_filename_length) , optional, intent(in) :: base_filename

    if (present(base_filename)) then
      this%base_filename = base_filename
    else
      this%base_filename = ''
    endif

  end subroutine BeTRSimulationSetBaseFilename
end module BeTRSimulation
