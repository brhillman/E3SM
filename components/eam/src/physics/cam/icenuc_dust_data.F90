module icenuc_dust_data

!-----------------------------------------------------------------------
! Purpose:
!
! Read and time-interpolate prescribed dust number concentration fields
! used only for ice nucleation. Input data should be provided as a
! tracer_data-compatible time series on a lon/lat grid, with values in
! #/cm3 on model levels or on levels that tracer_data can vertically
! interpolate to the model grid.
!-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile, trcdata_init, advance_trcdata, get_fld_data
  use ppgrid,       only : pcols, pver, begchunk, endchunk
  use physics_types,only : physics_state
  use physics_buffer, only : physics_buffer_desc
  use cam_abortutils, only : endrun
  use cam_history,    only : addfld, outfld
  use cam_logfile,    only : iulog

  implicit none
  private
  save

  type(trfld), pointer :: fields(:) => null()
  type(trfile)         :: file

  logical, public :: icenuc_dust_data_on = .false.
  character(len=256) :: icenuc_dust_filename = ' '
  character(len=32)  :: icenuc_dust_type = 'CYCLICAL'
  integer            :: icenuc_dust_cycle_yr = 0
  character(len=32)  :: icenuc_dust_varname = 'dst_num'

  integer :: number_flds = 0

  public :: icenuc_dust_data_readnl
  public :: icenuc_dust_data_init
  public :: icenuc_dust_data_advance
  public :: get_icenuc_dust_data

contains

  subroutine icenuc_dust_data_readnl(nlfile)

    use namelist_utils, only : find_group_name
    use units,          only : getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile

    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'icenuc_dust_data_readnl'

    namelist /icenuc_dust_data_nl/ icenuc_dust_data_on, icenuc_dust_filename, &
                                   icenuc_dust_type, icenuc_dust_cycle_yr,    &
                                   icenuc_dust_varname

    if (masterproc) then
       unitn = getunit()
       open(unitn, file=trim(nlfile), status='old')
       call find_group_name(unitn, 'icenuc_dust_data_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, icenuc_dust_data_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(icenuc_dust_data_on, 1, mpilog, 0, mpicom)
    call mpibcast(icenuc_dust_filename, len(icenuc_dust_filename), mpichar, 0, mpicom)
    call mpibcast(icenuc_dust_type, len(icenuc_dust_type), mpichar, 0, mpicom)
    call mpibcast(icenuc_dust_cycle_yr, 1, mpiint, 0, mpicom)
    call mpibcast(icenuc_dust_varname, len(icenuc_dust_varname), mpichar, 0, mpicom)
#endif

  end subroutine icenuc_dust_data_readnl

  subroutine icenuc_dust_data_init()

    character(len=32) :: specifier(1)

    if (.not. icenuc_dust_data_on) return

    if (len_trim(icenuc_dust_filename) == 0) then
       call endrun('icenuc_dust_data_init: icenuc_dust_filename must be set when icenuc_dust_data_on is true')
    end if

    if (masterproc) then
       write(iulog,*) 'icenuc_dust_data_init: reading prescribed ice nucleation dust from ', trim(icenuc_dust_filename)
    end if

    specifier(1) = trim(icenuc_dust_varname)

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf = .false.

    call trcdata_init(specifier, icenuc_dust_filename, ' ', ' ', fields, file, &
                      .false., icenuc_dust_cycle_yr, 0, 0, icenuc_dust_type)

    number_flds = 0
    if (associated(fields)) number_flds = size(fields)

    if (number_flds < 1) then
       call endrun('icenuc_dust_data_init: no prescribed dust fields were initialized')
    end if

    call addfld('INndustx', (/ 'lev' /), 'A', '1/m3', &
                'Injected dust number concentration used for ice nucleation')

  end subroutine icenuc_dust_data_init

  subroutine icenuc_dust_data_advance(state)

    type(physics_state), intent(in) :: state(begchunk:endchunk)

    if (.not. icenuc_dust_data_on) return

    call advance_trcdata(fields, file, state)

  end subroutine icenuc_dust_data_advance

  subroutine get_icenuc_dust_data(lchnk, ncol, pbuf, extra_dst_num)

    integer, intent(in) :: lchnk, ncol
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out) :: extra_dst_num(pcols,pver)

    extra_dst_num = 0._r8

    if (.not. icenuc_dust_data_on) return
    if (number_flds < 1) return

    call get_fld_data(fields, trim(icenuc_dust_varname), extra_dst_num, ncol, lchnk, pbuf)
    call outfld('INndustx', extra_dst_num*1.e6_r8, pcols, lchnk)

  end subroutine get_icenuc_dust_data

end module icenuc_dust_data
