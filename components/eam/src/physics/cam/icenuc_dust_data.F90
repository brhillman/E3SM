module icenuc_dust_data

!-----------------------------------------------------------------------
! Purpose:
! Optional prescribed 3-D number concentration field for additional dust
! used only by ice nucleation pathways.
!
! Expected input dataset:
!  - netCDF variable named ICENUC_DUST_NUM
!  - 12 monthly samples using either:
!      * zonal layout: ICENUC_DUST_NUM(lon=1,lev,lat,time) with lat(lat), lev(lev)
!      * full 3-D layout on the CAM physics grid: ICENUC_DUST_NUM(ncol,lev,time)
!        with PS(ncol,time) and hybi(lev+1)
!  - lon/lat files with lon>1 are not supported by this path; use ncol for full
!    horizontal structure
!  - required time variables: date(time)
!  - optional time variable: datesec(time)
!  - lev values are read as pressure levels in mb and converted internally
!  - the first date entry must be in January and the last in December
!-----------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: begchunk, endchunk, pcols, pver
use cam_abortutils, only: endrun
use cam_logfile,    only: iulog
use physics_types,  only: physics_state
use boundarydata,   only: boundarydata_type, boundarydata_init, boundarydata_update, &
                          boundarydata_vert_interp
use mpishorthand

implicit none
private
save

public :: &
   icenuc_dust_data_readnl,        &
   icenuc_dust_data_register,      &
   icenuc_dust_data_init,          &
   icenuc_dust_data_timestep_init

logical, public    :: icenuc_dust_data_on = .false.
character(len=256) :: bndtvi_dust_num = ' '
real(r8)           :: icenuc_dust_num_max = 1.0e12_r8 ! #/m3 sanity limit for prescribed input data

integer :: icenuc_dust_num_idx = -1

type(boundarydata_type) :: icenuc_dust_num_data
character(len=15), parameter, dimension(1) :: nc_name = (/'ICENUC_DUST_NUM'/)

!================================================================================================
contains
!================================================================================================

subroutine icenuc_dust_data_readnl(nlfile)

   use namelist_utils, only: find_group_name
   use units,          only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile

   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'icenuc_dust_data_readnl'

   namelist /icenuc_dust_data_nl/ icenuc_dust_data_on, bndtvi_dust_num, icenuc_dust_num_max

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
   call mpibcast(bndtvi_dust_num, len(bndtvi_dust_num), mpichar, 0, mpicom)
   call mpibcast(icenuc_dust_num_max, 1, mpir8, 0, mpicom)
#endif

   if (icenuc_dust_num_max <= 0._r8) then
      call endrun(subname // ':: icenuc_dust_num_max must be > 0')
   end if

end subroutine icenuc_dust_data_readnl

!================================================================================================

subroutine icenuc_dust_data_register()

   use physics_buffer, only: pbuf_add_field, dtype_r8

   if (icenuc_dust_data_on) then
      call pbuf_add_field('ICENUC_DUST_NUM', 'physpkg', dtype_r8, (/pcols,pver/), icenuc_dust_num_idx)
   end if

end subroutine icenuc_dust_data_register

!================================================================================================

subroutine icenuc_dust_data_init(phys_state)

   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

   if (.not. icenuc_dust_data_on) return

   if (len_trim(bndtvi_dust_num) == 0) then
      call endrun('icenuc_dust_data_init:: bndtvi_dust_num must be set when icenuc_dust_data_on=.true.')
   end if

   call boundarydata_init(bndtvi_dust_num, phys_state, nc_name, 1, icenuc_dust_num_data, 3)

   if (masterproc) then
      write(iulog,*) 'icenuc_dust_data_init: Initializing prescribed ice-nucleation dust number data'
      write(iulog,*) 'icenuc_dust_data_init: dataset is: ', trim(bndtvi_dust_num)
      write(iulog,*) 'icenuc_dust_data_init: max allowed value (#/m3): ', icenuc_dust_num_max
   end if

end subroutine icenuc_dust_data_init

!================================================================================================

subroutine icenuc_dust_data_timestep_init(pbuf2d, phys_state)

   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   real(r8), pointer :: tmpptr(:,:)
   integer :: lchnk

   if (.not. icenuc_dust_data_on) return

   call boundarydata_update(phys_state, icenuc_dust_num_data)

   do lchnk = begchunk, endchunk
      call pbuf_get_field(pbuf_get_chunk(pbuf2d, lchnk), icenuc_dust_num_idx, tmpptr)
      call icenuc_dust_data_get_num(phys_state(lchnk), tmpptr)
   end do

end subroutine icenuc_dust_data_timestep_init

!================================================================================================

subroutine icenuc_dust_data_get_num(state, q)

   type(physics_state), intent(in) :: state
   real(r8),            intent(inout) :: q(:,:)

   integer :: lchnk, i, k
   real(r8) :: dust_num_in(pcols,icenuc_dust_num_data%levsiz)
   character(len=*), parameter :: subname = 'icenuc_dust_data_get_num'

   lchnk = state%lchnk

   dust_num_in = 0._r8
   do k = 1, icenuc_dust_num_data%levsiz
      do i = 1, state%ncol
         if (icenuc_dust_num_data%isncol) then
            dust_num_in(i,k) = icenuc_dust_num_data%datainst(i,k,lchnk,1)
         else
            dust_num_in(i,k) = icenuc_dust_num_data%datainst(state%latmapback(i),k,lchnk,1)
         end if
      end do
   end do

   call boundarydata_vert_interp(lchnk, state%ncol, icenuc_dust_num_data%levsiz, &
                                 1, icenuc_dust_num_data%pin, state%pmid, dust_num_in, q)

   do k = 1, pver
      do i = 1, state%ncol
         if (q(i,k) < 0._r8) then
            if (q(i,k) > -1.e-12_r8) then ! allow tiny negative roundoff after interpolation
               q(i,k) = 0._r8
            else
               write(iulog,*) subname//': negative ICENUC_DUST_NUM value at (i,k)=', i, k, q(i,k)
               call endrun(subname//': negative prescribed dust number not allowed')
            end if
         end if

         if (q(i,k) > icenuc_dust_num_max) then
            write(iulog,*) subname//': ICENUC_DUST_NUM exceeds max at (i,k)=', i, k, q(i,k), &
                          ' max=', icenuc_dust_num_max
            call endrun(subname//': prescribed dust number exceeds icenuc_dust_num_max')
         end if
      end do
   end do

end subroutine icenuc_dust_data_get_num

!================================================================================================

end module icenuc_dust_data
