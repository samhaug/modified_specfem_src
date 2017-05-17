!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
! HMM
!
! generic heterogeneous mantle model
! Modified by Carlos Chaves at the University of Michigan - 10/02/2015
! Updated by Carlos Chaves at the University of Michigan - 05/01/2016
! Now it is working for N_R != N_THETA != N_PHI
!--------------------------------------------------------------------------------------------------

module model_heterogen_mantle_par

    ! heterogen_mantle_model_constants
    integer, parameter :: N_R = 201,N_THETA = 201,N_PHI = 201

    ! model array
    double precision,dimension(:),allocatable :: HMM_rho_in

end module model_heterogen_mantle_par

!
!--------------------------------------------------------------------------------------------------
!

subroutine model_heterogen_mntl_broadcast(myrank)

    ! standard routine to setup model

    use constants
    use model_heterogen_mantle_par

    implicit none

    integer :: myrank

    ! local parameters
    integer :: ier

    ! allocates model array
    allocate(HMM_rho_in(N_R*N_THETA*N_PHI),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating HMM array')

    ! master process reads in model
    if (myrank == 0) then
        write(IMAIN,*) 'Reading in model_heterogen_mantle.'
        call flush_IMAIN()

        call read_heterogen_mantle_model()

        write(IMAIN,*) 'model_heterogen_mantle is read in.'
        call flush_IMAIN()
    endif

    ! broadcast the information read on the master to the nodes
    call bcast_all_dp(HMM_rho_in,N_R*N_THETA*N_PHI)

    if (myrank == 0) then
        write(IMAIN,*) 'model_heterogen_mantle is broadcast.'
        write(IMAIN,*) 'First value in HMM:',HMM_rho_in(1)
        write(IMAIN,*) 'Last value in HMM:',HMM_rho_in(N_R*N_THETA*N_PHI)
        call flush_IMAIN()
    endif

end subroutine model_heterogen_mntl_broadcast

!
!-------------------------------------------------------------------------------------------------
!
subroutine read_heterogen_mantle_model()

    use constants
    use model_heterogen_mantle_par

    implicit none

    ! local parameters
    integer :: i,j,ier

    ! open heterogen.dat
    open(unit=IIN,file='./DATA/heterogen/heterogen.dat',access='direct',&
form='formatted',recl=20,status='old',action='read',iostat=ier)
    if (ier /= 0 ) call exit_MPI(0,'Error opening model file heterogen.dat')

    j = N_R*N_THETA*N_PHI

    do i = 1,j
        read(IIN,rec=i,fmt='(F20.15)') HMM_rho_in(i)
    enddo

    close(IIN)

end subroutine read_heterogen_mantle_model

!
!-------------------------------------------------------------------------------------------------
!

subroutine model_heterogen_mantle(radius,theta,phi,dvs,dvp,drho)

    use constants
    use model_heterogen_mantle_par

    implicit none

    integer :: pos
    double precision :: radius,theta,phi            ! input coordinates
    double precision :: rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8
    double precision :: dvs, dvp, drho
    double precision :: r_inner,r_outer ! lower and upper domain bounds for r
    double precision :: xd, yd, zd
    double precision :: dr, dlon, dlat, r, lat, lon, latmin
    double precision :: latmax, lonmin, lonmax, rmax, rmin,lon0,lon1,lat0,lat1

    r_inner = 3.480d6  !lower bound for heterogeneity zone
    ! NOTE: r_outer NEEDS TO BE (just) SMALLER THAN R_EARTH!!!!!!!!
    r_outer = R_EARTH-1.0  !6.300d6  !upper bound for heterogeneity zone (lower mantle: e.g. 4.500d6)

    lon0 = -10.0

    lon1 = 10.0

    lat0 = 80.0

    lat1 = 100.0

    dr=R_EARTH_KM/(real(N_R-1))

    !dlon=360.0/(real(N_PHI-1))

    !dlat=180.0/(real(N_THETA-1))

    dlon=20.0/(real(N_PHI-1))

    dlat=20.0/(real(N_THETA-1))

    r=radius*R_EARTH_KM

    radius = radius*R_EARTH

    lat=(RADIANS_TO_DEGREES)*theta !Colatitude, but just for indentification, I call it latitude. Colatitude, indeed is 90-lat;

    lon=(RADIANS_TO_DEGREES)*phi

    latmin=floor(lat/dlat)*dlat

    latmax=latmin+dlat

    lonmin=floor(lon/dlon)*dlon

    lonmax=lonmin+dlon

    rmax=floor(r/dr)*dr

    rmin=rmax+dr

    if (rmax*1000 < r_inner ) rmax=rmin

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((radius >= r_inner) .and. (radius <= r_outer) .and. (lat >= lat0) .and. &
    (lat <= lat1) .and. (lon >= lon0) .and. (lon <= lon1)) then

        xd=(lon-lonmin)/(dlon)
        yd=(lat-latmin)/(dlat)
        zd=(r-rmax)/(dr)

        ! rho1 at: x_low y_low z_low
        pos=(((int(((latmin)/dlat)+1.0)-1)*N_PHI)+int(((lonmin)/dlon)+1.0))+&
((int(((R_EARTH_KM-rmin)/dr)+1.0)-1)*N_THETA*N_PHI) !This part is new in the code!
        rho1 = HMM_rho_in(pos)

        ! rho2 at: x_high y_low z_low
        pos=(((int(((latmin)/dlat)+1.0)-1)*N_PHI)+int(((lonmax)/dlon)+1.0))+&
((int(((R_EARTH_KM-rmin)/dr)+1.0)-1)*N_THETA*N_PHI) !This part is new in the code!
        rho2 = HMM_rho_in(pos)

        ! rho3 at:x_low y_high z_low
        pos=(((int(((latmax)/dlat)+1.0)-1)*N_PHI)+int(((lonmin)/dlon)+1.0))+&
((int(((R_EARTH_KM-rmin)/dr)+1.0)-1)*N_THETA*N_PHI) !This part is new in the code!
        rho3 = HMM_rho_in(pos)

        ! rho4 at: x_high y_high z_low
        pos=(((int(((latmax)/dlat)+1.0)-1)*N_PHI)+int(((lonmax)/dlon)+1.0))+&
((int(((R_EARTH_KM-rmin)/dr)+1.0)-1)*N_THETA*N_PHI) !This part is new in the code!
        rho4 = HMM_rho_in(pos)

        ! rho5 at: x_low y_low z_high
        pos=(((int(((latmin)/dlat)+1.0)-1)*N_PHI)+int(((lonmin)/dlon)+1.0))+&
((int(((R_EARTH_KM-rmax)/dr)+1.0)-1)*N_THETA*N_PHI) !This part is new in the code!
        rho5 = HMM_rho_in(pos)

        ! rho6 at: x_high y_low z_high
        pos=(((int(((latmin)/dlat)+1.0)-1)*N_PHI)+int(((lonmax)/dlon)+1.0))+&
((int(((R_EARTH_KM-rmax)/dr)+1.0)-1)*N_THETA*N_PHI) !This part is new in the code!
        rho6 = HMM_rho_in(pos)

        ! rho7 at: x_low y_high z_high
        pos=(((int(((latmax)/dlat)+1.0)-1)*N_PHI)+int(((lonmin)/dlon)+1.0))+&
((int(((R_EARTH_KM-rmax)/dr)+1.0)-1)*N_THETA*N_PHI) !This part is new in the code!
        rho7 = HMM_rho_in(pos)

        ! rho8 at: x_high y_high z_high
        pos=(((int(((latmax)/dlat)+1.0)-1)*N_PHI)+int(((lonmax)/dlon)+1.0))+&
((int(((R_EARTH_KM-rmax)/dr)+1.0)-1)*N_THETA*N_PHI) !This part is new in the code!
        rho8 = HMM_rho_in(pos)

        drho = ((rho1*(1.-xd)+rho2*xd)*(1.-yd)+(rho3*(1.-xd)+rho4*xd)*yd)*(1.-zd)+&
((rho5*(1.-xd)+rho6*xd)*(1.-yd)+(rho7*(1.-xd)+rho8*xd)*yd)*zd

        dvp = (0.55/0.30)*drho !! I preserved the original values!
        dvs = (1.00/0.30)*drho !! I preserved the original values!

    else !outside of mantle domain
        drho = 0.
        dvs = 0.
        dvp = 0.
        endif

end subroutine model_heterogen_mantle
