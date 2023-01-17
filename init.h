!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     integer, parameter :: ns=8
!     integer, parameter :: nmax = 4900
!     integer, parameter :: ns=7
!     integer, parameter :: nmax = 1225
!     integer, parameter :: ns=6
!     integer, parameter :: nmax = 400
      integer, parameter :: ns=5
      integer, parameter :: nmax = 100
!     integer, parameter :: ns=4
!     integer, parameter :: nmax = 36
!     integer, parameter :: ns=3
!     integer, parameter :: nmax = 9
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer, parameter :: prozessoren = 36
    integer, parameter :: nmpara = 30
    integer, parameter :: ksteps = 100
    integer, parameter :: Iwmax = 2**10  
    integer, parameter :: Iwmaxreal = 2**10
    integer, parameter :: lattice_type = 5 !1: SC-3D, 2:FCC, 3: SC-2D, 4: bethe, 5: magnetic SC-2D, 6: magnetic SC-2D with more hopping
    logical, parameter :: symm = .false.
    logical, parameter :: gwcalc = .false.
    integer, parameter :: p = 0	!nominator in B field
    integer, parameter :: L = 3	!denominator in B field
