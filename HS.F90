program HS
!
! Jun INAGAKI 
! 2021.10.10      
!

  implicit none      
  
  character(10) :: outname
  real(8) :: r(521), ru(521), ru2(521), ru3(521), v(521)
  real(8) :: x(521)
  real(8) :: ee(24), wwnl(24), nnlz(24)
  real(8) :: tol, thresh
  real(8) :: deltax
  real(8) :: c, z, zzz, www
  real(8) :: twoz, twozzz

  integer :: istring, nfiles, ncards, nblock, LimMesh
  integer :: key, mesh, ipratt, maxit, nocopy, kut
  integer :: i, j, k, m, i1, i2
  integer :: ncspvs, ncores, nvales, ion, iz, twoion

! read input.
  open(11,file="in",status="old")

  read(11,'(a)') outname
  read(11,'(72a1)') 
  read(11,*) key, tol, thresh, mesh, ipratt, maxit, nocopy, kut
  read(11,*) LimMesh

!debug
  write(*,*) "outname ",outname
  write(*,*) "key ",key
  write(*,*) "tol ",tol
  write(*,*) "thresh ",thresh
  write(*,*) "mesh ",mesh
  write(*,*) "ipratt ",ipratt
  write(*,*) "maxit ",maxit
  write(*,*) "nocopy ",nocopy
  write(*,*) "kut ",kut
!end debug

  ncards = 90

  if (maxit <= 0) then
    maxit = 20
  end if

  nblock = mesh / 40

!debug
  write(*,*) "nblock ",nblock

! construct x mesh and r mesh

  i = 1
  x(i) = 0.0
  r(i) = 0.0
  deltax = 0.0025

  do j = 1, nblock
    do k = 1, 40
      i = i + 1
      x(i) = x(i-1) + deltax
    end do
    deltax = deltax + deltax
  end do

! assume key = 0  
  read(11,'(10f7.5)') (ru2(m),m=1,437,4)
  read(11,*) z, ncores, nvales, ion

  if (z .le. 0) then
    write(*,*) "error, z must be > 0."
  endif

  nfiles = nfiles + 1

    iz = z
    ncspvs = ncores + nvales
    c = 0.88534138 / (z**(1.0/3.0))
    twoion = ion + ion
    zzz = ion + 1
    twozzz = zzz + zzz

    do i = 2, mesh
      r(i) = c * x(i)
    end do

    read(11,*) (nnlz(i),wwnl(i),ee(i),i=1,ncspvs)

!debug
    write(*,*) nnlz

    www = 0.0
    do i = 1, ncspvs
      www = www + wwnl(i)
    end do
   
    if (abs(z + 1. - www - zzz) .ge. 0.001) then
      write(*,*) "z,www,zzz",z,www,zzz
      stop
    end if

! assume key = 0
    twoz = z + z

end program HS      
