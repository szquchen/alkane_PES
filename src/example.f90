program main
use pes_shell
implicit none

  integer::j
  real,allocatable::x(:),g(:),g0(:)
  real::pot, pot0
  character,allocatable::symb(:)

  !=== Initialize the PES
  call pes_init()

  !=== specify the number of carbon ===!
  !=== must be done before calling the PES ===!
  nc = 14
  nh = 2*nc + 2
  natm = nh + nc
  allocate(x(3*natm), symb(natm), g(3*natm), g0(3*natm))

  !== example, read in Cartesian coordinates x ==!
  !== then compute energy and gradient ==!
  !== x: dimension(1:3N), Cartesian coordinate in bohr ==!
  !== Atoms must be in order C C ... C H H ... H ==!
  open(21,file="test.xyz",status="old")
  open(22,file="test.pot",status="unknown")
  open(23,file="test.grd",status="unknown")

  read(21,*)
  read(21,*) pot0
  do j=1,natm
     read(21,*) symb(j), x(3*j-2:3*j), g0(3*j-2:3*j)
  end do
  close(21)
  x = x / auang ! convert to bohr

  call pot_grad(x, pot, g) ! output energy in hartree, g in hartree/bohr

  write(22,'(A, F15.8)') "Ab initio energy (hartree):", pot0
  write(22,'(A, F15.8)') "Predicted energy (hartree):", pot
  close(22)

  write(23,'(A)') "#  Ab inito   Predicted   diff. (hartree/bohr):"
  do j=1,3*natm
     write(23,'(3F15.8)') g0(j), g(j), g(j)-g0(j)
  end do
  close(23)

  deallocate(x, symb, g, g0)

end program
