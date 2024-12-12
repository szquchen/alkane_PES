module bemsa2

contains
  !==========================================!
  ! Evalute the basis                        !
  !   x: Morse variable of internuclear dist !
  !   p: dimension(nbasis), basis            !
  !==========================================!
  subroutine evp2(r, lambda, p)
    real,dimension(:),intent(out)::p
    real,intent(in)::r, lambda
    !:::::::::::::::::
    integer::i,np

    np = size(p) ! max. power
    p(1) = exp(-r / lambda)
    do i=2,np
       p(i) = p(i-1)*p(1)
    end do

    return
  end subroutine evp2

  !====================================================!
  ! Evalute the derivatives of the basis               !
  !   x1, x2: Cartesian of two atoms, in Ang           !
  !   r: internuclear distance in Ang                  !
  !   a2: Morse parameter                              !
  !   icomp: (int 1 to 6) the gradient component       !
  !   p: dimension(nbasis), the basis functions        !
  !   dp: dimension(nbasis), first derivative of basis !
  !====================================================!
  subroutine evdp2(x1, x2, r, a2, icomp, p, dp)
    real,intent(in)::r, a2, x1(3), x2(3)
    real,dimension(:),intent(out)::dp
    real,dimension(:),intent(in)::p
    integer,intent(in)::icomp
    !:::::::::::::::::::::::::
    integer::i,np

    np = size(p)
    if (icomp .eq. 1) then
       do i=1,np
          dp(i) = -real(i) * p(i) * (x1(1)-x2(1)) / (a2*r)
       end do
    elseif (icomp .eq. 2) then
       do i=1,np
          dp(i) = -real(i) * p(i) * (x1(2)-x2(2)) / (a2*r)
       end do
    elseif (icomp .eq. 3) then
       do i=1,np
          dp(i) = -real(i) * p(i) * (x1(3)-x2(3)) / (a2*r)
       end do
    elseif (icomp .eq. 4) then
       do i=1,np
          dp(i) = -real(i) * p(i) * (x2(1)-x1(1)) / (a2*r)
       end do
    elseif (icomp .eq. 5) then
       do i=1,np
          dp(i) = -real(i) * p(i) * (x2(2)-x1(2)) / (a2*r)
       end do
    elseif (icomp .eq. 6) then
       do i=1,np
          dp(i) = -real(i) * p(i) * (x2(3)-x1(3)) / (a2*r)
       end do
    end if

   return
 end subroutine evdp2

 !=====================================!
 ! Calculate one 2b energy             !
 !   r: internuclear distance          !
 !   lambda: Morse parameter           !
 !   coef: coefficients                !
 !   pot: potential energy in hartree  !
 !=====================================!
 function pot2b(r, lambda, coef)
   real,dimension(:),intent(in)::coef
   real,intent(in)::r, lambda
   real::pot2b
   !::::::::::::::::::
   real,dimension(size(coef))::p

   call evp2(r, lambda, p)
   pot2b = dot_product(coef, p)

   return
 end function pot2b

  !======================================!
  ! Calculate one 2b gradients           !
  !   x: dimension(6), Cartesian in bohr !
  !   r: internuclear distance           !
  !   coef: size (nbasis+1) coefs        !
  !   pot: potential energy in eV        !
  !======================================!
!  function grad2b(x1, x2, r, a2, coef)
!    real,dimension(:),intent(in)::coef
!    real,dimension(:),intent(in)::x1,x2
!    real,intent(in)::r,a2
!    real,dimension(6)::grad2b
!    !:::::::::::::::::::::::::::::
!    real,dimension(size(coef))::dp
!    integer::ncoef
!
!    do i=1,6
!       call evdp2(x1, x2, r, a2, i, dp)
!       grad2b(i) = dot_product(coef, dp)
!    end do
!
!    return
!  end function grad2b

end module
