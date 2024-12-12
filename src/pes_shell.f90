module pes_shell
use bemsa2
use bemsa3
use bemsa21
use bemsa4
use bemsa31
use bemsa22

  implicit none

  !=========================
  ! Setting parameters
  !=========================
  real,parameter::pi=3.14159265
  real,parameter::auang=0.5291772083
  real,parameter::aucm=219474.63
  real,parameter::c_mass= 12.0000000  !21874.66
  real,parameter::h_mass=  1.0078250  !1837.15
  real,parameter::d_mass=  2.0141018  !D
  real,parameter::o_mass= 15.9949146  !29156.95
  real,parameter::n_mass= 14.003074
  real,parameter::emass=1822.88848

  real,parameter::r2i= 8.5 / auang
  real,parameter::r2f= 9.5 / auang
  real,parameter::r3i= 6.5 / auang
  real,parameter::r3f= 7.5 / auang
  real,parameter::r4i= 4.5 / auang
  real,parameter::r4f= 5.5 / auang
  integer::nc,nh
  integer::natm

  integer,parameter::np2=10
  integer,parameter::np3=32 !32 23
  integer,parameter::np21=78 !78 55
  integer,parameter::np4=40 !77 40 18
  integer,parameter::np31=115 !234 115 49 
  integer,parameter::np22=174 !354 174 75
  integer,parameter::nm3=7 !7 7
  integer,parameter::nm21=4 !4 4
  integer,parameter::nm4=87 !111 87 62
  integer,parameter::nm31=32 !32 32 32
  integer,parameter::nm22=17 !17 17 17
  integer,parameter::nq3=6 !6 5
  integer,parameter::nq21=8 !8 7
  integer,parameter::nq4=12 !18 12 6
  integer,parameter::nq31=25 !41 25 15
  integer,parameter::nq22=41 !62 41 26
  integer,parameter::ncoef=3*np2 + 2*np3 + 2*np21 + 2*np4 + 2*np31 + np22

  real::p2(np2), dp2(np2)
  real::m3(nm3), m21(nm21)
  real::p3(np3), p21(np21)
  real::q3(nq3), q21(nq21)
  real::m4(nm4), m31(nm31), m22(nm22)
  real::p4(np4), p31(np31), p22(np22)
  real::q4(nq4), q31(nq31), q22(nq22)
  real,dimension(ncoef)::coef

  real,parameter::a_cc=2.5
  real,parameter::a_ch=2.5
  real,parameter::a_hh=2.5
  real,parameter::a_c3=1.8
  real,parameter::a_c2h=1.8
  real,parameter::a_h2c=1.8
  real,parameter::a_h3=1.8
  real,parameter::a_c4=1.2
  real,parameter::a_c3h=1.2
  real,parameter::a_c2h2=1.2
  real,parameter::a_h3c=1.2
  real,parameter::a_h4=1.2

contains
  subroutine pes_init()
    integer::i

    open(20,file="coeff.dat",status="old")
    do i=1,ncoef
       read(20,*) coef(i)
    end do
    close(20)

    return
  end subroutine pes_init

  !==================================================!
  ! switching functions for 2b,3b,4b
  !==================================================!
  subroutine f_switch(s,r,ri,rf)
    real,intent(out)::s
    real,intent(in)::r,ri,rf
    !::::::::::::::::::::
    real::ra,ra2,ra3

    ra=(r-ri)/(rf-ri)
    ra2=ra*ra
    ra3=ra2*ra

    s=10.0*ra3-15.0*ra*ra3+6.0*ra3*ra2
    s=1.0-s

  end subroutine f_switch

  !========================
  ! Calculating potential
  !========================
  subroutine getpot(x1d, pot)
    real,dimension(:),intent(in)::x1d
    real,intent(out)::pot
    !::::::::::::::::::::
    real,dimension(3,size(x1d)/3)::x
    real::rint(size(x1d)/3, size(x1d)/3),morse3(3),morse4(6)
    integer::j,k,l,m,istart,iend, natm
    real::vt,s,rmax
    integer::n4b

    natm = size(x1d)/3

    do j=1,natm
      x(:,j) = x1d(3*j-2:3*j)
    end do

    do j=1,natm-1
      do k=j+1,natm
        rint(j,k) = norm2(x(:,j)-x(:,k))
        rint(k,j) = rint(j,k)
      end do
    end do

    !==============================
    ! compute the predicted energy
    !==============================
    pot = -37.8519753684*real(nc) - 0.501257936933*real(nh)

    ! C-C 2-body
    istart=1
    iend = np2+istart-1
    do j=1,nc-1
      do k=j+1,nc
        if (rint(j,k) > r2f) cycle
        vt = pot2b(rint(j,k), a_cc, coef(istart:iend))
        if (rint(j,k) > r2i) then
          call f_switch(s, rint(j,k), r2i, r2f)
          pot = pot + vt*s
        else
          pot = pot + vt
        end if
      end do
    end do

    ! C-H 2-body
    istart=istart+np2
    iend = iend+np2
    do j=1,nc
      do k=nc+1,natm
        if (rint(j,k) > r2f) cycle
        vt = pot2b(rint(j,k), a_ch, coef(istart:iend))
        if (rint(j,k) > r2i) then
          call f_switch(s, rint(j,k), r2i, r2f)
          pot = pot + vt*s
        else
          pot = pot + vt
        end if
      end do
    end do

    ! H-H 2-body
    istart=istart+np2
    iend = iend+np2
    do j=nc+1,natm-1
      do k=j+1,natm
        if (rint(j,k) > r2f) cycle
        vt = pot2b(rint(j,k), a_hh, coef(istart:iend))
        if (rint(j,k) > r2i) then
          call f_switch(s, rint(j,k), r2i, r2f)
          pot = pot + vt*s
        else
          pot = pot + vt
        end if
      end do
    end do

    ! CCC 3b
    istart = istart+np2
    iend = iend+np3
    do j=1,nc-2
      do k=j+1,nc-1
        if (rint(j,k) > r3f) cycle
        morse3(1)=exp(-rint(j,k) / a_c3)
        do l=k+1,nc
          if (rint(j,l) > r3f .or. rint(k,l) > r3f) cycle
          morse3(2)=exp(-rint(j,l) / a_c3)
          morse3(3)=exp(-rint(k,l) / a_c3)
          vt = emsav3(morse3, coef(istart:iend))
          rmax = max(rint(j,k), rint(j,l), rint(k,l))
          if (rmax > r3i) then
            call f_switch(s, rmax, r3i, r3f)
            pot = pot + vt*s
          else
            pot = pot + vt
          end if
        end do
      end do
    end do

    ! C2H 3b
    istart = istart+np3
    iend = iend+np21
    do j=1,nc-1
      do k=j+1,nc
        if (rint(j,k) > r3f) cycle
        morse3(1)=exp(-rint(j,k) / a_c2h)
        do l=nc+1,natm
          if (rint(j,l) > r3f .or. rint(k,l) > r3f) cycle
          morse3(2)=exp(-rint(j,l) / a_c2h)
          morse3(3)=exp(-rint(k,l) / a_c2h)
          vt = emsav21(morse3, coef(istart:iend))
          rmax = max(rint(j,k), rint(j,l), rint(k,l))
          if (rmax > r3i) then
            call f_switch(s, rmax, r3i, r3f)
            pot = pot + vt*s
          else
            pot = pot + vt
          end if
        end do
      end do
    end do

    ! H2C 3b
    istart = istart+np21
    iend = iend+np21
    do j=nc+1,natm-1
      do k=j+1,natm
        if (rint(j,k) > r3f) cycle
        morse3(1)=exp(-rint(j,k) / a_h2c)
        do l=1,nc
          if (rint(j,l) > r3f .or. rint(k,l) > r3f) cycle
          morse3(2)=exp(-rint(j,l) / a_h2c)
          morse3(3)=exp(-rint(k,l) / a_h2c)
          vt = emsav21(morse3, coef(istart:iend))
          rmax = max(rint(j,k), rint(j,l), rint(k,l))
          if (rmax > r3i) then
            call f_switch(s, rmax, r3i, r3f)
            pot = pot + vt*s
          else
            pot = pot + vt
          end if
        end do
      end do
    end do

    ! HHH 3b
    istart = istart+np21
    iend = iend+np3
    do j=nc+1,natm-2
      do k=j+1,natm-1
        if (rint(j,k) > r3f) cycle
        morse3(1)=exp(-rint(j,k) / a_h3)
        do l=k+1,natm
          if (rint(j,l) > r3f .or. rint(k,l) > r3f) cycle
          morse3(2)=exp(-rint(j,l) / a_h3)
          morse3(3)=exp(-rint(k,l) / a_h3)
          vt = emsav3(morse3, coef(istart:iend))
          rmax = max(rint(j,k), rint(j,l), rint(k,l))
          if (rmax > r3i) then
            call f_switch(s, rmax, r3i, r3f)
            pot = pot + vt*s
          else
            pot = pot + vt
          end if
        end do
      end do
    end do

    ! CCCC 4-b
    n4b = 0
    istart = istart+np3
    iend = iend+np4
    do j=1,nc-3
      do k=j+1,nc-2
        if (rint(j,k) > r4f) cycle
        morse4(1)=exp(-rint(j,k) / a_c4)
        do l=k+1,nc-1
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          morse4(2)=exp(-rint(j,l) / a_c4)
          morse4(4)=exp(-rint(k,l) / a_c4)
          do m=l+1,nc
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            morse4(3)=exp(-rint(j,m) / a_c4)
            morse4(5)=exp(-rint(k,m) / a_c4)
            morse4(6)=exp(-rint(l,m) / a_c4)
            vt = emsav4(morse4, coef(istart:iend))
            n4b = n4b + 1
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
            else
              pot = pot + vt
            end if
          end do
        end do
      end do
    end do

    ! CCCH 4b
    istart = istart+np4
    iend = iend+np31
    do j=1,nc-2
      do k=j+1,nc-1
        if (rint(j,k) > r4f) cycle
        morse4(1)=exp(-rint(j,k) / a_c3h)
        do l=k+1,nc
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          morse4(2)=exp(-rint(j,l) / a_c3h)
          morse4(4)=exp(-rint(k,l) / a_c3h)
          do m=nc+1,natm
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            morse4(3)=exp(-rint(j,m) / a_c3h)
            morse4(5)=exp(-rint(k,m) / a_c3h)
            morse4(6)=exp(-rint(l,m) / a_c3h)
            vt = emsav31(morse4, coef(istart:iend))
            n4b = n4b + 1
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
            else
              pot = pot + vt
            end if
          end do
        end do
      end do
    end do

    ! CCHH 4b
    istart = istart+np31
    iend = iend+np22
    do j=1,nc-1
      do k=j+1,nc
        if (rint(j,k) > r4f) cycle
        morse4(1)=exp(-rint(j,k) / a_c2h2)
        do l=nc+1,natm-1
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          morse4(2)=exp(-rint(j,l) / a_c2h2)
          morse4(4)=exp(-rint(k,l) / a_c2h2)
          do m=l+1,natm
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            morse4(3)=exp(-rint(j,m) / a_c2h2)
            morse4(5)=exp(-rint(k,m) / a_c2h2)
            morse4(6)=exp(-rint(l,m) / a_c2h2)
            vt = emsav22(morse4, coef(istart:iend))
            n4b = n4b + 1
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
            else
              pot = pot + vt
            end if
          end do
        end do
      end do
    end do

    ! CHHH 4b
    istart = istart+np22
    iend = iend+np31
    do j=nc+1,natm-2
      do k=j+1,natm-1
        if (rint(j,k) > r4f) cycle
        morse4(1)=exp(-rint(j,k) / a_h3c)
        do l=k+1,natm
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          morse4(2)=exp(-rint(j,l) / a_h3c)
          morse4(4)=exp(-rint(k,l) / a_h3c)
          do m=1,nc
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            morse4(3)=exp(-rint(j,m) / a_h3c)
            morse4(5)=exp(-rint(k,m) / a_h3c)
            morse4(6)=exp(-rint(l,m) / a_h3c)
            vt = emsav31(morse4, coef(istart:iend))
            n4b = n4b + 1
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
            else
              pot = pot + vt
            end if
          end do
        end do
      end do
    end do

    istart = istart+np31
    iend = iend+np4
    do j=nc+1,natm-3
      do k=j+1,natm-2
        if (rint(j,k) > r4f) cycle
        morse4(1)=exp(-rint(j,k) / a_h4)
        do l=k+1,natm-1
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          morse4(2)=exp(-rint(j,l) / a_h4)
          morse4(4)=exp(-rint(k,l) / a_h4)
          do m=l+1,natm
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            morse4(3)=exp(-rint(j,m) / a_h4)
            morse4(5)=exp(-rint(k,m) / a_h4)
            morse4(6)=exp(-rint(l,m) / a_h4)
            vt = emsav4(morse4, coef(istart:iend))
            n4b = n4b + 1
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
            else
              pot = pot + vt
            end if
          end do
        end do
      end do
    end do
!    write(*,*) n4b

    return
  end subroutine getpot

  !=================================!
  ! Energy and Analytical gradients !
  !=================================!
  subroutine pot_grad(x1d, pot, grad1d)
    real,dimension(:),intent(in)::x1d
    real,intent(out)::pot
    real,dimension(:),intent(out)::grad1d
    !::::::::::::::::::::
    real,dimension(3,size(x1d)/3)::x
    real,dimension(3,size(grad1d)/3)::grad
    real::rint(size(x1d)/3, size(x1d)/3),morse3(3),morse4(6)
    real::rint3(3,3), rint4(4,4)
    real::xyz3(3,3), xyz4(4,3)
    integer::j,k,l,m,istart,iend,ig,natm
    real::vt,s,rmax,g2(6),g3(9),g4(12)

    natm = size(x1d)/3
    do j=1,natm
      x(:,j) = x1d(3*j-2:3*j)
    end do

    do j=1,natm-1
      do k=j+1,natm
        rint(j,k) = norm2(x(:,j)-x(:,k))
        rint(k,j) = rint(j,k)
      end do
    end do

    !==============================
    ! compute the predicted energy
    !==============================
    pot = -37.8519753684*real(nc) - 0.501257936933*real(nh)
    grad = 0.0

    ! C-C 2-body
    istart=1
    iend = np2+istart-1
    do j=1,nc-1
      do k=j+1,nc
        if (rint(j,k) > r2f) cycle
        call evp2(rint(j,k), a_cc, p2)
        vt = dot_product(coef(istart:iend), p2)
        do ig = 1,6
          call evdp2(x(:,j), x(:,k), rint(j,k), a_cc, ig, p2, dp2)
          g2(ig) = dot_product(coef(istart:iend), dp2)
        end do
        if (rint(j,k) > r2i) then
          call f_switch(s, rint(j,k), r2i, r2f)
          pot = pot + vt*s
          grad(:,j) = grad(:,j) + g2(1:3)*s
          grad(:,k) = grad(:,k) + g2(4:6)*s
        else
          pot = pot + vt
          grad(:,j) = grad(:,j) + g2(1:3)
          grad(:,k) = grad(:,k) + g2(4:6)
        end if
      end do
    end do

    ! C-H 2-body
    istart=istart+np2
    iend = iend+np2
    do j=1,nc
      do k=nc+1,natm
        if (rint(j,k) > r2f) cycle
        call evp2(rint(j,k), a_ch, p2)
        vt = dot_product(coef(istart:iend), p2)
        do ig =1,6
          call evdp2(x(:,j), x(:,k), rint(j,k), a_ch, ig, p2, dp2)
          g2(ig) = dot_product(coef(istart:iend), dp2)
        end do
        if (rint(j,k) > r2i) then
          call f_switch(s, rint(j,k), r2i, r2f)
          pot = pot + vt*s
          grad(:,j) = grad(:,j) + g2(1:3)*s
          grad(:,k) = grad(:,k) + g2(4:6)*s
        else
          pot = pot + vt
          grad(:,j) = grad(:,j) + g2(1:3)
          grad(:,k) = grad(:,k) + g2(4:6)
        end if
      end do
    end do

    ! H-H 2-body
    istart=istart+np2
    iend = iend+np2
    do j=nc+1,natm-1
      do k=j+1,natm
        if (rint(j,k) > r2f) cycle
        call evp2(rint(j,k), a_hh, p2)
        vt = dot_product(coef(istart:iend), p2)
        do ig = 1,6
          call evdp2(x(:,j), x(:,k), rint(j,k), a_hh, ig, p2, dp2)
          g2(ig) = dot_product(coef(istart:iend), dp2)
        end do
        if (rint(j,k) > r2i) then
          call f_switch(s, rint(j,k), r2i, r2f)
          pot = pot + vt*s
          grad(:,j) = grad(:,j) + g2(1:3)*s
          grad(:,k) = grad(:,k) + g2(4:6)*s
        else
          pot = pot + vt
          grad(:,j) = grad(:,j) + g2(1:3)
          grad(:,k) = grad(:,k) + g2(4:6)
        end if
      end do
    end do

    ! CCC 3b
    istart = istart+np2
    iend = iend+np3
    do j=1,nc-2
      xyz3(1,:) = x(:,j)
      do k=j+1,nc-1
        if (rint(j,k) > r3f) cycle
        rint3(1,2) = rint(j,k)
        rint3(2,1) = rint(j,k)
        xyz3(2,:) = x(:,k)
        morse3(1)=exp(-rint(j,k) / a_c3)
        do l=k+1,nc
          if (rint(j,l) > r3f .or. rint(k,l) > r3f) cycle
          rint3(1,3) = rint(j,l)
          rint3(3,1) = rint(j,l)
          rint3(2,3) = rint(k,l)
          rint3(3,2) = rint(k,l)
          xyz3(3,:) = x(:,l)
          morse3(2)=exp(-rint(j,l) / a_c3)
          morse3(3)=exp(-rint(k,l) / a_c3)
          call evmono3(morse3, m3)
          call evpoly3(m3, q3, p3)
          vt = dot_product(coef(istart:iend), p3)
          call deriv_rev3(coef(istart:iend),m3,q3,p3,xyz3,a_c3,rint3,g3)
          rmax = max(rint(j,k), rint(j,l), rint(k,l))
          if (rmax > r3i) then
            call f_switch(s, rmax, r3i, r3f)
            pot = pot + vt*s
            grad(:,j) = grad(:,j) + g3(1:3)*s
            grad(:,k) = grad(:,k) + g3(4:6)*s
            grad(:,l) = grad(:,l) + g3(7:9)*s
          else
            pot = pot + vt
            grad(:,j) = grad(:,j) + g3(1:3)
            grad(:,k) = grad(:,k) + g3(4:6)
            grad(:,l) = grad(:,l) + g3(7:9)
          end if
        end do
      end do
    end do

    ! C2H 3b
    istart = istart+np3
    iend = iend+np21
    do j=1,nc-1
      xyz3(1,:) = x(:,j)
      do k=j+1,nc
        if (rint(j,k) > r3f) cycle
        rint3(1,2) = rint(j,k)
        rint3(2,1) = rint(j,k)
        xyz3(2,:) = x(:,k)
        morse3(1)=exp(-rint(j,k) / a_c2h)
        do l=nc+1,natm
          if (rint(j,l) > r3f .or. rint(k,l) > r3f) cycle
          rint3(1,3) = rint(j,l)
          rint3(3,1) = rint(j,l)
          rint3(2,3) = rint(k,l)
          rint3(3,2) = rint(k,l)
          xyz3(3,:) = x(:,l)
          morse3(2)=exp(-rint(j,l) / a_c2h)
          morse3(3)=exp(-rint(k,l) / a_c2h)
          call evmono21(morse3, m21)
          call evpoly21(m21, q21, p21)
          vt = dot_product(coef(istart:iend), p21)
          call deriv_rev21(coef(istart:iend),m21,q21,p21,xyz3,a_c2h,rint3,g3)
          rmax = max(rint(j,k), rint(j,l), rint(k,l))
          if (rmax > r3i) then
            call f_switch(s, rmax, r3i, r3f)
            pot = pot + vt*s
            grad(:,j) = grad(:,j) + g3(1:3)*s
            grad(:,k) = grad(:,k) + g3(4:6)*s
            grad(:,l) = grad(:,l) + g3(7:9)*s
          else
            pot = pot + vt
            grad(:,j) = grad(:,j) + g3(1:3)
            grad(:,k) = grad(:,k) + g3(4:6)
            grad(:,l) = grad(:,l) + g3(7:9)
          end if
        end do
      end do
    end do

    ! H2C 3b
    istart = istart+np21
    iend = iend+np21
    do j=nc+1,natm-1
      xyz3(1,:) = x(:,j)
      do k=j+1,natm
        if (rint(j,k) > r3f) cycle
        rint3(1,2) = rint(j,k)
        rint3(2,1) = rint(j,k)
        xyz3(2,:) = x(:,k)
        morse3(1)=exp(-rint(j,k) / a_h2c)
        do l=1,nc
          if (rint(j,l) > r3f .or. rint(k,l) > r3f) cycle
          rint3(1,3) = rint(j,l)
          rint3(3,1) = rint(j,l)
          rint3(2,3) = rint(k,l)
          rint3(3,2) = rint(k,l)
          xyz3(3,:) = x(:,l)
          morse3(2)=exp(-rint(j,l) / a_h2c)
          morse3(3)=exp(-rint(k,l) / a_h2c)
          call evmono21(morse3, m21)
          call evpoly21(m21, q21, p21)
          vt = dot_product(coef(istart:iend), p21)
          call deriv_rev21(coef(istart:iend),m21,q21,p21,xyz3,a_h2c,rint3,g3)
          rmax = max(rint(j,k), rint(j,l), rint(k,l))
          if (rmax > r3i) then
            call f_switch(s, rmax, r3i, r3f)
            pot = pot + vt*s
            grad(:,j) = grad(:,j) + g3(1:3)*s
            grad(:,k) = grad(:,k) + g3(4:6)*s
            grad(:,l) = grad(:,l) + g3(7:9)*s
          else
            pot = pot + vt
            grad(:,j) = grad(:,j) + g3(1:3)
            grad(:,k) = grad(:,k) + g3(4:6)
            grad(:,l) = grad(:,l) + g3(7:9)
          end if
        end do
      end do
    end do

    ! HHH 3b
    istart = istart+np21
    iend = iend+np3
    do j=nc+1,natm-2
      xyz3(1,:) = x(:,j)
      do k=j+1,natm-1
        if (rint(j,k) > r3f) cycle
        rint3(1,2) = rint(j,k)
        rint3(2,1) = rint(j,k)
        xyz3(2,:) = x(:,k)
        morse3(1)=exp(-rint(j,k) / a_h3)
        do l=k+1,natm
          if (rint(j,l) > r3f .or. rint(k,l) > r3f) cycle
          rint3(1,3) = rint(j,l)
          rint3(3,1) = rint(j,l)
          rint3(2,3) = rint(k,l)
          rint3(3,2) = rint(k,l)
          xyz3(3,:) = x(:,l)
          morse3(2)=exp(-rint(j,l) / a_h3)
          morse3(3)=exp(-rint(k,l) / a_h3)
          call evmono3(morse3, m3)
          call evpoly3(m3, q3, p3)
          vt = dot_product(coef(istart:iend), p3)
          call deriv_rev3(coef(istart:iend),m3,q3,p3,xyz3,a_h3,rint3,g3)
          rmax = max(rint(j,k), rint(j,l), rint(k,l))
          if (rmax > r3i) then
            call f_switch(s, rmax, r3i, r3f)
            pot = pot + vt*s
            grad(:,j) = grad(:,j) + g3(1:3)*s
            grad(:,k) = grad(:,k) + g3(4:6)*s
            grad(:,l) = grad(:,l) + g3(7:9)*s
          else
            pot = pot + vt
            grad(:,j) = grad(:,j) + g3(1:3)
            grad(:,k) = grad(:,k) + g3(4:6)
            grad(:,l) = grad(:,l) + g3(7:9)
          end if
        end do
      end do
    end do

    ! CCCC 4-b
    istart = istart+np3
    iend = iend+np4
    do j=1,nc-3
      xyz4(1,:) = x(:,j)
      do k=j+1,nc-2
        if (rint(j,k) > r4f) cycle
        rint4(1,2) = rint(j,k)
        rint4(2,1) = rint(j,k)
        xyz4(2,:) = x(:,k)
        morse4(1)=exp(-rint(j,k) / a_c4)
        do l=k+1,nc-1
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          rint4(1,3) = rint(j,l)
          rint4(3,1) = rint(j,l)
          rint4(2,3) = rint(k,l)
          rint4(3,2) = rint(k,l)
          xyz4(3,:) = x(:,l)
          morse4(2)=exp(-rint(j,l) / a_c4)
          morse4(4)=exp(-rint(k,l) / a_c4)
          do m=l+1,nc
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            rint4(1,4) = rint(j,m)
            rint4(4,1) = rint(j,m)
            rint4(2,4) = rint(k,m)
            rint4(4,2) = rint(k,m)
            rint4(3,4) = rint(l,m)
            rint4(4,3) = rint(l,m)
            xyz4(4,:) = x(:,m)
            morse4(3)=exp(-rint(j,m) / a_c4)
            morse4(5)=exp(-rint(k,m) / a_c4)
            morse4(6)=exp(-rint(l,m) / a_c4)
            call evmono4(morse4, m4)
            call evpoly4(m4, q4, p4)
            vt = dot_product(coef(istart:iend), p4)
            call deriv_rev4(coef(istart:iend),m4,q4,p4,xyz4,a_c4,rint4,g4)
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
              grad(:,j) = grad(:,j) + g4(1:3)*s
              grad(:,k) = grad(:,k) + g4(4:6)*s
              grad(:,l) = grad(:,l) + g4(7:9)*s
              grad(:,m) = grad(:,m) + g4(10:12)*s
            else
              pot = pot + vt
              grad(:,j) = grad(:,j) + g4(1:3)
              grad(:,k) = grad(:,k) + g4(4:6)
              grad(:,l) = grad(:,l) + g4(7:9)
              grad(:,m) = grad(:,m) + g4(10:12)
            end if
          end do
        end do
      end do
    end do

    ! CCCH 4b
    istart = istart+np4
    iend = iend+np31
    do j=1,nc-2
      xyz4(1,:) = x(:,j)
      do k=j+1,nc-1
        if (rint(j,k) > r4f) cycle
        rint4(1,2) = rint(j,k)
        rint4(2,1) = rint(j,k)
        xyz4(2,:) = x(:,k)
        morse4(1)=exp(-rint(j,k) / a_c3h)
        do l=k+1,nc
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          rint4(1,3) = rint(j,l)
          rint4(3,1) = rint(j,l)
          rint4(2,3) = rint(k,l)
          rint4(3,2) = rint(k,l)
          xyz4(3,:) = x(:,l)
          morse4(2)=exp(-rint(j,l) / a_c3h)
          morse4(4)=exp(-rint(k,l) / a_c3h)
          do m=nc+1,natm
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            rint4(1,4) = rint(j,m)
            rint4(4,1) = rint(j,m)
            rint4(2,4) = rint(k,m)
            rint4(4,2) = rint(k,m)
            rint4(3,4) = rint(l,m)
            rint4(4,3) = rint(l,m)
            xyz4(4,:) = x(:,m)
            morse4(3)=exp(-rint(j,m) / a_c3h)
            morse4(5)=exp(-rint(k,m) / a_c3h)
            morse4(6)=exp(-rint(l,m) / a_c3h)
            call evmono31(morse4, m31)
            call evpoly31(m31, q31, p31)
            vt = dot_product(coef(istart:iend), p31)
            call deriv_rev31(coef(istart:iend),m31,q31,p31,xyz4,a_c3h,rint4,g4)
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
              grad(:,j) = grad(:,j) + g4(1:3)*s
              grad(:,k) = grad(:,k) + g4(4:6)*s
              grad(:,l) = grad(:,l) + g4(7:9)*s
              grad(:,m) = grad(:,m) + g4(10:12)*s
            else
              pot = pot + vt
              grad(:,j) = grad(:,j) + g4(1:3)
              grad(:,k) = grad(:,k) + g4(4:6)
              grad(:,l) = grad(:,l) + g4(7:9)
              grad(:,m) = grad(:,m) + g4(10:12)
            end if
          end do
        end do
      end do
    end do

    ! CCHH 4b
    istart = istart+np31
    iend = iend+np22
    do j=1,nc-1
      xyz4(1,:) = x(:,j)
      do k=j+1,nc
        if (rint(j,k) > r4f) cycle
        rint4(1,2) = rint(j,k)
        rint4(2,1) = rint(j,k)
        xyz4(2,:) = x(:,k)
        morse4(1)=exp(-rint(j,k) / a_c2h2)
        do l=nc+1,natm-1
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          rint4(1,3) = rint(j,l)
          rint4(3,1) = rint(j,l)
          rint4(2,3) = rint(k,l)
          rint4(3,2) = rint(k,l)
          xyz4(3,:) = x(:,l)
          morse4(2)=exp(-rint(j,l) / a_c2h2)
          morse4(4)=exp(-rint(k,l) / a_c2h2)
          do m=l+1,natm
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            rint4(1,4) = rint(j,m)
            rint4(4,1) = rint(j,m)
            rint4(2,4) = rint(k,m)
            rint4(4,2) = rint(k,m)
            rint4(3,4) = rint(l,m)
            rint4(4,3) = rint(l,m)
            xyz4(4,:) = x(:,m)
            morse4(3)=exp(-rint(j,m) / a_c2h2)
            morse4(5)=exp(-rint(k,m) / a_c2h2)
            morse4(6)=exp(-rint(l,m) / a_c2h2)
            call evmono22(morse4, m22)
            call evpoly22(m22, q22, p22)
            vt = dot_product(coef(istart:iend), p22)
            call deriv_rev22(coef(istart:iend),m22,q22,p22,xyz4,a_c2h2,rint4,g4)
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
              grad(:,j) = grad(:,j) + g4(1:3)*s
              grad(:,k) = grad(:,k) + g4(4:6)*s
              grad(:,l) = grad(:,l) + g4(7:9)*s
              grad(:,m) = grad(:,m) + g4(10:12)*s
            else
              pot = pot + vt
              grad(:,j) = grad(:,j) + g4(1:3)
              grad(:,k) = grad(:,k) + g4(4:6)
              grad(:,l) = grad(:,l) + g4(7:9)
              grad(:,m) = grad(:,m) + g4(10:12)
            end if
          end do
        end do
      end do
    end do

    ! CHHH 4b
    istart = istart+np22
    iend = iend+np31
    do j=nc+1,natm-2
      xyz4(1,:) = x(:,j)
      do k=j+1,natm-1
        if (rint(j,k) > r4f) cycle
        rint4(1,2) = rint(j,k)
        rint4(2,1) = rint(j,k)
        xyz4(2,:) = x(:,k)
        morse4(1)=exp(-rint(j,k) / a_h3c)
        do l=k+1,natm
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          rint4(1,3) = rint(j,l)
          rint4(3,1) = rint(j,l)
          rint4(2,3) = rint(k,l)
          rint4(3,2) = rint(k,l)
          xyz4(3,:) = x(:,l)
          morse4(2)=exp(-rint(j,l) / a_h3c)
          morse4(4)=exp(-rint(k,l) / a_h3c)
          do m=1,nc
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            rint4(1,4) = rint(j,m)
            rint4(4,1) = rint(j,m)
            rint4(2,4) = rint(k,m)
            rint4(4,2) = rint(k,m)
            rint4(3,4) = rint(l,m)
            rint4(4,3) = rint(l,m)
            xyz4(4,:) = x(:,m)
            morse4(3)=exp(-rint(j,m) / a_h3c)
            morse4(5)=exp(-rint(k,m) / a_h3c)
            morse4(6)=exp(-rint(l,m) / a_h3c)
            call evmono31(morse4, m31)
            call evpoly31(m31, q31, p31)
            vt = dot_product(coef(istart:iend), p31)
            call deriv_rev31(coef(istart:iend),m31,q31,p31,xyz4,a_h3c,rint4,g4)
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
              grad(:,j) = grad(:,j) + g4(1:3)*s
              grad(:,k) = grad(:,k) + g4(4:6)*s
              grad(:,l) = grad(:,l) + g4(7:9)*s
              grad(:,m) = grad(:,m) + g4(10:12)*s
            else
              pot = pot + vt
              grad(:,j) = grad(:,j) + g4(1:3)
              grad(:,k) = grad(:,k) + g4(4:6)
              grad(:,l) = grad(:,l) + g4(7:9)
              grad(:,m) = grad(:,m) + g4(10:12)
            end if
          end do
        end do
      end do
    end do

    ! HHHH 4b
    istart = istart+np31
    iend = iend+np4
    do j=nc+1,natm-3
      xyz4(1,:) = x(:,j)
      do k=j+1,natm-2
        if (rint(j,k) > r4f) cycle
        rint4(1,2) = rint(j,k)
        rint4(2,1) = rint(j,k)
        xyz4(2,:) = x(:,k)
        morse4(1)=exp(-rint(j,k) / a_h4)
        do l=k+1,natm-1
          if (rint(j,l) > r4f .or. rint(k,l) > r4f) cycle
          rint4(1,3) = rint(j,l)
          rint4(3,1) = rint(j,l)
          rint4(2,3) = rint(k,l)
          rint4(3,2) = rint(k,l)
          xyz4(3,:) = x(:,l)
          morse4(2)=exp(-rint(j,l) / a_h4)
          morse4(4)=exp(-rint(k,l) / a_h4)
          do m=l+1,natm
            if (rint(j,m)>r4f .or. rint(k,m)>r4f .or. rint(l,m)>r4f) cycle
            rint4(1,4) = rint(j,m)
            rint4(4,1) = rint(j,m)
            rint4(2,4) = rint(k,m)
            rint4(4,2) = rint(k,m)
            rint4(3,4) = rint(l,m)
            rint4(4,3) = rint(l,m)
            xyz4(4,:) = x(:,m)
            morse4(3)=exp(-rint(j,m) / a_h4)
            morse4(5)=exp(-rint(k,m) / a_h4)
            morse4(6)=exp(-rint(l,m) / a_h4)
            call evmono4(morse4, m4)
            call evpoly4(m4, q4, p4)
            vt = dot_product(coef(istart:iend), p4)
            call deriv_rev4(coef(istart:iend),m4,q4,p4,xyz4,a_h4,rint4,g4)
            rmax = max(rint(j,k), rint(j,l), rint(k,l), rint(j,m), rint(k,m), rint(l,m))
            if (rmax > r4i) then
              call f_switch(s, rmax, r4i, r4f)
              pot = pot + vt*s
              grad(:,j) = grad(:,j) + g4(1:3)*s
              grad(:,k) = grad(:,k) + g4(4:6)*s
              grad(:,l) = grad(:,l) + g4(7:9)*s
              grad(:,m) = grad(:,m) + g4(10:12)*s
            else
              pot = pot + vt
              grad(:,j) = grad(:,j) + g4(1:3)
              grad(:,k) = grad(:,k) + g4(4:6)
              grad(:,l) = grad(:,l) + g4(7:9)
              grad(:,m) = grad(:,m) + g4(10:12)
            end if
          end do
        end do
      end do
    end do

    do j=1,natm
      grad1d(3*j-2:3*j) = grad(:,j)
    end do

    return
  end subroutine pot_grad

  !=========================================!
  ! Numerical Hessian using analytical grad !
  !=========================================!
  subroutine hessian(x,H)
    real,dimension(:),intent(in)::x
    real,dimension(:,:),intent(inout)::H
    !:::::::::::::::::::::::::::
    real,dimension(1:size(x))::xt
    real,dimension(size(x),size(x))::gf,gb
    real,parameter::eps=0.005
    integer::ndim,i,j
    real::pot

    ndim = size(x)

    do i=1,ndim
      xt = x
      xt(i) = x(i) + eps
      call pot_grad(xt, pot, gf(:,i))

      xt = x
      xt(i) = x(i) - eps
      call pot_grad(xt, pot, gb(:,i))
    end do

    do i=1,ndim-1
      do j=i+1,ndim
        H(i,j) = (gf(i,j)-gb(i,j))/(4.0*eps) + (gf(j,i)-gb(j,i))/(4.0*eps)
        H(j,i) = H(i,j)
      end do
    end do

    do i=1,ndim
      H(i,i) = (gf(i,i)-gb(i,i))/(2.0*eps)
    end do

    return
  end subroutine hessian

  !======================!
  ! Numerical gradients  !
  !======================!
  subroutine numerical_grad(x, grad)
    real,dimension(:),intent(out)::grad
    real,dimension(:),intent(in)::x
    !::::::::::::::::::::::::::::
    real,dimension(size(x))::xt
    real::eps=0.001
    integer::i,j,k,natm
    real::vf,vb

    natm = size(x)/3
    do i=1,natm*3
      ! forward
      xt = x
      xt(i) = xt(i) + eps
      call getpot(xt, vf)

      ! backward
      xt = x
      xt(i) = xt(i) - eps
      call getpot(xt, vb)

      grad(i) = (vf - vb) / eps / 2.0
    end do

    return
  end subroutine numerical_grad

end module
