module bemsa21
  implicit none

contains
  function emsav21(x,c) result(v)
    real,dimension(1:3)::x
    real,dimension(1:78)::p
    real,dimension(1:78)::c
    real::v
    
    call bemsav21(x,p)
    v = dot_product(p,c)
    
    return
  end function emsav21
  
  subroutine bemsav21(x,p)
    real,dimension(1:3)::x
    real,dimension(1:4)::m
    real,dimension(1:78)::p
    real,dimension(1:8)::q
    ! ::::::::::::::::::::

    call evmono21(x,m)
    call evpoly21(m,q,p)
    
    return
  end subroutine bemsav21
  
  subroutine evmono21(x,m)
    real,dimension(1:3)::x
    real,dimension(1:4)::m
    !::::::::::::::::::::
  
    m(1) = x(3)
    m(2) = x(2)
    m(3) = x(1)
    m(4) = m(1)*m(2)

    return
  end subroutine evmono21

  subroutine evpoly21(m,q,p)
    real,dimension(1:4)::m
    real,dimension(1:78)::p
    real,dimension(1:8)::q
    !::::::::::::::::::::

    q(1) = m(1) + m(2)
    q(2) = m(3)
    p(1) = m(4)
    p(2) = q(2)*q(1)
    q(3) = q(1)*q(1) - p(1) - p(1)
    p(3) = q(2)*p(1)
    p(4) = p(1)*q(1)
    p(5) = q(2)*q(3)
    p(6) = q(2)*p(2)
    q(4) = q(1)*q(3) - p(4)
    p(7) = q(2)*p(4)
    p(8) = q(2)*p(3)
    p(9) = p(1)*p(1)
    p(10) = p(1)*q(3)
    p(11) = q(2)*q(4)
    p(12) = q(2)*p(5)
    p(13) = q(2)*p(6)
    q(5) = q(1)*q(4) - p(10)
    p(14) = q(2)*p(9)
    p(15) = q(2)*p(10)
    p(16) = q(2)*p(7)
    p(17) = q(2)*p(8)
    p(18) = p(1)*p(4)
    p(19) = p(1)*q(4)
    p(20) = q(2)*q(5)
    p(21) = q(2)*p(11)
    p(22) = q(2)*p(12)
    p(23) = q(2)*p(13)
    q(6) = q(1)*q(5) - p(19)
    p(24) = q(2)*p(18)
    p(25) = q(2)*p(19)
    p(26) = q(2)*p(14)
    p(27) = q(2)*p(15)
    p(28) = q(2)*p(16)
    p(29) = q(2)*p(17)
    p(30) = p(1)*p(9)
    p(31) = p(1)*p(10)
    p(32) = p(1)*q(5)
    p(33) = q(2)*q(6)
    p(34) = q(2)*p(20)
    p(35) = q(2)*p(21)
    p(36) = q(2)*p(22)
    p(37) = q(2)*p(23)
    q(7) = q(1)*q(6) - p(32)
    p(38) = q(2)*p(30)
    p(39) = q(2)*p(31)
    p(40) = q(2)*p(32)
    p(41) = q(2)*p(24)
    p(42) = q(2)*p(25)
    p(43) = q(2)*p(26)
    p(44) = q(2)*p(27)
    p(45) = q(2)*p(28)
    p(46) = q(2)*p(29)
    p(47) = p(1)*p(18)
    p(48) = p(1)*p(19)
    p(49) = p(1)*q(6)
    p(50) = q(2)*q(7)
    p(51) = q(2)*p(33)
    p(52) = q(2)*p(34)
    p(53) = q(2)*p(35)
    p(54) = q(2)*p(36)
    p(55) = q(2)*p(37)
    q(8) = q(1)*q(7) - p(49)
    p(56) = q(2)*p(47)
    p(57) = q(2)*p(48)
    p(58) = q(2)*p(49)
    p(59) = q(2)*p(38)
    p(60) = q(2)*p(39)
    p(61) = q(2)*p(40)
    p(62) = q(2)*p(41)
    p(63) = q(2)*p(42)
    p(64) = q(2)*p(43)
    p(65) = q(2)*p(44)
    p(66) = q(2)*p(45)
    p(67) = q(2)*p(46)
    p(68) = p(1)*p(30)
    p(69) = p(1)*p(31)
    p(70) = p(1)*p(32)
    p(71) = p(1)*q(7)
    p(72) = q(2)*q(8)
    p(73) = q(2)*p(50)
    p(74) = q(2)*p(51)
    p(75) = q(2)*p(52)
    p(76) = q(2)*p(53)
    p(77) = q(2)*p(54)
    p(78) = q(2)*p(55)

    return
  end subroutine evpoly21

  subroutine devmono21(dm,m,flag,xyz,a,r)
    real,dimension(1:4)::dm
    real,dimension(1:4)::m
    real::xyz(3,3),a,r(3,3)
    integer::flag
    !::::::::::::::::::::

    dm(1) =  - m(1)/a*drdx21(flag,3,xyz,r)
    dm(2) =  - m(2)/a*drdx21(flag,2,xyz,r)
    dm(3) =  - m(3)/a*drdx21(flag,1,xyz,r)
    dm(4) = dm(1)*m(2) + m(1)*dm(2)

    return
  end subroutine devmono21

  subroutine devpoly21(dm,q,p,dp)
    real,dimension(1:4)::dm
    real,dimension(1:78)::p,dp
    real,dimension(1:8)::q,dq
    !::::::::::::::::::::

    dq(1) = dm(1) + dm(2)
    dq(2) = dm(3)
    dp(1) = dm(4)
    dp(2) = dq(2)*q(1) + q(2)*dq(1)
    dq(3) = dq(1)*q(1) + q(1)*dq(1) - dp(1) - dp(1)
    dp(3) = dq(2)*p(1) + q(2)*dp(1)
    dp(4) = dp(1)*q(1) + p(1)*dq(1)
    dp(5) = dq(2)*q(3) + q(2)*dq(3)
    dp(6) = dq(2)*p(2) + q(2)*dp(2)
    dq(4) = dq(1)*q(3) + q(1)*dq(3) - dp(4)
    dp(7) = dq(2)*p(4) + q(2)*dp(4)
    dp(8) = dq(2)*p(3) + q(2)*dp(3)
    dp(9) = dp(1)*p(1) + p(1)*dp(1)
    dp(10) = dp(1)*q(3) + p(1)*dq(3)
    dp(11) = dq(2)*q(4) + q(2)*dq(4)
    dp(12) = dq(2)*p(5) + q(2)*dp(5)
    dp(13) = dq(2)*p(6) + q(2)*dp(6)
    dq(5) = dq(1)*q(4) + q(1)*dq(4) - dp(10)
    dp(14) = dq(2)*p(9) + q(2)*dp(9)
    dp(15) = dq(2)*p(10) + q(2)*dp(10)
    dp(16) = dq(2)*p(7) + q(2)*dp(7)
    dp(17) = dq(2)*p(8) + q(2)*dp(8)
    dp(18) = dp(1)*p(4) + p(1)*dp(4)
    dp(19) = dp(1)*q(4) + p(1)*dq(4)
    dp(20) = dq(2)*q(5) + q(2)*dq(5)
    dp(21) = dq(2)*p(11) + q(2)*dp(11)
    dp(22) = dq(2)*p(12) + q(2)*dp(12)
    dp(23) = dq(2)*p(13) + q(2)*dp(13)
    dq(6) = dq(1)*q(5) + q(1)*dq(5) - dp(19)
    dp(24) = dq(2)*p(18) + q(2)*dp(18)
    dp(25) = dq(2)*p(19) + q(2)*dp(19)
    dp(26) = dq(2)*p(14) + q(2)*dp(14)
    dp(27) = dq(2)*p(15) + q(2)*dp(15)
    dp(28) = dq(2)*p(16) + q(2)*dp(16)
    dp(29) = dq(2)*p(17) + q(2)*dp(17)
    dp(30) = dp(1)*p(9) + p(1)*dp(9)
    dp(31) = dp(1)*p(10) + p(1)*dp(10)
    dp(32) = dp(1)*q(5) + p(1)*dq(5)
    dp(33) = dq(2)*q(6) + q(2)*dq(6)
    dp(34) = dq(2)*p(20) + q(2)*dp(20)
    dp(35) = dq(2)*p(21) + q(2)*dp(21)
    dp(36) = dq(2)*p(22) + q(2)*dp(22)
    dp(37) = dq(2)*p(23) + q(2)*dp(23)
    dq(7) = dq(1)*q(6) + q(1)*dq(6) - dp(32)
    dp(38) = dq(2)*p(30) + q(2)*dp(30)
    dp(39) = dq(2)*p(31) + q(2)*dp(31)
    dp(40) = dq(2)*p(32) + q(2)*dp(32)
    dp(41) = dq(2)*p(24) + q(2)*dp(24)
    dp(42) = dq(2)*p(25) + q(2)*dp(25)
    dp(43) = dq(2)*p(26) + q(2)*dp(26)
    dp(44) = dq(2)*p(27) + q(2)*dp(27)
    dp(45) = dq(2)*p(28) + q(2)*dp(28)
    dp(46) = dq(2)*p(29) + q(2)*dp(29)
    dp(47) = dp(1)*p(18) + p(1)*dp(18)
    dp(48) = dp(1)*p(19) + p(1)*dp(19)
    dp(49) = dp(1)*q(6) + p(1)*dq(6)
    dp(50) = dq(2)*q(7) + q(2)*dq(7)
    dp(51) = dq(2)*p(33) + q(2)*dp(33)
    dp(52) = dq(2)*p(34) + q(2)*dp(34)
    dp(53) = dq(2)*p(35) + q(2)*dp(35)
    dp(54) = dq(2)*p(36) + q(2)*dp(36)
    dp(55) = dq(2)*p(37) + q(2)*dp(37)
    dq(8) = dq(1)*q(7) + q(1)*dq(7) - dp(49)
    dp(56) = dq(2)*p(47) + q(2)*dp(47)
    dp(57) = dq(2)*p(48) + q(2)*dp(48)
    dp(58) = dq(2)*p(49) + q(2)*dp(49)
    dp(59) = dq(2)*p(38) + q(2)*dp(38)
    dp(60) = dq(2)*p(39) + q(2)*dp(39)
    dp(61) = dq(2)*p(40) + q(2)*dp(40)
    dp(62) = dq(2)*p(41) + q(2)*dp(41)
    dp(63) = dq(2)*p(42) + q(2)*dp(42)
    dp(64) = dq(2)*p(43) + q(2)*dp(43)
    dp(65) = dq(2)*p(44) + q(2)*dp(44)
    dp(66) = dq(2)*p(45) + q(2)*dp(45)
    dp(67) = dq(2)*p(46) + q(2)*dp(46)
    dp(68) = dp(1)*p(30) + p(1)*dp(30)
    dp(69) = dp(1)*p(31) + p(1)*dp(31)
    dp(70) = dp(1)*p(32) + p(1)*dp(32)
    dp(71) = dp(1)*q(7) + p(1)*dq(7)
    dp(72) = dq(2)*q(8) + q(2)*dq(8)
    dp(73) = dq(2)*p(50) + q(2)*dp(50)
    dp(74) = dq(2)*p(51) + q(2)*dp(51)
    dp(75) = dq(2)*p(52) + q(2)*dp(52)
    dp(76) = dq(2)*p(53) + q(2)*dp(53)
    dp(77) = dq(2)*p(54) + q(2)*dp(54)
    dp(78) = dq(2)*p(55) + q(2)*dp(55)

    return
  end subroutine devpoly21

  function drdx21(flag,xindex,xyz,r)
    integer i,j,flag,xindex,xyzind,matom,m
    real::xyz(3,3),r(3,3)
    real::drdx21

    if (xindex.eq.1) then 
       i = 1
       j = 2
    elseif (xindex.eq.2) then 
       i = 1
       j = 3
    elseif (xindex.eq.3) then 
       i = 2
       j = 3
    endif

    m=flag
    matom=INT((dble(m)-0.00001d0)/3.d0)+1
    xyzind=MOD(m-1,3)+1

    drdx21 = 0.d0
    if (matom.eq.i.or.matom.eq.j) then
       drdx21=(xyz(i,xyzind)-xyz(j,xyzind))/r(i,j)
       if (matom.eq.j) then
          drdx21 = -drdx21
       endif
    endif

    return
  end function

  subroutine deriv_rev21(c,m,q,p,xyz,a,r,xxp)
    real::c(1:78),m(1:4),p(1:78)
    real::xyz(3,3),r(3,3),a
    !::::::::::::::::::::
    real::pp(1:78),mp(1:4),xxp(1:9)
    real::qp(1:8),q(1:8)

    qp(:)=0.d0
    pp(:)=0.d0
    mp(:)=0.d0
    xxp(:)=0.d0
    pp(78)=c(78)
    pp(77)=c(77)
    pp(76)=c(76)
    pp(75)=c(75)
    pp(74)=c(74)
    pp(73)=c(73)
    pp(72)=c(72)
    pp(71)=c(71)
    pp(70)=c(70)
    pp(69)=c(69)
    pp(68)=c(68)
    pp(67)=c(67)
    pp(66)=c(66)
    pp(65)=c(65)
    pp(64)=c(64)
    pp(63)=c(63)
    pp(62)=c(62)
    pp(61)=c(61)
    pp(60)=c(60)
    pp(59)=c(59)
    pp(58)=c(58)
    pp(57)=c(57)
    pp(56)=c(56)
    qp(8)=pp(72)*q(2)
    pp(55)=c(55)+pp(78)*q(2)
    pp(54)=c(54)+pp(77)*q(2)
    pp(53)=c(53)+pp(76)*q(2)
    pp(52)=c(52)+pp(75)*q(2)
    pp(51)=c(51)+pp(74)*q(2)
    pp(50)=c(50)+pp(73)*q(2)
    pp(49)=c(49)+pp(58)*q(2)+qp(8)*(-1.d0)
    pp(48)=c(48)+pp(57)*q(2)
    pp(47)=c(47)+pp(56)*q(2)
    pp(46)=c(46)+pp(67)*q(2)
    pp(45)=c(45)+pp(66)*q(2)
    pp(44)=c(44)+pp(65)*q(2)
    pp(43)=c(43)+pp(64)*q(2)
    pp(42)=c(42)+pp(63)*q(2)
    pp(41)=c(41)+pp(62)*q(2)
    pp(40)=c(40)+pp(61)*q(2)
    pp(39)=c(39)+pp(60)*q(2)
    pp(38)=c(38)+pp(59)*q(2)
    qp(7)=pp(71)*p(1)+qp(8)*q(1)+pp(50)*q(2)
    pp(37)=c(37)+pp(55)*q(2)
    pp(36)=c(36)+pp(54)*q(2)
    pp(35)=c(35)+pp(53)*q(2)
    pp(34)=c(34)+pp(52)*q(2)
    pp(33)=c(33)+pp(51)*q(2)
    pp(32)=c(32)+pp(70)*p(1)+pp(40)*q(2)+qp(7)*(-1.d0)
    pp(31)=c(31)+pp(69)*p(1)+pp(39)*q(2)
    pp(30)=c(30)+pp(68)*p(1)+pp(38)*q(2)
    pp(29)=c(29)+pp(46)*q(2)
    pp(28)=c(28)+pp(45)*q(2)
    pp(27)=c(27)+pp(44)*q(2)
    pp(26)=c(26)+pp(43)*q(2)
    pp(25)=c(25)+pp(42)*q(2)
    pp(24)=c(24)+pp(41)*q(2)
    qp(6)=pp(49)*p(1)+qp(7)*q(1)+pp(33)*q(2)
    pp(23)=c(23)+pp(37)*q(2)
    pp(22)=c(22)+pp(36)*q(2)
    pp(21)=c(21)+pp(35)*q(2)
    pp(20)=c(20)+pp(34)*q(2)
    pp(19)=c(19)+pp(48)*p(1)+pp(25)*q(2)+qp(6)*(-1.d0)
    pp(18)=c(18)+pp(47)*p(1)+pp(24)*q(2)
    pp(17)=c(17)+pp(29)*q(2)
    pp(16)=c(16)+pp(28)*q(2)
    pp(15)=c(15)+pp(27)*q(2)
    pp(14)=c(14)+pp(26)*q(2)
    qp(5)=pp(32)*p(1)+qp(6)*q(1)+pp(20)*q(2)
    pp(13)=c(13)+pp(23)*q(2)
    pp(12)=c(12)+pp(22)*q(2)
    pp(11)=c(11)+pp(21)*q(2)
    pp(10)=c(10)+pp(31)*p(1)+pp(15)*q(2)+qp(5)*(-1.d0)
    pp(9)=c(9)+pp(30)*p(1)+pp(14)*q(2)
    pp(8)=c(8)+pp(17)*q(2)
    pp(7)=c(7)+pp(16)*q(2)
    qp(4)=pp(19)*p(1)+qp(5)*q(1)+pp(11)*q(2)
    pp(6)=c(6)+pp(13)*q(2)
    pp(5)=c(5)+pp(12)*q(2)
    pp(4)=c(4)+pp(18)*p(1)+pp(7)*q(2)+qp(4)*(-1.d0)
    pp(3)=c(3)+pp(8)*q(2)
    qp(3)=pp(10)*p(1)+qp(4)*q(1)+pp(5)*q(2)
    pp(2)=c(2)+pp(6)*q(2)
    pp(1)=c(1)+pp(71)*q(7)+pp(70)*p(32)+pp(69)*p(31)+pp(68)*p(30)+pp(49)*q(6)+&
     pp(48)*p(19)+pp(47)*p(18)+pp(32)*q(5)+pp(31)*p(10)+pp(30)*p(9)+pp(19)*q(4)+&
    pp(18)*p(4)+pp(10)*q(3)+pp(9)*2*p(1)+pp(4)*q(1)+pp(3)*q(2)+qp(3)*(-2.d0)
    qp(2)=pp(78)*p(55)+pp(77)*p(54)+pp(76)*p(53)+pp(75)*p(52)+pp(74)*p(51)+&
     pp(73)*p(50)+pp(72)*q(8)+pp(67)*p(46)+pp(66)*p(45)+pp(65)*p(44)+pp(64)*p(43)+&
     pp(63)*p(42)+pp(62)*p(41)+pp(61)*p(40)+pp(60)*p(39)+pp(59)*p(38)+pp(58)*p(49)+&
     pp(57)*p(48)+pp(56)*p(47)+pp(55)*p(37)+pp(54)*p(36)+pp(53)*p(35)+pp(52)*p(34)+&
     pp(51)*p(33)+pp(50)*q(7)+pp(46)*p(29)+pp(45)*p(28)+pp(44)*p(27)+pp(43)*p(26)+&
     pp(42)*p(25)+pp(41)*p(24)+pp(40)*p(32)+pp(39)*p(31)+pp(38)*p(30)+pp(37)*p(23)+&
     pp(36)*p(22)+pp(35)*p(21)+pp(34)*p(20)+pp(33)*q(6)+pp(29)*p(17)+pp(28)*p(16)+&
     pp(27)*p(15)+pp(26)*p(14)+pp(25)*p(19)+pp(24)*p(18)+pp(23)*p(13)+pp(22)*p(12)+&
     pp(21)*p(11)+pp(20)*q(5)+pp(17)*p(8)+pp(16)*p(7)+pp(15)*p(10)+pp(14)*p(9)+&
     pp(13)*p(6)+pp(12)*p(5)+pp(11)*q(4)+pp(8)*p(3)+pp(7)*p(4)+pp(6)*p(2)+pp(5)*q(3)+&
    pp(3)*p(1)+pp(2)*q(1)
    qp(1)=qp(8)*q(7)+qp(7)*q(6)+qp(6)*q(5)+qp(5)*q(4)+qp(4)*q(3)+pp(4)*p(1)+&
    qp(3)*2*q(1)+pp(2)*q(2)
    mp(4)=pp(1)
    mp(3)=qp(2)
    mp(2)=qp(1)+mp(4)*m(1)
    mp(1)=qp(1)+mp(4)*m(2)
    xxp(1)=mp(3)*(-m(3)/a)*drdx21(1,1,xyz,r)+mp(2)*(-m(2)/a)*drdx21(1,2,xyz,r)
    xxp(2)=mp(3)*(-m(3)/a)*drdx21(2,1,xyz,r)+mp(2)*(-m(2)/a)*drdx21(2,2,xyz,r)
    xxp(3)=mp(3)*(-m(3)/a)*drdx21(3,1,xyz,r)+mp(2)*(-m(2)/a)*drdx21(3,2,xyz,r)
    xxp(4)=mp(3)*(-m(3)/a)*drdx21(4,1,xyz,r)+mp(1)*(-m(1)/a)*drdx21(4,3,xyz,r)
    xxp(5)=mp(3)*(-m(3)/a)*drdx21(5,1,xyz,r)+mp(1)*(-m(1)/a)*drdx21(5,3,xyz,r)
    xxp(6)=mp(3)*(-m(3)/a)*drdx21(6,1,xyz,r)+mp(1)*(-m(1)/a)*drdx21(6,3,xyz,r)
    xxp(7)=mp(2)*(-m(2)/a)*drdx21(7,2,xyz,r)+mp(1)*(-m(1)/a)*drdx21(7,3,xyz,r)
    xxp(8)=mp(2)*(-m(2)/a)*drdx21(8,2,xyz,r)+mp(1)*(-m(1)/a)*drdx21(8,3,xyz,r)
    xxp(9)=mp(2)*(-m(2)/a)*drdx21(9,2,xyz,r)+mp(1)*(-m(1)/a)*drdx21(9,3,xyz,r)

    return
  end subroutine deriv_rev21

end module bemsa21
