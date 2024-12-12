module bemsa31
  implicit none

contains
  function emsav31(x,c) result(v)
    real,dimension(1:6)::x
    real,dimension(1:115)::p
    real,dimension(1:115)::c
    real::v
    
    call bemsav31(x,p)
    v = dot_product(p,c)
    
    return
  end function emsav31
  
  subroutine bemsav31(x,p)
    real,dimension(1:6)::x
    real,dimension(1:32)::m
    real,dimension(1:115)::p
    real,dimension(1:25)::q
    ! ::::::::::::::::::::

    call evmono31(x,m)
    call evpoly31(m,q,p)
    
    return
  end subroutine bemsav31
  
  subroutine evmono31(x,m)
    real,dimension(1:6)::x
    real,dimension(1:32)::m
    !::::::::::::::::::::
  
    m(1) = x(6)
    m(2) = x(5)
    m(3) = x(3)
    m(4) = x(4)
    m(5) = x(2)
    m(6) = x(1)
    m(7) = m(1)*m(2)
    m(8) = m(1)*m(3)
    m(9) = m(2)*m(3)
    m(10) = m(3)*m(4)
    m(11) = m(2)*m(5)
    m(12) = m(1)*m(6)
    m(13) = m(4)*m(5)
    m(14) = m(4)*m(6)
    m(15) = m(5)*m(6)
    m(16) = m(1)*m(9)
    m(17) = m(1)*m(10)
    m(18) = m(2)*m(10)
    m(19) = m(1)*m(11)
    m(20) = m(3)*m(11)
    m(21) = m(2)*m(12)
    m(22) = m(3)*m(12)
    m(23) = m(2)*m(13)
    m(24) = m(3)*m(13)
    m(25) = m(1)*m(14)
    m(26) = m(3)*m(14)
    m(27) = m(1)*m(15)
    m(28) = m(2)*m(15)
    m(29) = m(4)*m(15)
    m(30) = m(2)*m(24)
    m(31) = m(1)*m(26)
    m(32) = m(1)*m(28)

    return
  end subroutine evmono31
  
  subroutine evpoly31(m,q,p)
    real,dimension(1:32)::m
    real,dimension(1:115)::p
    real,dimension(1:25)::q
    !::::::::::::::::::::

    q(1) = m(1) + m(2) + m(3)
    q(2) = m(4) + m(5) + m(6)
    q(3) = m(7) + m(8) + m(9)
    q(4) = m(10) + m(11) + m(12)
    q(5) = q(1)*q(2) - q(4)
    q(6) = m(13) + m(14) + m(15)
    q(7) = q(1)*q(1) - q(3) - q(3)
    q(8) = q(2)*q(2) - q(6) - q(6)
    p(1) = m(16)
    p(2) = m(17) + m(18) + m(19) + m(20) + m(21) + m(22)
    q(9) = q(2)*q(3) - p(2)
    p(3) = m(23) + m(24) + m(25) + m(26) + m(27) + m(28)
    q(10) = m(29)
    p(4) = q(1)*q(6) - p(3)
    q(11) = q(1)*q(3) - p(1) - p(1) - p(1)
    q(12) = q(1)*q(4) - p(2)
    q(13) = q(2)*q(7) - q(12)
    q(14) = q(2)*q(4) - p(3)
    q(15) = q(1)*q(8) - q(14)
    q(16) = q(2)*q(6) - q(10) - q(10) - q(10)
    q(17) = q(1)*q(7) - q(11)
    q(18) = q(2)*q(8) - q(16)
    p(5) = p(1)*q(2)
    p(6) = m(30) + m(31) + m(32)
    p(7) = q(3)*q(6) - p(6)
    p(8) = q(10)*q(1)
    p(9) = p(1)*q(1)
    p(10) = q(3)*q(4) - p(5)
    p(11) = q(1)*p(2) - p(5) - p(10) - p(5)
    p(12) = q(1)*p(3) - p(7) - p(6) - p(6)
    p(13) = q(1)*p(4) - p(7)
    p(14) = q(4)*q(5) - p(7) - p(12)
    q(19) = q(2)*q(9) - p(7)
    p(15) = q(4)*q(6) - p(8)
    p(16) = q(2)*p(3) - p(8) - p(15) - p(8)
    q(20) = q(10)*q(2)
    p(17) = q(2)*p(4) - p(8)
    q(21) = q(3)*q(3) - p(9) - p(9)
    q(22) = q(1)*q(12) - p(10)
    q(23) = q(1)*q(14) - p(14)
    q(24) = q(6)*q(6) - q(20) - q(20)
    q(25) = q(2)*q(14) - p(15)
    p(18) = p(1)*q(6)
    p(19) = q(10)*q(3)
    p(20) = p(1)*q(4)
    p(21) = p(1)*q(5)
    p(22) = q(1)*p(6) - p(18)
    p(23) = q(3)*p(3) - p(18) - p(22) - p(18)
    p(24) = q(3)*p(4) - p(18)
    p(25) = q(10)*q(7)
    p(26) = p(1)*q(8)
    p(27) = q(2)*p(6) - p(19)
    p(28) = q(10)*q(4)
    p(29) = q(4)*p(4) - p(25)
    p(30) = q(6)*q(9) - p(19)
    p(31) = q(10)*q(5)
    p(32) = p(1)*q(3)
    p(33) = p(1)*q(7)
    p(34) = q(3)*q(12) - p(20)
    p(35) = q(4)*q(11) - p(21) - p(34)
    p(36) = q(1)*p(11) - p(21) - p(35)
    p(37) = q(1)*p(12) - p(23) - p(22)
    p(38) = q(1)*p(13) - p(24)
    p(39) = q(3)*q(14) - p(26)
    p(40) = q(2)*p(11) - p(24) - p(22)
    p(41) = q(6)*q(12) - p(25)
    p(42) = q(1)*p(16) - p(30) - p(27)
    p(43) = q(2)*p(13) - p(25)
    p(44) = q(6)*p(3) - p(31) - p(28) - p(28)
    p(45) = q(1)*q(24) - p(44)
    p(46) = q(2)*p(14) - p(29) - p(27)
    p(47) = q(6)*q(14) - p(28)
    p(48) = q(2)*p(16) - p(31) - p(44)
    p(49) = q(8)*p(4) - p(28)
    p(50) = p(1)*q(10)
    p(51) = p(1)*p(3)
    p(52) = p(1)*p(4)
    p(53) = q(10)*q(11)
    p(54) = p(1)*q(16)
    p(55) = q(10)*p(2)
    p(56) = q(10)*q(9)
    p(57) = p(1)*q(12)
    p(58) = p(1)*p(2)
    p(59) = p(1)*q(9)
    p(60) = p(1)*q(13)
    p(61) = q(3)*p(6) - p(51)
    p(62) = q(7)*p(6) - p(52)
    p(63) = q(3)*p(12) - p(51) - p(62)
    p(64) = q(6)*q(21) - p(61)
    p(65) = q(3)*p(13) - p(52)
    p(66) = q(10)*q(17)
    p(67) = p(1)*q(14)
    p(68) = p(1)*q(15)
    p(69) = q(4)*p(6) - p(50) - p(50) - p(50)
    p(70) = q(1)*p(27) - p(54) - p(69)
    p(71) = q(10)*q(12)
    p(72) = p(4)*q(12) - p(66)
    p(73) = q(2)*p(23) - p(53) - p(72)
    p(74) = q(4)*p(13) - p(66)
    p(75) = q(1)*p(30) - p(54) - p(73)
    p(76) = q(10)*q(13)
    p(77) = q(6)*p(6) - p(55)
    p(78) = q(10)*p(3)
    p(79) = q(3)*q(24) - p(77)
    p(80) = q(10)*p(4)
    p(81) = p(1)*q(18)
    p(82) = q(8)*p(6) - p(56)
    p(83) = q(10)*q(14)
    p(84) = p(4)*q(14) - p(71)
    p(85) = q(6)*q(19) - p(56)
    p(86) = q(10)*q(15)
    p(87) = p(1)*p(1)
    p(88) = p(1)*q(11)
    p(89) = p(1)*q(17)
    p(90) = q(3)*q(22) - p(57)
    p(91) = q(4)*q(21) - p(59)
    p(92) = q(1)*p(35) - p(58) - p(91)
    p(93) = q(1)*p(36) - p(60) - p(92)
    p(94) = q(1)*p(37) - p(63) - p(62)
    p(95) = q(1)*p(38) - p(65)
    p(96) = q(3)*q(23) - p(67)
    p(97) = q(11)*q(14) - p(68) - p(96)
    p(98) = q(1)*p(40) - p(68) - p(97)
    p(99) = q(6)*q(22) - p(66)
    p(100) = q(1)*p(42) - p(73) - p(70)
    p(101) = q(2)*p(38) - p(66)
    p(102) = q(1)*p(44) - p(79) - p(77) - p(77)
    p(103) = q(1)*p(45) - p(79)
    p(104) = q(3)*q(25) - p(81)
    p(105) = q(2)*p(40) - p(74) - p(70)
    p(106) = q(6)*q(23) - p(71)
    p(107) = q(1)*p(48) - p(85) - p(82)
    p(108) = q(8)*p(13) - p(71)
    p(109) = q(4)*q(24) - p(80)
    p(110) = q(2)*p(44) - p(78) - p(109)
    p(111) = q(2)*p(45) - p(80)
    p(112) = q(2)*p(46) - p(84) - p(82)
    p(113) = q(6)*q(25) - p(83)
    p(114) = q(2)*p(48) - p(86) - p(110)
    p(115) = p(4)*q(18) - p(83)

    return
  end subroutine evpoly31

  subroutine devmono31(dm,m,flag,xyz,a,r)
    real,dimension(1:32)::dm
    real,dimension(1:32)::m
    real::xyz(4,3),a,r(4,4)
    integer::flag
    !::::::::::::::::::::

    dm(1) =  - m(1)/a*drdx31(flag,6,xyz,r)
    dm(2) =  - m(2)/a*drdx31(flag,5,xyz,r)
    dm(3) =  - m(3)/a*drdx31(flag,3,xyz,r)
    dm(4) =  - m(4)/a*drdx31(flag,4,xyz,r)
    dm(5) =  - m(5)/a*drdx31(flag,2,xyz,r)
    dm(6) =  - m(6)/a*drdx31(flag,1,xyz,r)
    dm(7) = dm(1)*m(2) + m(1)*dm(2)
    dm(8) = dm(1)*m(3) + m(1)*dm(3)
    dm(9) = dm(2)*m(3) + m(2)*dm(3)
    dm(10) = dm(3)*m(4) + m(3)*dm(4)
    dm(11) = dm(2)*m(5) + m(2)*dm(5)
    dm(12) = dm(1)*m(6) + m(1)*dm(6)
    dm(13) = dm(4)*m(5) + m(4)*dm(5)
    dm(14) = dm(4)*m(6) + m(4)*dm(6)
    dm(15) = dm(5)*m(6) + m(5)*dm(6)
    dm(16) = dm(1)*m(9) + m(1)*dm(9)
    dm(17) = dm(1)*m(10) + m(1)*dm(10)
    dm(18) = dm(2)*m(10) + m(2)*dm(10)
    dm(19) = dm(1)*m(11) + m(1)*dm(11)
    dm(20) = dm(3)*m(11) + m(3)*dm(11)
    dm(21) = dm(2)*m(12) + m(2)*dm(12)
    dm(22) = dm(3)*m(12) + m(3)*dm(12)
    dm(23) = dm(2)*m(13) + m(2)*dm(13)
    dm(24) = dm(3)*m(13) + m(3)*dm(13)
    dm(25) = dm(1)*m(14) + m(1)*dm(14)
    dm(26) = dm(3)*m(14) + m(3)*dm(14)
    dm(27) = dm(1)*m(15) + m(1)*dm(15)
    dm(28) = dm(2)*m(15) + m(2)*dm(15)
    dm(29) = dm(4)*m(15) + m(4)*dm(15)
    dm(30) = dm(2)*m(24) + m(2)*dm(24)
    dm(31) = dm(1)*m(26) + m(1)*dm(26)
    dm(32) = dm(1)*m(28) + m(1)*dm(28)

    return
  end subroutine devmono31

  subroutine devpoly31(dm,q,p,dp)
    real,dimension(1:32)::dm
    real,dimension(1:115)::p,dp
    real,dimension(1:25)::q,dq
    !::::::::::::::::::::

    dq(1) = dm(1) + dm(2) + dm(3)
    dq(2) = dm(4) + dm(5) + dm(6)
    dq(3) = dm(7) + dm(8) + dm(9)
    dq(4) = dm(10) + dm(11) + dm(12)
    dq(5) = dq(1)*q(2) + q(1)*dq(2) - dq(4)
    dq(6) = dm(13) + dm(14) + dm(15)
    dq(7) = dq(1)*q(1) + q(1)*dq(1) - dq(3) - dq(3)
    dq(8) = dq(2)*q(2) + q(2)*dq(2) - dq(6) - dq(6)
    dp(1) = dm(16)
    dp(2) = dm(17) + dm(18) + dm(19) + dm(20) + dm(21) + dm(22)
    dq(9) = dq(2)*q(3) + q(2)*dq(3) - dp(2)
    dp(3) = dm(23) + dm(24) + dm(25) + dm(26) + dm(27) + dm(28)
    dq(10) = dm(29)
    dp(4) = dq(1)*q(6) + q(1)*dq(6) - dp(3)
    dq(11) = dq(1)*q(3) + q(1)*dq(3) - dp(1) - dp(1) - dp(1)
    dq(12) = dq(1)*q(4) + q(1)*dq(4) - dp(2)
    dq(13) = dq(2)*q(7) + q(2)*dq(7) - dq(12)
    dq(14) = dq(2)*q(4) + q(2)*dq(4) - dp(3)
    dq(15) = dq(1)*q(8) + q(1)*dq(8) - dq(14)
    dq(16) = dq(2)*q(6) + q(2)*dq(6) - dq(10) - dq(10) - dq(10)
    dq(17) = dq(1)*q(7) + q(1)*dq(7) - dq(11)
    dq(18) = dq(2)*q(8) + q(2)*dq(8) - dq(16)
    dp(5) = dp(1)*q(2) + p(1)*dq(2)
    dp(6) = dm(30) + dm(31) + dm(32)
    dp(7) = dq(3)*q(6) + q(3)*dq(6) - dp(6)
    dp(8) = dq(10)*q(1) + q(10)*dq(1)
    dp(9) = dp(1)*q(1) + p(1)*dq(1)
    dp(10) = dq(3)*q(4) + q(3)*dq(4) - dp(5)
    dp(11) = dq(1)*p(2) + q(1)*dp(2) - dp(5) - dp(10) - dp(5)
    dp(12) = dq(1)*p(3) + q(1)*dp(3) - dp(7) - dp(6) - dp(6)
    dp(13) = dq(1)*p(4) + q(1)*dp(4) - dp(7)
    dp(14) = dq(4)*q(5) + q(4)*dq(5) - dp(7) - dp(12)
    dq(19) = dq(2)*q(9) + q(2)*dq(9) - dp(7)
    dp(15) = dq(4)*q(6) + q(4)*dq(6) - dp(8)
    dp(16) = dq(2)*p(3) + q(2)*dp(3) - dp(8) - dp(15) - dp(8)
    dq(20) = dq(10)*q(2) + q(10)*dq(2)
    dp(17) = dq(2)*p(4) + q(2)*dp(4) - dp(8)
    dq(21) = dq(3)*q(3) + q(3)*dq(3) - dp(9) - dp(9)
    dq(22) = dq(1)*q(12) + q(1)*dq(12) - dp(10)
    dq(23) = dq(1)*q(14) + q(1)*dq(14) - dp(14)
    dq(24) = dq(6)*q(6) + q(6)*dq(6) - dq(20) - dq(20)
    dq(25) = dq(2)*q(14) + q(2)*dq(14) - dp(15)
    dp(18) = dp(1)*q(6) + p(1)*dq(6)
    dp(19) = dq(10)*q(3) + q(10)*dq(3)
    dp(20) = dp(1)*q(4) + p(1)*dq(4)
    dp(21) = dp(1)*q(5) + p(1)*dq(5)
    dp(22) = dq(1)*p(6) + q(1)*dp(6) - dp(18)
    dp(23) = dq(3)*p(3) + q(3)*dp(3) - dp(18) - dp(22) - dp(18)
    dp(24) = dq(3)*p(4) + q(3)*dp(4) - dp(18)
    dp(25) = dq(10)*q(7) + q(10)*dq(7)
    dp(26) = dp(1)*q(8) + p(1)*dq(8)
    dp(27) = dq(2)*p(6) + q(2)*dp(6) - dp(19)
    dp(28) = dq(10)*q(4) + q(10)*dq(4)
    dp(29) = dq(4)*p(4) + q(4)*dp(4) - dp(25)
    dp(30) = dq(6)*q(9) + q(6)*dq(9) - dp(19)
    dp(31) = dq(10)*q(5) + q(10)*dq(5)
    dp(32) = dp(1)*q(3) + p(1)*dq(3)
    dp(33) = dp(1)*q(7) + p(1)*dq(7)
    dp(34) = dq(3)*q(12) + q(3)*dq(12) - dp(20)
    dp(35) = dq(4)*q(11) + q(4)*dq(11) - dp(21) - dp(34)
    dp(36) = dq(1)*p(11) + q(1)*dp(11) - dp(21) - dp(35)
    dp(37) = dq(1)*p(12) + q(1)*dp(12) - dp(23) - dp(22)
    dp(38) = dq(1)*p(13) + q(1)*dp(13) - dp(24)
    dp(39) = dq(3)*q(14) + q(3)*dq(14) - dp(26)
    dp(40) = dq(2)*p(11) + q(2)*dp(11) - dp(24) - dp(22)
    dp(41) = dq(6)*q(12) + q(6)*dq(12) - dp(25)
    dp(42) = dq(1)*p(16) + q(1)*dp(16) - dp(30) - dp(27)
    dp(43) = dq(2)*p(13) + q(2)*dp(13) - dp(25)
    dp(44) = dq(6)*p(3) + q(6)*dp(3) - dp(31) - dp(28) - dp(28)
    dp(45) = dq(1)*q(24) + q(1)*dq(24) - dp(44)
    dp(46) = dq(2)*p(14) + q(2)*dp(14) - dp(29) - dp(27)
    dp(47) = dq(6)*q(14) + q(6)*dq(14) - dp(28)
    dp(48) = dq(2)*p(16) + q(2)*dp(16) - dp(31) - dp(44)
    dp(49) = dq(8)*p(4) + q(8)*dp(4) - dp(28)
    dp(50) = dp(1)*q(10) + p(1)*dq(10)
    dp(51) = dp(1)*p(3) + p(1)*dp(3)
    dp(52) = dp(1)*p(4) + p(1)*dp(4)
    dp(53) = dq(10)*q(11) + q(10)*dq(11)
    dp(54) = dp(1)*q(16) + p(1)*dq(16)
    dp(55) = dq(10)*p(2) + q(10)*dp(2)
    dp(56) = dq(10)*q(9) + q(10)*dq(9)
    dp(57) = dp(1)*q(12) + p(1)*dq(12)
    dp(58) = dp(1)*p(2) + p(1)*dp(2)
    dp(59) = dp(1)*q(9) + p(1)*dq(9)
    dp(60) = dp(1)*q(13) + p(1)*dq(13)
    dp(61) = dq(3)*p(6) + q(3)*dp(6) - dp(51)
    dp(62) = dq(7)*p(6) + q(7)*dp(6) - dp(52)
    dp(63) = dq(3)*p(12) + q(3)*dp(12) - dp(51) - dp(62)
    dp(64) = dq(6)*q(21) + q(6)*dq(21) - dp(61)
    dp(65) = dq(3)*p(13) + q(3)*dp(13) - dp(52)
    dp(66) = dq(10)*q(17) + q(10)*dq(17)
    dp(67) = dp(1)*q(14) + p(1)*dq(14)
    dp(68) = dp(1)*q(15) + p(1)*dq(15)
    dp(69) = dq(4)*p(6) + q(4)*dp(6) - dp(50) - dp(50) - dp(50)
    dp(70) = dq(1)*p(27) + q(1)*dp(27) - dp(54) - dp(69)
    dp(71) = dq(10)*q(12) + q(10)*dq(12)
    dp(72) = dp(4)*q(12) + p(4)*dq(12) - dp(66)
    dp(73) = dq(2)*p(23) + q(2)*dp(23) - dp(53) - dp(72)
    dp(74) = dq(4)*p(13) + q(4)*dp(13) - dp(66)
    dp(75) = dq(1)*p(30) + q(1)*dp(30) - dp(54) - dp(73)
    dp(76) = dq(10)*q(13) + q(10)*dq(13)
    dp(77) = dq(6)*p(6) + q(6)*dp(6) - dp(55)
    dp(78) = dq(10)*p(3) + q(10)*dp(3)
    dp(79) = dq(3)*q(24) + q(3)*dq(24) - dp(77)
    dp(80) = dq(10)*p(4) + q(10)*dp(4)
    dp(81) = dp(1)*q(18) + p(1)*dq(18)
    dp(82) = dq(8)*p(6) + q(8)*dp(6) - dp(56)
    dp(83) = dq(10)*q(14) + q(10)*dq(14)
    dp(84) = dp(4)*q(14) + p(4)*dq(14) - dp(71)
    dp(85) = dq(6)*q(19) + q(6)*dq(19) - dp(56)
    dp(86) = dq(10)*q(15) + q(10)*dq(15)
    dp(87) = dp(1)*p(1) + p(1)*dp(1)
    dp(88) = dp(1)*q(11) + p(1)*dq(11)
    dp(89) = dp(1)*q(17) + p(1)*dq(17)
    dp(90) = dq(3)*q(22) + q(3)*dq(22) - dp(57)
    dp(91) = dq(4)*q(21) + q(4)*dq(21) - dp(59)
    dp(92) = dq(1)*p(35) + q(1)*dp(35) - dp(58) - dp(91)
    dp(93) = dq(1)*p(36) + q(1)*dp(36) - dp(60) - dp(92)
    dp(94) = dq(1)*p(37) + q(1)*dp(37) - dp(63) - dp(62)
    dp(95) = dq(1)*p(38) + q(1)*dp(38) - dp(65)
    dp(96) = dq(3)*q(23) + q(3)*dq(23) - dp(67)
    dp(97) = dq(11)*q(14) + q(11)*dq(14) - dp(68) - dp(96)
    dp(98) = dq(1)*p(40) + q(1)*dp(40) - dp(68) - dp(97)
    dp(99) = dq(6)*q(22) + q(6)*dq(22) - dp(66)
    dp(100) = dq(1)*p(42) + q(1)*dp(42) - dp(73) - dp(70)
    dp(101) = dq(2)*p(38) + q(2)*dp(38) - dp(66)
    dp(102) = dq(1)*p(44) + q(1)*dp(44) - dp(79) - dp(77) - dp(77)
    dp(103) = dq(1)*p(45) + q(1)*dp(45) - dp(79)
    dp(104) = dq(3)*q(25) + q(3)*dq(25) - dp(81)
    dp(105) = dq(2)*p(40) + q(2)*dp(40) - dp(74) - dp(70)
    dp(106) = dq(6)*q(23) + q(6)*dq(23) - dp(71)
    dp(107) = dq(1)*p(48) + q(1)*dp(48) - dp(85) - dp(82)
    dp(108) = dq(8)*p(13) + q(8)*dp(13) - dp(71)
    dp(109) = dq(4)*q(24) + q(4)*dq(24) - dp(80)
    dp(110) = dq(2)*p(44) + q(2)*dp(44) - dp(78) - dp(109)
    dp(111) = dq(2)*p(45) + q(2)*dp(45) - dp(80)
    dp(112) = dq(2)*p(46) + q(2)*dp(46) - dp(84) - dp(82)
    dp(113) = dq(6)*q(25) + q(6)*dq(25) - dp(83)
    dp(114) = dq(2)*p(48) + q(2)*dp(48) - dp(86) - dp(110)
    dp(115) = dp(4)*q(18) + p(4)*dq(18) - dp(83)

    return
  end subroutine devpoly31

  function drdx31(flag,xindex,xyz,r)
    integer i,j,flag,xindex,xyzind,matom,m
    real::xyz(4,3),r(4,4),drdx31

    if (xindex.eq.1) then 
       i = 1
       j = 2
    elseif (xindex.eq.2) then 
       i = 1
       j = 3
    elseif (xindex.eq.3) then 
       i = 1
       j = 4
    elseif (xindex.eq.4) then 
       i = 2
       j = 3
    elseif (xindex.eq.5) then 
       i = 2
       j = 4
    elseif (xindex.eq.6) then 
       i = 3
       j = 4
    endif

    m=flag
    matom=INT((dble(m)-0.00001d0)/3.d0)+1
    xyzind=MOD(m-1,3)+1

    drdx31 = 0.d0
    if (matom.eq.i.or.matom.eq.j) then
       drdx31=(xyz(i,xyzind)-xyz(j,xyzind))/r(i,j)
       if (matom.eq.j) then
          drdx31 = -drdx31
       endif
    endif

    return
  end function

  subroutine deriv_rev31(c,m,q,p,xyz,a,r,xxp)
    real::c(1:115),m(1:32),p(1:115)
    real::xyz(4,3),r(4,4),a
    !::::::::::::::::::::
    real::pp(1:115),mp(1:32),xxp(1:12)
    real::qp(1:25),q(1:25)

    qp(:)=0.d0
    pp(:)=0.d0
    mp(:)=0.d0
    xxp(:)=0.d0
    pp(115)=c(115)
    pp(114)=c(114)
    pp(113)=c(113)
    pp(112)=c(112)
    pp(111)=c(111)
    pp(110)=c(110)+pp(114)*(-1.d0)
    pp(109)=c(109)+pp(110)*(-1.d0)
    pp(108)=c(108)
    pp(107)=c(107)
    pp(106)=c(106)
    pp(105)=c(105)
    pp(104)=c(104)
    pp(103)=c(103)
    pp(102)=c(102)
    pp(101)=c(101)
    pp(100)=c(100)
    pp(99)=c(99)
    pp(98)=c(98)
    pp(97)=c(97)+pp(98)*(-1.d0)
    pp(96)=c(96)+pp(97)*(-1.d0)
    pp(95)=c(95)
    pp(94)=c(94)
    pp(93)=c(93)
    pp(92)=c(92)+pp(93)*(-1.d0)
    pp(91)=c(91)+pp(92)*(-1.d0)
    pp(90)=c(90)
    pp(89)=c(89)
    pp(88)=c(88)
    pp(87)=c(87)
    pp(86)=c(86)+pp(114)*(-1.d0)
    pp(85)=c(85)+pp(107)*(-1.d0)
    pp(84)=c(84)+pp(112)*(-1.d0)
    pp(83)=c(83)+pp(115)*(-1.d0)+pp(113)*(-1.d0)
    pp(82)=c(82)+pp(112)*(-1.d0)+pp(107)*(-1.d0)
    pp(81)=c(81)+pp(104)*(-1.d0)
    pp(80)=c(80)+pp(111)*(-1.d0)+pp(109)*(-1.d0)
    pp(79)=c(79)+pp(103)*(-1.d0)+pp(102)*(-1.d0)
    pp(78)=c(78)+pp(110)*(-1.d0)
    pp(77)=c(77)+pp(102)*(-2.d0)+pp(79)*(-1.d0)
    pp(76)=c(76)
    pp(75)=c(75)
    pp(74)=c(74)+pp(105)*(-1.d0)
    pp(73)=c(73)+pp(100)*(-1.d0)+pp(75)*(-1.d0)
    pp(72)=c(72)+pp(73)*(-1.d0)
    pp(71)=c(71)+pp(108)*(-1.d0)+pp(106)*(-1.d0)+pp(84)*(-1.d0)
    pp(70)=c(70)+pp(105)*(-1.d0)+pp(100)*(-1.d0)
    pp(69)=c(69)+pp(70)*(-1.d0)
    pp(68)=c(68)+pp(98)*(-1.d0)+pp(97)*(-1.d0)
    pp(67)=c(67)+pp(96)*(-1.d0)
    pp(66)=c(66)+pp(101)*(-1.d0)+pp(99)*(-1.d0)+pp(74)*(-1.d0)+pp(72)*(-1.d0)
    pp(65)=c(65)+pp(95)*(-1.d0)
    pp(64)=c(64)
    pp(63)=c(63)+pp(94)*(-1.d0)
    pp(62)=c(62)+pp(94)*(-1.d0)+pp(63)*(-1.d0)
    pp(61)=c(61)+pp(64)*(-1.d0)
    pp(60)=c(60)+pp(93)*(-1.d0)
    pp(59)=c(59)+pp(91)*(-1.d0)
    pp(58)=c(58)+pp(92)*(-1.d0)
    pp(57)=c(57)+pp(90)*(-1.d0)
    pp(56)=c(56)+pp(85)*(-1.d0)+pp(82)*(-1.d0)
    pp(55)=c(55)+pp(77)*(-1.d0)
    pp(54)=c(54)+pp(75)*(-1.d0)+pp(70)*(-1.d0)
    pp(53)=c(53)+pp(73)*(-1.d0)
    pp(52)=c(52)+pp(65)*(-1.d0)+pp(62)*(-1.d0)
    pp(51)=c(51)+pp(63)*(-1.d0)+pp(61)*(-1.d0)
    pp(50)=c(50)+pp(69)*(-3.d0)
    pp(49)=c(49)
    pp(48)=c(48)+pp(114)*q(2)+pp(107)*q(1)
    pp(47)=c(47)
    pp(46)=c(46)+pp(112)*q(2)
    pp(45)=c(45)+pp(111)*q(2)+pp(103)*q(1)
    pp(44)=c(44)+pp(110)*q(2)+pp(102)*q(1)+pp(48)*(-1.d0)+pp(45)*(-1.d0)
    pp(43)=c(43)
    pp(42)=c(42)+pp(100)*q(1)
    pp(41)=c(41)
    pp(40)=c(40)+pp(105)*q(2)+pp(98)*q(1)
    pp(39)=c(39)
    pp(38)=c(38)+pp(101)*q(2)+pp(95)*q(1)
    pp(37)=c(37)+pp(94)*q(1)
    pp(36)=c(36)+pp(93)*q(1)
    pp(35)=c(35)+pp(92)*q(1)+pp(36)*(-1.d0)
    pp(34)=c(34)+pp(35)*(-1.d0)
    pp(33)=c(33)
    pp(32)=c(32)
    pp(31)=c(31)+pp(48)*(-1.d0)+pp(44)*(-1.d0)
    pp(30)=c(30)+pp(75)*q(1)+pp(42)*(-1.d0)
    pp(29)=c(29)+pp(46)*(-1.d0)
    pp(28)=c(28)+pp(49)*(-1.d0)+pp(47)*(-1.d0)+pp(44)*(-2.d0)
    pp(27)=c(27)+pp(70)*q(1)+pp(46)*(-1.d0)+pp(42)*(-1.d0)
    pp(26)=c(26)+pp(39)*(-1.d0)
    pp(25)=c(25)+pp(43)*(-1.d0)+pp(41)*(-1.d0)+pp(29)*(-1.d0)
    pp(24)=c(24)+pp(40)*(-1.d0)+pp(38)*(-1.d0)
    pp(23)=c(23)+pp(73)*q(2)+pp(37)*(-1.d0)
    pp(22)=c(22)+pp(40)*(-1.d0)+pp(37)*(-1.d0)+pp(23)*(-1.d0)
    pp(21)=c(21)+pp(36)*(-1.d0)+pp(35)*(-1.d0)
    pp(20)=c(20)+pp(34)*(-1.d0)
    pp(19)=c(19)+pp(30)*(-1.d0)+pp(27)*(-1.d0)
    pp(18)=c(18)+pp(24)*(-1.d0)+pp(23)*(-2.d0)+pp(22)*(-1.d0)
    qp(25)=pp(113)*q(6)+pp(104)*q(3)
    qp(24)=pp(109)*q(4)+pp(79)*q(3)+pp(45)*q(1)
    qp(23)=pp(106)*q(6)+pp(96)*q(3)
    qp(22)=pp(99)*q(6)+pp(90)*q(3)
    qp(21)=pp(91)*q(4)+pp(64)*q(6)
    pp(17)=c(17)
    qp(20)=qp(24)*(-2.d0)
    pp(16)=c(16)+pp(48)*q(2)+pp(42)*q(1)
    pp(15)=c(15)+qp(25)*(-1.d0)+pp(16)*(-1.d0)
    qp(19)=pp(85)*q(6)
    pp(14)=c(14)+pp(46)*q(2)+qp(23)*(-1.d0)
    pp(13)=c(13)+pp(108)*q(8)+pp(74)*q(4)+pp(65)*q(3)+pp(43)*q(2)+pp(38)*q(1)
    pp(12)=c(12)+pp(63)*q(3)+pp(37)*q(1)+pp(14)*(-1.d0)
    pp(11)=c(11)+pp(40)*q(2)+pp(36)*q(1)
    pp(10)=c(10)+qp(22)*(-1.d0)+pp(11)*(-1.d0)
    pp(9)=c(9)+qp(21)*(-2.d0)
    pp(8)=c(8)+pp(17)*(-1.d0)+pp(16)*(-2.d0)+pp(15)*(-1.d0)
    pp(7)=c(7)+qp(19)*(-1.d0)+pp(14)*(-1.d0)+pp(13)*(-1.d0)+pp(12)*(-1.d0)
    pp(6)=c(6)+pp(82)*q(8)+pp(77)*q(6)+pp(69)*q(4)+pp(62)*q(7)+pp(61)*q(3)+&
    pp(27)*q(2)+pp(22)*q(1)+pp(12)*(-2.d0)+pp(7)*(-1.d0)
    pp(5)=c(5)+pp(11)*(-2.d0)+pp(10)*(-1.d0)
    qp(18)=pp(115)*p(4)+pp(81)*p(1)
    qp(17)=pp(89)*p(1)+pp(66)*q(10)
    qp(16)=pp(54)*p(1)+qp(18)*(-1.d0)
    qp(15)=pp(86)*q(10)+pp(68)*p(1)
    qp(14)=pp(97)*q(11)+pp(84)*p(4)+pp(83)*q(10)+pp(67)*p(1)+pp(47)*q(6)+&
    pp(39)*q(3)+qp(25)*q(2)+qp(23)*q(1)+qp(15)*(-1.d0)
    qp(13)=pp(76)*q(10)+pp(60)*p(1)
    qp(12)=pp(72)*p(4)+pp(71)*q(10)+pp(57)*p(1)+pp(41)*q(6)+pp(34)*q(3)+qp(22)*q(1)+&
    qp(13)*(-1.d0)
    qp(11)=pp(97)*q(14)+pp(88)*p(1)+pp(53)*q(10)+pp(35)*q(4)+qp(17)*(-1.d0)
    pp(4)=c(4)+pp(115)*q(18)+pp(84)*q(14)+pp(80)*q(10)+pp(72)*q(12)+pp(52)*p(1)+&
    pp(49)*q(8)+pp(29)*q(4)+pp(24)*q(3)+pp(17)*q(2)+pp(13)*q(1)
    qp(10)=pp(86)*q(15)+pp(83)*q(14)+pp(80)*p(4)+pp(78)*p(3)+pp(76)*q(13)+&
     pp(71)*q(12)+pp(66)*q(17)+pp(56)*q(9)+pp(55)*p(2)+pp(53)*q(11)+pp(50)*p(1)+&
    pp(31)*q(5)+pp(28)*q(4)+pp(25)*q(7)+pp(19)*q(3)+qp(20)*q(2)+pp(8)*q(1)+qp(16)*(-3.d0)
    pp(3)=c(3)+pp(78)*q(10)+pp(51)*p(1)+pp(44)*q(6)+pp(23)*q(3)+pp(16)*q(2)+&
    pp(12)*q(1)+qp(14)*(-1.d0)+pp(4)*(-1.d0)
    qp(9)=pp(59)*p(1)+pp(56)*q(10)+pp(30)*q(6)+qp(19)*q(2)
    pp(2)=c(2)+pp(58)*p(1)+pp(55)*q(10)+pp(11)*q(1)+qp(12)*(-1.d0)+qp(9)*(-1.d0)
    pp(1)=c(1)+pp(89)*q(17)+pp(88)*q(11)+pp(87)*2*p(1)+pp(81)*q(18)+pp(68)*q(15)+&
     pp(67)*q(14)+pp(60)*q(13)+pp(59)*q(9)+pp(58)*p(2)+pp(57)*q(12)+pp(54)*q(16)+&
     pp(52)*p(4)+pp(51)*p(3)+pp(50)*q(10)+pp(33)*q(7)+pp(32)*q(3)+pp(26)*q(8)+&
    pp(21)*q(5)+pp(20)*q(4)+pp(18)*q(6)+pp(9)*q(1)+pp(5)*q(2)+qp(11)*(-3.d0)
    qp(8)=pp(108)*p(13)+pp(82)*p(6)+pp(49)*p(4)+pp(26)*p(1)+qp(18)*q(2)+qp(15)*q(1)
    qp(7)=pp(62)*p(6)+pp(33)*p(1)+pp(25)*q(10)+qp(17)*q(1)+qp(13)*q(2)
    qp(6)=pp(113)*q(25)+pp(106)*q(23)+pp(99)*q(22)+pp(85)*q(19)+pp(77)*p(6)+&
     pp(64)*q(21)+pp(47)*q(14)+pp(44)*p(3)+pp(41)*q(12)+pp(30)*q(9)+pp(18)*p(1)+&
    qp(24)*2*q(6)+pp(15)*q(4)+pp(7)*q(3)+qp(16)*q(2)+pp(4)*q(1)+qp(8)*(-2.d0)
    qp(5)=pp(31)*q(10)+pp(21)*p(1)+pp(14)*q(4)
    qp(4)=pp(109)*q(24)+pp(91)*q(21)+pp(74)*p(13)+pp(69)*p(6)+pp(35)*q(11)+&
     pp(29)*p(4)+pp(28)*q(10)+pp(20)*p(1)+pp(15)*q(6)+pp(14)*q(5)+pp(10)*q(3)+&
    qp(14)*q(2)+qp(12)*q(1)+qp(5)*(-1.d0)
    qp(3)=pp(104)*q(25)+pp(96)*q(23)+pp(90)*q(22)+pp(79)*q(24)+pp(65)*p(13)+&
     pp(63)*p(12)+pp(61)*p(6)+pp(39)*q(14)+pp(34)*q(12)+pp(32)*p(1)+pp(24)*p(4)+&
     pp(23)*p(3)+pp(19)*q(10)+qp(21)*2*q(3)+pp(10)*q(4)+pp(7)*q(6)+qp(11)*q(1)+&
    qp(9)*q(2)+qp(7)*(-2.d0)
    qp(2)=pp(114)*p(48)+pp(112)*p(46)+pp(111)*p(45)+pp(110)*p(44)+pp(105)*p(40)+&
     pp(101)*p(38)+pp(73)*p(23)+pp(48)*p(16)+pp(46)*p(14)+pp(43)*p(13)+pp(40)*p(11)+&
     pp(27)*p(6)+qp(25)*q(14)+pp(17)*p(4)+qp(20)*q(10)+pp(16)*p(3)+qp(19)*q(9)+&
     pp(5)*p(1)+qp(18)*q(8)+qp(16)*q(6)+qp(14)*q(4)+qp(13)*q(7)+qp(9)*q(3)+qp(8)*2*q(2)+&
    qp(5)*q(1)
    qp(1)=pp(107)*p(48)+pp(103)*p(45)+pp(102)*p(44)+pp(100)*p(42)+pp(98)*p(40)+&
     pp(95)*p(38)+pp(94)*p(37)+pp(93)*p(36)+pp(92)*p(35)+pp(75)*p(30)+pp(70)*p(27)+&
     pp(45)*q(24)+pp(42)*p(16)+pp(38)*p(13)+pp(37)*p(12)+pp(36)*p(11)+pp(22)*p(6)+&
     qp(23)*q(14)+qp(22)*q(12)+pp(13)*p(4)+pp(12)*p(3)+pp(11)*p(2)+pp(9)*p(1)+&
     pp(8)*q(10)+qp(17)*q(7)+qp(15)*q(8)+qp(12)*q(4)+qp(11)*q(3)+pp(4)*q(6)+&
    qp(7)*2*q(1)+qp(5)*q(2)
    mp(32)=pp(6)
    mp(31)=pp(6)
    mp(30)=pp(6)
    mp(29)=qp(10)
    mp(28)=pp(3)+mp(32)*m(1)
    mp(27)=pp(3)
    mp(26)=pp(3)+mp(31)*m(1)
    mp(25)=pp(3)
    mp(24)=pp(3)+mp(30)*m(2)
    mp(23)=pp(3)
    mp(22)=pp(2)
    mp(21)=pp(2)
    mp(20)=pp(2)
    mp(19)=pp(2)
    mp(18)=pp(2)
    mp(17)=pp(2)
    mp(16)=pp(1)
    mp(15)=qp(6)+mp(29)*m(4)+mp(28)*m(2)+mp(27)*m(1)
    mp(14)=qp(6)+mp(26)*m(3)+mp(25)*m(1)
    mp(13)=qp(6)+mp(24)*m(3)+mp(23)*m(2)
    mp(12)=qp(4)+mp(22)*m(3)+mp(21)*m(2)
    mp(11)=qp(4)+mp(20)*m(3)+mp(19)*m(1)
    mp(10)=qp(4)+mp(18)*m(2)+mp(17)*m(1)
    mp(9)=qp(3)+mp(16)*m(1)
    mp(8)=qp(3)
    mp(7)=qp(3)
    mp(6)=qp(2)+mp(15)*m(5)+mp(14)*m(4)+mp(12)*m(1)
    mp(5)=qp(2)+mp(15)*m(6)+mp(13)*m(4)+mp(11)*m(2)
    mp(4)=qp(2)+mp(29)*m(15)+mp(14)*m(6)+mp(13)*m(5)+mp(10)*m(3)
    mp(3)=qp(1)+mp(26)*m(14)+mp(24)*m(13)+mp(22)*m(12)+mp(20)*m(11)+mp(10)*m(4)+&
    mp(9)*m(2)+mp(8)*m(1)
    mp(2)=qp(1)+mp(30)*m(24)+mp(28)*m(15)+mp(23)*m(13)+mp(21)*m(12)+mp(18)*m(10)+&
    mp(11)*m(5)+mp(9)*m(3)+mp(7)*m(1)
    mp(1)=qp(1)+mp(32)*m(28)+mp(31)*m(26)+mp(27)*m(15)+mp(25)*m(14)+mp(19)*m(11)+&
    mp(17)*m(10)+mp(16)*m(9)+mp(12)*m(6)+mp(8)*m(3)+mp(7)*m(2)
    xxp(1)=mp(6)*(-m(6)/a)*drdx31(1,1,xyz,r)+mp(5)*(-m(5)/a)*drdx31(1,2,xyz,r)+&
    mp(3)*(-m(3)/a)*drdx31(1,3,xyz,r)
    xxp(2)=mp(6)*(-m(6)/a)*drdx31(2,1,xyz,r)+mp(5)*(-m(5)/a)*drdx31(2,2,xyz,r)+&
    mp(3)*(-m(3)/a)*drdx31(2,3,xyz,r)
    xxp(3)=mp(6)*(-m(6)/a)*drdx31(3,1,xyz,r)+mp(5)*(-m(5)/a)*drdx31(3,2,xyz,r)+&
    mp(3)*(-m(3)/a)*drdx31(3,3,xyz,r)
    xxp(4)=mp(6)*(-m(6)/a)*drdx31(4,1,xyz,r)+mp(4)*(-m(4)/a)*drdx31(4,4,xyz,r)+&
    mp(2)*(-m(2)/a)*drdx31(4,5,xyz,r)
    xxp(5)=mp(6)*(-m(6)/a)*drdx31(5,1,xyz,r)+mp(4)*(-m(4)/a)*drdx31(5,4,xyz,r)+&
    mp(2)*(-m(2)/a)*drdx31(5,5,xyz,r)
    xxp(6)=mp(6)*(-m(6)/a)*drdx31(6,1,xyz,r)+mp(4)*(-m(4)/a)*drdx31(6,4,xyz,r)+&
    mp(2)*(-m(2)/a)*drdx31(6,5,xyz,r)
    xxp(7)=mp(5)*(-m(5)/a)*drdx31(7,2,xyz,r)+mp(4)*(-m(4)/a)*drdx31(7,4,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx31(7,6,xyz,r)
    xxp(8)=mp(5)*(-m(5)/a)*drdx31(8,2,xyz,r)+mp(4)*(-m(4)/a)*drdx31(8,4,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx31(8,6,xyz,r)
    xxp(9)=mp(5)*(-m(5)/a)*drdx31(9,2,xyz,r)+mp(4)*(-m(4)/a)*drdx31(9,4,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx31(9,6,xyz,r)
    xxp(10)=mp(3)*(-m(3)/a)*drdx31(10,3,xyz,r)+mp(2)*(-m(2)/a)*drdx31(10,5,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx31(10,6,xyz,r)
    xxp(11)=mp(3)*(-m(3)/a)*drdx31(11,3,xyz,r)+mp(2)*(-m(2)/a)*drdx31(11,5,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx31(11,6,xyz,r)
    xxp(12)=mp(3)*(-m(3)/a)*drdx31(12,3,xyz,r)+mp(2)*(-m(2)/a)*drdx31(12,5,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx31(12,6,xyz,r)

    return
  end subroutine deriv_rev31

end module bemsa31
