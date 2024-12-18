module bemsa4
  implicit none

contains
  function emsav4(x,c) result(v)
    real,dimension(1:6)::x
    real,dimension(1:40)::p
    real,dimension(1:40)::c
    real::v
    
    call bemsav4(x,p)
    v = dot_product(p,c)
    
    return
  end function emsav4
  
  subroutine bemsav4(x,p)
    real,dimension(1:6)::x
    real,dimension(1:87)::m
    real,dimension(1:40)::p
    real,dimension(1:12)::q
    ! ::::::::::::::::::::

    call evmono4(x,m)
    call evpoly4(m,q,p)
    
    return
  end subroutine bemsav4
  
  subroutine evmono4(x,m)
    real,dimension(1:6)::x
    real,dimension(1:87)::m
    !::::::::::::::::::::
  
    m(1) = x(6)
    m(2) = x(5)
    m(3) = x(4)
    m(4) = x(3)
    m(5) = x(2)
    m(6) = x(1)
    m(7) = m(3)*m(4)
    m(8) = m(2)*m(5)
    m(9) = m(1)*m(6)
    m(10) = m(1)*m(2)
    m(11) = m(1)*m(3)
    m(12) = m(2)*m(3)
    m(13) = m(1)*m(4)
    m(14) = m(2)*m(4)
    m(15) = m(1)*m(5)
    m(16) = m(3)*m(5)
    m(17) = m(4)*m(5)
    m(18) = m(2)*m(6)
    m(19) = m(3)*m(6)
    m(20) = m(4)*m(6)
    m(21) = m(5)*m(6)
    m(22) = m(1)*m(7)
    m(23) = m(2)*m(7)
    m(24) = m(1)*m(8)
    m(25) = m(2)*m(16)
    m(26) = m(2)*m(17)
    m(27) = m(3)*m(17)
    m(28) = m(1)*m(18)
    m(29) = m(1)*m(19)
    m(30) = m(1)*m(20)
    m(31) = m(3)*m(20)
    m(32) = m(1)*m(21)
    m(33) = m(2)*m(21)
    m(34) = m(1)*m(12)
    m(35) = m(1)*m(17)
    m(36) = m(2)*m(20)
    m(37) = m(3)*m(21)
    m(38) = m(1)*m(14)
    m(39) = m(1)*m(16)
    m(40) = m(2)*m(19)
    m(41) = m(4)*m(21)
    m(42) = m(2)*m(27)
    m(43) = m(1)*m(31)
    m(44) = m(1)*m(33)
    m(45) = m(1)*m(23)
    m(46) = m(1)*m(25)
    m(47) = m(1)*m(26)
    m(48) = m(1)*m(27)
    m(49) = m(1)*m(40)
    m(50) = m(1)*m(36)
    m(51) = m(2)*m(31)
    m(52) = m(1)*m(37)
    m(53) = m(2)*m(37)
    m(54) = m(1)*m(41)
    m(55) = m(2)*m(41)
    m(56) = m(3)*m(41)
    m(57) = m(1)*m(42)
    m(58) = m(1)*m(51)
    m(59) = m(1)*m(53)
    m(60) = m(1)*m(55)
    m(61) = m(1)*m(56)
    m(62) = m(2)*m(56)
    m(63) = m(1)*m(62)
    m(64) = m(2)*m(57)
    m(65) = m(3)*m(57)
    m(66) = m(4)*m(57)
    m(67) = m(5)*m(57)
    m(68) = m(1)*m(58)
    m(69) = m(3)*m(58)
    m(70) = m(4)*m(58)
    m(71) = m(1)*m(59)
    m(72) = m(2)*m(59)
    m(73) = m(1)*m(60)
    m(74) = m(2)*m(60)
    m(75) = m(1)*m(61)
    m(76) = m(2)*m(62)
    m(77) = m(3)*m(61)
    m(78) = m(3)*m(62)
    m(79) = m(4)*m(61)
    m(80) = m(4)*m(62)
    m(81) = m(5)*m(59)
    m(82) = m(5)*m(60)
    m(83) = m(5)*m(62)
    m(84) = m(6)*m(58)
    m(85) = m(6)*m(59)
    m(86) = m(6)*m(60)
    m(87) = m(6)*m(61)

    return
  end subroutine evmono4

  subroutine evpoly4(m,q,p)
    real,dimension(1:87)::m
    real,dimension(1:40)::p
    real,dimension(1:12)::q
    !::::::::::::::::::::

    q(1) = m(1) + m(2) + m(3) + m(4) + m(5) + m(6)
    q(2) = m(7) + m(8) + m(9)
    q(3) = m(10) + m(11) + m(12) + m(13) + m(14) + m(15) + m(16) +  &
         m(17) + m(18) + m(19) + m(20) + m(21)
    q(4) = q(1)*q(1) - q(3) - q(2) - q(3) - q(2)
    p(1) = m(22) + m(23) + m(24) + m(25) + m(26) + m(27) + m(28) +  &
         m(29) + m(30) + m(31) + m(32) + m(33)
    q(5) = m(34) + m(35) + m(36) + m(37)
    p(2) = m(38) + m(39) + m(40) + m(41)
    q(6) = q(1)*q(2) - p(1)
    q(7) = q(1)*q(3) - q(5) - p(2) - p(1) - q(5) - p(2) - p(1) -  &
         q(5) - p(2)
    q(8) = q(1)*q(4) - q(7) - q(6)
    p(3) = m(42) + m(43) + m(44)
    p(4) = m(45) + m(46) + m(47) + m(48) + m(49) + m(50) + m(51) +  &
         m(52) + m(53) + m(54) + m(55) + m(56)
    p(5) = q(2)*q(3) - p(4)
    p(6) = q(1)*p(1) - p(4) - p(3) - p(5) - p(4) - p(3) - p(3) - p(3)
    q(9) = q(1)*q(5) - p(4)
    p(7) = q(1)*p(2) - p(4)
    q(10) = q(2)*q(2) - p(3) - p(3)
    q(11) = q(3)*q(3) - p(4) - p(3) - q(9) - p(7) - p(6) - p(4) -  &
         p(3) - q(9) - p(7) - p(6) - p(4) - p(3) - p(4) - p(3)
    q(12) = q(3)*q(4) - q(9) - p(7) - p(5)
    p(8) = m(57) + m(58) + m(59) + m(60) + m(61) + m(62)
    p(9) = q(1)*p(3) - p(8)
    p(10) = q(2)*q(5)
    p(11) = q(2)*p(2)
    p(12) = q(1)*p(4) - p(8) - p(10) - p(11) - p(8) - p(8) - p(8)
    p(13) = q(2)*p(1) - p(8) - p(9) - p(8)
    p(14) = q(3)*p(1) - p(8) - p(12) - p(10) - p(11) - p(9) - p(8) -  &
         p(10) - p(11) - p(9) - p(8) - p(8)
    p(15) = q(3)*p(2) - p(8) - p(12) - p(8)
    p(16) = q(2)*q(7) - p(12) - p(14)
    p(17) = q(1)*p(6) - p(12) - p(9) - p(14)
    p(18) = q(4)*p(2) - p(10)
    p(19) = m(63)
    p(20) = m(64) + m(65) + m(66) + m(67) + m(68) + m(69) + m(70) +  &
         m(71) + m(72) + m(73) + m(74) + m(75) + m(76) + m(77) +  &
         m(78) + m(79) + m(80) + m(81) + m(82) + m(83) + m(84) +  &
         m(85) + m(86) + m(87)
    p(21) = q(1)*p(8) - p(19) - p(20) - p(19) - p(19) - p(19) -  &
         p(19) - p(19)
    p(22) = q(2)*p(3) - p(19) - p(19) - p(19)
    p(23) = q(2)*p(4) - p(20)
    p(24) = q(3)*p(3) - p(20)
    p(25) = p(1)*q(5) - p(20)
    p(26) = p(1)*p(2) - p(20)
    p(27) = q(5)*p(2) - p(19) - p(19) - p(19) - p(19)
    p(28) = q(4)*p(3) - p(21)
    p(29) = q(2)*q(9) - p(25)
    p(30) = q(2)*p(7) - p(26)
    p(31) = q(4)*p(4) - p(20) - p(29) - p(30)
    p(32) = q(2)*p(6) - p(21) - p(28) - p(21)
    p(33) = p(2)*p(2) - p(21) - p(21)
    p(34) = q(3)*q(10) - p(23)
    p(35) = q(2)*q(11) - p(27)
    p(36) = q(3)*p(6) - p(20) - p(31) - p(25) - p(26) - p(24) - p(24)
    p(37) = p(2)*q(7) - p(20) - p(31) - p(25)
    p(38) = q(2)*q(12) - p(31) - p(36)
    p(39) = q(1)*p(17) - p(31) - p(28) - p(36)
    p(40) = p(2)*q(8) - p(29)

    return
  end subroutine evpoly4

  subroutine devmono4(dm,m,flag,xyz,a,r)
    real,dimension(1:87)::dm
    real,dimension(1:87)::m
    real::xyz(4,3),a,r(4,4)
    integer::flag
    !::::::::::::::::::::

    dm(1) =  - m(1)/a*drdx4(flag,6,xyz,r)
    dm(2) =  - m(2)/a*drdx4(flag,5,xyz,r)
    dm(3) =  - m(3)/a*drdx4(flag,4,xyz,r)
    dm(4) =  - m(4)/a*drdx4(flag,3,xyz,r)
    dm(5) =  - m(5)/a*drdx4(flag,2,xyz,r)
    dm(6) =  - m(6)/a*drdx4(flag,1,xyz,r)
    dm(7) = dm(3)*m(4) + m(3)*dm(4)
    dm(8) = dm(2)*m(5) + m(2)*dm(5)
    dm(9) = dm(1)*m(6) + m(1)*dm(6)
    dm(10) = dm(1)*m(2) + m(1)*dm(2)
    dm(11) = dm(1)*m(3) + m(1)*dm(3)
    dm(12) = dm(2)*m(3) + m(2)*dm(3)
    dm(13) = dm(1)*m(4) + m(1)*dm(4)
    dm(14) = dm(2)*m(4) + m(2)*dm(4)
    dm(15) = dm(1)*m(5) + m(1)*dm(5)
    dm(16) = dm(3)*m(5) + m(3)*dm(5)
    dm(17) = dm(4)*m(5) + m(4)*dm(5)
    dm(18) = dm(2)*m(6) + m(2)*dm(6)
    dm(19) = dm(3)*m(6) + m(3)*dm(6)
    dm(20) = dm(4)*m(6) + m(4)*dm(6)
    dm(21) = dm(5)*m(6) + m(5)*dm(6)
    dm(22) = dm(1)*m(7) + m(1)*dm(7)
    dm(23) = dm(2)*m(7) + m(2)*dm(7)
    dm(24) = dm(1)*m(8) + m(1)*dm(8)
    dm(25) = dm(2)*m(16) + m(2)*dm(16)
    dm(26) = dm(2)*m(17) + m(2)*dm(17)
    dm(27) = dm(3)*m(17) + m(3)*dm(17)
    dm(28) = dm(1)*m(18) + m(1)*dm(18)
    dm(29) = dm(1)*m(19) + m(1)*dm(19)
    dm(30) = dm(1)*m(20) + m(1)*dm(20)
    dm(31) = dm(3)*m(20) + m(3)*dm(20)
    dm(32) = dm(1)*m(21) + m(1)*dm(21)
    dm(33) = dm(2)*m(21) + m(2)*dm(21)
    dm(34) = dm(1)*m(12) + m(1)*dm(12)
    dm(35) = dm(1)*m(17) + m(1)*dm(17)
    dm(36) = dm(2)*m(20) + m(2)*dm(20)
    dm(37) = dm(3)*m(21) + m(3)*dm(21)
    dm(38) = dm(1)*m(14) + m(1)*dm(14)
    dm(39) = dm(1)*m(16) + m(1)*dm(16)
    dm(40) = dm(2)*m(19) + m(2)*dm(19)
    dm(41) = dm(4)*m(21) + m(4)*dm(21)
    dm(42) = dm(2)*m(27) + m(2)*dm(27)
    dm(43) = dm(1)*m(31) + m(1)*dm(31)
    dm(44) = dm(1)*m(33) + m(1)*dm(33)
    dm(45) = dm(1)*m(23) + m(1)*dm(23)
    dm(46) = dm(1)*m(25) + m(1)*dm(25)
    dm(47) = dm(1)*m(26) + m(1)*dm(26)
    dm(48) = dm(1)*m(27) + m(1)*dm(27)
    dm(49) = dm(1)*m(40) + m(1)*dm(40)
    dm(50) = dm(1)*m(36) + m(1)*dm(36)
    dm(51) = dm(2)*m(31) + m(2)*dm(31)
    dm(52) = dm(1)*m(37) + m(1)*dm(37)
    dm(53) = dm(2)*m(37) + m(2)*dm(37)
    dm(54) = dm(1)*m(41) + m(1)*dm(41)
    dm(55) = dm(2)*m(41) + m(2)*dm(41)
    dm(56) = dm(3)*m(41) + m(3)*dm(41)
    dm(57) = dm(1)*m(42) + m(1)*dm(42)
    dm(58) = dm(1)*m(51) + m(1)*dm(51)
    dm(59) = dm(1)*m(53) + m(1)*dm(53)
    dm(60) = dm(1)*m(55) + m(1)*dm(55)
    dm(61) = dm(1)*m(56) + m(1)*dm(56)
    dm(62) = dm(2)*m(56) + m(2)*dm(56)
    dm(63) = dm(1)*m(62) + m(1)*dm(62)
    dm(64) = dm(2)*m(57) + m(2)*dm(57)
    dm(65) = dm(3)*m(57) + m(3)*dm(57)
    dm(66) = dm(4)*m(57) + m(4)*dm(57)
    dm(67) = dm(5)*m(57) + m(5)*dm(57)
    dm(68) = dm(1)*m(58) + m(1)*dm(58)
    dm(69) = dm(3)*m(58) + m(3)*dm(58)
    dm(70) = dm(4)*m(58) + m(4)*dm(58)
    dm(71) = dm(1)*m(59) + m(1)*dm(59)
    dm(72) = dm(2)*m(59) + m(2)*dm(59)
    dm(73) = dm(1)*m(60) + m(1)*dm(60)
    dm(74) = dm(2)*m(60) + m(2)*dm(60)
    dm(75) = dm(1)*m(61) + m(1)*dm(61)
    dm(76) = dm(2)*m(62) + m(2)*dm(62)
    dm(77) = dm(3)*m(61) + m(3)*dm(61)
    dm(78) = dm(3)*m(62) + m(3)*dm(62)
    dm(79) = dm(4)*m(61) + m(4)*dm(61)
    dm(80) = dm(4)*m(62) + m(4)*dm(62)
    dm(81) = dm(5)*m(59) + m(5)*dm(59)
    dm(82) = dm(5)*m(60) + m(5)*dm(60)
    dm(83) = dm(5)*m(62) + m(5)*dm(62)
    dm(84) = dm(6)*m(58) + m(6)*dm(58)
    dm(85) = dm(6)*m(59) + m(6)*dm(59)
    dm(86) = dm(6)*m(60) + m(6)*dm(60)
    dm(87) = dm(6)*m(61) + m(6)*dm(61)

    return
  end subroutine devmono4

  subroutine devpoly4(dm,q,p,dp)
    real,dimension(1:87)::dm
    real,dimension(1:40)::p,dp
    real,dimension(1:12)::q,dq
    !::::::::::::::::::::

    dq(1) = dm(1) + dm(2) + dm(3) + dm(4) + dm(5) + dm(6)
    dq(2) = dm(7) + dm(8) + dm(9)
    dq(3) = dm(10) + dm(11) + dm(12) + dm(13) + dm(14) + dm(15) +  &
         dm(16) + dm(17) + dm(18) + dm(19) + dm(20) + dm(21)
    dq(4) = dq(1)*q(1) + q(1)*dq(1) - dq(3) - dq(2) - dq(3) - dq(2)
    dp(1) = dm(22) + dm(23) + dm(24) + dm(25) + dm(26) + dm(27) +  &
         dm(28) + dm(29) + dm(30) + dm(31) + dm(32) + dm(33)
    dq(5) = dm(34) + dm(35) + dm(36) + dm(37)
    dp(2) = dm(38) + dm(39) + dm(40) + dm(41)
    dq(6) = dq(1)*q(2) + q(1)*dq(2) - dp(1)
    dq(7) = dq(1)*q(3) + q(1)*dq(3) - dq(5) - dp(2) - dp(1) - dq(5) -  &
         dp(2) - dp(1) - dq(5) - dp(2)
    dq(8) = dq(1)*q(4) + q(1)*dq(4) - dq(7) - dq(6)
    dp(3) = dm(42) + dm(43) + dm(44)
    dp(4) = dm(45) + dm(46) + dm(47) + dm(48) + dm(49) + dm(50) +  &
         dm(51) + dm(52) + dm(53) + dm(54) + dm(55) + dm(56)
    dp(5) = dq(2)*q(3) + q(2)*dq(3) - dp(4)
    dp(6) = dq(1)*p(1) + q(1)*dp(1) - dp(4) - dp(3) - dp(5) - dp(4) -  &
         dp(3) - dp(3) - dp(3)
    dq(9) = dq(1)*q(5) + q(1)*dq(5) - dp(4)
    dp(7) = dq(1)*p(2) + q(1)*dp(2) - dp(4)
    dq(10) = dq(2)*q(2) + q(2)*dq(2) - dp(3) - dp(3)
    dq(11) = dq(3)*q(3) + q(3)*dq(3) - dp(4) - dp(3) - dq(9) - dp(7) -  &
         dp(6) - dp(4) - dp(3) - dq(9) - dp(7) - dp(6) -  &
         dp(4) - dp(3) - dp(4) - dp(3)
    dq(12) = dq(3)*q(4) + q(3)*dq(4) - dq(9) - dp(7) - dp(5)
    dp(8) = dm(57) + dm(58) + dm(59) + dm(60) + dm(61) + dm(62)
    dp(9) = dq(1)*p(3) + q(1)*dp(3) - dp(8)
    dp(10) = dq(2)*q(5) + q(2)*dq(5)
    dp(11) = dq(2)*p(2) + q(2)*dp(2)
    dp(12) = dq(1)*p(4) + q(1)*dp(4) - dp(8) - dp(10) - dp(11) -  &
         dp(8) - dp(8) - dp(8)
    dp(13) = dq(2)*p(1) + q(2)*dp(1) - dp(8) - dp(9) - dp(8)
    dp(14) = dq(3)*p(1) + q(3)*dp(1) - dp(8) - dp(12) - dp(10) -  &
         dp(11) - dp(9) - dp(8) - dp(10) - dp(11) - dp(9) -  &
         dp(8) - dp(8)
    dp(15) = dq(3)*p(2) + q(3)*dp(2) - dp(8) - dp(12) - dp(8)
    dp(16) = dq(2)*q(7) + q(2)*dq(7) - dp(12) - dp(14)
    dp(17) = dq(1)*p(6) + q(1)*dp(6) - dp(12) - dp(9) - dp(14)
    dp(18) = dq(4)*p(2) + q(4)*dp(2) - dp(10)
    dp(19) = dm(63)
    dp(20) = dm(64) + dm(65) + dm(66) + dm(67) + dm(68) + dm(69) +  &
         dm(70) + dm(71) + dm(72) + dm(73) + dm(74) +  &
         dm(75) + dm(76) + dm(77) + dm(78) + dm(79) +  &
         dm(80) + dm(81) + dm(82) + dm(83) + dm(84) +  &
         dm(85) + dm(86) + dm(87)
    dp(21) = dq(1)*p(8) + q(1)*dp(8) - dp(19) - dp(20) - dp(19) -  &
         dp(19) - dp(19) - dp(19) - dp(19)
    dp(22) = dq(2)*p(3) + q(2)*dp(3) - dp(19) - dp(19) - dp(19)
    dp(23) = dq(2)*p(4) + q(2)*dp(4) - dp(20)
    dp(24) = dq(3)*p(3) + q(3)*dp(3) - dp(20)
    dp(25) = dp(1)*q(5) + p(1)*dq(5) - dp(20)
    dp(26) = dp(1)*p(2) + p(1)*dp(2) - dp(20)
    dp(27) = dq(5)*p(2) + q(5)*dp(2) - dp(19) - dp(19) - dp(19) -  &
         dp(19)
    dp(28) = dq(4)*p(3) + q(4)*dp(3) - dp(21)
    dp(29) = dq(2)*q(9) + q(2)*dq(9) - dp(25)
    dp(30) = dq(2)*p(7) + q(2)*dp(7) - dp(26)
    dp(31) = dq(4)*p(4) + q(4)*dp(4) - dp(20) - dp(29) - dp(30)
    dp(32) = dq(2)*p(6) + q(2)*dp(6) - dp(21) - dp(28) - dp(21)
    dp(33) = dp(2)*p(2) + p(2)*dp(2) - dp(21) - dp(21)
    dp(34) = dq(3)*q(10) + q(3)*dq(10) - dp(23)
    dp(35) = dq(2)*q(11) + q(2)*dq(11) - dp(27)
    dp(36) = dq(3)*p(6) + q(3)*dp(6) - dp(20) - dp(31) - dp(25) -  &
         dp(26) - dp(24) - dp(24)
    dp(37) = dp(2)*q(7) + p(2)*dq(7) - dp(20) - dp(31) - dp(25)
    dp(38) = dq(2)*q(12) + q(2)*dq(12) - dp(31) - dp(36)
    dp(39) = dq(1)*p(17) + q(1)*dp(17) - dp(31) - dp(28) - dp(36)
    dp(40) = dp(2)*q(8) + p(2)*dq(8) - dp(29)

    return
  end subroutine devpoly4

  function drdx4(flag,xindex,xyz,r)
    integer i,j,flag,xindex,xyzind,matom,m
    real::xyz(4,3),r(4,4),drdx4

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

    drdx4 = 0.d0
    if (matom.eq.i.or.matom.eq.j) then
       drdx4=(xyz(i,xyzind)-xyz(j,xyzind))/r(i,j)
       if (matom.eq.j) then
          drdx4 = -drdx4
       endif
    endif

    return
  end function

  subroutine deriv_rev4(c,m,q,p,xyz,a,r,xxp)
    real::c(1:40),m(1:87),p(1:40)
    real::xyz(4,3),r(4,4),a
    !::::::::::::::::::::
    real::pp(1:40),mp(1:87),xxp(1:12)
    real::qp(1:12),q(1:12)

    qp(:)=0.d0
    pp(:)=0.d0
    mp(:)=0.d0
    xxp(:)=0.d0
    pp(40)=c(40)
    pp(39)=c(39)
    pp(38)=c(38)
    pp(37)=c(37)
    pp(36)=c(36)+pp(39)*(-1.d0)+pp(38)*(-1.d0)
    pp(35)=c(35)
    pp(34)=c(34)
    pp(33)=c(33)
    pp(32)=c(32)
    pp(31)=c(31)+pp(39)*(-1.d0)+pp(38)*(-1.d0)+pp(37)*(-1.d0)+pp(36)*(-1.d0)
    pp(30)=c(30)+pp(31)*(-1.d0)
    pp(29)=c(29)+pp(40)*(-1.d0)+pp(31)*(-1.d0)
    pp(28)=c(28)+pp(39)*(-1.d0)+pp(32)*(-1.d0)
    pp(27)=c(27)+pp(35)*(-1.d0)
    pp(26)=c(26)+pp(36)*(-1.d0)+pp(30)*(-1.d0)
    pp(25)=c(25)+pp(37)*(-1.d0)+pp(36)*(-1.d0)+pp(29)*(-1.d0)
    pp(24)=c(24)+pp(36)*(-2.d0)
    pp(23)=c(23)+pp(34)*(-1.d0)
    pp(22)=c(22)
    pp(21)=c(21)+pp(33)*(-2.d0)+pp(32)*(-2.d0)+pp(28)*(-1.d0)
    pp(20)=c(20)+pp(37)*(-1.d0)+pp(36)*(-1.d0)+pp(31)*(-1.d0)+pp(26)*(-1.d0)+&
    pp(25)*(-1.d0)+pp(24)*(-1.d0)+pp(23)*(-1.d0)+pp(21)*(-1.d0)
    pp(19)=c(19)+pp(27)*(-4.d0)+pp(22)*(-3.d0)+pp(21)*(-6.d0)
    pp(18)=c(18)
    pp(17)=c(17)+pp(39)*q(1)
    pp(16)=c(16)
    pp(15)=c(15)
    pp(14)=c(14)+pp(17)*(-1.d0)+pp(16)*(-1.d0)
    pp(13)=c(13)
    pp(12)=c(12)+pp(17)*(-1.d0)+pp(16)*(-1.d0)+pp(15)*(-1.d0)+pp(14)*(-1.d0)
    pp(11)=c(11)+pp(14)*(-2.d0)+pp(12)*(-1.d0)
    pp(10)=c(10)+pp(18)*(-1.d0)+pp(14)*(-2.d0)+pp(12)*(-1.d0)
    pp(9)=c(9)+pp(17)*(-1.d0)+pp(14)*(-2.d0)+pp(13)*(-1.d0)
    pp(8)=c(8)+pp(21)*q(1)+pp(15)*(-2.d0)+pp(14)*(-4.d0)+pp(13)*(-2.d0)+&
    pp(12)*(-4.d0)+pp(9)*(-1.d0)
    qp(12)=pp(38)*q(2)
    qp(11)=pp(35)*q(2)
    qp(10)=pp(34)*q(3)
    pp(7)=c(7)+pp(30)*q(2)+qp(12)*(-1.d0)+qp(11)*(-2.d0)
    qp(9)=pp(29)*q(2)+qp(12)*(-1.d0)+qp(11)*(-2.d0)
    pp(6)=c(6)+pp(36)*q(3)+pp(32)*q(2)+pp(17)*q(1)+qp(11)*(-2.d0)
    pp(5)=c(5)+qp(12)*(-1.d0)+pp(6)*(-1.d0)
    pp(4)=c(4)+pp(31)*q(4)+pp(23)*q(2)+pp(12)*q(1)+qp(11)*(-4.d0)+pp(7)*(-1.d0)+&
    qp(9)*(-1.d0)+pp(6)*(-2.d0)+pp(5)*(-1.d0)
    pp(3)=c(3)+pp(28)*q(4)+pp(24)*q(3)+pp(22)*q(2)+pp(9)*q(1)+qp(11)*(-4.d0)+&
    qp(10)*(-2.d0)+pp(6)*(-4.d0)
    qp(8)=pp(40)*p(2)
    qp(7)=pp(37)*p(2)+pp(16)*q(2)+qp(8)*(-1.d0)
    qp(6)=qp(8)*(-1.d0)
    pp(2)=c(2)+pp(40)*q(8)+pp(37)*q(7)+pp(33)*2*p(2)+pp(27)*q(5)+pp(26)*p(1)+&
    pp(18)*q(4)+pp(15)*q(3)+pp(11)*q(2)+pp(7)*q(1)+qp(7)*(-3.d0)
    qp(5)=pp(27)*p(2)+pp(25)*p(1)+pp(10)*q(2)+qp(9)*q(1)+qp(7)*(-3.d0)
    pp(1)=c(1)+pp(26)*p(2)+pp(25)*q(5)+pp(14)*q(3)+pp(13)*q(2)+pp(6)*q(1)+&
    qp(7)*(-2.d0)+qp(6)*(-1.d0)
    qp(4)=pp(31)*p(4)+pp(28)*p(3)+pp(18)*p(2)+qp(12)*q(3)+qp(8)*q(1)
    qp(3)=pp(36)*p(6)+pp(34)*q(10)+pp(24)*p(3)+pp(15)*p(2)+pp(14)*p(1)+qp(12)*q(4)+&
    qp(11)*2*q(3)+pp(5)*q(2)+qp(7)*q(1)+qp(4)*(-2.d0)
    qp(2)=pp(38)*q(12)+pp(35)*q(11)+pp(32)*p(6)+pp(30)*p(7)+pp(29)*q(9)+pp(23)*p(4)+&
     pp(22)*p(3)+pp(16)*q(7)+pp(13)*p(1)+pp(11)*p(2)+pp(10)*q(5)+qp(10)*2*q(2)+&
    pp(5)*q(3)+qp(6)*q(1)+qp(4)*(-2.d0)
    qp(1)=pp(39)*p(17)+pp(21)*p(8)+pp(17)*p(6)+pp(12)*p(4)+pp(9)*p(3)+pp(7)*p(2)+&
    qp(9)*q(5)+pp(6)*p(1)+qp(8)*q(4)+qp(7)*q(3)+qp(6)*q(2)+qp(4)*2*q(1)
    mp(87)=pp(20)
    mp(86)=pp(20)
    mp(85)=pp(20)
    mp(84)=pp(20)
    mp(83)=pp(20)
    mp(82)=pp(20)
    mp(81)=pp(20)
    mp(80)=pp(20)
    mp(79)=pp(20)
    mp(78)=pp(20)
    mp(77)=pp(20)
    mp(76)=pp(20)
    mp(75)=pp(20)
    mp(74)=pp(20)
    mp(73)=pp(20)
    mp(72)=pp(20)
    mp(71)=pp(20)
    mp(70)=pp(20)
    mp(69)=pp(20)
    mp(68)=pp(20)
    mp(67)=pp(20)
    mp(66)=pp(20)
    mp(65)=pp(20)
    mp(64)=pp(20)
    mp(63)=pp(19)
    mp(62)=pp(8)+mp(83)*m(5)+mp(80)*m(4)+mp(78)*m(3)+mp(76)*m(2)+mp(63)*m(1)
    mp(61)=pp(8)+mp(87)*m(6)+mp(79)*m(4)+mp(77)*m(3)+mp(75)*m(1)
    mp(60)=pp(8)+mp(86)*m(6)+mp(82)*m(5)+mp(74)*m(2)+mp(73)*m(1)
    mp(59)=pp(8)+mp(85)*m(6)+mp(81)*m(5)+mp(72)*m(2)+mp(71)*m(1)
    mp(58)=pp(8)+mp(84)*m(6)+mp(70)*m(4)+mp(69)*m(3)+mp(68)*m(1)
    mp(57)=pp(8)+mp(67)*m(5)+mp(66)*m(4)+mp(65)*m(3)+mp(64)*m(2)
    mp(56)=pp(4)+mp(62)*m(2)+mp(61)*m(1)
    mp(55)=pp(4)+mp(60)*m(1)
    mp(54)=pp(4)
    mp(53)=pp(4)+mp(59)*m(1)
    mp(52)=pp(4)
    mp(51)=pp(4)+mp(58)*m(1)
    mp(50)=pp(4)
    mp(49)=pp(4)
    mp(48)=pp(4)
    mp(47)=pp(4)
    mp(46)=pp(4)
    mp(45)=pp(4)
    mp(44)=pp(3)
    mp(43)=pp(3)
    mp(42)=pp(3)+mp(57)*m(1)
    mp(41)=pp(2)+mp(56)*m(3)+mp(55)*m(2)+mp(54)*m(1)
    mp(40)=pp(2)+mp(49)*m(1)
    mp(39)=pp(2)
    mp(38)=pp(2)
    mp(37)=qp(5)+mp(53)*m(2)+mp(52)*m(1)
    mp(36)=qp(5)+mp(50)*m(1)
    mp(35)=qp(5)
    mp(34)=qp(5)
    mp(33)=pp(1)+mp(44)*m(1)
    mp(32)=pp(1)
    mp(31)=pp(1)+mp(51)*m(2)+mp(43)*m(1)
    mp(30)=pp(1)
    mp(29)=pp(1)
    mp(28)=pp(1)
    mp(27)=pp(1)+mp(48)*m(1)+mp(42)*m(2)
    mp(26)=pp(1)+mp(47)*m(1)
    mp(25)=pp(1)+mp(46)*m(1)
    mp(24)=pp(1)
    mp(23)=pp(1)+mp(45)*m(1)
    mp(22)=pp(1)
    mp(21)=qp(3)+mp(41)*m(4)+mp(37)*m(3)+mp(33)*m(2)+mp(32)*m(1)
    mp(20)=qp(3)+mp(36)*m(2)+mp(31)*m(3)+mp(30)*m(1)
    mp(19)=qp(3)+mp(40)*m(2)+mp(29)*m(1)
    mp(18)=qp(3)+mp(28)*m(1)
    mp(17)=qp(3)+mp(35)*m(1)+mp(27)*m(3)+mp(26)*m(2)
    mp(16)=qp(3)+mp(39)*m(1)+mp(25)*m(2)
    mp(15)=qp(3)
    mp(14)=qp(3)+mp(38)*m(1)
    mp(13)=qp(3)
    mp(12)=qp(3)+mp(34)*m(1)
    mp(11)=qp(3)
    mp(10)=qp(3)
    mp(9)=qp(2)
    mp(8)=qp(2)+mp(24)*m(1)
    mp(7)=qp(2)+mp(23)*m(2)+mp(22)*m(1)
    mp(6)=qp(1)+mp(87)*m(61)+mp(86)*m(60)+mp(85)*m(59)+mp(84)*m(58)+mp(21)*m(5)+&
    mp(20)*m(4)+mp(19)*m(3)+mp(18)*m(2)+mp(9)*m(1)
    mp(5)=qp(1)+mp(83)*m(62)+mp(82)*m(60)+mp(81)*m(59)+mp(67)*m(57)+mp(21)*m(6)+&
    mp(17)*m(4)+mp(16)*m(3)+mp(15)*m(1)+mp(8)*m(2)
    mp(4)=qp(1)+mp(80)*m(62)+mp(79)*m(61)+mp(70)*m(58)+mp(66)*m(57)+mp(41)*m(21)+&
    mp(20)*m(6)+mp(17)*m(5)+mp(14)*m(2)+mp(13)*m(1)+mp(7)*m(3)
    mp(3)=qp(1)+mp(78)*m(62)+mp(77)*m(61)+mp(69)*m(58)+mp(65)*m(57)+mp(56)*m(41)+&
     mp(37)*m(21)+mp(31)*m(20)+mp(27)*m(17)+mp(19)*m(6)+mp(16)*m(5)+mp(12)*m(2)+&
    mp(11)*m(1)+mp(7)*m(4)
    mp(2)=qp(1)+mp(76)*m(62)+mp(74)*m(60)+mp(72)*m(59)+mp(64)*m(57)+mp(62)*m(56)+&
     mp(55)*m(41)+mp(53)*m(37)+mp(51)*m(31)+mp(42)*m(27)+mp(40)*m(19)+mp(36)*m(20)+&
     mp(33)*m(21)+mp(26)*m(17)+mp(25)*m(16)+mp(23)*m(7)+mp(18)*m(6)+mp(14)*m(4)+&
    mp(12)*m(3)+mp(10)*m(1)+mp(8)*m(5)
    mp(1)=qp(1)+mp(75)*m(61)+mp(73)*m(60)+mp(71)*m(59)+mp(68)*m(58)+mp(63)*m(62)+&
     mp(61)*m(56)+mp(60)*m(55)+mp(59)*m(53)+mp(58)*m(51)+mp(57)*m(42)+mp(54)*m(41)+&
     mp(52)*m(37)+mp(50)*m(36)+mp(49)*m(40)+mp(48)*m(27)+mp(47)*m(26)+mp(46)*m(25)+&
     mp(45)*m(23)+mp(44)*m(33)+mp(43)*m(31)+mp(39)*m(16)+mp(38)*m(14)+mp(35)*m(17)+&
     mp(34)*m(12)+mp(32)*m(21)+mp(30)*m(20)+mp(29)*m(19)+mp(28)*m(18)+mp(24)*m(8)+&
    mp(22)*m(7)+mp(15)*m(5)+mp(13)*m(4)+mp(11)*m(3)+mp(10)*m(2)+mp(9)*m(6)
    xxp(1)=mp(6)*(-m(6)/a)*drdx4(1,1,xyz,r)+mp(5)*(-m(5)/a)*drdx4(1,2,xyz,r)+&
    mp(4)*(-m(4)/a)*drdx4(1,3,xyz,r)
    xxp(2)=mp(6)*(-m(6)/a)*drdx4(2,1,xyz,r)+mp(5)*(-m(5)/a)*drdx4(2,2,xyz,r)+&
    mp(4)*(-m(4)/a)*drdx4(2,3,xyz,r)
    xxp(3)=mp(6)*(-m(6)/a)*drdx4(3,1,xyz,r)+mp(5)*(-m(5)/a)*drdx4(3,2,xyz,r)+&
    mp(4)*(-m(4)/a)*drdx4(3,3,xyz,r)
    xxp(4)=mp(6)*(-m(6)/a)*drdx4(4,1,xyz,r)+mp(3)*(-m(3)/a)*drdx4(4,4,xyz,r)+&
    mp(2)*(-m(2)/a)*drdx4(4,5,xyz,r)
    xxp(5)=mp(6)*(-m(6)/a)*drdx4(5,1,xyz,r)+mp(3)*(-m(3)/a)*drdx4(5,4,xyz,r)+&
    mp(2)*(-m(2)/a)*drdx4(5,5,xyz,r)
    xxp(6)=mp(6)*(-m(6)/a)*drdx4(6,1,xyz,r)+mp(3)*(-m(3)/a)*drdx4(6,4,xyz,r)+&
    mp(2)*(-m(2)/a)*drdx4(6,5,xyz,r)
    xxp(7)=mp(5)*(-m(5)/a)*drdx4(7,2,xyz,r)+mp(3)*(-m(3)/a)*drdx4(7,4,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx4(7,6,xyz,r)
    xxp(8)=mp(5)*(-m(5)/a)*drdx4(8,2,xyz,r)+mp(3)*(-m(3)/a)*drdx4(8,4,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx4(8,6,xyz,r)
    xxp(9)=mp(5)*(-m(5)/a)*drdx4(9,2,xyz,r)+mp(3)*(-m(3)/a)*drdx4(9,4,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx4(9,6,xyz,r)
    xxp(10)=mp(4)*(-m(4)/a)*drdx4(10,3,xyz,r)+mp(2)*(-m(2)/a)*drdx4(10,5,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx4(10,6,xyz,r)
    xxp(11)=mp(4)*(-m(4)/a)*drdx4(11,3,xyz,r)+mp(2)*(-m(2)/a)*drdx4(11,5,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx4(11,6,xyz,r)
    xxp(12)=mp(4)*(-m(4)/a)*drdx4(12,3,xyz,r)+mp(2)*(-m(2)/a)*drdx4(12,5,xyz,r)+&
    mp(1)*(-m(1)/a)*drdx4(12,6,xyz,r)

    return
  end subroutine deriv_rev4

end module bemsa4
