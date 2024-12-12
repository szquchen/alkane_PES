module bemsa3
  implicit none

contains
  function emsav3(x,c) result(v)
    real,dimension(1:3)::x
    real,dimension(1:32)::p
    real,dimension(1:32)::c
    real::v
    
    call bemsav3(x,p)
    v = dot_product(p,c)
    
    return
  end function emsav3
  
  subroutine bemsav3(x,p)
    real,dimension(1:3)::x
    real,dimension(1:7)::m
    real,dimension(1:32)::p
    real,dimension(1:6)::q
    ! ::::::::::::::::::::

    call evmono3(x,m)
    call evpoly3(m,q,p)
    
    return
  end subroutine bemsav3
  
  subroutine evmono3(x,m)
    real,dimension(1:3)::x
    real,dimension(1:7)::m
    !::::::::::::::::::::
  
    m(1) = x(3)
    m(2) = x(2)
    m(3) = x(1)
    m(4) = m(1)*m(2)
    m(5) = m(1)*m(3)
    m(6) = m(2)*m(3)
    m(7) = m(1)*m(6)

    return
  end subroutine evmono3

  subroutine evpoly3(m,q,p)
    real,dimension(1:7)::m
    real,dimension(1:32)::p
    real,dimension(1:6)::q
    !::::::::::::::::::::

    q(1) = m(1) + m(2) + m(3)
    p(1) = m(4) + m(5) + m(6)
    q(2) = q(1)*q(1) - p(1) - p(1)
    p(2) = m(7)
    p(3) = q(1)*p(1) - p(2) - p(2) - p(2)
    q(3) = q(1)*q(2) - p(3)
    p(4) = p(2)*q(1)
    p(5) = p(1)*p(1) - p(4) - p(4)
    p(6) = p(1)*q(2) - p(4)
    q(4) = q(1)*q(3) - p(6)
    p(7) = p(2)*p(1)
    p(8) = p(2)*q(2)
    p(9) = q(1)*p(5) - p(7)
    p(10) = p(1)*q(3) - p(8)
    q(5) = q(1)*q(4) - p(10)
    p(11) = p(2)*p(2)
    p(12) = p(2)*p(3)
    p(13) = p(2)*q(3)
    p(14) = p(1)*p(5) - p(12)
    p(15) = q(1)*p(9) - p(12) - p(14) - p(14)
    p(16) = p(1)*q(4) - p(13)
    q(6) = q(1)*q(5) - p(16)
    p(17) = p(2)*p(4)
    p(18) = p(2)*p(5)
    p(19) = p(2)*p(6)
    p(20) = p(2)*q(4)
    p(21) = q(1)*p(14) - p(18)
    p(22) = q(3)*p(5) - p(17)
    p(23) = p(1)*q(5) - p(20)
    p(24) = p(2)*p(7)
    p(25) = p(2)*p(8)
    p(26) = p(2)*p(9)
    p(27) = p(2)*p(10)
    p(28) = p(2)*q(5)
    p(29) = p(1)*p(14) - p(26)
    p(30) = q(2)*p(14) - p(24)
    p(31) = p(5)*q(4) - p(25)
    p(32) = p(1)*q(6) - p(28)

    return
  end subroutine evpoly3

  subroutine devmono3(dm,m,flag,xyz,a,r)
    real,dimension(1:7)::dm
    real,dimension(1:7)::m
    real::xyz(3,3),a,r(3,3)
    integer::flag
    !::::::::::::::::::::

    dm(1) =  - m(1)/a*drdx3(flag,3,xyz,r)
    dm(2) =  - m(2)/a*drdx3(flag,2,xyz,r)
    dm(3) =  - m(3)/a*drdx3(flag,1,xyz,r)
    dm(4) = dm(1)*m(2) + m(1)*dm(2)
    dm(5) = dm(1)*m(3) + m(1)*dm(3)
    dm(6) = dm(2)*m(3) + m(2)*dm(3)
    dm(7) = dm(1)*m(6) + m(1)*dm(6)

    return
  end subroutine devmono3

  subroutine devpoly3(dm,q,p,dp)
    real,dimension(1:7)::dm
    real,dimension(1:32)::p,dp
    real,dimension(1:6)::q,dq
    !::::::::::::::::::::

    dq(1) = dm(1) + dm(2) + dm(3)
    dp(1) = dm(4) + dm(5) + dm(6)
    dq(2) = dq(1)*q(1) + q(1)*dq(1) - dp(1) - dp(1)
    dp(2) = dm(7)
    dp(3) = dq(1)*p(1) + q(1)*dp(1) - dp(2) - dp(2) - dp(2)
    dq(3) = dq(1)*q(2) + q(1)*dq(2) - dp(3)
    dp(4) = dp(2)*q(1) + p(2)*dq(1)
    dp(5) = dp(1)*p(1) + p(1)*dp(1) - dp(4) - dp(4)
    dp(6) = dp(1)*q(2) + p(1)*dq(2) - dp(4)
    dq(4) = dq(1)*q(3) + q(1)*dq(3) - dp(6)
    dp(7) = dp(2)*p(1) + p(2)*dp(1)
    dp(8) = dp(2)*q(2) + p(2)*dq(2)
    dp(9) = dq(1)*p(5) + q(1)*dp(5) - dp(7)
    dp(10) = dp(1)*q(3) + p(1)*dq(3) - dp(8)
    dq(5) = dq(1)*q(4) + q(1)*dq(4) - dp(10)
    dp(11) = dp(2)*p(2) + p(2)*dp(2)
    dp(12) = dp(2)*p(3) + p(2)*dp(3)
    dp(13) = dp(2)*q(3) + p(2)*dq(3)
    dp(14) = dp(1)*p(5) + p(1)*dp(5) - dp(12)
    dp(15) = dq(1)*p(9) + q(1)*dp(9) - dp(12) - dp(14) - dp(14)
    dp(16) = dp(1)*q(4) + p(1)*dq(4) - dp(13)
    dq(6) = dq(1)*q(5) + q(1)*dq(5) - dp(16)
    dp(17) = dp(2)*p(4) + p(2)*dp(4)
    dp(18) = dp(2)*p(5) + p(2)*dp(5)
    dp(19) = dp(2)*p(6) + p(2)*dp(6)
    dp(20) = dp(2)*q(4) + p(2)*dq(4)
    dp(21) = dq(1)*p(14) + q(1)*dp(14) - dp(18)
    dp(22) = dq(3)*p(5) + q(3)*dp(5) - dp(17)
    dp(23) = dp(1)*q(5) + p(1)*dq(5) - dp(20)
    dp(24) = dp(2)*p(7) + p(2)*dp(7)
    dp(25) = dp(2)*p(8) + p(2)*dp(8)
    dp(26) = dp(2)*p(9) + p(2)*dp(9)
    dp(27) = dp(2)*p(10) + p(2)*dp(10)
    dp(28) = dp(2)*q(5) + p(2)*dq(5)
    dp(29) = dp(1)*p(14) + p(1)*dp(14) - dp(26)
    dp(30) = dq(2)*p(14) + q(2)*dp(14) - dp(24)
    dp(31) = dp(5)*q(4) + p(5)*dq(4) - dp(25)
    dp(32) = dp(1)*q(6) + p(1)*dq(6) - dp(28)

    return
  end subroutine devpoly3

  function drdx3(flag,xindex,xyz,r)
    integer i,j,flag,xindex,xyzind,matom,m
    real::xyz(3,3),r(3,3),drdx3

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

    drdx3 = 0.d0
    if (matom.eq.i.or.matom.eq.j) then
       drdx3=(xyz(i,xyzind)-xyz(j,xyzind))/r(i,j)
       if (matom.eq.j) then
          drdx3 = -drdx3
       endif
    endif

    return
  end function

  subroutine deriv_rev3(c,m,q,p,xyz,a,r,xxp)
    real::c(1:32),m(1:7),p(1:32)
    real::xyz(3,3),r(3,3),a
    !::::::::::::::::::::
    real::pp(1:32),mp(1:7),xxp(1:9)
    real::qp(1:6),q(1:6)

    qp(:)=0.d0
    pp(:)=0.d0
    mp(:)=0.d0
    xxp(:)=0.d0
    pp(32)=c(32)
    pp(31)=c(31)
    pp(30)=c(30)
    pp(29)=c(29)
    pp(28)=c(28)+pp(32)*(-1.d0)
    pp(27)=c(27)
    pp(26)=c(26)+pp(29)*(-1.d0)
    pp(25)=c(25)+pp(31)*(-1.d0)
    pp(24)=c(24)+pp(30)*(-1.d0)
    pp(23)=c(23)
    pp(22)=c(22)
    pp(21)=c(21)
    pp(20)=c(20)+pp(23)*(-1.d0)
    pp(19)=c(19)
    pp(18)=c(18)+pp(21)*(-1.d0)
    pp(17)=c(17)+pp(22)*(-1.d0)
    qp(6)=pp(32)*p(1)
    pp(16)=c(16)+qp(6)*(-1.d0)
    pp(15)=c(15)
    pp(14)=c(14)+pp(30)*q(2)+pp(29)*p(1)+pp(21)*q(1)+pp(15)*(-2.d0)
    pp(13)=c(13)+pp(16)*(-1.d0)
    pp(12)=c(12)+pp(15)*(-1.d0)+pp(14)*(-1.d0)
    pp(11)=c(11)
    qp(5)=pp(28)*p(2)+pp(23)*p(1)+qp(6)*q(1)
    pp(10)=c(10)+pp(27)*p(2)+qp(5)*(-1.d0)
    pp(9)=c(9)+pp(26)*p(2)+pp(15)*q(1)
    pp(8)=c(8)+pp(25)*p(2)+pp(10)*(-1.d0)
    pp(7)=c(7)+pp(24)*p(2)+pp(9)*(-1.d0)
    qp(4)=pp(31)*p(5)+pp(20)*p(2)+pp(16)*p(1)+qp(5)*q(1)
    pp(6)=c(6)+pp(19)*p(2)+qp(4)*(-1.d0)
    pp(5)=c(5)+pp(31)*q(4)+pp(22)*q(3)+pp(18)*p(2)+pp(14)*p(1)+pp(9)*q(1)
    pp(4)=c(4)+pp(17)*p(2)+pp(6)*(-1.d0)+pp(5)*(-2.d0)
    qp(3)=pp(22)*p(5)+pp(13)*p(2)+pp(10)*p(1)+qp(4)*q(1)
    pp(3)=c(3)+pp(12)*p(2)+qp(3)*(-1.d0)
    pp(2)=c(2)+pp(28)*q(5)+pp(27)*p(10)+pp(26)*p(9)+pp(25)*p(8)+pp(24)*p(7)+&
     pp(20)*q(4)+pp(19)*p(6)+pp(18)*p(5)+pp(17)*p(4)+pp(13)*q(3)+pp(12)*p(3)+&
    pp(11)*2*p(2)+pp(8)*q(2)+pp(7)*p(1)+pp(4)*q(1)+pp(3)*(-3.d0)
    qp(2)=pp(30)*p(14)+pp(8)*p(2)+pp(6)*p(1)+qp(3)*q(1)
    pp(1)=c(1)+pp(32)*q(6)+pp(29)*p(14)+pp(23)*q(5)+pp(16)*q(4)+pp(14)*p(5)+&
    pp(10)*q(3)+pp(7)*p(2)+pp(6)*q(2)+pp(5)*2*p(1)+pp(3)*q(1)+qp(2)*(-2.d0)
    qp(1)=pp(21)*p(14)+qp(6)*q(5)+pp(15)*p(9)+qp(5)*q(4)+pp(9)*p(5)+qp(4)*q(3)+&
    pp(4)*p(2)+qp(3)*q(2)+pp(3)*p(1)+qp(2)*2*q(1)
    mp(7)=pp(2)
    mp(6)=pp(1)+mp(7)*m(1)
    mp(5)=pp(1)
    mp(4)=pp(1)
    mp(3)=qp(1)+mp(6)*m(2)+mp(5)*m(1)
    mp(2)=qp(1)+mp(6)*m(3)+mp(4)*m(1)
    mp(1)=qp(1)+mp(7)*m(6)+mp(5)*m(3)+mp(4)*m(2)
    xxp(1)=mp(3)*(-m(3)/a)*drdx3(1,1,xyz,r)+mp(2)*(-m(2)/a)*drdx3(1,2,xyz,r)
    xxp(2)=mp(3)*(-m(3)/a)*drdx3(2,1,xyz,r)+mp(2)*(-m(2)/a)*drdx3(2,2,xyz,r)
    xxp(3)=mp(3)*(-m(3)/a)*drdx3(3,1,xyz,r)+mp(2)*(-m(2)/a)*drdx3(3,2,xyz,r)
    xxp(4)=mp(3)*(-m(3)/a)*drdx3(4,1,xyz,r)+mp(1)*(-m(1)/a)*drdx3(4,3,xyz,r)
    xxp(5)=mp(3)*(-m(3)/a)*drdx3(5,1,xyz,r)+mp(1)*(-m(1)/a)*drdx3(5,3,xyz,r)
    xxp(6)=mp(3)*(-m(3)/a)*drdx3(6,1,xyz,r)+mp(1)*(-m(1)/a)*drdx3(6,3,xyz,r)
    xxp(7)=mp(2)*(-m(2)/a)*drdx3(7,2,xyz,r)+mp(1)*(-m(1)/a)*drdx3(7,3,xyz,r)
    xxp(8)=mp(2)*(-m(2)/a)*drdx3(8,2,xyz,r)+mp(1)*(-m(1)/a)*drdx3(8,3,xyz,r)
    xxp(9)=mp(2)*(-m(2)/a)*drdx3(9,2,xyz,r)+mp(1)*(-m(1)/a)*drdx3(9,3,xyz,r)

    return
  end subroutine deriv_rev3

end module bemsa3
