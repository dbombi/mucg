!--------------------------------------------------------------------------    
    module parsubfuns
    implicit none
! Parameters
    integer, parameter :: sz = 1000 ! Size of calculation array
    real(8), parameter :: R = 8.3144621d+0 ! Gas cosnstant J/molK
    real(8), parameter :: mt=1793d+0 !
    real(8), parameter :: t0=273.15d0 !
    real(8), parameter :: w1=48570d+0 !
    real(8), parameter :: xa=1d-3 ! 0.001d0 
    real(8), parameter :: h=38575d+0
    real(8), parameter :: s=13.48d+0
    real(8), parameter :: strain=50d0
! Functions   
    contains
    subroutine composition(cin)
! Subroutine reads composition of 11 elements, calculates base (Fe) and prints to screen    
    implicit none
    real(8), intent(out), dimension(12) :: cin
!
    print 2 !C
    call rread(cin(2)) !C
    if (cin(2)<0.1D-03 .or. cin(2)>2.0D+00) call bound(cin(2),0.1D-03,2.0D+00) !C
    print 3 !Si
    call rread(cin(3)) !Si
    if (cin(3)<0.0D+00 .or. cin(3)>2.5D+00) call bound(cin(3),0.0D+00,2.5D+00) !Si
    print 4 !Mn
    call rread(cin(4)) !Mn
    if (cin(4)<0.0D+00 .or. cin(4)>5.5D+00) call bound(cin(4),0.0D+00,5.5D+00) !Mn
    print 5 !Ni
    call rread(cin(5)) !Ni
    if (cin(5)<0.0D+00 .or. cin(5)>4.5D+00) call bound(cin(5),0.0D+00,4.5D+00) !Ni
    print 6 !Mo
    call rread(cin(6)) !Mo
    if (cin(6)<0.0D+00 .or. cin(6)>1.5D+00) call bound(cin(6),0.0D+00,1.5D+00) !Mo
    print 7 !Cr
    call rread(cin(7)) !Cr
    if (cin(7)<0.0D+00 .or. cin(7)>3.5D+00) call bound(cin(7),0.0D+00,3.5D+00) !Cr
    print 8 !V
    call rread(cin(8)) !V
    if (cin(8)<0.0D+00 .or. cin(8)>1.5D+00) call bound(cin(8),0.0D+00,1.5D+00) !V
    print 9 !Co
    call rread(cin(9)) !Co
    if (cin(9)<0.0D+00 .or. cin(8)>4.0D+00) call bound(cin(9),0.0D+00,4.0D+00) !Co
    print 10 !Cu
    call rread(cin(10)) !Cu
    if (cin(10)<0.0D+00 .or. cin(10)>4.0D+00) call bound(cin(10),0.0D+00,4.0D+00) !Cu
    print 11 !Al
    call rread(cin(11)) !Al
    if (cin(11)<0.0D+00 .or. cin(11)>2.0D+00) call bound(cin(11),0.0D+00,2.0D+00) !Al
    print 12 !W
    call rread(cin(12)) !W
    if (cin(12)<0.0D+00 .or. cin(12)>4.0D+00) call bound(cin(12),0.0D+00,4.0D+00) !W
! calculate composition of base element Fe
    cin(1)=100d00 - sum(cin(2:12))
!
! Write compositional data to output    
    print 13
    print 14
    write(*,15) cin(2:7)
    write(*,16) cin(8:12),cin(1)
    print 13
!    
2   format (6x, 'Carbon, between 1e-4 and 2.0 wt % ?')
3   format (6x,' Silicon, between 0.0 and 2.5 wt % ?')
4   format (6x,' Manganese, between 0.0 and 5.5 wt % ?')
5   format (6x,' Nickel, between 0.0 and 4.5 wt % ?')
6   format (6x,' Molybdenum, between 0.0 and 1.5 wt % ?')
7   format (6x,' Chromium, between 0.0 and 3.5 wt % ?')
8   format (6x,' Vanadium, between 0.0 and 1.5 wt % ?')
9   format (6x,' Cobalt, between 0.0 and 4.0 wt % ?')
10  format (6x,' Copper, between 0.0 and 4.0 wt % ?')
11  format (6x,' Aluminium, between 0.0 and 2.0 wt % ?')
12  format (6x,' Tungsten, between 0.0 and 4.0 wt % ?')
13  format ('----------------------------------------------------------------------------------------------'/)
14  format (6x, 'Entered chemical composition in wt %')
15  format (8H      C=,f8.4, 6H   Si=,f8.4, 6H   Mn=,f8.4, 6H   Ni=,f8.4, 6H   Mo=,f8.4, 6H   Cr=,f8.4)
16  format (8H      V=,f8.4, 6H   Co=,f8.4, 6H   Cu=,f8.4, 6H   Al=,f8.4, 6H   W= ,f8.4, 6H   Fe=,f8.4/)
    end subroutine composition
!
!-----------------------------------------------------------------------
!
    subroutine mech(mechd,mechsta,gs1)
! Soubroutine reads, grain size and deformation properties
! Default grain size is 40 um
    implicit none
    integer :: stre, stra
    real(8), intent(out) :: mechd,mechsta,gs1
    real(8) :: gs, stress, prestr, s11, s33
!
    print 1
    print 2
    read *, gs
    print 3
    print 4
    read *, stre
    select case(stre)
    case (1) ! Under uniaxial tension
        print 5
        read *, stress
        mechd = 0.86d0*stress
    case (2) ! Under uniaxial compression
        print 6
        read *, stress
        mechd = 0.58d0*stress
    case default ! No stress
        stress=0d0
        mechd=0d0
    end select
    print 7
    print 8
    read *, stra
    select case(stra)
    case (1) ! Tensile pre-strain
        print 9
        read *, prestr
        mechsta=(8.0D10*2.52D-10/(8*3.141593*0.73)*((2.0D13+2.0D14*prestr)**0.5&
            &-(2.0D13)**0.5))/7.81D6*55.8
        s11=exp(prestr)
        s33=s11**(-0.5)
    case (2) ! Compressive  pre-strain
        print 10
        read *, prestr
        mechsta=(8.0D10*2.52D-10/(8*3.141593*0.73)*((2.0D13+2.0D14*prestr)**0.5&
            &-(2.0D13)**0.5))/7.81D6*55.8
        s11=exp(-prestr)
        s33=s11**(-0.5)
    case default ! No pre-strain
        prestr=0d0
        mechsta=0d0
        s11 = 1
        s33 = 1
    end select
    print 11
    if (gs==0) then
        gs=40d0
        write(*,12) stress,prestr
    else
        write(*,13) gs,stress,prestr
    endif
    print 11
    gs1=(3*S11*S33*(2*3**0.5+1))/(S11+3*(S11*(1+2*S11**2*S33**4)**0.5+(S11**2+2*S33**2)**0.5)&
        &+S33*(2*(1+S11**4*S33**2))**0.5)*GS
!
1   format (6x,'Prior austenite grain size in micro meter?')
2   format (6x,'For default value enter 0')  
3   format (6x, 'Stress conditions:')
4   format (6x, '1 - Under uniaxial tension, 2 - Under uniaxial compression, Other - No stress')
5   format (6x,'Uniaxial tension stress in MPa?')
6   format (6x,'Uniaxial compression stress in MPa?')
7   format (6x, 'Strain conditions:')
8   format (6x, '1 - Tensile pres-strain, 2 - Compressive pre-strain, Other - No pre-strain')
9   format (6x, 'Tensile pre-strain in mm/mm?')
10  format (6x, 'Compressive pre-strain in mm/mm?')
11  format ('----------------------------------------------------------------------------------------------'/)
12  format (6x,'Austenite GS= Default (40 micron), Applied stress=',f6.0,' MPa Pre-strain=',f6.2' mm/mm'/)
13  format (6x,'Austenite GS=',f6.0,' Applied stress=',f6.0,' MPa  Pre-strain=',f6.2' mm/mm'/)
    end subroutine mech
!
!--------------------------------------------------------------------------    
!

    
!
!-----------------------------------------------------------------------
!
    subroutine rread(A)
! Subroutine to read a real number in a way which to some extent traps typing errors
!
    implicit none
    real(8), intent(inout) :: A
    integer :: i=0
!
    do
        i=i+1
        read(*,*,err=11) A
        exit
11      if (i<3) then
            print 1
            cycle
        else
            print 2
            stop
        endif
    enddo
!    
1   format (6x,' Incorrect Input. Try again'/)
2   format (1H ,6x,'Too many attempts.  Program terminated.')
    end subroutine rread  
!
!-----------------------------------------------------------------------
!
    subroutine bound(C,L,H)
! Subroutine checks low and high limit and offers correction when value out of bounds
    implicit none
    integer :: i=0
    integer :: loop=1
    real(8), intent(inout) :: C 
    real(8), intent(in) :: L, H
    
    do 
        i=i+1
        print 1, L,H
        print 2 
        call rread(C)
        if (C>=L .and. C<=H) then 
            return
        else if (i<=loop) then
            cycle
        else
            print 3
            stop
        endif
    enddo
!
    return
1   format (12X,' Value out of bounds'/  14X,' The limits are ',F8.4,' to ', F8.3/)
2   format (1H ,'Input new value: ')
3   format ( 1H ,'Too many attempts. Program terminated.')
    end subroutine bound
!
!-----------------------------------------------------------------------      
! Function for logarithm activity of Fe in austenite
    function afeg(xeq,aj)
    implicit none
    real(8), intent(in) :: xeq, aj
    real(8) :: afeg
    real(8) :: deq
!
    deq=sqrt(1-2*(1+2*aj)*xeq+(1+8*aj)*(xeq**2))
    afeg=5d0*log((1-xeq)/(1-2*xeq))+6d0*log((1-2*AJ+(4*aj-1)*xeq-deq)/(2*aj*(2*xeq-1)))
    end function afeg
!
!-----------------------------------------------------------------------      
! Function for differential of logarithm activity of Fe in austenite
    function dafeg(xeq,aj)
    implicit none
    real(8),  intent(in) :: xeq, aj
    real(8) :: dafeg
    real(8) :: eteq, eteq2,deq
!
    deq=sqrt(1-2*(1+2*aj)*xeq+(1+8*aj)*(xeq**2))
    eteq=5d0*((1/(xeq-1))+2/(1-2*xeq))
    eteq2=6d0*((4*aj-1-(deq/2)*(-2-4*aj+2*xeq+16*xeq*aj))/(1-2*aj+(4*aj-1)*xeq-deq))+6*(4*aj/(2*aj*(2*xeq-1)))
    dafeg=eteq+eteq2
    end function dafeg
!
!-----------------------------------------------------------------------      
! Function for logarithm activity of C in austenite
    function cg(x,temp,W,R)
    implicit none
    real(8), intent(in) :: x, temp, W, R
    real(8) :: cg
    real(8) :: aj,dg,eg
!
    aj=1d0-exp(-W/(R*temp))
    if (x <= 1.0d-10)  then
        cg=log(1.0d-10)
    else
        dg=sqrt(1d0-2d0*(1d0+2d0*aj)*x+(1d0+8d0*aj)*(x**2))
        eg=5d0*log((1d0-2d0*x)/x)+6d0*W/(R*temp)+(38575d0-13.48*temp)/(R*temp)
        cg=eg+6d0*log((dg-1d0+3d0*x)/(dg+1d0-3d0*x))
    endif
    end function cg
!
!-----------------------------------------------------------------------      
! Function for differential (in respect to x) of logarithm activity of C in austenite
    function dcg(x,temp,W,R)
    implicit none
    real(8), intent(in) :: x, temp, W, R
    real(8) :: dcg
    real(8) :: aj,dg,ddg
!
    aj=1d0-exp(-W/(R*temp))
    dg=sqrt(1d0-2d0*(1d0+2d0*aj)*x+(1d0+8d0*aj)*(x**2))
    ddg=(dg/2)*(-2d0-4d0*aj+2d0*x+16d0*aj*x)
    dcg=-((10d0/(1d0-2d0*x))+(5d0/x))+6d0*((ddg+3d0)/(dg-1d0+3d0*x)-(ddg-3d0)/(dg+1d0-3d0*x))
    end function dcg
!
!-----------------------------------------------------------------------      
!
    function energy(temp,t10,t20,stored)
    implicit none
    real(8), intent(in) :: temp,t10,t20,stored
    real(8) :: energy
    real(8) :: f,t7,t8
!
    t7=temp-100d0*t20
    t8=t7-1140d0
    if (t7>=940) then
        f=t8*(-3.58434d-9*T8 + 2.70013d-6) - 1.04923d-3
        f=t8*(F*T8 + 0.26557D+0) - 8.88909d+0
    end if
    if (t7<300) f=1.38*T7-1499.0
    if (t7<700 .and. T7>=300)  f=1.65786*t7-1581d0
    if (t7<940 .and. T7>=700)  f=1.30089*t7-1331d0
    energy=(141d0*T10 + f)*4.187
! Decresase stability of austenite
    energy=energy-stored
    end function energy
!
!-----------------------------------------------------------------------
!
    function fto1(H,S,X,temp,W,W1,H1,S1,F,AJ,AJ1,R)
    implicit none
    real(8), intent(in) :: AJ,AJ1,F,H,H1,R,S,S1,temp,W,W1,X
    real(8) :: fto1
    real(8) :: D,D11,T60,zener,zen1,zen2
!
    d=sqrt(1-2*(1+2*AJ)*X+(1+8*AJ)*X*X)
    d11=sqrt(9-6*X*(2*AJ1+3)+(9+16*AJ1)*X*X)
    T60 = temp*(1-X)/(28080*X)
    if (t60>1d0) then
        zener=0.0
    else if (t60<0.25) then
        zen2=1d0
        zen1=3.295d0
        zener=4.187*((((zen2*X)**2)*(-50898.56)/(1-X))+zen1*temp*X*0.6623741)
    else
        zen1=0.2307+42.7974*T60-233.8631*(T60**2)+645.4485*(T60**3)&
            &-954.3995*(T60**4)+711.8095*(T60**5)-211.5136*(T60**6)
        zen2=-2.6702+45.6337*T60-225.3965*(T60**2)+567.7112*(T60**3)&
            &-771.6466*(T60**4)+538.1778*(T60**5)-151.3818*(T60**6)
        zener=4.187*((((zen2*X)**2)*(-50898.56)/(1-X))+zen1*temp*X*0.6623741)
    endif
    fto1=X*R*temp*log(X*X)+X*(H1-(H)-(S1-(S))*temp+4*W1-6*W)-R*temp*(1-X)*log((1-X)**4)&
        &+5*R*temp*(1-2*X)*log(1-2*X)-R*temp*X*log(((D-1+3*X)/(D+1-3*X))**6)&
        &-R*temp*(1-X)*log(((1-2*AJ+(4*AJ-1)*X-D)/(2*AJ*(2*X-1)))**6)&
        &+3*R*temp*X*log(3-4*X)+R*temp*X*log(((D11-3+5*X)/(D11+3-5*X))**4)+(1-X)*F+zener
    end function fto1
!
!-----------------------------------------------------------------------      
! Function WSFUN, Gn curve for Widmanstatten ferrite nucleation, 
! Acta Metall. 1981, J/mol, based on FPRO rather than GMAAX
    function wsfun(ktemp,t0)
    implicit none
    real(8), intent(in) :: ktemp,t0
    real(8) :: wsfun
!
    wsfun=3.251745029D+00*(ktemp-t0)-2183.014D+0
! When equal to or greater than zero, no problem with nucleation
    end function wsfun
!
!-----------------------------------------------------------------------      
    subroutine gmaax(A1,A,W1,F,R,temp,X,AFE,H1,S1)
    implicit none
    real(8), intent(out) :: a1, x
    real(8), intent(in) :: a, w1, f, r, temp, afe, h1, s1
    real(8) :: AJ1,D1,B1,B2,B3,A1FE,TEST,DA1,DA2,ERR,DA1FE
! CALCULATION OF THE OPTIMUM NUCLEUS C CONTENT AND ACTIVITY OF C IN
! FERRITE NUCLEUS, AVER=1.2
! Estimate a value for X, the optimum carbon concentration of the
! ferrite nucleus, in mole fraction
!
    x=6.3998D-07*temp-3.027D-04
    if (x <= 0.0D+00) x=1.0D-13
    if (x > 0.32D-03) x=0.32D-03
! Estimation complete
!
    aj1=1.0D+00-exp(-W1/(R*temp))
    do
        d1=sqrt(9.0D+00-6.0D+00*X*(2.0D+00*AJ1+3.0D+00)+(9.0D+00+16.0D+00*AJ1)*X**2)
        b1=((D1-3.0D+00+5.0D+00*X)/(D1+3.0D+00-5.0D+00*X))
        b2=((3.0D+00-4.0D+00*X)/X)
        b3=(H1-S1*temp + 4.0D+00*W1)/(R*temp)
        a1=b3+log(b1**4*b2**3)
!
        if (x > 1d-08) then
            a1fe=log(1e0-x)
        else
            a1fe=x
        endif
!
        test=F+R*temp*(A1FE-AFE) - R*temp*(A1-A)
!
! Newton iteration
!
        if (abs(TEST) > 10e0) then
            da1=(3.0D+00*X/(3.0D+00-4.0D+00*X))*((4.0D+00*X-3.0D+00)/(X**2)-4.0D+00/X)
            da2=(0.5D+00/D1)*(-12.0D+00*AJ1-18.0D+00+18.0D+00*X+32.0D+00*AJ1*X)
            da2=4.0D+00*(((DA2+5.0D+00)/(D1-3.0D+00+5.0D+00*X))-((DA2-5.0D+00)/(D1+3.0D+00-5.0D+00*X)))
            da1=da1+da2
            da1fe=1.0D+00/(X-1.0D+00)
            err=TEST/(R*temp*(DA1FE-DA1))
            if (err>x) err=0.3D+00*X
            x=abs(X-ERR)
        else
            exit
        endif
    enddo
! End of iteration
    end subroutine gmaax
!
!-----------------------------------------------------------------------      
!
    subroutine omega(c,w,xbar,t10,t20,j)
    implicit none
    integer, intent(in) :: j
    integer :: iu
    real(8), intent(inout), dimension(12) :: c
    real(8), intent(out) :: w, xbar, t10,  t20
    real(8), dimension(12) :: p, y
    real(8) :: b1, b2, b3
! Subroutine to calculate the carbon carbon interaction energy in austenite,
! as a function of alloy composition.  Based on .mucg18, in j/mol.
    b1=0d0
    b2=0d0
    b3=0d0
    if (j==1) then
        b1=sum(c)
    else
        c(1)=100d00 - sum(c(2:12))
        c(1)=c(1)/55.84D+00 !Fe
        c(2)=c(2)/12.0115D+00 !C
        c(3)=c(3)/28.09D+00 !Si
        c(4)=c(4)/54.94D+00 !Mn
        c(5)=c(5)/58.71D+00 !Ni
        c(6)=c(6)/95.94D+00 !Mo
        c(7)=c(7)/52.0D+00 !Cr
        c(8)=c(8)/50.94D+00 !V
        c(9)=c(9)/58.94D+00 !Co
        c(10)=c(10)/63.54D+00 !Cu
        c(11)=c(11)/26.98D+00 !Al
        c(12)=c(12)/183.85D+00 !W
        b1=sum(c)
    endif
    do iu=3,12
        y(iu)=c(iu)/c(1)
    enddo
    do iu=1,12
        c(iu)=c(iu)/b1
    enddo
    xbar=c(1)
    xbar=(int(10000.0D+00*xbar))/10000
    t10=y(3)*(-3)+y(4)*2+y(5)*12+y(6)*(-9)+y(7)*(-1)+y(8)*(-12)+y(9)*3.5+y(10)*7+y(11)*(-7)+y(12)*(-9)
    t20=-3*y(3)-37.5*y(4)-6*y(5)-26*y(6)-19*y(7)-44*y(8)+19.5*y(9)-4.5*y(10)+8*y(11)-26*y(12)
!
! Polynomials representing carbon-carbon interaction energy (J/mol) in
! austenite as a function of the molecular fraction of individual solutes.
!
    p(3)=-2.4233D+07+6.9547D+07*c(3)
    p(3)=+3.864D+06+p(3)*c(3)
    p(3)=+45802.87+(-280061.63+P(3)*C(3))*C(3) 
    p(3)=2013.0341+(763.8167+P(3)*C(3))*C(3)
    p(4)=2.0119D+06+(3.1716D+07-1.3885D+08*C(4))*C(4)
    p(4)=+6287.52+(-21647.96+P(4)*C(4))*C(4)
    p(4)=2012.067+(-1764.095+P(4)*C(4))*C(4)
    p(5)=-2.4968D+07+(1.8838D+08-5.5531D+08*C(5))*C(5)
    p(5)=-54915.32+(1.6216D+06+P(5)*C(5))*C(5)
    p(5)=2006.8017+(2330.2424+P(5)*C(5))*C(5)
    p(6)=-1.3306D+07+(8.411D+07-2.0826D+08*C(6))*C(6)
    p(6)=-37906.61+(1.0328D+06+P(6)*C(6))*C(6)
    p(6)=2006.834+(-2997.314+P(6)*C(6))*C(6)
    p(7)=+8.5676D+06+(-6.7482D+07+2.0837D+08*C(7))*C(7)
    p(7)=+33657.8+(-566827.83+P(7)*C(7))*C(7)
    p(7)=2012.367+(-9224.2655+P(7)*C(7))*C(7)
    p(8)=+5411.7566+(250118.1085-4.1676D+06*C(8))*C(8)
    p(8)=2011.9996+(-6247.9118+P(8)*C(8))*C(8)
    p(9)=(8427.00D0+5986D0*c(9))/4.187D0
    p(10)=2011.0D0
    p(11)=2011.0D0
    p(12)=2006.834D0-2997.314D0*c(12)-37906.61D0*c(12)**2+1.0328D+06*c(12)**3&
        &-1.3306D+07*c(12)**4+8.411D+07*c(12)**5-2.0826D+08*c(12)**6
    do iu=3,12
        b3=b3+p(iu)*y(iu)
        b2=b2+y(iu)
    enddo
    if (b2==0e0) then
        w=8054.0
    else
        w=(b3/b2)*4.187
    endif
    end subroutine omega
!
!-----------------------------------------------------------------------      
!
    subroutine tttt(temp,gmax,sheart,difft,R)
    implicit none
    real(8),  intent(in) :: temp, gmax, R
    real(8),  intent(out) :: sheart,difft
!
    sheart=exp((0.2432D+06/(R*temp))-0.135D+03+20.0*log(temp)-5*log(abs(GMAX)))
    difft=exp((0.6031D+06/(R*temp))-0.1905D+03+20.0*log(temp)-4.0D+00*log(abs(GMAX)))
    end subroutine tttt
!
!-----------------------------------------------------------------------      
! WSTINE
    subroutine wstine(j,dt4,dfpro,ws,ws1,sz,t0)
    implicit none
    integer, intent(in) :: j, sz
    integer :: i,k
    real(8), intent(in), dimension(sz) :: dt4, dfpro
    real(8), intent(in) :: t0
    real(8), intent(out) :: ws,ws1
    real(8) :: v14
    k=0 
    do i=1,j
        v14=wsfun((dt4(i)+t0),t0)
! Stored energy for growth of Wid ferrite = 50 J/mol
        if(dfpro(i)>-50.0D+00) exit
! WS1 is the lower intersection of Gn and Fpro curves
! WS is the upper intersection of Gn and Fpro curves
        if(dfpro(i)<=v14 .and. k==0) ws1=dt4(i)
        if(dfpro(i)>=v14) k=1
        if(dfpro(i)<=v14) ws=dt4(i)
    enddo
    end subroutine wstine
!
!-----------------------------------------------------------------------      
! To calculate the carbon concentration (mole fraction) at the T-zero phase
! boundary at a specified temperature.
    subroutine axto(H,S,H1,S1,XTO,temp,W,W1,F,AJ,AJ1,STORE,JN,J,R)
    implicit none
    integer, intent(in) :: jn
    integer, intent(out) :: j
    real(8), intent(in) :: h,s,h1,s1,temp,w,w1,f,aj,aj1,store,R
    real(8), intent(out) :: xto
    real(8) :: da,dfto,g9
!!    real(8), external :: fto1, g91 
! XTO is the carbon concentration at the T-zero boundary, T is the
! absolute temperature, STORE is the stored energy in Joules per mole
! H, S, H1, S1, W, W1, F, AJ, AJ1 are thermodynamic quantities defined in 
! subroutines OMEGA, ENERGY and main program
!
! Initialize J
    j=0
! Calculate the driving force for the transformation of  austenite to ferrite
! without any change in chemical composition, allowing for Zener ordering
    do
        j=j+1
        dfto=fto1(H,S,XTO,temp,W,W1,H1,S1,F,AJ,AJ1,R)+STORE
! Solve for XTO using Newton's iterative method
        da=abs(dfto)
        if (da<=10.0D+00) exit
!
! Obtain differential of FTO1 with respect to XTO
        g9=g91(xto,temp,W,W1,H1,S1,F,H,S,AJ,AJ1,R)
! Modify original guess of XTO, assuming convergence occurs within nine iterations
        if (j>jn) exit
        xto=xto-dfto/g9
        if(xto<=1.0e-12) then
            xto=0.0000D+00
            exit
        endif
    enddo
    end subroutine axto    
!
!-----------------------------------------------------------------------      
!
    function g91(xto,temp,W,W1,H1,S1,F,H,S,AJ,AJ1,R)
! To calculate the differential of FTO1 with respect to XTO. See subroutine 
! FTO1 for description of parameters. The differentiation is used in order 
! to apply the Newton Iteration.
    implicit none
    real(8), intent(in) :: xto, temp, w, w1, s, s1, f, h, h1, aj, aj1, R
    real(8) :: g91
    real(8) :: dt6,dzen1,dzen2,dzen3,dzen6,dzen7,dzen8,fd,fd1,g1
    real(8) :: v1,v2,v3,v4,v5,v6,v7,v8,v9
!
    fd=sqrt(1e0-2e0*(1e0+2e0*AJ)*xto+(1e0+8e0*AJ)*xto**2)
    fd1=sqrt(9e0-6e0*xto*(2e0*AJ1+3D0)+(9e0+16e0*AJ1)*xto**2)
    dt6=temp*(1e0-xto)/(28080e0*xto)
!
    if (dt6 > 1e0) then
        dzen3=0.0e0
    else if (dt6 < 1e0) then
        dzen2=1e0
        dzen1=3.295e0
        dzen3=4.187e0*((((dzen2*xto)**2)*(-50898.56e0)/(1e0-xto))+&
            &dzen1*temp*xto*(0.6623741e0))
    else
        dzen1=0.2307d0+42.7974d0*dt6-233.8631d0*(dt6**2)+645.4485d0*(dt6**3)&
            &-954.3995d0*(dt6**4)+711.8095d0*(dt6**5)-211.5136d0*(dt6**6)
        dzen2=-2.6702d0 +45.6337d0*dt6-225.3965d0*(dt6**2)+567.7112d0*(dt6**3)&
            &-771.6466d0*(dt6**4)+538.1778d0*(dt6**5)-151.3818d0*(dt6**6)
        dzen3=4.187e0*((((dzen2*xto)**2)*(-50898.56e0)/(1e0-xto))&
            &+dzen1*temp*xto*(0.6623741e0))
    endif
!
    v1=fd-1d0+3d0*xto
    v2=fd+1d0-3d0*xto
    v3=1d0-2d0*aj+(4d0*aj-1d0)*xto-fd
    v4=2d0*aj*(2d0*xto-1d0)
    v5=fd1-3d0+5d0*xto
    v6=fd1+3d0-5d0*xto
    v7=(xto-1d0-2d0*aj+8d0*xto*aj)/fd
    v8=(9d0*xto-9d0-6d0*aj1+16d0*aj1*xto)/fd1
    v9=h1-(h)-(s1-(s))*temp-6d0*w+4d0*w1
    g1=2e0+log(xto**2)+4d0+log((1e0-xto)**4)-10e0-log((1d0-2d0*xto)**10)-log((v1/v2)**6)&
        &-6d0*((xto/v1)*(v7+3d0+v1*(3d0-v7)/v2))+log((v3/v4)**6)&
        &-6d0*((1d0-xto)/v3)*(4d0*aj*(1d0-(v3/v4))-v7)+3d0*(log(3d0-4d0*xto)&
        &-4d0*xto/(3d0-4d0*xto))+log((v5/v6)**4) + 4d0*(xto/v5)*(v8+5d0+(v5/v6)*(5d0-v8))
!
    if (dt6>1e0) then
        dzen8=0.0e0
    else
        dzen6=(-3.3948d0+13.6112d0*dt6-13.4376d0*(dt6**2))*temp/(28080d0*((1d0-xto)**2))
        dzen7=(-3.3118d0+15.7462d0*dt6-23.2449d0*(dt6**2))*temp/(28080d0*((1d0-xto)**2))
        dzen8=50898d0*((dzen2*xto)**2)/((1d0-xto)**2)+(-50898d0*(2e0*dzen2*dzen6*(xto**2)&
            &+2d0*xto*(dzen2**2))/(1d0-xto))+dzen1*temp*0.6623741d0+dzen7*temp*xto*0.6623741d0

    endif
    g91=v9+R*temp*g1+dzen8-F
    end function g91
!
!-----------------------------------------------------------------------      
!
    subroutine analy(J8,J9,const,slope,corr,X,Y,sz)
    implicit none
    integer, intent(in) :: j8,j9,sz
    integer :: i,i1
    real(8), intent(in), dimension(sz) :: x,y
    real(8), intent(out) :: const, slope, corr
    real(8) :: ax,ax2,ay,ay2,axy
!
    ax=0e0
    ay=0e0
    ax2=0e0
    ay2=0e0
    axy=0e0
    i1=j8-j9
    do i=1, i1
        ax=ax+x(i)
        ay=ay+y(i)
        axy=axy+x(i)*y(i)
        ax2=ax2+x(i)*x(i)
        ay2=ay2+y(i)*y(i)
    enddo
    const=(ay*ax2-ax*axy)/(i1*ax2-ax*ax)
    slope=((i1*axy-ax*ay)/(i1*ax2-ax*ax))
    corr=(i1*axy-ax*ay)/(sqrt((i1*ax2-ax*ax)*(i1*ay2-ay*ay)))
    end subroutine analy
!
!-----------------------------------------------------------------------      
! Caluclate MS and BS temperature
    subroutine msbstemp(j,x,y,x1,mechd,mechsta,sz,ms,bs)
    implicit none
    integer, intent(in) :: j, sz
    integer :: i,k,l
    real(8), intent(in), dimension(sz) :: x,y
    real(8), intent(in) :: x1, mechd, mechsta
    real(8), intent(out) :: ms, bs
    real(8) :: bainstr, marstr
    real(8) :: const, slope, corr
    logical :: bain, mart
    bain=.false.
    mart=.false.
    const=1.0e+00
    corr=1.0e+00
    slope=1.0e+00
    bainstr=-400e0
    marstr=-1120d0-10568d0*x1+94.1d0+mechd-mechsta
!!    marstr=-1086.9d0+2941.8d0*x1-9.6016d5*x1**2+2.5382d7*x1**3-1.8488d8*x1**4+mechd-mechsta
!
    do i=1,j
        if(y(i)<marstr) then
            k=i
            mart=.true.
        endif
        if(y(i)<bainstr) then
            l=i
            bain=.true.
        endif
    enddo
    if ((bain .neqv. .true.) .or. (mart .neqv. .true.)) then
        call analy(j,(int(j/2)),const,slope,corr,x,y,sz)
    endif
    if (bain .eqv. .true.) then
        bs=x(l-1)+(x(l)-x(l-1))*((y(l-1)-bainstr)/(y(l-1)-y(l)))
    else
        bs=(bainstr-const)/slope
    endif
    if (mart .eqv. .true.) then
        ms=x(k-1)+(x(k)-x(k-1))*((y(k-1)-marstr)/(y(k-1)-y(k)))
    else
        ms=(marstr-const)/slope
        !if (ms<-273.15) ms=-273.15
    endif
    end subroutine msbstemp
!
!-----------------------------------------------------------------------      
! Main subroutine
    subroutine mucg(c,ws,bs,ms,stored,mechd,mechsta,gs1,step,prn,stctemp,stxeq,stxeq2)
    implicit none
! Variables
    integer, intent(in) :: step
    integer :: i, j, j98, j99, jn
    logical, intent(in) :: prn
    logical :: dat
    real(8), intent(inout), dimension(12) :: c
    real(8), intent(out), dimension(sz) :: stctemp, stxeq, stxeq2
    real(8), dimension(sz) :: bdif, bsh, btem, dfpro, ddfto, dt4, dxq
    real(8), intent(in) :: mechd, mechsta, gs1, stored
    real(8), intent(out) :: ws, bs, ms
    real(8) :: xto, xto400, x44, x1, x, xeq, xeq2
    real(8) :: temp, ctemp
    real(8) :: h1, s1, t10, t20,  v14, w, ws1, f, f44
    real(8) :: a, a1, a44, aeq, afe, afeq, afe44, aj, aj1
    real(8) :: ms0, da44, dafe44, df441, difft
    real(8) :: eteq, eteq2, fpro, fproa, fson, fto, gmax
    real(8) :: sheart, teq, teq2, xms, xbar
! Inicialization of variables
    dat = .false. !Don't write to file
!Write results to file
    if (dat) then
! Open file for TTT diagram
    open(unit=2,file='ttt.txt')
!    write(2,*) "CTEMP, SHEART, DIFFT"
    endif
!
    call omega(c,w,xbar,t10,t20,0)
    x1=c(2)
    if (prn) then
        print 10
        write(*,11) c(2:7)
        write(*,12) c(8:12),c(1)
        write(*,13) x1,t10,t20,w
        print 14
        print 15
    endif
    j=0
    fto=-1e0
    xeq=0.2d+00 !Paraequilibrium carbon concentration of austenite in mole fraction
    xeq2=0.15d+00 !Same as xeq but with 50 J/mol  of strain energy in the ferrite in order to do calculations for Widmanstatten ferrite
    xto=0.07d0 !T-zero carbon concentration in mole fraction
    xto400=0.06d0 !Same as xto, but with 400 J/mol of stored energy in the ferrite, to allow calculations on bainite
    x44=0.1d+00
    jn=1000
!
    do i=473,1173,step !loop 27 (Main loop)
        j=j+1
        xto=xto+1d-4
        xto400=xto400+1d-4
        x44=0.3*xeq
        temp=dble(i)
        ctemp=temp-dble(int(t0))
        if (temp <= 1000d0) then
            h1=111918.0d+00
            s1=51.44d+00
        else
            h1=105525.0d+00
            s1=45.34521d+00
        endif
        f=energy(temp,t10,t20,stored)
        aj=1d0-exp(-w/(r*temp))
        aj1=1d+0-exp(-w1/(r*temp))
        do !8
            teq=R*temp*afeg(xeq,aj)-f
            if (abs(teq)<1d0) exit
            eteq=dafeg(xeq,aj)*R*temp
            xeq=xeq-teq/eteq
        enddo
        do !10
            teq2=R*temp*afeg(xeq2,aj)-F-strain
            if (abs(teq2)<1e0) exit
            eteq2=dafeg(xeq2,aj)*R*temp
            xeq2=xeq2-teq2/eteq2
        enddo
        aeq=cg(xeq,temp,w,R)
        afeq=afeg(xeq,aj)
        if (xeq < (x1+0.0004)) exit
        a=cg(x1,temp,w,R)
        afe=afeg(x1,aj)
        fson = R*temp*(xa*(aeq-a)+(1d0-xa)*(afeq-afe))
        fpro = R*temp*(x1*(aeq-a)+(1d0-x1)*(afeq-afe))
! V14 is the free energy needed to nucleate Widmanstatten ferrite.
! If V14=0, goto 18.
        v14=wsfun(temp,t0)
        if (v14>=0d0) then
            x44 = 99999999999.99999
        else
            do !14
                a44=cg(x44,temp,W,R)
                afe44=afeg(x44,aj)
                f44=R*temp*(x44*(aeq-a44)+(1d0-x44)*(afeq-afe44))
                f44=f44-v14
                if (abs(f44)>=10d0) then
                    da44=dcg(x44,temp,W,R)
                    dafe44=dafeg(x44,aj)
                    df441=R*temp*(aeq-x44*da44-afeq+x44*dafe44-dafe44-a44+afe44)
                    x44=x44-f44/df441
                else
                    exit
                endif
            enddo
        endif
        fproa=fpro*((xeq-xa)/(xeq-x1))
        call gmaax(a1,a,W1,F,R,temp,x,afe,h1,s1)
        gmax=R*temp*(a1-a)
        if (fto>=0d0) then
            fto=0d0
        else
            fto=fto1(H,S,X1,temp,W,W1,H1,S1,F,AJ,AJ1,R)
        endif
        call axto(H,S,H1,S1,XTO,temp,W,W1,F,AJ,AJ1,0.0D+00,JN,J98,R)
        call axto(H,S,H1,S1,XTO400,temp,W,W1,F,AJ,AJ1,400.0D+00,JN,J99,R)
        call tttt(temp,GMAX,SHEART,DIFFT,R)
        if (x44==0d0) sheart=1d+20
        if (prn) then
            write(*,16) FPRO,FPROA,GMAX,CTEMP,X,FSON,XEQ,XEQ2,FTO,XTO,X44,XTO400,SHEART,DIFFT
        endif
        if (dat) then
            write(2,16) FPRO,FPROA,GMAX,CTEMP,X,FSON,XEQ,XEQ2,FTO,XTO,X44,XTO400,SHEART,DIFFT
        endif
!
! Output for structure
        stxeq(j)=xeq
        stxeq2(j)=xeq2
        stctemp(j)=ctemp      
! Plot 
!        write(2,"(f8.2,2f10.4)") ctemp, log(SHEART), log(DIFFT)
        bsh(j)=sheart
        bdif(j)=difft
        btem(j)=ctemp
!       WRITE(*,1090)j,BSH(j),BDIF(j),BTEM(j)
        dxq(j)=xeq
        ddfto(j)=fto
        dt4(j)=ctemp
        dfpro(j)=fpro     
    enddo !loop 27
    j=j-1 !j8
    if (prn) then
        print 17
        print 14
        print 18
    endif
    call wstine(j,dt4,dfpro,ws,ws1,sz,t0)
!
! Modification made February 1991 to allow for large carbon concentrations
! Based on Metal Science paper on thermodynamics of martensitic transformations
! in plain carbon steels
    xms=x1
    if (xms>0.0594d0) xms=0.0594d0
    call msbstemp(j,dt4,ddfto,xms,mechd,mechsta,sz,ms,bs)
    ms0=1d0/0.2689*log(1/(4/3*3.141593*(40*1.0D-3/2)**3)*(exp(-log(0.99)/0.05)-1)+1)+ms
    ms=ms0-1d0/0.2689*log(1/(4/3*3.141593*(gs1*1.0D-3/2)**3)*(exp(-log(0.99)/0.05)-1)+1)
!   
    if (prn) write(*,19) ws1,ws
    if (ws<bs) then
        bs=ws
        if (prn) write(*,20) bs
    else
        if (prn) write(*,21) bs
    endif
    if (prn) then
    write(*,22) ms
    print 14
    endif
!
    if (dat) then
        write(2,*) "ws1, ws, bs, ms"
        write(2,"(4f8.2)") ws1, ws, bs, ms
        close(2)
    endif
!
! Format statements
10  format (6x, 'Entered chemical composition in at. %')
11  format (8H      C=,f8.4, 6H   Si=,f8.4, 6H   Mn=,f8.4, 6H   Ni=,f8.4, 6H   Mo=,f8.4, 6H   Cr=,f8.4)
12  format (8H      V=,f8.4, 6H   Co=,f8.4, 6H   Cu=,f8.4, 6H   Al=,f8.4, 6H   W= ,f8.4, 6H   Fe=,f8.4/)
13   format (21H      Carbon content=,f10.5,6H  T10=,f10.6,6H  T20=,f10.6,10H   WGamma=,f7.0/)
14  format ('----------------------------------------------------------------------------------------------')
15  format ('   FPRO   FPROA    GMAX   CTEMP',12H   X NUCLEUS,6H  FSON,9H      XEQ,10H    XEQ50 ,&
     &6H   FTO,9H      XTO,11H      X44  ,8H XTO400 ,13H    SHEART   ,10H  DIFFT   )
16 format (1H ,F8.1,2F8.1,F7.1,D10.2,F8.1,2F9.4,F8.1,F8.5,1H ,2F8.5,1H ,E11.2,E11.2)
17 format (4H    )
18 format ('    ***** FTO versus TEMPERATURE *****   '/)
19 format (' Widmanstatten start temperature range ',F7.1,' C & ',F7.1,' C')
20 format (' Nucleation limited bainite start temperature=',F7.1,' C')
21 format (' Growth limited bainite start temperature=',F7.1,' C')
22 format (' Martensite start temperature=',F7.1,' C')
    end subroutine mucg
!
    end module parsubfuns
!--------------------------------------------------------------------------
