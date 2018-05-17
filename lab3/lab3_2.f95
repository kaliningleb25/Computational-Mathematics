program Lab3
    integer, parameter :: NEQN = 2
    REAL :: T=0.0,Y(NEQN),TOUT =0.0,RELERR=0.1E-07,ABSERR=0.00
    REAL :: TFINAL=0.1,TPRINT=0.05,WORK(27)
    INTEGER :: IWORK(5),IFLAG

    Y(1)=-3
    Y(2)=1
    IFLAG=1
    
    write(*,*) "=========================================================================="
    IFLAG=1
    write(*,*) "RKF45 с шагом 0,05"
    write(*,*)
    
    write(*, *) "      T            OutVal(1)         OutVal(2)"
    do i = 1, 21
      call RKF45(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK,IWORK)
      write(*, "(f10.4, f16.4, f16.4)") T, Y(1), Y(2)
      tout = tout + tprint
    end do
  
    write(*,*) "=========================================================================="
    write(*,*) "Метод Эйлера-Коши с шагом 0,05"
    write(*,*)
    call Euler(0.05)

    IFLAG=1

    write(*,*) "=========================================================================="
    write(*,*) "Метод Эйлера-Коши с шагом 0,0001"
    write(*,*)
    call Euler(0.0001)
    
contains

    subroutine Euler(h)
    REAL    :: K1(NEQN),K2(NEQN),YN(NEQN),Ytmp(NEQN),Y(NEQN)
    REAL    :: t=0.0, finish = 1.1, i
    integer :: counter=0
      T=0.0
      i = -0.00001
      Y(1) = -3.0
      Y(2) = 1.0

      write(*, *) "      T            OutVal(1)         OutVal(2)"

      do while (T <= finish)
157     Ytmp=Y
        call F(t, Ytmp, YN)
        K1 = h*YN
      
        Ytmp = Y + K1
        call F(t+h,Ytmp,YN)
        K2= h*YN
      
        Ytmp = Y + (K2)
        YN = Y +1.0/2.0*(K1+K2)

        Y =YN

        if (T > i) then
          write(*, "(f10.4, f16.4, f16.4)") t, Y(1), Y(2)
          i = i+ 0.1
        end if 
        
        T=h+T
        counter = counter +1
      end do
      
end subroutine Euler

subroutine F(t, Y, YP)
    Real t, Y(NEQN), YP(NEQN)

    YP(1) = -73*Y(1)-210*Y(2)+log(t*t+1)
    YP(2) = Y(1)+exp(-t)+t*t+1

end subroutine F


 SUBROUTINE RKF45(F,NEQN,Y,T,TOUT,RELERR,ABSERR, &
                      & IFLAG,WORK,IWORK)

      EXTERNAL F
      INTEGER NEQN,IFLAG,IWORK(5)
      REAL Y(NEQN),T,TOUT,RELERR,ABSERR,WORK(1)

      REAL K1,K2,K3,K4,K5,K6,K1M

      K1M=NEQN+1
      K1=K1M+1
      K2=K1+NEQN
      K3=K2+NEQN
      K4=K3+NEQN
      K5=K4+NEQN
      K6=K5+NEQN

      CALL RKFS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,&
              & WORK(1),WORK(K1M),WORK(K1),WORK(K2),&
              & WORK(K3),WORK(K4),WORK(K5),WORK(K6),&
              & WORK(K6+1),IWORK(1),IWORK(2),&
              & IWORK(3),IWORK(4),IWORK(5))
      RETURN
      END


      SUBROUTINE RKFS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,&
     &                YP,H,F1,F2,F3,F4,F5,SAVRE,SAVAE,&
     &                NFE,KOP,INIT,JFLAG,KFLAG)
      LOGICAL HFAILD,OUTPUT
      INTEGER NEQN,IFLAG,NFE,KOP,INIT,JFLAG,KFLAG
      REAL Y(NEQN),T,TOUT,RELERR,ABSERR,H,YP(NEQN),&
     &    F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN),&
     &     SAVRE,SAVAE
      EXTERNAL F
      REAL A,AE,DT,EE,EEOET,ESTTOL,ET,HMIN,REMIN,&
     &     RER,S,SCALE,TOL,TOLN,U26,EPSP1,EPS,YPK
      save U26, EPS
      INTEGER K,MAXNFE,MFLAG
      REAL AMAX1,AMIN1
      DATA REMIN/1.E-12/
      DATA MAXNFE/3000/
      IF(NEQN.LT.1)GO TO 10
      IF((RELERR.LT.0.0).OR.(ABSERR.LT.0.0))GO TO 10
      MFLAG=IABS(IFLAG)
      IF((MFLAG.EQ.0).OR.(MFLAG.GT.8))GO TO 10
      IF(MFLAG.NE.1)GO TO 20
      EPS=1.0
    5 EPS=EPS/2.0
      EPSP1=EPS+1.
      IF(EPSP1.GT.1.)GO TO 5
      U26=26.*EPS
      GO TO 50
   10 IFLAG=8
      RETURN
   20 IF((T.EQ.TOUT).AND.(KFLAG.NE.3))GO TO 10
      IF(MFLAG.NE.2)GO TO 25
      IF((KFLAG.EQ.3).OR.(INIT.EQ.0))GO TO 45
      IF(KFLAG.EQ.4)GO TO 40
      IF((KFLAG.EQ.5).AND.(ABSERR.EQ.0.0))GO TO 30
      IF((KFLAG.EQ.6).AND.(RELERR.LE.SAVRE).AND.&
     &(ABSERR.LE.SAVAE))GO TO 30
      GO TO 50
   25 IF(IFLAG.EQ.3)GO TO 45
      IF(IFLAG.EQ.4)GO TO 40
      IF((IFLAG.EQ.5).AND.(ABSERR.GT.0.0))GO TO 45
   30 PRINT 35
   35 FORMAT(/20X)
      STOP
   40 NFE=0
      IF(MFLAG.EQ.2)GO TO 50
   45 IFLAG=JFLAG
      IF(KFLAG.EQ.3)MFLAG=IABS(IFLAG)
   50 JFLAG=IFLAG
      KFLAG=0
      SAVRE=RELERR
      SAVAE=ABSERR
      RER=2.*EPS+REMIN
      IF(RELERR.GE.RER)GO TO 55
      RELERR=RER
      IFLAG=3
      KFLAG=3
      RETURN
   55 DT=TOUT-T
      IF(MFLAG.EQ.1)GO TO 60
      IF(INIT.EQ.0)GO TO 65
      GO TO 80
   60 INIT=0
      KOP=0
      A=T
      CALL F(A,Y,YP)
      NFE=1
      IF(T.NE.TOUT)GO TO 65
      IFLAG=2
      RETURN
   65 INIT=1
      H=ABS(DT)
      TOLN=0.
      DO 70 K=1,NEQN
      TOL=RELERR*ABS(Y(K))+ABSERR
      IF(TOL.LE.0)GO TO 70
      TOLN=TOL
      YPK=ABS(YP(K))
      IF(YPK*H**5.GT.TOL)H=(TOL/YPK)**0.2
   70 CONTINUE
      IF(TOLN.LE.0.0)H=0.0
      H=AMAX1(H,U26*AMAX1(ABS(T),ABS(DT)))
      JFLAG=ISIGN(2,IFLAG)
   80 H=SIGN(H,DT)
      IF(ABS(H).GE.2.0*ABS(DT))KOP=KOP+1
      IF(KOP.NE.100)GO TO 85
      KOP=0
      IFLAG=7
      RETURN
   85 IF(ABS(DT).GT.U26*ABS(T))GO TO 95
      DO 90 K=1,NEQN
   90 Y(K)=Y(K)+DT*YP(K)
      A=TOUT
      CALL F(A,Y,YP)
      NFE=NFE+1
      GO TO 300
   95 OUTPUT=.FALSE.
      SCALE=2./RELERR
      AE=SCALE*ABSERR
  100 HFAILD=.FALSE.
      HMIN=U26*ABS(T)
      DT=TOUT-T
      IF(ABS(DT).GE.2.*ABS(H))GO TO 200
      IF(ABS(DT).GT.ABS(H))GO TO 150
      OUTPUT=.TRUE.
      H=DT
      GO TO 200
  150 H=0.5*DT
  200 IF(NFE.LE.MAXNFE)GO TO 220
      IFLAG=4
      KFLAG=4
      RETURN
  220 CALL FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,F1)
      NFE=NFE+5
      EEOET=0.
      DO 250 K=1,NEQN
      ET=ABS(Y(K))+ABS(F1(K))+AE
      IF(ET.GT.0.)GO TO 240
      IFLAG=5
      KFLAG=5
      RETURN
  240 EE=ABS((-2090.*YP(K)+(21970.*F3(K)-15048.*F4(K)))&
     &  +(22528.*F2(K)-27360.*F5(K)))
  250 EEOET=AMAX1(EEOET,EE/ET)
      ESTTOL=ABS(H)*EEOET*SCALE/752400.
      IF(ESTTOL.LE.1.0)GO TO 260
      HFAILD=.TRUE.
      OUTPUT=.FALSE.
      S=0.1
      IF(ESTTOL.LT.59049.)S=0.9/ESTTOL**0.2
      H=S*H
      IF(ABS(H).GT.HMIN)GO TO 200
      IFLAG=6
      KFLAG=6
      RETURN
  260 T=T+H
      DO 270 K=1,NEQN
  270 Y(K)=F1(K)
      A=T
      CALL F(A,Y,YP)
      NFE=NFE+1
      S=5.
      IF(ESTTOL.GT.1.889568E-4)S=0.9/ESTTOL**0.2
      IF(HFAILD)S=AMIN1(S,1.0)
      H=SIGN(AMAX1(S*ABS(H),HMIN),H)
      IF(OUTPUT)GO TO 300
      IF(IFLAG.GT.0)GO TO 100
      IFLAG=-2
      RETURN
  300 T=TOUT
      IFLAG=2
      RETURN
      END

SUBROUTINE FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,S)
      INTEGER NEQN
      REAL Y(NEQN),YP(NEQN),F1(NEQN),F2(NEQN),&
      &    F3(NEQN),F4(NEQN),F5(NEQN),S(NEQN)
      REAL CH, T, H
      INTEGER K
      CH=H/4.0
      DO 221 K=1,NEQN
  221 F5(K)=Y(K)+CH*YP(K)
      CALL F(T+CH,F5,F1)
      CH=3.0*H/32.0
      DO 222 K=1,NEQN
  222 F5(K)=Y(K)+CH*(YP(K)+3.0*F1(K))
      CALL F(T+3.0*H/8.0,F5,F2)
      CH=H/2197.0
      DO 223 K=1,NEQN
  223 F5(K)=Y(K)+CH*(1932.0*YP(K)+(7296.0*F2(K)-7200.0*F1(K)))
      CALL F(T+12.0*H/13.0,F5,F3)
      CH=H/4104.0
      DO 224 K=1,NEQN
  224 F5(K)=Y(K)+CH*((8341.0*YP(K)-845.0*F3(K))+ &
      &     (29440.0*F2(K)-32832.0*F1(K)))
      CALL F(T+H,F5,F4)
      CH=H/20520.0
      DO 225 K=1,NEQN
  225 F1(K)=Y(K)+CH*((-6080.0*YP(K)+(9295.0*F3(K)-&
     &      5643.0*F4(K)))+(41040.0*F1(K)-28352.0*F2(K)))
      CALL F(T+H/2.0,F1,F5)
      CH=H/7618050.0
      DO 230 K=1,NEQN
  230 S(K)=Y(K)+CH*((902880.0*YP(K)+(3855735.0*F3(K)-&
     &     1371249.0*F4(K)))+(3953664.0*F2(K)+277020.0*F5(K)))
      RETURN
      END
end program Lab3