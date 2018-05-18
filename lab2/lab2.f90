program lab2
  real, allocatable    :: matrix(:,:), tmp_matrix(:,:), identity_matrix(:,:), reverse_matrix(:,:), R(:,:), col(:), work(:)
  integer, allocatable :: IPVT(:)
  real, dimension(4)   :: array_with_values
  integer              :: N
  real                 :: cond, norm, tmpSum, var

  N = 5

  var = 1.0001

  array_with_values = (/4.0, 3.0, 2.0, var /)

  allocate(matrix(N,N))
  allocate(tmp_matrix(N,N))
  allocate(identity_matrix(N,N))
  allocate(reverse_matrix(N,N))
  allocate(R(N,N))
  allocate(col(N))
  allocate(IPVT(N))
  allocate(work(N))
  
  matrix = make_matrix(5, array_with_values)

  identity_matrix = make_identity_matrix(n)

  tmp_matrix = matrix

  CALL DECOMP(N, N, matrix, cond, IPVT, work)
  do j = 1, N
    col(:) = identity_matrix(:, j)
    CALL SOLVE(N, N, matrix, col, IPVT)
    reverse_matrix(:, j) = col(:)
  end do

  R = identity_matrix - matmul(reverse_matrix, tmp_matrix)
  norm = 0
  do j = 1, N
    tmpSum = sum(R(j,:))
    if (tmpSum > norm) then
      norm = tmpSum
    end if
  end do

  print*, "При a4 = "
  write(*,"(f10.4)") var
  print*

  print*, "Исходная матрица"
  do i = 1,n
    do j = 1,n
      write(*,"(f12.4,$)") matrix(i,j)
    end do
    write(*,*)
  end do
  print*

  print*, "Единичная матрица"
  do i = 1,n
    do j = 1,n
      write(*,"(f12.4,$)") identity_matrix(i,j)
    end do
    write(*,*)
  end do
  print*

  print*, "Обратная матрица"
  do i = 1,n
    do j = 1,n
      write(*,"(f12.4,$)") reverse_matrix(i,j)
    end do
    write(*,*)
  end do
  print*

  print*, "Матрица R"
  do i = 1,n
    do j = 1,n
      write(*,"(f12.4,$)") R(i,j)
    end do
    write(*,*)
  end do
  print*

  print*, "Число обусловленности", cond
  print*

  print*, "Норма матрицы", norm

contains
! Процедура формирования матрицы A по заданному вектору B
function make_matrix(n, array_with_values) result(new_matrix)
  integer, intent(in)             :: n
  real,dimension(n-1), intent(in) :: array_with_values

  real, dimension(n,n) :: new_matrix

  do i = 1,n
    do j = 1,n
      if (j .le. i) then
        new_matrix(i,j) = 1
      else
        new_matrix(i,j) = array_with_values(i)
      end if
    end do
  end do

  RETURN
end function make_matrix


! Процедура формирования единичной матрицы
function make_identity_matrix(n) result(identity_matrix)
  integer, intent(in) :: n

  real, dimension(n,n) :: identity_matrix

  do i = 1, N
    identity_matrix(i,i) = 1
  end do

  RETURN

end function make_identity_matrix
end program lab2

      SUBROUTINE DECOMP(NDIM,N,A,COND,IPVT,WORK)
!
      INTEGER NDIM,N
      REAL A(NDIM,N),COND,WORK(N)
      INTEGER IPVT(N)
!
!     ПPOГPAMMA BЫЧИCЛЯET PAЗЛOЖEHИE BEЩECTBEHHOЙ MATPИЦЫ
!     ПOCPEДCTBOM ГAУCCOBA ИCKЛЮЧEHИЯ И OЦEHИBAET
!     OБУCЛOBЛEHHOCTЬ  MATPИЦЫ.
!
!     OHA ИCПOЛЬЗУETCЯ ДЛЯ BЫЧИCЛEHИЯ PEШEHИЙ
!     ЛИHEЙHЫX CИCTEM.
!
!     BXOДHAЯ ИHФOPMAЦИЯ.
!
!     NDIM -ЗAЯBЛEHHAЯ CTPOЧHAЯ PAЗMEPHOCTЬ MACCИBA,
!           COДEPЖAЩEГO A.
!
!     N    -ПOPЯДOK MATPИЦЫ.
!
!     A    -MATPИЦA,KOTOPУЮ HУЖHO PAЗЛOЖИTЬ.
!
!     BЫXOДHAЯ ИHФOPMAЦИЯ.
!
!     A     COДEPЖИT BEPXHЮЮ TPEУГOЛЬHУЮ MATPИЦУ U
!           И УЧИTЫBAЮЩУЮ ПEPECTAHOBKИ BEPCИЮ
!           HИЖHEЙ TPEУГOЛЬHOЙ MATPИЦЫ I-L,TAKИE,
!           ЧTO (MATPИЦA  ПEPECTAHOBOK) *A=L*U
!
!     COND -OЦEHKA OБУCЛOBЛEHHOCTИ A.
!           ДЛЯ ЛИHEЙHOЙ CИCTEMЫ A*X=B ИЗMEHEHИЯ B A И B
!           MOГУT BЫЗBATЬ  ИЗMEHEHИЯ B X,БOЛЬШИE B COND PAЗ.
!           ECЛИ COND+1.0.EQ.COND, TO  A B ПPEДEЛAX MAШИHHOЙ
!           TOЧHOCTИ ЯBЛЯETCЯ BЫPOЖДEHHOЙ MATPИЦEЙ. COND
!           ПOЛAГAETCЯ PABHЫM 1.0E+32,ECЛИ OБHAPУЖEHA TOЧHAЯ
!           BЫPOЖДEHHOCTЬ.
!
!     IPVT -BEKTOP BEДУЩИX ЭЛEMEHTOB.
!           IPVT(K)=ИHДEKC K-Й BEДУЩEЙ CTPOKИ
!           IPVT(N)=(-1)**(ЧИCЛO ПEPECTAHOBOK)
!
!     PAБOЧEE ПOЛE. BEKTOP WORK ДOЛЖEH БЫTЬ OПИCAH И
!             BKЛЮЧEH B BЫЗOB. EГO BXOДHOE COДEPЖAHИE OБЫЧHO
!             HE ДAET BAЖHOЙ ИHФOPMAЦИИ.
!
!     OПPEДEЛИTEЛЬ MATPИЦЫ A MOЖET БЫTЬ ПOЛУЧEH HA BЫXOДE
!     ПO ФOPMУЛE:
!          DET(A)=IPVT(N)*A(1,1)*A(2,2)*...*A(N,N).
!
      REAL EK,T,ANORM,YNORM,ZNORM
      INTEGER NM1,I,J,K,KP1,KB,KM1,M
!
      IPVT(N)=1
      IF(N.EQ.1)GO TO 80
      NM1=N-1
!
!     BЫЧИCЛИTЬ 1-HOPMУ MATPИЦЫ A
!
      ANORM=0.0
      DO 10 J=1,N
        T=0.0
        DO 5 I=1,N
          T=T+ABS(A(I,J))
    5   CONTINUE
        IF(T.GT.ANORM) ANORM=T
   10 CONTINUE
!
!     ГAУCCOBO ИCKЛЮЧEHИE ! ЧACTИЧHЫM BЫБOPOM
!     BEДУЩEГO ЭЛEMEHTA
!
      DO 35 K=1,NM1
        KP1=K+1
!
!       HAЙTИ BEДУЩИЙ ЭЛEMEHT
!
        M=K
        DO 15 I=KP1,N
          IF(ABS(A(I,K)).GT.ABS(A(M,K))) M=I
   15   CONTINUE
        IPVT(K)=M
        IF(M.NE.K)IPVT(N)=-IPVT(N)
        T=A(M,K)
        A(M,K)=A(K,K)
        A(K,K)=T
!
!       ПPOПУCTИTЬ ЭTOT ШAГ,ECЛИ BEДУЩИЙ ЭЛEMEHT PABEH HУЛЮ
!
        IF(T.EQ.0.0)GO TO 35
!
!       BЫЧИCЛИTЬ MHOЖИTEЛИ
!
        DO 20 I=KP1,N
          A(I,K)=-A(I,K)/T
   20   CONTINUE
!
!       ПEPECTABЛЯTЬ И ИCKЛЮЧATЬ ПO CTOЛБЦAM
!
        DO 30 J=KP1,N
          T=A(M,J)
          A(M,J)=A(K,J)
          A(K,J)=T
          IF(T.EQ.0.0)GO TO 30
          DO 25 I=KP1,N
            A(I,J)=A(I,J)+A(I,K)*T
   25     CONTINUE
   30   CONTINUE
   35 CONTINUE
!
!     COND=(1-HOPMA MATPИЦЫ A)*(OЦEHKA ДЛЯ 1-HOPMЫ MATPИЦЫ,
!     OБPATHOЙ K A)
!     OЦEHKA ПOЛУЧAETCЯ ПOCPEДCTBOM OДHOГO ШAГA METOДA
!     OБPATHЫX ИTEPAЦИЙ ДЛЯ HAИMEHЬШEГO CИHГУЛЯPHOГO
!     BEKTOPA. ЭTO TPEБУET PEШEHИЯ ДBУX CИCTEM УPABHEHИЙ,
!     (TPAHCПOHИPOBAHHAЯ ДЛЯ A) *Y=E И A*Z=Y, ГДE E-BEKTOP
!     ИЗ +1 И -1, BЫБPAHHЫЙ TAK, ЧTOБЫ MAKCИMИЗИPOBATЬ
!     BEЛИЧИHУ Y.
!     ESTIMATE=(1-HOPMA Z)/(1-HOPMA Y)
!
!     PEШИTЬ CИCTEMУ (TPAHCПOHИPOBAHHAЯ ДЛЯ A)*Y=E
!
      DO 50 K=1,N
        T=0.0
        IF(K.EQ.1)GO TO 45
        KM1=K-1
        DO 40 I=1,KM1
          T=T+A(I,K)*WORK(I)
   40   CONTINUE
   45   EK=1.0
        IF(T.LT.0.0)EK=-1.0
        IF(A(K,K).EQ.0.0)GO TO 90
        WORK(K)=-(EK+T)/A(K,K)
   50 CONTINUE
      DO 60 KB=1,NM1
        K=N-KB
        T=WORK(K)
        KP1=K+1
        DO 55 I=KP1,N
          T=T+A(I,K)*WORK(I)
   55   CONTINUE
        WORK(K)=T
        M=IPVT(K)
        IF(M.EQ.K)GO TO 60
        T=WORK(M)
        WORK(M)=WORK(K)
        WORK(K)=T
   60 CONTINUE
!
      YNORM=0.0
      DO 65 I=1,N
        YNORM=YNORM+ABS(WORK(I))
   65 CONTINUE
!
!     PEШИTЬ CИCTEMУ A*Z=Y
!
      CALL SOLVE(NDIM,N,A,WORK,IPVT)
!
      ZNORM=0.0
      DO 70 I=1,N
        ZNORM=ZNORM+ABS(WORK(I))
   70 CONTINUE
!
!     OЦEHИTЬ OБУCЛOBЛEHHOCTЬ
!
      COND=ANORM*ZNORM/YNORM
      IF(COND.LT.1.0)COND=1.0
      RETURN
!
!     CЛУЧAЙ MATPИЦЫ 1*1
!
   80 COND=1.0
      IF(A(1,1).NE.0.0)RETURN
!
!     TOЧHAЯ BЫPOЖДEHHOCTЬ
!
   90 CONTINUE
      COND=1.0E+32
      RETURN
      END


SUBROUTINE SOLVE(NDIM,N,A,B,IPVT)
!
      INTEGER NDIM,N,IPVT(N)
      REAL A(NDIM,N),B(N)
!
!     PEШEHИE ЛИHEЙHOЙ CИCTEMЫ A*X=B
!     ПOДПPOГPAMMY HE CЛEДYET ИCПOЛЬЗOBATЬ,
!     ECЛИ DECOMP OБHAPYЖИЛ BЫPOЖДEHHOCTЬ
!
!     BXOДHAЯ ИHФOPMAЦИЯ.
!
!     NDIM -ЗAЯBЛEHHAЯ CTPOЧHAЯ PAЗMEPHOCTЬ
!           MACCИBA,COДEPЖAЩEГO A.
!     N    -ПOPЯДOK MATPИЦЫ.
!     A    -ФAKTOPИЗOBAHHAЯ MATPИЦA,ПOЛYЧEHHAЯ ИЗ DECOMP
!     B    -BEKTOP ПPABЫX ЧACTEЙ.
!     IPVT -BEKTOP BEДYЩИX ЭЛEMEHTOB,ПOЛYЧEHHЫЙ ИЗ DECOMP
!
!     BЫXOДHAЯ ИHФOPMAЦИЯ.
!
!     B    =BEKTOP PEШEHИЯ X.
!
      INTEGER KB,KM1,NM1,KP1,I,K,M
      REAL T
!
!     ПPЯMOЙ XOД
!
      IF(N.EQ.1) GO TO 50
      NM1=N-1
      DO 20 K=1,NM1
        KP1=K+1
        M=IPVT(K)
        T=B(M)
        B(M)=B(K)
        B(K)=T
        DO 10 I=KP1,N
          B(I)=B(I)+A(I,K)*T
   10   CONTINUE
   20 CONTINUE
!
!     OБPATHAЯ ПOДCTAHOBKA
!
      DO 40 KB=1,NM1
        KM1=N-KB
        K=KM1+1
        B(K)=B(K)/A(K,K)
        T=-B(K)
        DO 30 I=1,KM1
          B(I)=B(I)+A(I,K)*T
   30   CONTINUE
   40 CONTINUE
   50 B(1)=B(1)/A(1,1)
      RETURN
      END


