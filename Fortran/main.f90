
program main
  implicit none
  character*16::filename
  character*3::num
  integer(8),  parameter :: Nr_max =50, Nz_max =100, Nr_min=Nr_max/9, Nz_min=Nz_max/9
  integer(8):: Nr,Nz,Nr_i,Nr_f,Nr_inc, nume, Nri,iexp,GCI_bol,i
  integer(8),dimension(3):: n_gci
  real(8):: Ta,Lr,Lz,Tinf,Erro_p_max, Erro_p_med
  real(8),dimension(Nr_min*Nz_min,3):: T_gci
  real(8),dimension(3):: h_gci
  real(8),dimension(Nr_min*Nz_min):: Tpas,error95,Rl_pas,Zl_pas

  Ta=37.d0
  Tinf=20.d0

  Nr_i=3
  Nr_f=3
  Nr_inc=1

  Lr=1.d0!4.d-2
  Lz=1.d0!20.d-2

  GCI_bol=0
  IF (GCI_bol==0) THEN
    iexp=-1
    Nr_i=100
    Nr_f=Nr_i
    Nr_inc=1
    WRITE(unit=num,fmt='(I2.2)')Nr+1
    filename="converg3.txt"
    OPEN(unit=10,file=filename,status="replace")
    do Nr=Nr_i,Nr_f,Nr_inc
      !Nz=(Lz/Lr)*Nr
      Nz=Nr
      CALL diffusion(Nr,Nz,Lr,Lz,Ta,Tinf,Erro_p_max,Erro_p_med,iexp,Nr_min,Nz_min,Tpas,Rl_pas,Zl_pas)
      PRINT *, "foi"
      PRINT *, "Erro", Erro_p_max
      WRITE(unit=10,fmt=*) Nr, Erro_p_max , Erro_p_med
    end do
    close(unit=11)

  ELSE

    do iexp=3,1,-1
      Nr=Nr_max/(3.d0**(iexp-1.d0))
      Nz=Nz_max/(3.d0**(iexp-1.d0))
      CALL diffusion(Nr,Nz,Lr,Lz,Ta,Tinf,Erro_p_max,Erro_p_med,iexp,Nr_min,Nz_min,Tpas,Rl_pas,Zl_pas)
      T_gci(:,iexp)=Tpas
      n_gci(iexp)=Nr*Nz
      h_gci(iexp)=(Lr*Lz/(dble(n_gci(iexp))))**0.5d0
      print *,"Nr",Nr,"Lr",Lr,"Nz",Nz,"Lz",Lz,"dr",Lr/dble(Nr),"dz",Lz/dble(Nz)
      print*, "h:",h_gci(iexp)

    end do

    CALL GCI(Nr_min, Nz_min,T_gci ,n_gci, h_gci,error95)


    filename="Temp_gci.txt"
    OPEN(unit=11,file=filename,status="replace")
    DO i=1,Nr_min*Nz_min
      WRITE(unit=11,fmt='(f10.6,A,f10.7,A,f9.7,A,f9.7)')T_gci(i,1), ";", error95(i), ";", Rl_pas(i),";", Zl_pas(i)
    END DO
    CLOSE(unit=11)

    filename="Temp_gci_sec.txt"
    OPEN(unit=11,file=filename,status="replace")
    DO i=1,Nr_min*Nz_min,Nz_min
      WRITE(unit=11,fmt='(f10.6,A,f10.7,A,f9.7,A,f9.7)')T_gci(i,1), ";", error95(i), ";", Rl_pas(i),";", Zl_pas(i)
    END DO
    CLOSE(unit=11)

  END IF

  close(unit=11)

end program main

SUBROUTINE diffusion(Nr,Nz,Lr,Lz,Ta,Tinf,Erro_p_max,Erro_p_med,iexp,Nr_min,Nz_min,Tpas,Rl_pas,Zl_pas)
  implicit none
  character*15::filename
  character(12) :: steste
  integer(8)::Nr, Nz,iexp,Nr_min,Nz_min
  integer(8)::i,j,integ, c,wr,raz_pel,raz_mus,raz_os, raz_total
  real(8)::L,Lr,Lz,Ta,Tb,dr,dz,Cb,Rhop,Wp,h,Tinf,Ateste,Erro_p_max,Erro_p_med, Ts,T0, sum,K_pel,K_os,K_mus,W_pel, W_mus,W_os, NLin
  real(8):: K_tumor, W_tumor
  real(8)::T3,T2,T1,R3,R2,R1,qp,qp1,qp2,qp3, qpsurface
  real(8)::C1,C2,C3,C4,Wc, sourterm, P,P2
  real(8),  parameter :: PI  = 4 * atan (1.0_8)
  real(8)::AreaT, Area1,Area2, Area3
  real(8),dimension(Nr,Nz)::R,Z,K, Erro_mat,W
  real(8),dimension(Nr*Nz,Nr*Nz)::A
  real(8),dimension(Nr*Nz)::Kl,Rl,Zl,Wl, B
  real(8),dimension(Nr*Nz)::Area,T,Tn
  real(8),dimension(Nr*Nz)::T_a,Erro,Erro_p,Tf
  real(8)::tempo1,tempo2,tempo3, diftempo, diftempo2
  real(8),dimension(Nr_min*Nz_min)::Tpas,Rl_pas,Zl_pas
  real(8)::ZE,ZW,RN,RS,kls,wps,kp,T_inf,MMS

  A=0.d0
  B=0.d0
  R=0.d0
  Z=0.d0
  K=0.d0
  W=0.d0
  T=0.d0

  MMS=1.d0
  NLin=1.d0

  Rhop=1100.d0
  Cb=3600.d0
  h=10.d0

  dr=Lr/dble(Nr)		!Distancia entre os elementos no sentido de r
  dz=Lz/dble(Nz)		!Distancia entre os elementos no sentido de z

  K_pel=0.5d0;

  print *, "---------------------------------------------------------------------------------------------------"
  print * , Lr , Lz, Nr, Nz
  !read(*,*)
  do i=Nr,1,-1		!Looping para atribuir uma posição em r a cada elemento
    R(i,1:Nz)=dr*(0.5d0+dble(Nr-i))
  enddo

  do i=1,Nz		!Looping para atribuir uma posição em z a cada elemento
    Z(1:Nr,i)=dz*(0.5d0+dble(i-1))
  enddo

  K_pel=0.293d0;
  K_mus=0.449d0;
  K_os=0.36d0;

  W_pel=2.22d-5;
  W_mus=7.03d-4;
  W_os=6.57d-4;

  K_tumor=0.449d0;
  W_tumor= W_mus*10.d0;

  raz_pel=1
  raz_mus=4
  raz_os=3
  raz_total=raz_pel+raz_mus+raz_os

  do i=1,Nr*(raz_pel)/raz_total
    do j=1,Nz
      K(i,j)=K_pel
      W(i,j)=W_pel
    end do
  end do

  do i=Nr*(raz_pel)/raz_total+1,Nr*(raz_pel+raz_mus)/raz_total
    do j=1,Nz
      K(i,j)=K_mus
      W(i,j)=W_mus
    end do
  end do

  do i=Nr*(raz_pel+raz_mus)/raz_total+1,Nr
    do j=1,Nz
      K(i,j)=K_os
      W(i,j)=W_os
    end do
  end do

  print *, "dr:",dr,"dz",dz

  do i=26,31
    do j=28,29
      K(i,j)=K_tumor
      W(i,j)=W_tumor
    end do
  end do

  do i=1,Nr
    do j=1,Nz
      Rl(Nz*(i-1)+j)=R(i,j)
      Zl(Nz*(i-1)+j)=Z(i,j)
      Kl(Nz*(i-1)+j)=K(i,j)
      Wl(Nz*(i-1)+j)=W(i,j)
    enddo
  enddo

  IF (MMS==1) THEN
    K=0.5d0;
    Kl=0.5d0;
    Wc=5d-4;
    W=Wc
    Wl=Wc
    K_pel=0.5d0;
    C1=4*PI/(Lr)
    C2=PI/2.d0
    C3=4*PI/(Lz)
    C4=PI/2.d0
    P=Tinf

    do i=1,Nr
      do j=1,Nz
        c=(i-1)*Nz+j
        P2=-dsin(C3*Z(i,j)+C4)
        Tf(Nz*(i-1)+j)=P+P2+dsin(C1*R(i,j)+C2)*dsin(C3*Z(i,j)+C4)
        !print*, dcos(C1*Lr+C2),"ok"
        !read(*,*)
      enddo
    enddo
    filename="Temp_f.txt"
    OPEN(unit=11,file=filename,status="replace")
    DO i=1,Nr*Nz
      WRITE(unit=11,fmt='(f10.6,A,f9.7,A,f9.7)')Tf(i), ";", Rl(i), ";", Zl(i)
    END DO
    CLOSE(unit=11)
  END IF


  do i=2,Nr-1            !Looping dos elementos internos
    do j=2,Nz-1
      c=(i-1)*Nz+j
      A(c,c-Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-Nz)))*(Rl(c)+Rl(c-Nz))/2.d0*dz/dr !North
      A(c,c+Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+Nz)))*(Rl(c)+Rl(c+Nz))/2.d0*dz/dr !South
      A(c,c+1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+1)))*Rl(c)*dr/dz !Eact
      A(c,c-1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-1)))*Rl(c)*dr/dz !West
      A(c,c)=-1.d0*(A(c,c-Nz)+A(c,c+Nz)+A(c,c+1)+A(c,c-1)+Rhop*cb*Wl(c)*Rl(c)*dr*dz*NLin)
      B(c) = -Rhop*cb*Wl(c)*Rl(c)*dr*dz*Ta
    enddo
  enddo

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!Condições de contorno!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Fronteira Esquerda!
    do i=2,Nr-1
      do j=1,1
        c=(i-1)*Nz+j
        A(c,c-Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-Nz)))*(Rl(c)+Rl(c-Nz))/2.d0*dz/dr !Nzorth
        A(c,c+Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+Nz)))*(Rl(c)+Rl(c+Nz))/2.d0*dz/dr !South
        A(c,c+1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+1)))*Rl(c)*dr/dz !East
        A(c,c-1)=0.d0                                         !West
        A(c,c)=-1.d0*(A(c,c-Nz)+A(c,c+Nz)+A(c,c+1)+Rhop*Cb*Wl(c)*Rl(c)*dr*dz*NLin)
        B(c) = -Rhop*Cb*Wl(c)*Rl(c)*dr*dz*Ta
      enddo
    enddo

    !Fronteira Direita!
    do i=2,Nr-1
      do j=Nz,Nz
        c=(i-1)*Nz+j
        A(c,c-Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-Nz)))*(Rl(c)+Rl(c-Nz))/2.d0*dz/dr !Nzorth
        A(c,c+Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+Nz)))*(Rl(c)+Rl(c+Nz))/2.d0*dz/dr !South
        A(c,c+1)=0.d0                                         !East
        A(c,c-1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-1)))*Rl(c)*dr/dz !West
        A(c,c)=-1.d0*(A(c,c-Nz)+A(c,c+Nz)+A(c,c-1)+Rhop*Cb*Wl(c)*Rl(c)*dr*dz*NLin)
        B(c) = -Rhop*Cb*Wl(c)*Rl(c)*dr*dz*Ta
      enddo
    enddo

    !Fronteira Inferior!
    do i=Nr,Nr
      do j=2,Nz-1
        c=(i-1)*Nz+j
        A(c,c-Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-Nz)))*(Rl(c)+Rl(c-Nz))/2.d0*dz/dr !North
        A(c,c+1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+1)))*Rl(c)*dr/dz !East
        A(c,c-1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-1)))*Rl(c)*dr/dz !West
        A(c,c)=-1.d0*(A(c,c-Nz)+A(c,c+1)+A(c,c-1)+Rhop*Cb*Wl(c)*Rl(c)*dr*dz*NLin)
        B(c) = -Rhop*cb*Wl(c)*Rl(c)*dr*dz*Ta
      enddo
    enddo

    !Fronteira Superior!
    do i=1,1
      do j=2,Nz-1
        c=(i-1)*Nz+j
        A(c,c+Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+Nz)))*(Rl(c)+Rl(c+Nz))/2.d0*dz/dr !South
        A(c,c+1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+1)))*Rl(c)*dr/dz !Eact
        A(c,c-1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-1)))*Rl(c)*dr/dz !West
        A(c,c)=-1.d0*(A(c,c+Nz)+A(c,c+1)+A(c,c-1)+Rhop*Cb*Wl(c)*Rl(c)*dr*dz*NLin+h/(1.d0+h*dr/Kl(c))*(Rl(c)+dr/2.d0)*dz)
        B(c) = -Rhop*Cb*Wl(c)*Rl(c)*dr*dz*Ta-h/(1.d0+h*dr/Kl(c))*(Rl(c)+dr/2.d0)*dz*Tinf
      enddo
    enddo

    !Canto N-W!
    i=1
    j=1
    c=(i-1)*Nz+j
    A(c,c+Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+Nz)))*(Rl(c)+Rl(c+Nz))/2.d0*dz/dr !South
    A(c,c+1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+1)))*Rl(c)*dr/dz !East
    A(c,c)=-1.d0*(A(c,c+Nz)+A(c,c+1)+Rhop*Cb*Wl(c)*Rl(c)*dr*dz*NLin+h/(1.d0+h*dr/Kl(c))*(Rl(c)+dr/2.d0)*dz)
    B(c) = -Rhop*Cb*Wl(c)*Rl(c)*dr*dz*Ta-h/(1+h*dr/Kl(c))*(Rl(c)+dr/2.d0)*dz*Tinf

    !Canto N-E!
    i=1
    j=Nz
    c=(i-1)*Nz+j
    A(c,c+Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+Nz)))*(Rl(c)+Rl(c+Nz))/2.d0*dz/dr !South
    A(c,c-1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-1)))*Rl(c)*dr/dz !West
    A(c,c)=-1.d0*(A(c,c+Nz)+A(c,c-1)+Rhop*Cb*Wl(c)*Rl(c)*dr*dz*NLin+h/(1.d0+h*dr/Kl(c))*(Rl(c)+dr/2.d0)*dz)
    B(c) = -Rhop*Cb*Wl(c)*Rl(c)*dr*dz*Ta-h/(1+h*dr/Kl(c))*(Rl(c)+dr/2.d0)*dz*Tinf

    !Canto S-W!
    i=Nr
    j=1
    c=(i-1)*Nz+j
    A(c,c-Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-Nz)))*(Rl(c)+Rl(c-Nz))/2.d0*dz/dr !Nzorth
    A(c,c+1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c+1)))*Rl(c)*dr/dz !East
    A(c,c)=-1.d0*(A(c,c-Nz)+A(c,c+1)+Rhop*Cb*Wl(c)*Rl(c)*dr*dz*NLin)
    B(c) =- Rhop*Cb*Wl(c)*Rl(c)*dr*dz*Ta

    !Canto S-E!
    i=Nr
    j=Nz
    c=(i-1)*Nz+j
    A(c,c-Nz)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-Nz)))*(Rl(c)+Rl(c-Nz))/2.d0*dz/dr!   Nzorth
    A(c,c-1)=(2.d0/(1.d0/Kl(c)+1.d0/Kl(c-1)))*Rl(c)*dr/dz                   !West
    A(c,c)=-1.d0*(A(c,c-Nz)+A(c,c-1)+Rhop*Cb*Wl(c)*Rl(c)*dr*dz*NLin)
    B(c) = -Rhop*Cb*Wl(c)*Rl(c)*dr*dz*Ta

    kp=K_pel
    T_inf=Tinf

    IF (MMS==1) THEN
      do i=1,Nr
        do j=1,Nz
          c=(i-1)*Nz+j
          kls=Kl(c)
          ZE=Z(i,j)+dz/2.d0
          ZW=Z(i,j)-dz/2.d0
          RN=R(i,j)+dr/2.d0
          RS=R(i,j)-dr/2.d0
          wps=Wl(c)

 sourterm=-C1*RN*kls*cos(C1*RN + C2)*cos(C3*ZE + C4)/C3 + C1*RS*kls*cos(C1*RS+ C2)*cos(C3*ZE + C4)/C3 &
 - 1.0d0/2.0d0*C3*RN**2*kls*cos(C3*ZE + C4) + (1.0d0/2.0d0)*C3*RS**2*kls*cos(C3*ZE + C4) &
 - 1.0d0/2.0d0*Cb*RN**2*Rhop*T_inf*wps*ZE + (1.0d0/2.0d0)*Cb*RN**2*Rhop*Ta*wps*ZE &
 + (1.0d0/2.0d0)*Cb*RS**2*Rhop*T_inf*wps*ZE - 1.0d0/2.0d0*Cb*RS**2*Rhop*Ta*wps*ZE &
 - 1.0d0/2.0d0*Cb*RN**2*Rhop*wps*cos(C3*ZE + C4)/C3 +(1.0d0/2.0d0)*Cb*RS**2*Rhop*wps*cos(C3*ZE + C4)/C3 &
 - C3*RN*kls*cos(C1*RN + C2)*cos(C3*ZE + C4)/C1 + C3*RS*kls*cos(C1*RS + C2)*cos(C3*ZE + C4)/C1 &
 - Cb*RN*Rhop*wps*cos(C1*RN + C2)*cos(C3*ZE + C4)/(C1*C3) + Cb*RS*Rhop*wps*cos(C1*RS + C2)*cos(C3*ZE + C4)/(C1*C3) &
 + C3*kls*sin(C1*RN + C2)*cos(C3*ZE + C4)/C1**2 - C3*kls*sin(C1*RS + C2)*cos(C3*ZE + C4)/C1**2 &
 + Cb*Rhop*wps*sin(C1*RN + C2)*cos(C3*ZE + C4)/(C1**2*C3) - Cb*Rhop*wps*sin(C1*RS + C2)*cos(C3*ZE + C4)/(C1**2*C3)&
 +C1*RN*kls*cos(C1*RN + C2)*cos(C3*ZW + C4)/C3 - C1*RS*kls*cos(C1*RS+ C2)*cos(C3*ZW + C4)/C3 &
 + (1.0d0/2.0d0)*C3*RN**2*kls*cos(C3*ZW +C4) - 1.0d0/2.0d0*C3*RS**2*kls*cos(C3*ZW + C4) &
 + (1.0d0/2.0d0)*Cb*RN**2*Rhop*T_inf*wps*ZW - 1.0d0/2.0d0*Cb*RN**2*Rhop*Ta*wps*ZW &
 - 1.0d0/2.0d0*Cb*RS**2*Rhop*T_inf*wps*ZW + (1.0d0/2.0d0)*Cb*RS**2*Rhop*Ta*wps*ZW &
 + (1.0d0/2.0d0)*Cb*RN**2*Rhop*wps*cos(C3*ZW + C4)/C3- 1.0d0/2.0d0*Cb*RS**2*Rhop*wps*cos(C3*ZW + C4)/C3 &
 + C3*RN*kls*cos(C1*RN + C2)*cos(C3*ZW + C4)/C1 - C3*RS*kls*cos(C1*RS + C2)*cos(C3*ZW + C4)/C1 &
 + Cb*RN*Rhop*wps*cos(C1*RN + C2)*cos(C3*ZW + C4)/(C1*C3) - Cb*RS*Rhop*wps*cos(C1*RS + C2)*cos(C3*ZW + C4)/(C1*C3)&
 - C3*kls*sin(C1*RN + C2)*cos(C3*ZW + C4)/C1**2 + C3*kls*sin(C1*RS + C2)*cos(C3*ZW + C4)/C1**2 &
 - Cb*Rhop*wps*sin(C1*RN + C2)*cos(C3*ZW + C4)/(C1**2*C3) + Cb*Rhop*wps*sin(C1*RS + C2)*cos(C3*ZW + C4)/(C1**2*C3)

          B(c) = B(c)+sourterm
          sourterm=0
        enddo
      enddo
    END IF

    tempo2=time()
    CALL sor(T,A,B,Nr*Nz)
    tempo3=time()
    diftempo=tempo2-tempo1
    diftempo2=tempo3-tempo2
    PRINT * , "Tempo:", diftempo2

    IF (iexp==-1)  THEN
      IF (NLin==0.d0) THEN
        CALL analitico(Nr,Nz,Lr,raz_pel,raz_mus,raz_total,Rhop,Cb,h,W_pel,W_mus,W_os,K_pel,K_mus,K_os,Wl,Rl,Tinf,Ta,T&
                                                              ,T_a,Erro,Erro_p,Erro_p_max,Erro_p_med)
      END IF
    END IF

    IF (iexp==-1) THEN
      Tpas=0
    ELSE IF ( iexp==3 ) THEN
      Tpas=T
      Rl_pas=Rl
      Zl_pas=Zl
    ELSE IF ( iexp==2 ) THEN
      DO i=1,Nr_min
        DO j=1,Nz_min
          Tpas((i-1)*Nr_min+j)=T(3*3*Nz_min*(i-1)+3*Nz_min+2+(j-1)*3)
          !print*, (i-1)*Nr_min+j , "->",3*3*Nr_min*(i-1)+3*Nz_min+2+(j-1)*3
        END DO
      END DO
    ELSE
      DO i=1,Nr_min
        DO j=1,Nz_min
          Tpas((i-1)*Nr_min+j)=T(9*9*Nz_min*(i-1)+9*Nz_min*4+5+(j-1)*9)
          !print*, (i-1)*Nr_min+j , "->",9*9*Nz_min*(i-1)+9*Nz_min*4+5+(j-1)*9
        END DO
      END DO
    END IF

    filename="Temp_a.txt"
    OPEN(unit=11,file=filename,status="replace")
    DO i=20,Nr*Nz,Nz
      WRITE(unit=11,fmt='(f11.6,A,f11.6,A,f11.6,A,f11.6)')T(i), ";", T_a(i), ";", Erro_p(i), ";" , Rl(i)
    END DO
    CLOSE(unit=11)

    filename="Temp_2d.txt"
    OPEN(unit=11,file=filename,status="replace")
    DO i=1,Nr*Nz
      WRITE(unit=11,fmt='(f10.6,A,f9.7,A,f9.7)')abs( T(i)-Tf(i)) ,";", Rl(i), ";", Zl(i)
    END DO
    CLOSE(unit=11)

    filename="Temp_2d_r.txt"
    OPEN(unit=11,file=filename,status="replace")
    DO i=1,Nr
      WRITE(unit=11,fmt='(f10.6)') R(i,1)
    END DO
    CLOSE(unit=11)

    filename="Temp_2d_z.txt"
    OPEN(unit=11,file=filename,status="replace")
    DO j=1,Nz
      WRITE(unit=11,fmt='(f10.6)') Z(1,j)
    END DO
    CLOSE(unit=11)

    filename="Temp_2d_su.txt"
    OPEN(unit=11,file=filename,status="replace")
    DO i=1,Nz
      WRITE(unit=11,fmt='(2f10.6)') T(i), Zl(i)
    END DO
    CLOSE(unit=11)

END SUBROUTINE diffusion

SUBROUTINE analitico(Nr,Nz,Lr,raz_pel,raz_mus,raz_total,Rhop,Cb,h,W_pel,W_mus,W_os,K_pel,K_mus,K_os,Wl,Rl,Tinf,Ta,T&
                                                          &,T_a,Erro,Erro_p,Erro_p_max,Erro_p_med)
  IMPLICIT NONE
  integer(8)::Nr, Nz
  integer(8)::i,j, c,wr,raz_pel,raz_mus,raz_os, raz_total
  real(8)::L,Lr,Ta,Cb,Rhop,Wp,h,Tinf,Erro_p_max,Erro_p_med,sum,K_pel,K_os,K_mus,W_pel, W_mus,W_os
  real(8)::T3,T2,T1,R3,R2,R1,qp1,qp2,qp3, qsurface
  real(8),  parameter :: PI  = 4 * atan (1.0_8)
  real(8)::AreaT, Area1,Area2, Area3
  real(8),dimension(Nr*Nz)::Rl,Wl
  real(8),dimension(Nr*Nz)::T
  real(8),dimension(Nr*Nz)::T_a,Erro,Erro_p,Tn

  qp3=Rhop*Cb*W_pel*Ta
  qp2=Rhop*Cb*W_mus*Ta
  qp1=Rhop*Cb*W_os*Ta
  R3=Lr
  R2=Lr*(1.d0-dble(raz_pel)/dble(raz_total))!R(Nr*(raz_pel)/raz_total+1,1)
  R1=Lr*(1.d0-dble(raz_pel+raz_mus)/dble(raz_total))!R(Nr*(raz_pel+raz_mus)/raz_total+1,1)
  AreaT=PI*(R3**2.d0)
  Area3=PI*(R3+R2)*(R3-R2)
  Area2=PI*(R2+R1)*(R2-R1)
  Area1=PI*(R1**2)
  qsurface=(qp1*Area1+qp2*Area2+qp3*Area3)/AreaT

  T3=Tinf+qsurface*R3/(2.d0*h)

  T2=T3+qp2*(R2**2)*dlog(R3/R2)/(2*K_pel)+(R1**2)*dlog(R3/R2)*(qp1-qp2)/(2*K_pel)-qp3*R3*R2*dlog(R3/R2)/(4*K_pel)&
                                                                              &+qp3*(R3**3)*(1-(R2/R3)**2)/K_pel
  T1=T2+(qp1-qp2)*(dlog(R2/R1))*(R1**2)/(2.d0*K_mus)+qp2*(R2**2)*(1-(R1/R2)**2.d0)/(4.d0*K_mus)

  sum=0.d0
  DO i=1,Nr*Nz
    IF (i<(Nr*(raz_pel)/raz_total)*Nz ) THEN
      T_a(i)=T3+qp3*(R3**2)*(1-(Rl(i)/R3)**2.d0)/(4.d0*K_pel)-(qp3*(R3**2)*(1-(R2/R3)**2.d0)/(4.d0*K_pel)+T3-T2)*&
                                                                              &(dlog(R3/Rl(i)))/dlog(R3/R2)
    ELSE IF( i<(Nr*(raz_pel+raz_mus)/raz_total)*Nz ) THEN
      T_a(i)=T2+qp2*(R2**2)*(1-(Rl(i)/R2)**2.d0)/(4.d0*K_mus)-(qp2*(R2**2)*(1-(R1/R2)**2.d0)/(4.d0*K_mus)+T2-T1)*&
                                                                              &(dlog(R2/Rl(i)))/dlog(R2/R1)
    ELSE
      T_a(i)=T1+qp1*(R1**2)*(1-(Rl(i)/R1)**2.d0)/(4.d0*K_os)
    END IF

   Erro(i)=abs( T(i)-T_a(i))
   Erro_p(i)=Erro(i)/abs(T_a(i))
   sum=sum+Erro_p(i)
  END DO

  Erro_p_max=maxval(Erro_p)
  Erro_p_med=sum/(Nr*Nz)
    print*, "T3:",T3,"T2:",T2,"T1:",T1

END SUBROUTINE analitico

SUBROUTINE sor(vx,A,vb,N)
  IMPLICIT NONE
  INTEGER(8)::i,j,N
  REAL(8), PARAMETER::pi=4.D0*datan(1.D0)
  REAL(8),DIMENSION(N)::vb,vx,vxnew
  REAL(8),DIMENSION(N,N)::A
  INTEGER(8)::contador
  REAL(8)::tol,eps,wopt,h

  h=1.d0/(dble(N)+1.d0)
  wopt=2.d0/(1.d0+dsin(pi*h))
  print *, "wopt:",wopt
  wopt=wopt*0.98
  wopt=1.d0
  vxnew=0.d0
  vx=0.d0
  tol=1.d-3
  eps=1.d0
  contador=0

  DO WHILE(eps.gt.tol)
    DO i = 1,N        !Varia a linha da matriz
      vxnew(i)=vb(i)
      DO j = 1,i-1
        vxnew(i)=vxnew(i)-A(i,j)*vxnew(j)
      ENDDO
      DO j = i+1, N
        vxnew(i)=vxnew(i)-A(i,j)*vx(j)
      ENDDO
      vxnew(i)=wopt*vxnew(i)/A(i,i)+(1.d0-wopt)*vx(i)
    ENDDO
    eps=DABS(MAXVAL(vxnew-vx))
    vx=vxnew
    contador=contador+1
  ENDDO
  PRINT *, "contador:", contador

END SUBROUTINE sor

SUBROUTINE GCI(Nr_min, Nz_min,T_gci ,n_gci, h_gci,error95)
  integer(8):: Nr_min, Nz_min ,i
  integer(8),dimension(3):: n_gci
  real(8),dimension(Nr_min*Nz_min,3):: T_gci
  real(8),dimension(Nr_min*Nz_min):: p, ea, error95
  real(8),dimension(3):: h_gci
  real(8),dimension(2):: ref
  real(8):: FS

  FS=1.25d0

  ref(1)=h_gci(2)/h_gci(1)  !ref: 1 -> 2
  ref(2)=h_gci(3)/h_gci(2)
    !ref: 1 -> 2

  DO i=1,Nr_min*Nz_min
    p(i)=(1/log(ref(1)))*(log(abs((T_gci(i,3)-T_gci(i,2))/(T_gci(i,2)-T_gci(i,1)))))
    ea(i)=abs(T_gci(i,1)-T_gci(i,2))
    error95(i)=FS*ea(i)/(ref(1)-1.d0)
  END DO
END SUBROUTINE
