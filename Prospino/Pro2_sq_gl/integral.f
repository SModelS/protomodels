      REAL*8 FUNCTION SL1B0A(MS2,MG2,MT2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2
      INTEGER INUM
      IF (INUM.EQ.1) SL1B0A = -LOG(MG2/MS2)
C***  FINITE PART OF B0(K1,MG,MG)
      IF (INUM.EQ.2) SL1B0A = -LOG(MT2/MS2)
C***  FINITE PART OF B0(K1,MT,MT)
      IF (INUM.EQ.3) SL1B0A = 1.D0 -MG2/(MG2 -MS2)*LOG(MG2/MS2)
C***  FINITE PART OF B0(K1,MG,MS)
ctp      print*, " SL1B0A ",SL1B0A
      RETURN
      END

      REAL*8 FUNCTION SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,
     +     KBETAT,KXT,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,S,BETA,XS,A,B,X1,X2
      COMPLEX*16 KBETAG,KXG,KBETAT,KXT
      INTEGER INUM
      IF (INUM.EQ.1) SL1B0B = 2.D0 - LOG(S/MS2)
C***  FINITE PART OF B0(K1+K2,0,0)
      IF (INUM.EQ.2) SL1B0B = 2.D0 + BETA*LOG(XS)
C***  FINITE PART OF B0(K1+K2,MS,MS)
      IF (INUM.EQ.3) SL1B0B =  
     +     2.D0 -LOG(MG2/MS2) +REAL(KBETAG*LOG(-KXG))
C***  FINITE PART OF B0(K1+K2,MG,MG)
      IF (INUM.EQ.4) SL1B0B =  
     +     2.D0 -LOG(MT2/MS2) +REAL(KBETAT*LOG(-KXT))
C***  FINITE PART OF B0(K1+K2,MT,MT)
      IF (INUM.EQ.5) THEN
         A = 1.D0 +(MS2 -MG2)/S
         B = MS2/S
         X1 = A/2.D0 + SQRT(A**2/4.D0 - B)
         X2 = A/2.D0 - SQRT(A**2/4.D0 - B)
         SL1B0B = 2.D0 -LOG(S/MS2) - (1.D0 -X1)*LOG(ABS(1-X1)) 
     +        -X1*LOG(ABS(-X1)) - (1.D0 -X2)*LOG(ABS(1 -X2))
     +        -X2*LOG(ABS(-X2))
C***  FINITE PART OF B0(K1+K2,MS,MG)
      END IF
ctp      print*, " SL1B0B ",SL1B0B
      RETURN
      END

      REAL*8 FUNCTION SL1B0C(MS2,MG2,MT2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2
      COMPLEX*16 KDEL,KA,KB,KX1,KX2
      INTEGER INUM
      IF (INUM.EQ.1) SL1B0C = 2.D0
C***  FINITE PART OF B0(-P2,0,MS)
      IF (INUM.EQ.2)  SL1B0C =  
     +   2.D0 -LOG(MG2/MS2) -(1.D0-MG2/MS2)*LOG(ABS(1.D0-MS2/MG2) )
C***  FINITE PART OF B0(-P2,0,MG)
      IF (INUM.EQ.3) THEN
         KDEL = (0.D0,1.D-15)
         KA = DCMPLX(1.D0 +(MG2 -MT2)/MS2)
         KB = MG2/MS2 -KDEL
         KX1 = KA/2.D0 + SQRT(KA**2/4.D0 - KB)
         KX2 = KA/2.D0 - SQRT(KA**2/4.D0 - KB)
         SL1B0C = REAL(2.D0 - (1.D0 -KX1)*LOG(1.D0-KX1) 
     +        -KX1*LOG(-KX1) -(1.D0 -KX2)*LOG(1.D0 -KX2) 
     +        -KX2*LOG(-KX2) )
C***  FINITE PART OF B0(-P2,MG,MT)
      END IF
ctp      print*, " SL1B0C ",SL1B0C
      RETURN
      END


      REAL*8 FUNCTION SL1B0D(MS2,MG2,MT2,T,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,T,A,B,X1,X2
      INTEGER INUM
      IF (INUM.EQ.1) SL1B0D = 2.D0 -(T-MS2)/T*LOG((MS2-T)/MS2)
C***  FINITE PART OF B0(K2-P2,0,MS)
      IF (INUM.EQ.2) SL1B0D = 
     +     2.D0 -LOG(MG2/MS2) -(T-MG2)/T*LOG((MG2-T)/MG2)
C***  FINITE PART OF B0(K2-P2,0,MG)      
      IF (INUM.EQ.3) THEN
C***  FINITE PART OF B0(K2-P2,MS,MT)
         A = 1.D0 +(MS2 -MT2)/T
         B = MS2/T
         X1 = A/2.D0 + SQRT(A**2/4.D0 - B)
         X2 = A/2.D0 - SQRT(A**2/4.D0 - B)
         SL1B0D = 2.D0 -LOG(-T/MS2) - (1.D0 -X1)*LOG(ABS(1-X1)) 
     +        -X1*LOG(ABS(-X1)) - (1.D0 -X2)*LOG(ABS(1 -X2))
     +        -X2*LOG(ABS(-X2))
      END IF
      IF (INUM.EQ.4) THEN
C***  FINITE PART OF B0(K2-P2,MG,MT)
         A = 1.D0 +(MG2 -MT2)/T
         B = MG2/T
         X1 = A/2.D0 + SQRT(A**2/4.D0 - B)
         X2 = A/2.D0 - SQRT(A**2/4.D0 - B)
         SL1B0D = 2.D0 -LOG(-T/MS2) - (1.D0 -X1)*LOG(ABS(1-X1)) 
     +        -X1*LOG(ABS(-X1)) - (1.D0 -X2)*LOG(ABS(1 -X2))
     +        -X2*LOG(ABS(-X2))
      END IF
ctp      print*, " SL1B0D ",SL1B0D
      RETURN
      END

      REAL*8 FUNCTION SL1B0E(MS2,MG2,MT2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2
      INTEGER INUM
      COMPLEX*16 KDEL, KA, KB, KX1, KX2
      IF (INUM.EQ.1) SL1B0E = 2.D0 -LOG(MG2/MS2)
C***  FINITE PART OF B0(KG = MG,0,MG)
      IF (INUM.EQ.2) SL1B0E = 
     +     2.D0 - (1.D0 -MS2/MG2)*LOG(ABS(1.D0 -MG2/MS2))
C***  FINITE PART OF B0(KG = MG,0,MS)
      IF (INUM.EQ.3) THEN
         KDEL = (0.D0,1.D-15)
         KA = 1.D0 +(MS2 -MT2)/MG2
         KB = MS2/MG2 -KDEL
         KX1 = KA/2.D0 + SQRT(KA**2/4.D0 - KB)
         KX2 = KA/2.D0 - SQRT(KA**2/4.D0 - KB)
         SL1B0E = REAL(2.D0 -LOG(MG2/MS2) -(1.D0-KX1)*LOG(1.D0-KX1) 
     +        -KX1*LOG(-KX1) -(1.D0 -KX2)*LOG(1.D0 -KX2) 
     +        -KX2*LOG(-KX2))
C***  FINITE PART OF B0(KG = MG,MT,MS)
      END IF
ctp      print*, " SL1B0E ",SL1B0E
      RETURN
      END

      REAL*8 FUNCTION SL1BP(MS2,MG2,MT2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,M2
      COMPLEX*16 KMG2,KMS2,KMT2,KA,KB,KX1,KX2,KF1
      INTEGER INUM
      IF (INUM.EQ.1) SL1BP = -1.D0/MS2
C***  FINITE PART OF B0P(-P2,0,MS)
      IF (INUM.EQ.2) SL1BP  = 
     +      1.D0/MS2*( -1.D0 -MG2/MS2*LOG(ABS(1.D0 -MS2/MG2)))
C***  FINITE PART OF B0P(-P2,0,MG)
      IF (INUM.EQ.3) THEN
C***  FINITE PART OF B0P(K1,MS,MG) 
         M2 = MG2 -MS2         
         SL1BP  = -1.D0/2.D0/M2 +MG2/M2**2 - MG2*MS2/M2**3*LOG(MG2/MS2)
      END IF
      IF (INUM.EQ.4) SL1BP = 1.D0/6.D0/MS2
C***  FINITE PART OF B0P(K1,MS,MS) 
      IF (INUM.EQ.5) THEN
C***  FINITE PART OF B0P(-P2,MG,MT) 
C***  HERE WE USE THE WIDTHS OF THE PARTICLES
         KMS2 = MS2 -(0.D0,1.D0)*MS2*0.25D0*(MS2-MG2)**2/MS2**2
         KMG2 = MG2 -(0.D0,1.D0)*SQRT(MG2)*0.2D0
         KMT2 = MT2 -(0.D0,1.D0)*SQRT(MT2)*1.5D0
         KA = 1.D0 +(KMG2 -KMT2)/KMS2
         KB = KMG2/KMS2 
         KX1 = KA/2.D0 + SQRT(KA**2/4.D0 - KB)
         KX2 = KA/2.D0 - SQRT(KA**2/4.D0 - KB)
         KF1 = - 1.D0 + (KX1 -KX1**2)/(KX1-KX2)*LOG((KX1-1.D0)/KX1)
     +        - (KX2 -KX2**2)/(KX1-KX2)*LOG((KX2-1.D0)/KX2) 
         SL1BP = REAL(KF1/KMS2)
      END IF
      IF (INUM.EQ.6) SL1BP = (-1.D0 -1.D0/2.D0*LOG(MS2/MG2))/MG2
C***  FINITE PART OF B0P(KG=MG,0,MG)
      IF (INUM.EQ.7)      
     +     SL1BP = 1.D0/MG2*( -1.D0 -MS2/MG2*LOG(ABS(1.D0 -MG2/MS2)))
C***  FINITE PART OF B0P(KG=MG,0,MS)
ctp      print*, " SL1BP ",SL1BP
      RETURN
      END

      REAL*8 FUNCTION SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,S,XS,ZETA2,SPENCE,XG,XLAM,X1,X2
      COMPLEX*16 KXG,KXT,KSPENC
      INTEGER INUM
      IF (INUM.EQ.1) SL1C0A = 
     +     1.D0/S*( 0.5D0*LOG(S/MS2)**2 - 3.5D0*ZETA2)         
C***  FINITE PART OF C0(K1,K2,0,0,0)
      IF (INUM.EQ.2) SL1C0A = 
     +     1.D0/S*(0.5D0*LOG(XS)**2 - 3.D0*ZETA2)
C***  FINITE PART OF C0(K1,K2,MS,MS,MS)
      IF (INUM.EQ.3) SL1C0A = 1.D0/S * 
     +     REAL( LOG(-KXG)**2 )/2.D0
C***  FINITE PART OF C0(K1,K2,MG,MG,MG)
      IF (INUM.EQ.4) SL1C0A = 1.D0/S * 
     +     REAL( LOG(-KXT)**2 )/2.D0
C***  FINITE PART OF C0(K1,K2,MT,MT,MT)
      IF (INUM.EQ.5) SL1C0A = 1.D0/S * (
     +     + LOG(XS)**2 -2.D0*SPENCE(1.D0 -MG2/MS2) -6.D0*ZETA2
     +     + SPENCE(1.D0 + MG2/MS2*XS) +SPENCE(1.D0 +MG2/MS2/XS) )
C***  FINITE PART OF C0(K1,K2,MS,MG,MS)
      IF (INUM.EQ.6) THEN
         IF (MG2.GE.S/4) THEN
            SL1C0A = (REAL( LOG(-KXG)**2 +KSPENC(1.D0 +MS2/MG2*KXG) 
     +           + KSPENC(1.D0 +MS2/MG2/KXG))
     +           -2.D0*SPENCE(1.D0 -MS2/MG2)) /S
         ELSE
            XG = REAL(KXG)
            SL1C0A = (LOG(XG)**2 -6.D0*ZETA2 +SPENCE(1.D0+MS2/MG2*XG) 
     +           +SPENCE(1.D0+MS2/MG2/XG) 
     +           -2.D0*SPENCE(1.D0 -MS2/MG2)) /S
         END IF
C***  FINITE PART OF C0(K1,K2,MG,MS,MG)         
      END IF
      IF (INUM.EQ.7) THEN
         XLAM = SQRT((S -MS2 -MG2)**2 -4*MS2*MG2)
         X1 = (1.D0 + MG2/S - MS2/S)/2.D0 + XLAM/S/2.D0
         X2 = (1.D0 + MG2/S - MS2/S)/2.D0 - XLAM/S/2.D0
         SL1C0A = 1.D0/S*(    SPENCE((MG2-MS2)/MG2) 
     +        -SPENCE(1.D0/X1) -SPENCE(1.D0/X2) ) 
C***  FINITE PART OF C0(K1,K2,MG,MG,MS)
      END IF
      IF (INUM.EQ.8) THEN
         XLAM = SQRT((S -MS2 -MG2)**2 -4*MS2*MG2)
         X1 = (1.D0 + MS2/S - MG2/S)/2.D0 + XLAM/S/2.D0
         X2 = (1.D0 + MS2/S - MG2/S)/2.D0 - XLAM/S/2.D0
         SL1C0A = 1.D0/S*(    SPENCE((MS2-MG2)/MS2) 
     +        -SPENCE(1.D0/X1) -SPENCE(1.D0/X2) ) 
C***  FINITE PART OF C0(K1,K2,MS,MS,MG)
      END IF
ctp      print*, " SL1C0A ",SL1C0A
      RETURN
      END


      REAL*8 FUNCTION SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,S,SB,XS,ZETA2,SPENCE,C0GENS
      REAL*8 BETA,X1S,X2S,XWG,X1G,X2G,XG
      COMPLEX*16 KXWG,KX1G,KX2G,KXG,KSPENC
      INTEGER INUM
      IF (INUM.EQ.1) SL1C0B = 1.D0/SB *(
     +      - 2.D0*SPENCE(XS) - 2.D0*LOG(XS)*LOG(1-XS)
     +     + LOG(XS)**2/2.D0 - 4.D0*ZETA2 )
C***  FINITE PART OF C0(P1,-P2,MS,0,MS)

      IF (INUM.EQ.2) THEN
         BETA = SQRT( 1.D0 - 4.D0*MS2/S)
         X1S = 1.D0/2.D0*( 1.D0 + BETA)
         X2S = 1.D0/2.D0*( 1.D0 - BETA)
         KXWG = SQRT( DCMPLX(1.D0 - 4.D0*MG2/S)) 
         IF (MG2.LE.S/4.D0) THEN
            XWG = REAL(KXWG)
            X1G = 1.D0/2.D0*( 1.D0 + XWG)
            X2G = 1.D0/2.D0*( 1.D0 - XWG)
            XG = X2G/X1G
            SL1C0B = 2.D0/SB*( SPENCE(X2S/(X1G -X1S))
     +           - SPENCE(X1S/(X1S-X1G)) +SPENCE(X2S/(X2G-X1S)) 
     +           - SPENCE(X1S/(X1S-X2G)))
         ELSE
            KX1G = 1.D0/2.D0*( 1.D0 + KXWG)
            KX2G = 1.D0/2.D0*( 1.D0 - KXWG)
            KXG = KX2G/KX1G
            SL1C0B = REAL(2.D0/SB*( KSPENC(X2S/(KX1G -X1S))
     +           - KSPENC(X1S/(X1S-KX1G)) +KSPENC(X2S/(KX2G-X1S))    
     +           - KSPENC(X1S/(X1S-KX2G))))
         END IF
C***  FINITE PART OF C0(-P1,-P2,MG,0,MG)
      END IF

      IF (INUM.EQ.3) 
     +     SL1C0B= 1.D0/SB*(2.D0*SPENCE(-XS) + 0.5D0*LOG(XS)**2 +ZETA2)
C***  FINITE PART OF C0(-P1,-P2,0,MS,0)

      IF (INUM.EQ.4) SL1C0B = 1.D0/SB *(
     +     -2.D0*SPENCE(1.D0 -(1.D0 -MG2/MS2)/(1.D0 +XS))
     +     +2.D0*SPENCE(1.D0 -(1.D0 -MG2/MS2)/(1.D0 +1.D0/XS))
     +     -SPENCE(1.D0 +(MG2/MS2)*XS)
     +     +SPENCE(1.D0 +(MG2/MS2)/XS)
     +     + LOG(XS)**2 -2.D0*LOG(XS)*LOG(1.D0 +XS) )
C***  FINITE PART OF C0(-P1,-P2,0,MG,0)
      IF (INUM.EQ.5) SL1C0B = 
     +     C0GENS(S,MS2,MG2,MT2)
C***  FINITE PART OF C0(-P1,-P2,MT,MG,MT)
      IF (INUM.EQ.6) SL1C0B = 
     +     C0GENS(S,MS2,MT2,MG2)
C***  FINITE PART OF C0(-P1,-P2,MG,MT,MG)
ctp      print*, " SL1C0B ",SL1C0B
      RETURN
      END


      REAL*8 FUNCTION SL1C0C(MS2,MG2,MT2,T,ZETA2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,T,T1,SPENCE,ZETA2,A,B,C0GENT
      INTEGER INUM
      T1 = T -MS2
      IF (INUM.EQ.1) SL1C0C = 1.D0/T1*(ZETA2 - SPENCE(T/MS2))
C***  FINITE PART OF C0(K2,-P2,MS,MS,0)
      IF (INUM.EQ.2) SL1C0C = 1.D0/T1*( SPENCE(MS2/MG2) -SPENCE(T/MG2))
C***  FINITE PART OF C0(K2,-P2,MG,MG,0)
      IF (INUM.EQ.3) SL1C0C = 
     +     1.D0/T1*(+ 0.25D0*ZETA2 + LOG(-T1/MS2)**2 +SPENCE(T/MS2))
C***  FINITE PART OF C0(K2,-P2,0,0,MS)
      IF (INUM.EQ.4) THEN
C***  FINITE PART OF C0(K2,-P2,0,0,MG)
         A = LOG(ABS((MG2-MS2)/(MG2 -T)))
         IF (MG2.GE.MS2) THEN
            B = A*LOG((MG2-MS2)*(MG2 -T)/MG2/MS2)
         ELSE
            B = A*LOG((MS2-MG2)*(MG2 -T)/MG2/MS2) - 6.D0*ZETA2 
         END IF
         SL1C0C = 1.D0/T1*(-B + SPENCE(T/MG2) - SPENCE(MS2/MG2) )
      END IF
      IF (INUM.EQ.5) SL1C0C = 
     +     1.D0/T1*(ZETA2 -SPENCE(T/MG2) -SPENCE(1.D0-MS2/MG2) 
     +     +LOG(MG2/MS2)*LOG(1.D0 -T/MG2) )
C***  FINITE PART OF C0(K2,-P2,MG,MS,0)
      IF (INUM.EQ.6) SL1C0C = 
     +     1.D0/T1*(ZETA2 - SPENCE(T/MS2) + SPENCE(1.D0 +(MG2-MS2)/T1)
     +     -SPENCE(1.D0 +MS2/MG2*(MG2-MS2)/T1) -0.5D0*LOG(MG2/MS2)**2
     +     +LOG(MG2/MS2)*LOG(ABS((MS2-MG2)/T1)))
C***  FINITE PART OF C0(K2,-P2,MS,MG,0)
      IF (INUM.EQ.7) SL1C0C = 
     +     C0GENT(T,MS2,MG2,MT2)
C***  FINITE PART OF C0(K2,-P2,MT,MT,MG)
      IF (INUM.EQ.8) SL1C0C = 
     +     C0GENT(T,MS2,MT2,MG2)
C***  FINITE PART OF C0(K2,-P2,MG,MG,MT)

ctp      print*, " SL1C0C ",SL1C0C
      RETURN
      END


      REAL*8 FUNCTION SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,M2,S,T,U,SB,XS,SPENCE,ZETA2,A,T1,TG,U1,UG
      REAL*8 MG,MT,RD0REG
      COMPLEX*16 KDEL,KY0,KOMS,KBETAG,KXGG,KROOTX,KROOTXB,KXB
      COMPLEX*16 KY0TERM,KTP,KTM,KTERM1,KTERM2,KTERM3,KTERM4,KSPENC
      INTEGER INUM
      M2 = MG2 -MS2
      T1 = T -MS2
      TG = T -MG2
      U1 = U -MS2
      UG = U -MG2
      IF (INUM.EQ.1) SL1D0 = 
     +     1.D0/S/T1*(2.D0*LOG(S/MS2)*LOG(-T1/MS2) -4.D0*ZETA2)
C***  FINITE PART OF D0(K1,K2,-P2,0,0,0,MS)

      IF (INUM.EQ.2) SL1D0 = 1.D0/SB/T1 *(
     +     -2.D0*LOG(XS)*LOG(1.D0-XS) +2.D0*LOG(XS)*LOG(1.D0+XS)
     +     -2.D0*LOG(XS)*LOG(-T1/MS2) -2.D0*SPENCE(XS)
     +     +2.D0*SPENCE(-XS) -3.D0*ZETA2 )
C***  FINITE PART OF D0(K1,K2,-P2,MS,MS,MS,0)

      IF (INUM.EQ.3) SL1D0 = 1.D0/T1/U1 *(
     +     2.D0*LOG(-T1/MS2)*LOG(-U1/MS2) -7.D0/2.D0*ZETA2 )
C***  FINITE PART OF D0(K1,-P2,K2,MS,MS,0,0)

      IF (INUM.EQ.4) THEN 
C***  FINITE PART OF D0(K1,K2,-P2,0,0,0,MG)
         IF (MG2.GE.MS2) A = -2.D0*LOG(M2/MS2)*LOG(M2/MG2)
         IF (MG2.LT.MS2) A =
     +        -2.D0*LOG(-M2/MS2)*LOG(-M2/MG2) +12.D0*ZETA2
         SL1D0 = 1.D0/S/TG *( A 
     +    -4.D0*SPENCE(1.D0+M2/TG) -1.5D0*ZETA2 +0.5D0*LOG(S/MS2)**2
     +    -0.5D0*LOG(S/MG2)**2 +2.D0*LOG(S/MS2)*LOG(-TG/MG2)
     +    -SPENCE(1.D0 + M2**2/S/MG2) )
      END IF

      IF (INUM.EQ.5) THEN
         KDEL = (0.D0,1.D-12)
         KY0=(MG2-T-KDEL)/(MG2-MS2-KDEL)
         KOMS = S +KDEL
         KBETAG = SQRT(1.D0 -4.D0*MG2/KOMS)
         KXGG = (KBETAG -1.D0)/(KBETAG +1.D0)
         KROOTX = 1D0 - 8D0*MG2/KY0/(4D0*MG2-KOMS)
     +        + 4D0*MG2/KY0**2/(4D0*MG2-KOMS)
         KROOTXB = SQRT(KROOTX)
         KXB = ( - 8D0*MG2/KY0/(4D0*MG2-KOMS)
     +        + 4D0*MG2/KY0**2/(4D0*MG2-KOMS) )/(KROOTXB+1D0)**2
         KY0TERM = 2D0*KY0 - 1D0
         KTP = 1D0+SQRT(KY0TERM)*SQRT(-KXB)
         KTM = 1D0+SQRT(KY0TERM)/SQRT(-KXB)
         KTERM1 = 1D0+1D0/SQRT(KXGG)
         KTERM3 = 1D0+SQRT(KXGG)
         SL1D0 = -2D0*REAL( ( LOG(KXGG)*( LOG(KTM) - LOG(KTP) )
     +        + 2D0*KSPENC(1.D0 -KTERM3/KTP) 
     +        + 2D0*KSPENC(1.D0 -KTERM1/KTM)
     +        - 2D0*KSPENC(1.D0 -KTERM3/KTM) 
     +        - 2D0*KSPENC(1.D0 -KTERM1/KTP) )
     +        /S/(T-MG2)/KBETAG/KROOTXB )
C***  FINITE PART OF D0(K1,K2,-P2,MG,MG,MG,0)
      END IF

      IF (INUM.EQ.6) THEN
         IF (MG2.GE.MS2) THEN
            SL1D0 = -1.D0/(TG*UG -M2**2) *( 
     +           LOG(M2**2/MG2/MS2)*(LOG(-M2/TG) +LOG(-M2/UG))
     +           + LOG(-M2/TG)**2 +LOG(-M2/UG)**2
     +           +4.D0*SPENCE(1.D0 +M2/TG) + 4.D0*SPENCE(1.D0 +M2/UG)
     +           + 2.D0*SPENCE(1.D0 -TG*UG/M2**2) )
         ELSE
            SL1D0 = -1.D0/(TG*UG -M2**2) *( 
     +           LOG(M2**2/MG2/MS2)*(LOG(M2/TG) +LOG(M2/UG))
     +           + LOG(M2/TG)**2 +LOG(M2/UG)**2 -12.D0*ZETA2
     +           +4.D0*SPENCE(1.D0 +M2/TG) + 4.D0*SPENCE(1.D0 +M2/UG)
     +           + 2.D0*SPENCE(1.D0 -TG*UG/M2**2) )
         END IF
C***  FINITE PART OF D0(K1,-P2,K2,0,0,MG,MG)
      END IF
       
      IF (INUM.EQ.7) SL1D0 = 1.D0/SB/TG *(
     +     -(2.D0*SPENCE(XS) + SPENCE(1.D0 +XS*MG2/MS2)
     +     +SPENCE(1.D0 +XS*MS2/MG2) +ZETA2 +0.5D0*LOG(MG2/MS2)**2
     +     +2.D0*LOG(XS)*(LOG(1.D0 -XS) +0.5D0*LOG(TG**2/MS2/MG2))) )
C***  FINITE PART OF D0(K1,K2,-P2,MS,MG,MS,0)

      IF (INUM.EQ.8) SL1D0 = 1.D0/TG/U1 *(
     +     ( LOG(-U1/MS2)**2 +2.D0*LOG(-U1/MS2)*LOG(ABS(TG/M2))
     +     + SPENCE(1.D0 +U1/M2) +SPENCE(1.D0 +U1*MG2/MS2/M2)
     +     -2.D0*SPENCE(1.D0 + M2/TG) -0.75D0*ZETA2 ) )
C***     FINITE PART OF D0(K1,-P2,K2,MS,MG,0,0)

      IF (INUM.EQ.9) THEN
         KDEL = DCMPLX(0.D0,1.D-12)
         KY0 = -(T1 +KDEL)/(MG2 -MS2 -KDEL)
         KOMS = S +KDEL
         KBETAG = SQRT(1.D0 -4.D0*MG2/KOMS)
         KXGG = (KBETAG -1.D0)/(KBETAG +1.D0)
         KROOTXB = 1D0 - 4D0*(MS2+MG2)/KY0/(4D0*MG2-KOMS)
     +        + 4D0*MS2/KY0**2/(4D0*MG2-KOMS)
         KROOTXB = SQRT(KROOTXB)
         KXB = ( - 4D0*(MS2+MG2)/KY0/(4D0*MG2-KOMS)
     +        + 4D0*MS2/KY0**2/(4D0*MG2-KOMS) )/(KROOTXB+1D0)**2
         KY0TERM = KY0*(1D0+MG2/MS2) - 1D0
         KTP = 1D0+SQRT(KY0TERM)*SQRT(-KXB)
         KTM = 1D0+SQRT(KY0TERM)/SQRT(-KXB)
         KTERM1 = 1D0+SQRT(MG2/MS2)/SQRT(KXGG)
         KTERM2 = 1D0+SQRT(MS2/MG2)/SQRT(KXGG)
         KTERM3 = 1D0+SQRT(MG2/MS2)*SQRT(KXGG)
         KTERM4 = 1D0+SQRT(MS2/MG2)*SQRT(KXGG)
         SL1D0 = -2D0*REAL( ( LOG(KXGG)*( LOG(KTM) - LOG(KTP) )
     +        + KSPENC(1.D0 -KTERM3/KTP) + KSPENC(1.D0 -KTERM4/KTP)
     +        + KSPENC(1.D0 -KTERM1/KTM) + KSPENC(1.D0 -KTERM2/KTM)
     +        - KSPENC(1.D0 -KTERM3/KTM) - KSPENC(1.D0 -KTERM4/KTM)
     +        - KSPENC(1.D0 -KTERM1/KTP) - KSPENC(1.D0 -KTERM2/KTP) )
     +        /S/T1/KBETAG/KROOTXB )
C***  FINITE PART OF D0(K1,K2,-P2,MG,MS,MG,0)
      END IF
      IF (INUM.EQ.10) THEN
C***  FINITE PART OF D0(K1,K2,-P2,MT,MT,MT,MG)
         MT = SQRT(MT2)
         MG = SQRT(MG2)
         SL1D0 = RD0REG(MS2,0.D0,0.D0,MS2,S,T,MG,MT,MT,MT)
      END IF
      IF (INUM.EQ.11) THEN
C***  FINITE PART OF D0(K1,K2,-P2,MG,MG,MG,MT)
         MT = SQRT(MT2)
         MG = SQRT(MG2)
         SL1D0 = RD0REG(MS2,0.D0,0.D0,MS2,S,T,MT,MG,MG,MG)
      END IF
      IF (INUM.EQ.12) THEN
C***  FINITE PART OF D0(K1,-P2,K2,MG,MG,MT,MT)
         MT = SQRT(MT2)
         MG = SQRT(MG2)
         SL1D0 =  RD0REG(0.D0,MS2,0.D0,MS2,U,T,MG,MG,MT,MT)
      END IF

ctp      print*, " SL1D0 ",SL1D0
      RETURN
      END


      REAL*8 FUNCTION SL2C0B(MS2,MG2,MT2,S,ZETA2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,S,ZETA2,SPENCE
      REAL*8 MS,MG,BE,XSW,XSG,A1,XLAM,CCDUKE
      COMPLEX*16 KDEL,KLAM,KX1,KX2,KY1,KY2
      INTEGER INUM
      IF (INUM.EQ.1) THEN
C***  FINITE PART OF C0(-P1=MS,-P2=MG,MS,0,MG)
         MS = SQRT(MS2)
         MG = SQRT(MG2)
         BE = SQRT(1.D0 -4*MS*MG/(S -(MS-MG)**2))
         XSW = (BE -1.D0)/(BE +1.D0)
         A1 = LOG(ABS(XSW))*( LOG(MG/MS) -0.5D0*LOG(ABS(XSW)) 
     +        +2*LOG(1.D0-XSW**2) ) 
     +        + 2*ZETA2 + 0.5D0*LOG(MS/MG)**2 + SPENCE(XSW**2)
     +        + SPENCE(1.D0 -MS/MG*XSW) + SPENCE(1.D0 -MG/MS*XSW)
         SL2C0B = A1 * XSW/MS/MG/(1.D0-XSW**2) 
      END IF
      IF (INUM.EQ.2) THEN
C***  FINITE PART OF C0(-P1=MS,-P2=MG,MG,0,MS)
         KDEL = (0.D0,1.D-15)
         XLAM = SQRT((S -MS2 -MG2)**2 -4*MS2*MG2)
         KLAM = SQRT((S -MS2 -MG2)**2 -4*MS2*MG2 +KDEL)
         KX1 = (1.D0 + MS2/S - MG2/S)/2.D0 + KLAM/S/2.D0
         KX2 = (1.D0 + MS2/S - MG2/S)/2.D0 - KLAM/S/2.D0
         KY1 = (1.D0 + MG2/S - MS2/S)/2.D0 + KLAM/S/2.D0
         KY2 = (1.D0 + MG2/S - MS2/S)/2.D0 - KLAM/S/2.D0
         SL2C0B = 1.D0/XLAM*(
     +        +CCDUKE(-KY1,1.D0,MS2/S-MG2/S-KDEL,2*MG2/S-2*MS2/S) 
     +        -CCDUKE(-KY2,1.D0,MS2/S-MG2/S-KDEL,2*MG2/S-2*MS2/S) 
     +        -CCDUKE(-KY1,1.D0,-KX1,1.D0) 
     +        +CCDUKE(-KY2,1.D0,-KX1,1.D0) 
     +        -CCDUKE(-KY1,1.D0,-KX2,1.D0) 
     +        +CCDUKE(-KY2,1.D0,-KX2,1.D0)  )
      END IF
      IF (INUM.EQ.3) THEN
*     FINITE PART OF C0(-P1=MS,-P2=MG,0,MS,0)
         MS = SQRT(MS2)
         MG = SQRT(MG2)
         XLAM = SQRT((S -MS2 -MG2)**2 -4*MS2*MG2)
         BE = SQRT(1.D0 -4*MG*MS/(S-(MG-MS)**2))
         XSG = (1.D0 -BE)/(1.D0 +BE)
c         print*, " test 1 ",xlam,be,xsg
      SL2C0B = 1D0/XLAM*( 
     +        -SPENCE(MG2/MS2) +2*SPENCE(-MG*XSG/MS) 
     +        +2*SPENCE(1 -(1-MG2/MS2)/(1+XSG*MG/MS)) 
     +        +0.5D0*( LOG(1+MG*XSG/MS) -LOG(1 +MS*XSG/MG) 
     +                +LOG(MS*XSG/MG) )**2 )
c      print*, " test 2 ",SPENCE(MG2/MS2)
c      print*, " test 3 ",SPENCE(-MG*XSG/MS)
c      print*, " test 4 ",SPENCE(1 -(1-MG2/MS2))
c      print*, " test 5 ",LOG(1 +MG*XSG/MS),1 +MG*XSG/MS
c      print*, " test 6 ",LOG(1 +MS*XSG/MG),1 +MS*XSG/MG
c      print*, " test 7 ",LOG(MS*XSG/MG),MS*XSG/MG
      END IF
      IF (INUM.EQ.4) THEN
C***  FINITE PART OF C0(-P1,-P2,0,MG,0)
         MS = SQRT(MS2)
         MG = SQRT(MG2)
         XLAM = SQRT((S -MS2 -MG2)**2 -4*MS2*MG2)
         BE = SQRT(1.D0 -4*MG*MS/(S-(MG-MS)**2))
         XSG = (1.D0 -BE)/(1.D0 +BE)
      SL2C0B = 1D0/XLAM*( 
     +        -SPENCE(MS2/MG2) +2*SPENCE(-MS*XSG/MG) 
     +        +2*SPENCE(1 -(1-MS2/MG2)/(1+XSG*MS/MG)) 
     +        +0.5D0*( LOG(1+MS*XSG/MG) -LOG(1 +MG*XSG/MS) 
     +                +LOG(MG*XSG/MS) )**2 )
      END IF
ctp      print*, " SL2C0B ",SL2C0B,MS2,MG2,MT2,S,ZETA2,INUM
      RETURN
      END


      REAL*8 FUNCTION SL2C0C(MS2,MG2,MT2,T,ZETA2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,T,T1,SPENCE,ZETA2,A,B
      INTEGER INUM
      T1 = T - MS2
      IF (INUM.EQ.1) SL2C0C = 1.D0/T1*(ZETA2 - SPENCE(T/MS2))
C***  FINITE PART OF C0(K1,-P1=MS,MS,MS,0)
      IF (INUM.EQ.2) SL2C0C = 
     +     1.D0/T1*( SPENCE(MS2/MG2) -SPENCE(T/MG2))
C***  FINITE PART OF C0(K2,-P1=MS,MG,MG,0)
      IF (INUM.EQ.3) SL2C0C = 
     + 1.D0/T1*(+0.25D0*ZETA2 +LOG(-T1/MS2)**2 +SPENCE(T/MS2))
C***  FINITE PART OF C0(K1,-P1=MS,0,0,MS)
      IF (INUM.EQ.4) THEN
         A = LOG(ABS((MG2-MS2)/(MG2 -T)))
         IF (MG2.GE.MS2) THEN
            B = A*LOG((MG2-MS2)*(MG2 -T)/MG2/MS2)
         ELSE
            B = A*LOG((MS2-MG2)*(MG2 -T)/MG2/MS2) - 6.D0*ZETA2 
         END IF
         SL2C0C = 1.D0/T1*( -B  +SPENCE(T/MG2) -SPENCE(MS2/MG2))
C***  FINITE PART OF C0(K1,-P1=MS,0,0,MG)
      END IF
      IF (INUM.EQ.5) SL2C0C = 
     +     (ZETA2 -SPENCE(T/MG2) -SPENCE(1.D0-MS2/MG2) 
     +     +LOG(MG2/MS2)*LOG(1.D0 -T/MG2))/T1 
C***  FINITE PART OF C0(K1,-P1=MS,MG,MS,0)
      IF (INUM.EQ.6) SL2C0C = 
     +     (ZETA2 - SPENCE(T/MS2) + SPENCE(1.D0 +(MG2-MS2)/T1)
     +     -SPENCE(1.D0 +MS2/MG2*(MG2-MS2)/T1) -0.5D0*LOG(MG2/MS2)**2
     +     +LOG(MG2/MS2)*LOG(ABS((MS2-MG2)/T1)))/T1
C***  FINITE PART OF C0(K1,-P1,MS,MG,0)
ctp      print*, " SL2C0C ",SL2C0C
      RETURN
      END

      REAL*8 FUNCTION SL2C0D(MS2,MG2,MT2,T,ZETA2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,T,TG,SPENCE,ZETA2,A,B
      REAL*8 WG,W1,W2,W3,Y1,Y2,Z1,Z2
      COMPLEX*16 KA1,KA2,KB1,KB2,KSPENC
      INTEGER INUM
      TG = T -MG2
      IF (INUM.EQ.1) SL2C0D = 1.D0/TG*(ZETA2 - SPENCE(T/MG2))
C***  FINITE PART OF C0(K2,-P2=MG,MG,MG,0)
      IF (INUM.EQ.2) SL2C0D = 1.D0/TG*( SPENCE(MG2/MS2) -SPENCE(T/MS2))
C***  FINITE PART OF C0(K2,-P2=MG,MS,MS,0)
      IF (INUM.EQ.3) SL2C0D = 
     +     1.D0/TG*(+ 0.25D0*ZETA2 + LOG(-TG/MG2)**2 +SPENCE(T/MG2)
     +     -LOG(MS2/MG2)*LOG(-TG/MG2) +0.25D0*LOG(MS2/MG2)**2 )
C***  FINITE PART OF C0(K2,-P2=MG,0,0,MG)
      IF (INUM.EQ.4) THEN
         A = LOG(ABS((MS2-MG2)/(MS2 -T)))
         IF (MS2.GE.MG2) THEN
            B = A*LOG((MS2-MG2)*(MS2 -T)/MS2/MG2)
         ELSE
            B = A*LOG((MG2-MS2)*(MS2 -T)/MS2/MG2) - 6.D0*ZETA2 
         END IF
         SL2C0D = 1.D0/TG*
     +     ( -B +SPENCE(T/MS2) -SPENCE(MG2/MS2) +A*LOG(MS2/MG2))
C***  FINITE PART OF C0(K2,-P2=MG,0,0,MS)
      END IF
      IF (INUM.EQ.5) SL2C0D = 
     +     (ZETA2 -SPENCE(T/MS2) -SPENCE(1.D0-MG2/MS2) 
     +     +LOG(MS2/MG2)*LOG(1.D0 -T/MS2))/TG
C***  FINITE PART OF C0(K2,-P2=MG,MS,MG,0)

      IF (INUM.EQ.6) SL2C0D = 
     +     (ZETA2 - SPENCE(T/MG2) + SPENCE(1.D0 +(MS2-MG2)/TG)
     +     -SPENCE(1.D0 +MG2/MS2*(MS2-MG2)/TG) -0.5D0*LOG(MS2/MG2)**2
     +     +LOG(MS2/MG2)*LOG(ABS((MG2-MS2)/TG)))/TG
C***  FINITE PART OF C0(K2,-P2=MG,MG,MS,0)
      IF (INUM.EQ.7) THEN
         WG = - MG2/TG
         W1 = - MT2/TG
         W2 = - MS2/TG
         W3 = W2 -W1 -WG
         Y1 = (W3 +1.D0)/(WG -1.D0)
         Y2 = W1/(WG -1.D0)
         Z1 = W3/WG
         Z2 = W1/WG
         KA1 = - Y1/2.D0 + SQRT(DCMPLX(Y1**2/4.D0 - Y2))
         KA2 = - Y1/2.D0 - SQRT(DCMPLX(Y1**2/4.D0 - Y2))
         KB1 = - Z1/2.D0 + SQRT(DCMPLX(Z1**2/4.D0 - Z2))
         KB2 = - Z1/2.D0 - SQRT(DCMPLX(Z1**2/4.D0 - Z2))
         SL2C0D = 1/TG*REAL( - KSPENC(1.D0/KA1) - KSPENC(1.D0/KA2) 
     +        + KSPENC(1.D0/KB1) + KSPENC(1.D0/KB2) )
C***  FINITE PART OF C0(K2,-P2=MG,MT,MT,MS)
      END IF
      IF (INUM.EQ.8) THEN
         WG = - MG2/TG
         W1 = - MS2/TG
         W2 = - MT2/TG
         W3 = W2 -W1 -WG
         Y1 = (W3 +1.D0)/(WG -1.D0)
         Y2 = W1/(WG -1.D0)
         Z1 = W3/WG
         Z2 = W1/WG
         KA1 = - Y1/2.D0 + SQRT(DCMPLX(Y1**2/4.D0 - Y2))
         KA2 = - Y1/2.D0 - SQRT(DCMPLX(Y1**2/4.D0 - Y2))
         KB1 = - Z1/2.D0 + SQRT(DCMPLX(Z1**2/4.D0 - Z2))
         KB2 = - Z1/2.D0 - SQRT(DCMPLX(Z1**2/4.D0 - Z2))
         SL2C0D = 1/TG*REAL( - KSPENC(1.D0/KA1) - KSPENC(1.D0/KA2) 
     +             + KSPENC(1.D0/KB1) + KSPENC(1.D0/KB2) )
C***  FINITE PART OF C0(K2,-P2=MG,MS,MS,MT)
      END IF
ctp      print*, " SL2C0D ",SL2C0D
      RETURN
      END

      REAL*8 FUNCTION SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,INUM)
      IMPLICIT NONE
      REAL*8 MS2,MG2,MT2,M2,S,T,U,SPENCE,ZETA2,T1,TG,U1,UG
      REAL*8 MS,MG,CCDUKE,BE,XS
      COMPLEX*16 KDEL,KLAM,KY1,KY2,KHB,KHC,KZ1,KZ2
      INTEGER INUM
      KDEL = (0.D0,1.D-16)
      T1 = T -MS2
      TG = T -MG2
      U1 = U -MS2
      UG = U -MG2
      M2 = MG2 -MS2

      IF (INUM.EQ.1) SL2D0  = 1.D0/S/T1 * (
     +     -13.D0/4.D0*ZETA2
     +     -REAL(LOG((MS2 -MG2 +KDEL)/MS2)**2)
     +     +2.D0*LOG((MS2-T)/MS2)*LOG(S/MS2)
     +     -2*SPENCE((MG2-T)/(MS2-T)) )
C***  FINITE PART OF D0(K1,K2,-P2=MG,0,0,0,MS)
      IF (INUM.EQ.2) SL2D0  = 1.D0/S/TG * (
     +     +0.75D0*LOG(MS2/MG2)**2 
     +     + LOG(MS2/MG2)*LOG(ABS(MS2-MG2)/S)
     +     -2*LOG(MS2/MG2)*LOG((MG2-T)/MG2)
     +     -13.D0/4.D0*ZETA2
     +     -REAL(LOG((MG2 -MS2 +KDEL)/MG2)**2)
     +     +2.D0*LOG((MG2-T)/MG2)*LOG(S/MG2)
     +     -2*SPENCE((MS2-T)/(MG2-T)) )
      IF (INUM.EQ.3) THEN
         KLAM = SQRT((S -MG2 -MS2)**2 -4*MG2*MS2 +KDEL)
         KY1 = 0.5D0/S* (S +MG2 -MS2 + KLAM)
         KY2 = 0.5D0/S* (S +MG2 -MS2 - KLAM)
         KHB = (-S*T1 -M2*T1 -M2*S +KDEL)/S/T1
         KHC = (MS2*T1 -2*MS2*M2 -M2*U1 -KDEL)/S/T1
         KZ1 = 0.5D0*( -KHB + SQRT(KHB**2 -4*KHC))
         KZ2 = 0.5D0*( -KHB - SQRT(KHB**2 -4*KHC))
         SL2D0 = 1.D0/S/T1/REAL(KZ1 -KZ2) *(
     +        LOG(MG2/S)*REAL(LOG((1.D0 -KZ1)*KZ2/KZ1/(1.D0 -KZ2)))
     +        + CCDUKE(-KZ1,1D0,M2/S-KDEL,-2*M2/S)
     +        - CCDUKE(-KZ2,1D0,M2/S-KDEL,-2*M2/S)
     +        - CCDUKE(-KZ1,1D0,M2/S-KDEL,-T1/S)
     +        + CCDUKE(-KZ2,1D0,M2/S-KDEL,-T1/S)
     +        - CCDUKE(-KZ1,1D0,-KY1,1.D0)
     +        + CCDUKE(-KZ2,1D0,-KY1,1.D0)
     +        - CCDUKE(-KZ1,1D0,-KY2,1.D0)
     +        + CCDUKE(-KZ2,1D0,-KY2,1.D0) )
C***  FINITE PART OF D0(K1,K2,-P2=MG,MG,MG,MS,0)
      END IF
      IF (INUM.EQ.4) THEN
         KLAM = SQRT((S -MG2 -MS2)**2 -4*MG2*MS2 +KDEL)
         KY1 = 0.5D0/S* (S -MG2 +MS2 + KLAM)
         KY2 = 0.5D0/S* (S -MG2 +MS2 - KLAM)
         KHB = (-S*TG +M2*TG +M2*S +KDEL)/S/TG
         KHC = (MG2*TG +2*MG2*M2 +M2*UG -KDEL)/S/TG
         KZ1 = 0.5D0*( -KHB + SQRT(KHB**2 -4*KHC))
         KZ2 = 0.5D0*( -KHB - SQRT(KHB**2 -4*KHC))
         SL2D0 = 1.D0/S/TG/REAL(KZ1 -KZ2) * (
     +        LOG(MS2/S)*REAL(LOG((1.D0 -KZ1)*KZ2/KZ1/(1.D0 -KZ2)))
     +        + CCDUKE(-KZ1,1D0,-M2/S-KDEL,2*M2/S)
     +        - CCDUKE(-KZ2,1D0,-M2/S-KDEL,2*M2/S)
     +        - CCDUKE(-KZ1,1D0,-M2/S-KDEL,-TG/S)
     +        + CCDUKE(-KZ2,1D0,-M2/S-KDEL,-TG/S)
     +        - CCDUKE(-KZ1,1D0,-KY1,1.D0)
     +        + CCDUKE(-KZ2,1D0,-KY1,1.D0)
     +        - CCDUKE(-KZ1,1D0,-KY2,1.D0)
     +        + CCDUKE(-KZ2,1D0,-KY2,1.D0) )
C***  FINITE PART OF D0(K1,K2,-P2=MG,MG,MS,MS,0)
      END IF
      IF (INUM.EQ.5) THEN
         MS = SQRT(MS2)
         MG = SQRT(MG2)
         BE = SQRT(1.D0 -4*MS*MG/(S -(MS-MG)**2))
         XS = (BE -1.D0)/(BE +1.D0)
         SL2D0 = XS/MS/MG/(1.D0-XS**2)/(T -MS2)  * (
     +       LOG(ABS(XS))*( 2*LOG((MS2 -T)/MS2) +2*LOG(1.D0 -XS**2))
     +     - ZETA2 + 0.25D0*LOG(MS2/MG2)**2 + SPENCE(XS**2)
     +     + 2*SPENCE(1.D0 -MS/MG*XS) +2*SPENCE(1.D0 -MG/MS*XS)  )
C***  FINITE PART OF D0(K1,K2,-P2=MG,MS,MS,MG,0)
      END IF
      IF (INUM.EQ.6) THEN
         MS = SQRT(MS2)
         MG = SQRT(MG2)
         BE = SQRT(1.D0 -4*MS*MG/(S -(MG-MS)**2))
         XS = (BE -1.D0)/(BE +1.D0)
         SL2D0 =  XS/MS/MG/(1.D0-XS**2)/(T -MG2) * (
     +       LOG(ABS(XS))*( 2*LOG((MG2 -T)/MG2) +2*LOG(1.D0 -XS**2))
     +     -LOG(ABS(XS)) * LOG(MS2/MG2)
     +     - ZETA2 + 0.25D0*LOG(MG2/MS2)**2 + SPENCE(XS**2)
     +     + 2*SPENCE(1.D0 -MG/MS*XS) +2*SPENCE(1.D0 -MS/MG*XS) )
C***  FINITE PART OF D0(K1,K2,-P2=MG,MS,MG,MG,0)
      END IF

      IF (INUM.EQ.7) SL2D0 = 1.D0/(T1*UG +M2**2) * (
     +     0.5D0*LOG(MS2/MG2)*LOG(M2**2/T1/UG)
     +     +0.5D0*LOG(MS2/MG2)**2 
     +     + 0.5D0*LOG(MS2/MG2)*LOG(UG/T1)
     +     +2*LOG(-UG/MG2)*LOG(-T1/MS2) 
     +     -LOG(ABS(M2/MG2))**2 -LOG(ABS(M2/MS2))**2 +6*ZETA2
     +     -2*SPENCE(1.D0 +M2/UG) -2*SPENCE(1.D0 -M2/T1)
     +     -SPENCE(1.D0 -M2/UG) -SPENCE(1.D0 +M2/T1)
     +     -SPENCE(1.D0 -MG2/MS2*M2/UG) -SPENCE(1.D0 +MS2/MG2*M2/T1)
     +     +2*SPENCE(1.D0 +M2**2/T1/UG) )
C***  FINITE PART OF D0(K1,-P1=MS,K2,MS,MG,0,0)

      IF (INUM.EQ.8) SL2D0 = 1.D0/TG/U1 * (
     +     - 3.5D0*ZETA2 
     +     - 0.25D0*LOG(MS2/MG2)**2 
     +     + 2*LOG(-TG/SQRT(MG2)/SQRT(MS2))*LOG(-U1/MS2) )

C***  FINITE PART OF D0(K1,-P1=MS,K2,MG,MS,0,0)

      IF (INUM.EQ.9) SL2D0 = 1.D0/T1/U1 * (
     +     -0.75D0*ZETA2 + LOG(-T1/MS2)**2 
     +     +2*LOG(-T1/MS2)*LOG(-U1/ABS(MS2-MG2))
     +     +2*SPENCE((T -MG2)/(MS2-MG2)) -2*SPENCE((MG2-U)/(MS2-U)) )
C***  FINITE PART OF D0(K1,-P1=MS,K2,0,0,MS,MS)
      
      IF (INUM.EQ.10) SL2D0 = 1.D0/TG/UG * (
     +     +0.25D0*LOG(MS2/MG2)**2 
     +     + LOG(MS2/MG2) * LOG(ABS((MS2-MG2)*MG2/TG/UG))
     +     -0.75D0*ZETA2 + LOG(-TG/MG2)**2 
     +     +2*LOG(-TG/MG2)*LOG(-UG/ABS(MS2-MG2))
     +     +2*SPENCE((T -MS2)/(MG2-MS2)) -2*SPENCE((MS2-U)/(MG2-U)) )
C***  FINITE PART OF D0(K1,-P1=MS,K2,0,0,MG,MG)
ctp      print*, " SL2D0 ",SL2D0
      RETURN
      END

      REAL*8 FUNCTION C0GENS(S,MS2,MG2,MT2)
*     FINITE PART OF C0(-P1,-P2,MT,MG,MT)
*     FROM W. BEENAKKER, PHD THESIS
      REAL*8 A,B,C,D,E,F,BE,AL,S,MS2,MG2,MT2,X
      COMPLEX*16 Y,SUM,KSPENC
      INTEGER I,J
      DIMENSION X(1:3)
      DIMENSION Y(1:3,1:2)

      A = MS2
      B = S
      C = -S
      D = MT2 -MG2 -MS2
      E = 0.D0
      F = MG2
      BE = SQRT(1.D0 -4.D0*MS2/S)
      AL = (1.D0 -BE)/2.D0
      
      X(1) = - (D +E*AL)/(-S*BE) +AL
      X(2) = - (D +E*AL)/(1-AL)/(-S*BE)
      X(3) = + (D +E*AL)/(AL)/(-S*BE)
      Y(1,1) = (-C -E +SQRT(DCMPLX((C+E)**2 -4.D0*B*(A+D+F))))/2.D0/B
      Y(1,2) = (-C -E -SQRT(DCMPLX((C+E)**2 -4.D0*B*(A+D+F))))/2.D0/B
      Y(2,1) = (-D -E +
     +     SQRT(DCMPLX((D+E)**2 -4.D0*F*(A+B+C))))/2.D0/(A+B+C)
      Y(2,2) = (-D -E -
     +     SQRT(DCMPLX((D+E)**2 -4.D0*F*(A+B+C))))/2.D0/(A+B+C)
      Y(3,1) = (-D +SQRT(DCMPLX(D**2 -4.D0*A*F)))/2.D0/A
      Y(3,2) = (-D -SQRT(DCMPLX(D**2 -4.D0*A*F)))/2.D0/A
      
      SUM = DCMPLX(0.D0)
      DO 100, I =1,3
         DO 100 J=1,2
            SUM = SUM +(-1.D0)**I*
     +           ( + KSPENC(X(I)/(X(I) -Y(I,J))) 
     +             - KSPENC((X(I)-1.D0)/(X(I) -Y(I,J))) )

 100     CONTINUE

      C0GENS =  - REAL(SUM)/S/BE
      RETURN
      END


      DOUBLE PRECISION FUNCTION C0GENT(T,MS2,MG2,MT2)
*     * THIS IS C0(K2,-P2,MT,MT,MG)
      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      IMPLICIT COMPLEX*16 (K)
      REAL*8 T,T1,MS2,MG2,MT2,WG,W1,W2,W3,Y1,Y2
      COMPLEX*16 KA1,KA2,KB1,KB2,KSPENC
      T1 = T -MS2
      WG = - MS2/T1
      W1 = - MT2/T1
      W2 = - MG2/T1
      W3 = W2 -W1 -WG
      Y1 = (W3 +1.D0)/(WG -1.D0)
      Y2 = W1/(WG -1.D0)
      Z1 = W3/WG
      Z2 = W1/WG
      KA1 = - Y1/2.D0 + SQRT(DCMPLX(Y1**2/4.D0 - Y2))
      KA2 = - Y1/2.D0 - SQRT(DCMPLX(Y1**2/4.D0 - Y2))
      KB1 = - Z1/2.D0 + SQRT(DCMPLX(Z1**2/4.D0 - Z2))
      KB2 = - Z1/2.D0 - SQRT(DCMPLX(Z1**2/4.D0 - Z2))

      C0GENT = 1/T1*REAL( - KSPENC(1.D0/KA1) - KSPENC(1.D0/KA2) 
     +             + KSPENC(1.D0/KB1) + KSPENC(1.D0/KB2) )
      RETURN
      END

      REAL*8 FUNCTION RD0REG(P1,P2,P3,P4,P24,P13,M1,M2,M3,M4)
      REAL*8 P1,P2,P3,P4,P24,P13,M1,M2,M3,M4,MM2
      COMPLEX*16 D0REG
      MM2 = M2
      IF (P4.EQ.M1**2+M2**2) MM2 = M2 + 0.01D0
      IF (P4.EQ.M2**2+M3**2) MM2 = M2 + 0.01D0
      RD0REG = REAL(D0REG(P1,P2,P3,P4,P24,P13,M1,MM2,M3,M4))
      RETURN
      END




      REAL*8 FUNCTION A4P0P0()
      A4P0P0 = 2.D0
      RETURN
      END


      REAL*8 FUNCTION A4P1P2(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      IF((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X)).GT.0.D0)THEN
       A4P1P2 = 
     +     (2.D0*A*(BU**2+CU**2) -2.D0*B*AU*BU)/(AU**2-BU**2-CU**2)/X
     +     + B*(B*AU -A*BU)/X**(3.D0/2.D0) 
     +     *LOG((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X)))
      ELSE
       A4P1P2 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION A4P2P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      IF((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X)).GT.0.D0)THEN
      A4P2P1 = 
     +     2*B*(B*AU -A*BU)/(A**2 -B**2)/X
     +    + (A*(BU**2 +CU**2)-B*AU*BU)/X**(3.D0/2.D0)
     +    *LOG( (A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X)))
      ELSE
       A4P2P1 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION A4P1P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      IF((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X)).GT.0.D0)THEN
       A4P1P1 = 
     +  LOG((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X)))/SQRT(X) 
      ELSE
       A4P1P1 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION A4P1P0(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      IF((A +B)/(A-B).GT.0.D0)THEN
       A4P1P0 = 1.D0/B*LOG((A +B)/(A-B))
      ELSE
       A4P1P0 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION A4M1P0(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      A4M1P0 = 2.D0*A
      RETURN
      END


      REAL*8 FUNCTION A4M1P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X,Y
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      Y = BU**2 + CU**2
      IF((AU +SQRT(Y))/(AU -SQRT(Y)).GT.0.D0)THEN
       A4M1P1 = 
     +     2.D0*B*BU/Y
     +    + (A*Y -B*AU*BU)/Y**(3.D0/2.D0)
     +    * LOG((AU +SQRT(Y))/(AU -SQRT(Y)))
      ELSE
       A4M1P1 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION A4M2P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X,Y
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      Y = BU**2 + CU**2
      IF((AU +SQRT(Y))/(AU -SQRT(Y)).GT.0.D0)THEN
       A4M2P1 = 
     +     4.D0*A*B*BU/Y
     +    + B**2*AU*(CU**2 -2.D0*BU**2)/Y**2
     +    + ( (A*Y -B*AU*BU)**2/Y**(5.D0/2.D0)
     +    -B**2*CU**2*(AU**2 -BU**2 -CU**2)/2.D0/Y**(5.D0/2.D0) )
     +    * LOG((AU +SQRT(Y))/(AU -SQRT(Y)))
      ELSE
       A4M2P1 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION A4P2M2(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      IF((A + B)/(A-B).GT.0.D0)THEN
      A4P2M2 = 
     +     2.D0*(BU**2 -CU**2)/B**2
     +    + 2.D0*(B*AU -A*BU)**2/B**2/(A**2 -B**2)
     +    + (A*CU**2 +2.D0*BU*(B*AU -A*BU))/B**3
     +    * LOG((A + B)/(A-B))
      ELSE
       A4P2M2 = 0
      ENDIF
      RETURN
      END

      REAL*8 FUNCTION A4M1P2(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,Y
      Y = BU**2 +CU**2
      IF((AU +SQRT(Y))/(AU -SQRT(Y)).GT.0.D0)THEN
      A4M1P2 = 
     +     2.D0*(A*Y - B*AU*BU)/Y/(AU**2 -Y)
     +    + B*BU/Y**(3.D0/2.D0)*LOG((AU +SQRT(Y))/(AU -SQRT(Y)))
      ELSE
       A4M1P2 = 0
      ENDIF
      RETURN
      END

      REAL*8 FUNCTION A4P0P2(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      A4P0P2 = 2.D0/(AU**2 -BU**2 -CU**2 )
      RETURN
      END

      REAL*8 FUNCTION A4P2P2(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X,Y
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      Y = BU**2 +CU**2
      IF((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X)).GT.0.D0)THEN
      A4P2P2 = 
     +     2.D0*B**2/(A**2 -B**2)/X + 2.D0*Y/(AU**2 -Y)/X
     +     - 6.D0*B**2*CU**2/X**2
     +    +( B*BU/X**(3.D0/2.D0)
     +    + 3.D0*B*(B*AU -A*BU)*(A*Y -B*AU*BU)/X**(5.D0/2.D0))
     +    *LOG((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X)))
      ELSE
       A4P2P2 = 0
      ENDIF
      RETURN
      END

      REAL*8 FUNCTION AHP1P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      IF(2.D0*AU/(AU +BU).GT.0.D0)THEN
       AHP1P1 = -2.D0/A/(AU +BU)*LOG(2.D0*AU/(AU +BU)) 
      ELSE
       AHP1P1 = 0
      ENDIF
      RETURN
      END

      REAL*8 FUNCTION ABP2P0(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      ABP2P0 = -1.D0/A**2 
      RETURN
      END


      REAL*8 FUNCTION ABP2P2(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X,Y
      X = AU +BU
      Y = BU**2 +CU**2
      IF(X**2/(AU**2 -Y).GT.0.D0)THEN
      ABP2P2 = 1.D0/A**2/X**2 * (
     +     (3.D0*CU**2/X**2 + 2.D0*BU/X)
     +     * LOG(X**2/(AU**2 -Y))
     +     - 8.D0*CU**2/X**2 + 2.D0*Y/(AU**2 -Y) -1.D0 )
      ELSE
       ABP2P2 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION ABP1P2(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X,Y
      X = AU +BU
      Y = BU**2 +CU**2
      IF(X**2/(AU**2 -Y).GT.0.D0)THEN
      ABP1P2 =  1.D0/A/X**2
     +    * ( LOG(X**2/(AU**2 -Y)) +2.D0*(Y +AU*BU)/(AU**2 -Y) )
      ELSE
       ABP1P2 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION ABP2P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X,Y
      X = AU +BU
      Y = BU**2 +CU**2
      IF(X**2/(AU**2 -Y).GT.0.D0)THEN
      ABP2P1 = 1.D0/A**2/X
     +  * ((Y +AU*BU)/X**2*LOG(X**2/(AU**2-Y)) -2.D0*CU**2/X**2 -1.D0)
      ELSE
       ABP2P1 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION ABP1P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      IF((AU +BU)**2/(AU**2 -BU**2 -CU**2).GT.0.D0)THEN
       ABP1P1 = LOG((AU +BU)**2/(AU**2 -BU**2 -CU**2))/A/(AU +BU)
      ELSE
       ABP1P1 = 0
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION ABP1M1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      ABP1M1 = -2.D0*BU/A
      RETURN
      END

      REAL*8 FUNCTION ABP1M2(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      ABP1M2 = (CU**2 -4.D0*AU*BU -2.D0*BU**2)/A
      RETURN
      END


      REAL*8 FUNCTION C4P1P2(A,B,AU,BU,CU)
C***  THE REAL PART OF A COMPLEX INTEGRAL
      REAL*8 A,B,AU,BU,CU,X
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      C4P1P2 = 
     +     (2.D0*A*(BU**2+CU**2) -2.D0*B*AU*BU)/(AU**2-BU**2-CU**2)/X
     +     + B*(B*AU -A*BU)/X**(3.D0/2.D0) 
     +     *LOG(ABS((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X))))
      RETURN
      END


      REAL*8 FUNCTION C4M1P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X,Y
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      Y = BU**2 + CU**2
      C4M1P1 = 
     +     2.D0*B*BU/Y
     +    + (A*Y -B*AU*BU)/Y**(3.D0/2.D0)
     +    * LOG(ABS((AU +SQRT(Y))/(AU -SQRT(Y))))
      RETURN
      END

      REAL*8 FUNCTION C4M2P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X,Y
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      Y = BU**2 + CU**2
      C4M2P1 = 
     +     4.D0*A*B*BU/Y
     +    + B**2*AU*(CU**2 -2.D0*BU**2)/Y**2
     +    + ( (A*Y -B*AU*BU)**2/Y**(5.D0/2.D0)
     +    -B**2*CU**2*(AU**2 -BU**2 -CU**2)/2.D0/Y**(5.D0/2.D0) )
     +    * LOG(ABS((AU +SQRT(Y))/(AU -SQRT(Y)))) 
      RETURN
      END



      REAL*8 FUNCTION C4P1P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,X
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      C4P1P1 = 
     +  LOG(ABS((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X))))
     +     /SQRT(X) 
      RETURN
      END


      REAL*8 FUNCTION C4P1P0(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      C4P1P0 = 1.D0/B*LOG(ABS((A +B)/(A-B)))
      RETURN
      END

      REAL*8 FUNCTION CBP1P1(A,B,AU,BU,CU,DEL)
C***  HERE WE USE A WIDTH OF THE PARTICLE
      REAL*8 A,B,AU,BU,CU,DEL
      CBP1P1 = LOG(ABS((AU +BU)**2/(AU**2 -BU**2 -CU**2)))/A *
     +     (AU +BU)/((AU +BU)**2 +DEL)
      RETURN
      END

      REAL*8 FUNCTION K4P1P1(A,B,AU,BU,CU)
C***  THE IMAGINARY PART OF A COMPLEX INTEGRAL
      REAL*8 A,B,AU,BU,CU,X
      COMPLEX*16 KA,KDEL,KONE
      KDEL = (0.D0,1.D-12)
      KONE = (0.D0,-1.D0)
      KA = A - KDEL
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      K4P1P1 = REAL( KONE * 
     +     LOG((KA*AU -B*BU + SQRT(X))/(KA*AU -B*BU - SQRT(X)) )
     +     /SQRT(X) )
      RETURN
      END


      REAL*8 FUNCTION K4M1P1(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU,Y
      COMPLEX*16 KAU,KDEL,KONE
      KDEL = (0.D0,1.D-12)
      KONE = (0.D0,-1.D0)
      KAU = AU - KDEL
      Y = BU**2 + CU**2
      K4M1P1 = 
     +    + (A*Y -B*AU*BU)/Y**(3.D0/2.D0)
     +    * REAL( KONE * LOG((KAU +SQRT(Y))/(KAU -SQRT(Y))) )
      RETURN
      END


      REAL*8 FUNCTION K4P1P0(A,B,AU,BU,CU)
      REAL*8 A,B,AU,BU,CU
      COMPLEX*16 KA,KDEL,KONE
      KDEL = (0.D0,1.D-12)
      KONE = (0.D0,-1.D0)
      KA = A - KDEL
      K4P1P0 = 1.D0/B* REAL(KONE * LOG((KA +B)/(KA-B)))
      RETURN
      END





      REAL*8 FUNCTION CCDUKE(KA,B,KC,E)
*     THE INTEGRAL 3.13.1 FROM DEVOTO AND DUKE 
*     LOG[C + E Y]/(A + B Y),{Y,0,1}
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT COMPLEX*16 (K)
      A = REAL(KA)
      C = REAL(KC)
      CCDUKE = 1.D0/B * ( REAL(LOG(KC - A*E/B)*LOG((KA+B)/KA)) 
     +     -SPENCE(E*(A+B)/(A*E -B*C)) 
     +     +SPENCE(A*E/(A*E -B*C)) )

      RETURN
      END



      REAL*8 FUNCTION SPENCE(X)
*     -----------------------
      COMPLEX*16 KSPENC, Y, Z
      REAL*8 X,RSPENC
      IF (X.GE.1) THEN
      Y = DCMPLX(X)
      Z = KSPENC(Y)
      SPENCE = REAL(Z)
      ELSE
         SPENCE = RSPENC(X)
      ENDIF   
      END

      COMPLEX*16 FUNCTION KSPENC(X)
*     ----------------------------
* HANS KUIJF, 1988
* SPENCE(X) CALCS THE COMPLEX SPENCE-FUNCTION, THROUGH MAPPING ON
* THE AREA WHERE THERE IS A QUICKLY CONVERGENT SERIES.
      REAL*8 PI
      COMPLEX*16 X, SPEN
      PI=4.D0*ATAN(1.D0)
* MAP THE X ON THE UNIT CIRCLE.
* BUT SO THAT X IS NOT IN THE NEIGHBOURHOOD OF (1,0)
* ABS(Z)=-LOG(1D0-X) IS ALWAYS SMALLER THAN 1.10
* BUT (1.10)^19/(19!)*BERNOULLI(19)=2.7D-15
      IF (ABS(1.D0-X).LT.1.D-13) THEN
        KSPENC=PI*PI/6.D0
      ELSE IF (ABS(1.D0-X).LT.0.5D0) THEN
        KSPENC=PI*PI/6.D0-LOG(1.D0-X)*LOG(X)-SPEN(1.D0-X)
      ELSE IF (ABS(X).GT.1.D0) THEN
        KSPENC=-PI*PI/6.D0-0.5D0*LOG(-X)*LOG(-X)-SPEN(1.D0/X)
      ELSE
        KSPENC = SPEN(X)
      END IF
      END

      COMPLEX*16 FUNCTION SPEN(X)
*     ---------------------------
      COMPLEX*16 X,SUM,Z,Z2
      Z=-LOG(1D0-X )
      Z2=Z*Z
* HORNER'S RULE FOR THE POWERS Z^3 THROUGH Z^19
      SUM=43867.D0/798.D0
      SUM=SUM*Z2/342.D0-3617.D0/510.D0
      SUM=SUM*Z2/272.D0+7.D0/6.D0
      SUM=SUM*Z2/210.D0-691.D0/2730.D0
      SUM=SUM*Z2/156.D0+5.D0/66.D0
      SUM=SUM*Z2/110.D0-1.D0/30.D0
      SUM=SUM*Z2/ 72.D0+1.D0/42.D0
      SUM=SUM*Z2/ 42.D0-1.D0/30.D0
      SUM=SUM*Z2/ 20.D0+1.D0/6.D0
* THE FIRST THREE TERMS OF THE POWER SERIES
      SUM=Z2*Z*SUM/6.D0-0.25D0*Z2+Z
      SPEN=SUM
      END

      REAL*8 FUNCTION RSPENC(X)
*     ----------------------------
      REAL*8 PI, X, RSPEN
      PI=4.D0*ATAN(1.D0)
      IF (ABS(1.D0-X).LT.1.D-13) THEN
        RSPENC=PI*PI/6.D0
      ELSE IF (ABS(1.D0-X).LT.0.5D0) THEN
        RSPENC=PI*PI/6.D0-LOG(1.D0-X)*LOG(X)-RSPEN(1.D0-X)
      ELSE IF (ABS(X).GT.1.D0) THEN
        RSPENC=-PI*PI/6.D0-0.5D0*LOG(-X)*LOG(-X)-RSPEN(1.D0/X)
      ELSE
        RSPENC = RSPEN(X)
      END IF
      END

      REAL*8 FUNCTION RSPEN(X)
*     ---------------------------
      REAL*8 X,SUM,Z,Z2
      Z=-LOG(1D0-X )
      Z2=Z*Z
* HORNER'S RULE FOR THE POWERS Z^3 THROUGH Z^19
      SUM=43867.D0/798.D0
      SUM=SUM*Z2/342.D0-3617.D0/510.D0
      SUM=SUM*Z2/272.D0+7.D0/6.D0
      SUM=SUM*Z2/210.D0-691.D0/2730.D0
      SUM=SUM*Z2/156.D0+5.D0/66.D0
      SUM=SUM*Z2/110.D0-1.D0/30.D0
      SUM=SUM*Z2/ 72.D0+1.D0/42.D0
      SUM=SUM*Z2/ 42.D0-1.D0/30.D0
      SUM=SUM*Z2/ 20.D0+1.D0/6.D0
* THE FIRST THREE TERMS OF THE POWER SERIES
      SUM=Z2*Z*SUM/6.D0-0.25D0*Z2+Z
      RSPEN=SUM
      END


************************************************************************
*                                                                      *
*  GENERAL 4-POINT FUNCTION REGULAR CASE WITH 4 NON-ZERO MASSES        *
*                                                                      *
************************************************************************
* FUNCTIONS:                                                           *
* D0REG,CSPENC,CSPENH,CSPCOE,CDLN,  ETAE                               *
************************************************************************
      FUNCTION CDLN(CZ,EPS)
************************************************************************
*     COMPLEX LOGARITHM OF CZ + I*EPS                                  *
*----------------------------------------------------------------------*
*     09.01.90 ANSGAR DENNER                                           *
************************************************************************
      IMPLICIT   NONE
      REAL*8     EPS
      COMPLEX*16 CDLN,CZ
      REAL*8     PI,PI2_6
      COMMON /PI/     PI,PI2_6
 
      IF(AIMAG(CZ).EQ.0D0.AND.REAL(CZ).LE.0.D0)THEN
        IF (EPS.EQ.0D0) THEN
C$$$          WRITE(80,*) 'CDLN:  ARGUMENT ON CUT '
C$$$          WRITE(80,*) 'CDLN:  EPS = 0'
C$$$          WRITE(80,*) 'CDLN:  CZ  = ',CZ
        END IF
        CDLN=LOG(-CZ)+DCMPLX(0D0,PI)*SIGN(1D0,EPS)
      ELSE
        CDLN=LOG(CZ)
      END IF
      END
************************************************************************
      FUNCTION CSPENC(CZ,EPS)
************************************************************************
*       COMPLEX SPENCE FUNCTION  OF CZ + I*EPS                         *
*       CALCULATED BY MAPPING ON THE AREA WHERE THERE IS A QUICKLY     *
*       CONVERGENT SERIES                                              *
*----------------------------------------------------------------------*
*     08.01.90 ANSGAR DENNER                                           *
************************************************************************
      IMPLICIT   NONE
      REAL*8     EPS
      COMPLEX*16 CSPENC,CZ
      REAL*8     PI,PI2_6
      REAL*8     AZ,RZ,AZ1
      COMPLEX*16 CZ1,CSPENH,CDLN
      INTEGER    LTEST
 
      COMMON /PI/     PI,PI2_6
      COMMON /LTEST/  LTEST
 
C     PI      = 4D0*ATAN(1D0)
C     PI2_6   = PI*PI/6D0
C     PI2_6   = 1.64493406684822643D0
      CZ1     = 1D0-CZ
      AZ1     = ABS(CZ1)
      AZ      = ABS(CZ)
      RZ      = REAL(CZ)
 
      IF (AZ1.LT.1D-15) THEN
         CSPENC = PI2_6
      ELSE IF (RZ.LT.0.5D0) THEN
         IF (AZ.LT.1D0) THEN
            CSPENC = CSPENH(CZ,EPS)
         ELSE
            CSPENC = -PI2_6 - .5D0*CDLN(-CZ,-EPS)**2
     &                - CSPENH(1D0/CZ,-EPS)
         END IF
      ELSE
         IF (AZ1.LT.1D0) THEN
            CSPENC =  PI2_6 - CDLN(CZ,EPS)*CDLN(CZ1,-EPS)
     &                       - CSPENH(CZ1,-EPS)
         ELSE
            CSPENC = 2D0*PI2_6 + .5D0*CDLN(-CZ1,-EPS)**2
     &              - CDLN(CZ,EPS)*CDLN(CZ1,-EPS)
     &              + CSPENH(1D0/CZ1,EPS)
         END IF
      END IF
      END
************************************************************************
      FUNCTION CSPENH(CZ,EPS)
************************************************************************
*       COMPLEX SPENCE FUNCTION OF CZ + I*EPS                          *
*       IN CONVERGENCE REGION                                          *
*       CALCULATION OF BERNOULLI SERIES                                *
*----------------------------------------------------------------------*
*     09.01.90 ANSGAR DENNER                                           *
************************************************************************
      IMPLICIT   NONE
      COMPLEX*16 CSPENH,CDLN,CZ,X,X2
      REAL*8     EPS
      REAL*8     B(11)
ctp      data B(11)/
ctp     1   0.1666666666666666666666666667D0,
ctp     2  -0.0333333333333333333333333333D0,
ctp     3   0.0238095238095238095238095238D0,
ctp     4  -0.0333333333333333333333333333D0,
ctp     5   0.0757575757575757575757575758D0,
ctp     6  -0.2531135531135531135531135531D0,
ctp     7   1.1666666666666666666666666667D0,
ctp     8  -7.0921568627450980392156862745D0,
ctp     9  54.97117794486215538847117794486D0,
ctp     +  -529.124242424242424242424242424242D0,
ctp     1  6192.123188405797101449275362318D0  /
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     B(10)=-174611./330.
C     B(11)=854513./138.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...
C
      INTEGER    J
      REAL*8     FACTOR
      COMPLEX*16 POWER,TERM,CSP
 
      B(11)  =    854513D0/ 138D0
      B(10)  =  - 174611D0/ 330D0
      B(9)   =     43867D0/ 798D0
      B(8)   =  -   3617D0/ 510D0
      B(7)   =         7D0/   6D0
      B(6)   =  -    691D0/2730D0
      B(5)   =         5D0/  66D0
      B(4)   =  -      1D0/  30D0
      B(3)   =         1D0/  42D0
      B(2)   =  -      1D0/  30D0
      B(1)   =         1D0/   6D0
      X      =  -CDLN(1D0-CZ,-EPS)
C     WRITE(80,*)  'CSPENH'
      X2     =  X*X
      POWER  =  X
      FACTOR =  1D0
      CSPENH =  X - X2/4D0
      DO 10 J=2,22,2
         FACTOR = FACTOR / J / (J+1)
         POWER  = POWER * X2
         TERM   = B(J/2) * FACTOR * POWER
         CSP    = CSPENH + TERM
         IF (CSP.EQ.CSPENH) RETURN
         CSPENH = CSP
10    CONTINUE
C      WRITE(80,*) 'CSPENH CONVERGES BADLY  ',CZ,X
C      WRITE(80,*) 'CSPENH CONVERGES BADLY  ',CSP-TERM,CSPENH,TERM
      END
************************************************************************
      FUNCTION CSPCOE(Z1,Z2,I1,I2)
************************************************************************
*  COMPLEX SPENCE FUNCTION PLUS CONTINUATION TERMS                     *
*----------------------------------------------------------------------*
*  08.07.94 ANSGAR DENNER      LAST CHANGED  14.07.94                  *
************************************************************************
      IMPLICIT   NONE
      COMPLEX*16 CSPCOE,CSPENC,ETAE,CDLN,Z1,Z2,Z12,ETAA
      REAL*8     I1,I2
      REAL*8     PI,PI2_6
      COMMON /PI/     PI,PI2_6
 
      Z12 = Z1*Z2
ctp      IF(REAL(Z12).GT.0.5D0) THEN
      IF(REAL(Z12).GT.0.5D0) THEN
        CSPCOE = CSPENC(1D0-Z12,0D0)
        ETAA   = ETAE(Z1,Z2,I1,I2,-I1)
        IF(ETAA.NE.0D0) THEN
          CSPCOE = CSPCOE + ETAA*CDLN(1D0-Z12,-I1)
        END IF
      ELSE IF(ABS(Z12).LT.1D-4) THEN
        CSPCOE = PI2_6-CSPENC(Z12,0D0)
     &           + (CDLN(Z1,I1)+CDLN(Z2,I2))
     &                  *Z12*(1D0+Z12/2D0+Z12*Z12/3D0+Z12*Z12*Z12/4D0)
      ELSE
        CSPCOE = PI2_6-CSPENC(Z12,I1)
     &           - (CDLN(Z1,I1)+CDLN(Z2,I2))
     &                  *CDLN(1D0-Z12,-I1)
      END IF
 
      END
************************************************************************
      FUNCTION ETAE(C1,C2,I1,I2,I12)
************************************************************************
*     COMPLEX ETA-FUNCTION
*----------------------------------------------------------------------*
*     8.06.90    ANSGAR DENNER       LAST CHANGED   11.07.94
************************************************************************
      IMPLICIT     NONE
      COMPLEX*16 ETAE,C1,C2
      REAL*8     I1,I2,I12
      REAL*8     IM1,IM2,IM12
      REAL*8     PI,PI2_6
      COMMON /PI/     PI,PI2_6
 
      IM1    = AIMAG(C1)
      IM2    = AIMAG(C2)
      IM12   = AIMAG(C1*C2)
      IF(IM1 .EQ.0D0) IM1  = I1
      IF(IM2 .EQ.0D0) IM2  = I2
      IF(IM12.EQ.0D0) IM12 = I12
 
      IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN
          ETAE = DCMPLX(0D0,2D0*PI)
      ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN
          ETAE = DCMPLX(0D0,-2D0*PI)
      ELSE
          ETAE = DCMPLX(0D0)
          IF(IM1.EQ.0.AND.REAL(C1).LT.0D0 .OR.
     &       IM2.EQ.0.AND.REAL(C2).LT.0D0 .OR.
     &       IM12.EQ.0.AND.REAL(C1*C2).LT.0D0) THEN
C$$$             WRITE(80,*) ' ETA NOT DEFINED '
C$$$             WRITE(80,*) ' ETA:  C1  = ',C1
C$$$             WRITE(80,*) ' ETA:  C2  = ',C2
C$$$             WRITE(80,*) ' ETA:  C12 = ',C1*C2
          END IF
      END IF
      END
************************************************************************
      FUNCTION D0REG(Q12,Q23,Q34,Q14,Q24,Q13,M1,M2,M3,M4)
************************************************************************
*  SCALAR 4-POINT FUNCTION  FOR  Q13 < 0 (ONE R REAL POSITIVE)         *
*  REGULAR CASE                                                        *
*  IMAGINARY PART SOMETIMES WRONG                                      *
*  QIJ = (PI-PJ)**2   AND   PI IS MOMENTUM OF PROPAGATOR WITH MASS MI  *
*----------------------------------------------------------------------*
*  07.01.94 ANSGAR DENNER       LAST CHANGED 13.07.94 ANSGAR DENNER    *
************************************************************************
      IMPLICIT   NONE
      REAL*8     Q12,Q23,Q34,Q14,Q13,Q24,M1,M2,M3,M4
      REAL*8     K12,K13,K14,K23,K24,K34
      REAL*8     M12,M22,M32,M42
      REAL*8     IR12,IR14,IR23,IR24,IR34
      REAL*8     IX(2,4),IS(4)
      COMPLEX*16 R12,R13,R14,R23,R24,R34
      COMPLEX*16 A,B,C
      COMPLEX*16 X(2,4),S(4)
      COMPLEX*16 D0REG,CSPCOE,ETAE,CDLN
      INTEGER    K,J
      REAL*8     PI,PI2_6
      COMMON /PI/     PI,PI2_6
 
      PI=3.1415926535897932384626433832D0
      PI2_6 = 1.6449340668482264364724151666D0
 
      M12 = M1*M1
      M22 = M2*M2
      M32 = M3*M3
      M42 = M4*M4
      K12 = (M12+M22-Q12)/M1/M2
      K13 = (M12+M32-Q13)/M1/M3
      K14 = (M12+M42-Q14)/M1/M4
      K23 = (M22+M32-Q23)/M2/M3
      K24 = (M22+M42-Q24)/M2/M4
      K34 = (M32+M42-Q34)/M3/M4
      IF (K13.LT.2D0) THEN
       WRITE(*,*) ' D0REG: CASE NOT IMPLEMENTED'
       WRITE(*,*) ' K13 = ',K13
       WRITE(*,*) ' Q13 = ',Q13
       WRITE(*,*) ' M1  = ',M1
       WRITE(*,*) ' M3  = ',M3
      END IF
      R12 = K12/2D0*(1D0+SQRT(DCMPLX(1D0-4D0/K12**2)))
      R13 = K13/2D0*(1D0+SQRT(DCMPLX(1D0-4D0/K13**2)))
      R13 = 1D0/R13
      R14 = K14/2D0*(1D0+SQRT(DCMPLX(1D0-4D0/K14**2)))
      R23 = K23/2D0*(1D0+SQRT(DCMPLX(1D0-4D0/K23**2)))
      R24 = K24/2D0*(1D0+SQRT(DCMPLX(1D0-4D0/K24**2)))
      R24 = 1D0/R24
      R34 = K34/2D0*(1D0+SQRT(DCMPLX(1D0-4D0/K34**2)))
      A   =  K34/R24-K23 + (K12-K14/R24)*R13
      B   =  (1D0/R13-R13)*(1D0/R24-R24)+K12*K34-K14*K23
      C   =  K34*R24-K23 + (K12-K14*R24)/R13
      X(1,4) = (-B+SQRT(B*B-4D0*A*C))/2D0/A
      X(2,4) = (-B-SQRT(B*B-4D0*A*C))/2D0/A
      IF(ABS(X(1,4)).GT.ABS(X(2,4))) THEN
        X(2,4) = C/A/X(1,4)
      ELSE
        X(1,4) = C/A/X(2,4)
      END IF
 
      IF(K12.LT.-2D0) THEN
        IR12 = SIGN(1D1,1D0-ABS(R12))
      ELSE
        IR12 = 0D0
      END IF
      IF(K14.LT.-2D0) THEN
        IR14 = SIGN(1D1,1D0-ABS(R14))
      ELSE
        IR14 = 0D0
      END IF
      IF(K23.LT.-2D0) THEN
        IR23 = SIGN(1D1,1D0-ABS(R23))
      ELSE
        IR23 = 0D0
      END IF
      IF(K24.LT.-2D0) THEN
        IR24 = SIGN(1D1,1D0-ABS(R24))
      ELSE
        IR24 = 0D0
      END IF
      IF(K34.LT.-2D0) THEN
        IR34 = SIGN(1D1,1D0-ABS(R34))
      ELSE
        IR34 = 0D0
      END IF
      IF(REAL(X(1,4)).GT.0D0) THEN
         IX(1,4) = 1D0
      ELSE
         IX(1,4) = 0D0
      END IF
      IF(REAL(X(2,4)).GT.0D0) THEN
         IX(2,4) = -1D0
      ELSE
         IX(2,4) = 0D0
      END IF
 
      X(1,1) = X(1,4)/R24
      X(2,1) = X(2,4)/R24
      X(1,2) = X(1,4)/R24*R13
      X(2,2) = X(2,4)/R24*R13
      X(1,3) = X(1,4)*R13
      X(2,3) = X(2,4)*R13
      S(1)  = R12
      S(2)  = R23
      S(3)  = R34
      S(4)  = R14
 
      IS(1)  = IR12
      IS(2)  = IR23
      IS(3)  = IR34
      IS(4)  = IR14
      IF(REAL(X(1,1)).GT.0D0) THEN
         IX(1,1) = IX(1,4) + IR24
      ELSE
         IX(1,1) = 0D0
      END IF
      IF(REAL(X(2,1)).GT.0D0) THEN
         IX(2,1) = IX(2,4) + IR24
      ELSE
         IX(2,1) = 0D0
      END IF
      IX(1,3) = IX(1,4)
      IX(2,3) = IX(2,4)
      IX(1,2) = IX(1,1)
      IX(2,2) = IX(2,1)
 
      D0REG = DCMPLX(0D0)
      DO 20 K=1,2
        DO 10 J=1,4
        D0REG = D0REG + (-1)**(J+K) * (
     &           CSPCOE(-X(K,J),S(J),-IX(K,J),IS(J))
     &         + CSPCOE(-X(K,J),1D0/S(J),-IX(K,J),-IS(J)) )
10      CONTINUE
      D0REG = D0REG - (-1)**K *
     &            ETAE(-X(K,4),1D0/R24,-IX(K,4),-IR24,-IX(K,1)) *
     &   CDLN((1D0+K14*X(K,4)+X(K,4)**2)/(1D0+K34*X(K,3)+X(K,3)**2),
     &         -REAL(1D0+K34*X(K,3)+X(K,3)**2))
20    CONTINUE
      D0REG = D0REG/M1/M2/M3/M4/SQRT(B*B-4D0*A*C)

      END




c$$$************************************************************************
c$$$*                                                                      *
c$$$*     SUBROUTINE INTEG (CALLING VEGAS):                                *
c$$$*     ---------------------------------                                *
c$$$*                                                                      *
c$$$*     IDIM-DIMENSIONAL MONTE CARLO INTEGRATION OF FUNCTION "FUNCT"     *
c$$$*     1ST RUN (VEGAS) : IPOINT  CALLS, ITER  ITERATIONS                *
c$$$*     2ND RUN (VEGAS1): IPOINT1 CALLS, ITER1 ITERATIONS                *
c$$$*     RESULT "RES" IS RETURNED WITH DESIRED ACCURACY "ACC"             *
c$$$*                                                                      *
c$$$************************************************************************
c$$$
c$$$
c$$$
c$$$C----------------------------- MAIN PROGRAM ----------------------------
c$$$C
c$$$C     IMPLICIT REAL*8 (A-H,K-Z)
c$$$C
c$$$C     EXTERNAL FUNCT
c$$$C
c$$$C     COMMON/IOUT/IPRINT
c$$$C
c$$$C     ...
c$$$C
c$$$C     PARAMETERS FOR VEGAS
c$$$C     IDIM    =  ...
c$$$C     IPOINT  =  ...
c$$$C     ITER    =   10
c$$$C     IPOINT1 =  ...
c$$$C     ITER1   =    3
c$$$C     ACC     = 1.D-5
c$$$C
c$$$C     NO / LONG / SHORT PRINTOUT OF INTEGRALS: IPRINT = 0 / 1 / 10
c$$$C     IPRINT = 0
c$$$C
c$$$C     ...
c$$$C
c$$$C     CALL INTEG(FUNCT,IDIM,IPOINT,ITER,IPOINT1,ITER1,ACC,RES)
c$$$C
c$$$C     ...
c$$$C
c$$$C     STOP
c$$$C     END
c$$$C
c$$$C-----------------------------------------------------------------------
c$$$
c$$$
c$$$
c$$$
c$$$      SUBROUTINE INTEG(FXN,IDIM,IPOINT,ITER,IPOINT1,ITER1,ACC,RES)
c$$$
c$$$      IMPLICIT REAL*8 (A-H,K-Z)
c$$$
c$$$      COMMON/BVEG2/INDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
c$$$     1,D(50,10),DI(50,10),INXI(50,10)
c$$$      COMMON/RESULT/S1,S2,S3,S4
c$$$      COMMON/IOUT/IPRINT
c$$$
c$$$      EXTERNAL FXN
c$$$      IF (IPRINT.NE.0)
c$$$     .PRINT *,'------------------- CALL TO VEGAS -------------------'
c$$$      CALL VEGAS(FXN,ACC,IDIM,IPOINT,ITER,IPRINT,0)
c$$$      IF (IPRINT.NE.0)
c$$$     .PRINT *,'------------------- CALL TO VEGAS1 ------------------'
c$$$      CALL VEGAS1(FXN,ACC,IDIM,IPOINT1,ITER1,IPRINT,0)
c$$$      IF (IPRINT.NE.0) PRINT *,'        '
c$$$
c$$$      RES=S1
c$$$
c$$$      RETURN
c$$$      END
c$$$
c$$$
c$$$
c$$$
c$$$      SUBROUTINE VEGAS(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
c$$$      IMPLICIT REAL*8 (A-H,O-Z)
c$$$      REAL*4 UNIV
c$$$      COMMON/BVEG2/INDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
c$$$     1,D(50,10),DI(50,10),INXI(50,10)
c$$$C***  ADDED BY RH (FOR SOME FORTRAN COMPILERS): 
c$$$      COMMON/BVEG3/ACC
c$$$C***  END OF ADDITION.
c$$$      DIMENSION XIN(50),R(50),DX(10),IA(10),KG(10),DT(10)
c$$$      DIMENSION XL(10),XU(10),QRAN(10),X(10)
c$$$      COMMON/RESULT/S1,S2,S3,S4
c$$$       EXTERNAL FXN
c$$$      DATA XL,XU/10*0.D0,10*1.D0/
c$$$      DATA NDMX/50/,ALPH/1.5D0/,ONE/1.D0/,MDS/1/
c$$$      IPR=1
c$$$      IF(NPRN.GT.0)IPR=0
c$$$      INDO=1
c$$$      DO 1 J=1,NDIM
c$$$1     XI(1,J)=ONE
c$$$      ENTRY VEGAS1(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
c$$$      NOW=IGRAPH
c$$$C---  NEXT LINE TAKEN OUT ON 04/06/92 BY J. ZUNFT
c$$$C---  IF(IGRAPH.GT.0)CALL INPLOT(NOW,F1,W)
c$$$      IT=0
c$$$      SI=0.D0
c$$$      SI2=SI
c$$$      SWGT=SI
c$$$      SCHI=SI
c$$$      SCALLS=SI
c$$$      ENTRY VEGAS2(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
c$$$      ND=NDMX
c$$$      NG=1
c$$$      IF(MDS.EQ.0) GO TO 2
c$$$      NG=(NCALL*0.5D0)**(1.D0/NDIM)
c$$$      MDS=1
c$$$      IF((2*NG-NDMX).LT.0) GO TO 2
c$$$      MDS=-1
c$$$      NPG=NG/NDMX+1
c$$$      ND=NG/NPG
c$$$      NG=NPG*ND
c$$$2     K=NG**NDIM
c$$$      NPG=NCALL/K
c$$$      IF(NPG.LT.2)NPG=2
c$$$      CALLS=NPG*K
c$$$      DXG=ONE/NG
c$$$      DV2G=DXG**(2*NDIM)/NPG/NPG/(NPG-ONE)
c$$$      XND=ND
c$$$      NDM=ND-1
c$$$      DXG=DXG*XND
c$$$      XJAC=ONE
c$$$      DO 3 J=1,NDIM
c$$$      DX(J)=XU(J)-XL(J)
c$$$3     XJAC=XJAC*DX(J)
c$$$      IF(ND.EQ.INDO) GO TO 8
c$$$      RC=INDO/XND
c$$$      DO 7 J=1,NDIM
c$$$      K=0
c$$$      XN=0.D0
c$$$      DR=XN
c$$$      I=K
c$$$4     K=K+1
c$$$      DR=DR+ONE
c$$$      XO=XN
c$$$      XN=XI(K,J)
c$$$5     IF(RC.GT.DR) GO TO 4
c$$$      I=I+1
c$$$      DR=DR-RC
c$$$      XIN(I)=XN-(XN-XO)*DR
c$$$      IF(I.LT.NDM) GO TO 5
c$$$      DO 6  I=1,NDM
c$$$6     XI(I,J)=XIN(I)
c$$$7     XI(ND,J)=ONE
c$$$      INDO=ND
c$$$      ACC=BCC
c$$$8     IF(NPRN.NE.0.AND.NPRN.NE.10)PRINT 200,NDIM,CALLS,IT,ITMX
c$$$     1,ACC,MDS,ND
c$$$      IF(NPRN.EQ.10)PRINT 290,NDIM,CALLS,ITMX,ACC,MDS,ND
c$$$      ENTRY VEGAS3(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
c$$$9     IT=IT+1
c$$$      TI=0.D0
c$$$      TSI=TI
c$$$C---  NEXT LINE TAKEN OUT ON 04/06/92 BY J. ZUNFT
c$$$C---  IF(IGRAPH.GT.0)CALL REPLOT(NOW,F1,W)
c$$$      DO 10 J=1,NDIM
c$$$      KG(J)=1
c$$$      DO 10 I=1,ND
c$$$      INXI(I,J)=0
c$$$      D(I,J)=TI
c$$$10    DI(I,J)=TI
c$$$11    FB=0.D0
c$$$      F2B=FB
c$$$      K=0
c$$$12    K=K+1
c$$$      DO 121 J=1,NDIM
c$$$121   QRAN(J)=DBLE(UNIV())
c$$$      WGT=XJAC
c$$$      DO 15 J=1,NDIM
c$$$      XN=(KG(J)-QRAN(J))*DXG+ONE
c$$$      IA(J)=XN
c$$$      IAJ=IA(J)
c$$$      IAJ1=IAJ-1
c$$$      IF(IAJ.GT.1) GO TO 13
c$$$      XO=XI(IAJ,J)
c$$$      RC=(XN-IAJ)*XO
c$$$      GO TO 14
c$$$13    XO=XI(IAJ,J)-XI(IAJ1,J)
c$$$      RC=XI(IAJ1,J)+(XN-IAJ)*XO
c$$$14    X(J)=XL(J)+RC*DX(J)
c$$$15    WGT=WGT*XO*XND
c$$$      F=FXN(X)*WGT
c$$$      F1=F/CALLS
c$$$      W=WGT/CALLS
c$$$C---  NEXT LINE TAKEN OUT ON 04/06/92 BY J. ZUNFT
c$$$C---  IF(IGRAPH.GT.0)CALL XPLOT(NOW,F1,W)
c$$$      F2=F*F
c$$$      FB=FB+F
c$$$      F2B=F2B+F2
c$$$      DO 16 J=1,NDIM
c$$$      IAJ=IA(J)
c$$$      INXI(IAJ,J)=INXI(IAJ,J)+1
c$$$      DI(IAJ,J)=DI(IAJ,J)+F/CALLS
c$$$16    IF(MDS.GE.0)  D(IAJ,J)=D(IAJ,J)+F2
c$$$      IF(K.LT.NPG) GO TO 12
c$$$      F2B=F2B*NPG
c$$$      F2B=SQRT(F2B)
c$$$      F2B=(F2B-FB)*(F2B+FB)
c$$$      TI=TI+FB
c$$$      TSI=TSI+F2B
c$$$      IF(MDS.GE.0) GO TO 18
c$$$      DO 17 J=1,NDIM
c$$$      IAJ=IA(J)
c$$$17    D(IAJ,J)=D(IAJ,J)+F2B
c$$$18    K=NDIM
c$$$19    KG(K)=MOD(KG(K),NG)+1
c$$$      IF(KG(K).NE.1) GO TO 11
c$$$      K=K-1
c$$$      IF(K.GT.0) GO TO 19
c$$$      TI=TI/CALLS
c$$$      TSI=TSI*DV2G
c$$$      TI2=TI*TI
c$$$C--------- CHANGE BY J. ZUNFT, 05/06/92------------------
c$$$      IF (ABS(TSI).LT.1.D-35) TSI = SIGN(1.D-35,TSI)
c$$$C--------------------------------------------------------
c$$$      WGT=TI2/TSI
c$$$      SI=SI+TI*WGT
c$$$      SI2=SI2+TI2
c$$$      SWGT=SWGT+WGT
c$$$      SCHI=SCHI+TI2*WGT
c$$$      SCALLS=SCALLS+CALLS
c$$$C--------- CHANGE BY J. ZUNFT, 05/06/92------------------
c$$$      IF (ABS(SWGT).LT.1.D-35) SWGT = SIGN(1.D-35,SWGT)
c$$$      IF (ABS(SI2) .LT.1.D-35) SI2  = SIGN(1.D-35,SI2)
c$$$C--------------------------------------------------------
c$$$      AVGI=SI/SWGT
c$$$      SD=SWGT*IT/SI2
c$$$      CHI2A=0.
c$$$      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
c$$$      SD=ONE/SD
c$$$C--------- CHANGE BY J. ZUNFT, 04/14/92-----------------
c$$$      IF (SD.LT.0.D0) THEN
c$$$      PRINT *,' SD  = ',SD,' < 0: SIGN CHANGED'
c$$$      SD = ABS(SD)
c$$$      END IF
c$$$      IF (TSI.LT.0.D0) THEN
c$$$      PRINT *,' TSI = ',TSI,' < 0: SIGN CHANGED'
c$$$      TSI = ABS(TSI)
c$$$      END IF
c$$$C-------------------------------------------------------
c$$$      SD=SQRT(SD)
c$$$      IF(NPRN.EQ.0) GO TO 21
c$$$      TSI=SQRT(TSI)
c$$$      IF(NPRN.NE.10)PRINT 201,IPR,IT,TI,TSI,AVGI,SD,CHI2A
c$$$      IF(NPRN.EQ.10)PRINT 203,IT,TI,TSI,AVGI,SD,CHI2A
c$$$      IF(NPRN.GE.0) GO TO 21
c$$$      DO 20 J=1,NDIM
c$$$      PRINT 202,J
c$$$20    PRINT 204,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
c$$$C--------- CHANGE BY J. ZUNFT, 05/06/92------------------
c$$$21    IF (ABS(AVGI).LT.1.D-35) AVGI = SIGN(1.D-35,AVGI)
c$$$C--------------------------------------------------------
c$$$      IF(ABS(SD/AVGI).LE.ABS(ACC).OR.IT.GE.ITMX)NOW=2
c$$$      S1=AVGI
c$$$      S2=SD
c$$$      S3=TI
c$$$      S4=TSI
c$$$C---  NEXT LINE TAKEN OUT ON 04/06/92 BY J. ZUNFT
c$$$C---  IF(IGRAPH.GT.0)CALL PLOTIT(NOW,F1,W)
c$$$C      DO 23 J=1,NDIM
c$$$C      XO=D(1,J)
c$$$C      XN=D(2,J)
c$$$C      D(1,J)=(XO+XN)*0.5D0
c$$$C      DT(J)=D(1,J)
c$$$C      DO 22 I=2,NDM
c$$$C      D(I,J)=XO+XN
c$$$C      XO=XN
c$$$C      XN=D(I+1,J)
c$$$C      D(I,J)=(D(I,J)+XN)/3.D0
c$$$C22    DT(J)=DT(J)+D(I,J)
c$$$C      D(ND,J)=(XN+XO)*0.5D0
c$$$C23    DT(J)=DT(J)+D(ND,J)
c$$$C-----THIS PART OF THE VEGAS-ALGORITHM IS UNSTABLE
c$$$C-----IT SHOULD BE REPLACED BY
c$$$      DO 23 J=1,NDIM
c$$$      DT(J)=0.D0
c$$$      DO 23 I=1,ND
c$$$      IF(INXI(I,J).GT.0)D(I,J)=D(I,J)/INXI(I,J)
c$$$23    DT(J)=DT(J)+D(I,J)
c$$$      DO 28 J=1,NDIM
c$$$      RC=0.D0
c$$$      DO 24 I=1,ND
c$$$      R(I)=0.D0
c$$$C--------- CHANGE BY J. ZUNFT, 04/15/92 ---------
c$$$C---  IF(D(I,J).LE.0.D0)GO TO 24
c$$$      IF(D(I,J).EQ.0.D0.OR.DT(J)/D(I,J).LE.0.D0) GO TO 24
c$$$C------------------------------------------------
c$$$      XO=DT(J)/D(I,J)
c$$$C--------- CHANGE BY J. ZUNFT, 05/06/92 ---------
c$$$      IF (XO.EQ.1.D0) XO = 0.99999999999999999D0
c$$$C------------------------------------------------
c$$$      R(I)=((XO-ONE)/XO/LOG(XO))**ALPH
c$$$24    RC=RC+R(I)
c$$$      RC=RC/XND
c$$$      K=0
c$$$      XN=0.D0
c$$$      DR=XN
c$$$      I=K
c$$$25    K=K+1
c$$$      DR=DR+R(K)
c$$$      XO=XN
c$$$      XN=XI(K,J)
c$$$26    IF(RC.GT.DR) GO TO 25
c$$$      I=I+1
c$$$      DR=DR-RC
c$$$C----------- CHANGE BY J. ZUNFT, 04/14/92 --------------
c$$$C---  XIN(I)=XN-(XN-XO)*DR/R(K)
c$$$      IF (DR.EQ.0.D0) THEN
c$$$      XIN(I)=XN
c$$$      ELSE
c$$$      XIN(I)=XN-(XN-XO)*DR/R(K)
c$$$      END IF
c$$$C-------------------------------------------------------
c$$$      IF(I.LT.NDM) GO TO 26
c$$$      DO 27 I=1,NDM
c$$$27    XI(I,J)=XIN(I)
c$$$28    XI(ND,J)=ONE
c$$$      IF(IT.LT.ITMX.AND.ABS(ACC).LT.ABS(SD/AVGI))GO TO 9
c$$$200   FORMAT(35H0INPUT PARAMETERS FOR VEGAS   NDIM=,I3
c$$$     1,8H  NCALL=,F8.0/28X,5H  IT=,I5,8H  ITMX =,I5/28X
c$$$     2,6H  ACC=,G9.3/28X,6H  MDS=,I3,6H   ND=,I4//)
c$$$290   FORMAT(13H0VEGAS  NDIM=,I3,8H  NCALL=,F8.0,8H  ITMX =,I5
c$$$     1,6H  ACC=,G9.3,6H  MDS=,I3,6H   ND=,I4)
c$$$201   FORMAT(/I1,20HINTEGRATION BY VEGAS/13H0ITERATION NO,I3,
c$$$     114H.   INTEGRAL =,G14.8/20X,10HSTD DEV  =,G10.4/
c$$$     234H ACCUMULATED RESULTS.   INTEGRAL =,G14.8/
c$$$     324X,10HSTD DEV  =,G10.4 / 24X,18HCHI**2 PER ITN   =,G10.4)
c$$$202   FORMAT(14H0DATA FOR AXIS,I2 / 7X,1HX,7X,10H  DELT I  ,
c$$$     12X,11H CONVCE    ,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE
c$$$     2,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE     /)
c$$$204   FORMAT(1X,3G12.4,5X,3G12.4,5X,3G12.4)
c$$$203   FORMAT(1H ,I3,G20.8,G12.4,G20.8,G12.4,G12.4)
c$$$C---  NEXT 3 LINES TAKEN OUT ON 04/06/92 BY J. ZUNFT
c$$$C---  S1=AVGI
c$$$C---  S2=SD
c$$$C---  S3=CHI2A
c$$$      RETURN
c$$$      END
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$C-----------------------------------------------------------------------
c$$$C   INITIALIZING THE RANDOM NUMBER GENERATOR
c$$$C     IF OLD CONFIGURATION EXISTS THEN
c$$$C     CALL UREAD(11)
c$$$C     ELSE
c$$$C     CALL RSTART(12,34,56,78)
c$$$C     END IF
c$$$C-----------------------------------------------------------------------
c$$$
c$$$
c$$$C----------------------------------------------------------------------
c$$$C  A UNIVERSAL RANDOM NUMBER GENERATOR
c$$$
c$$$        FUNCTION UNIV()
c$$$        REAL U(97)
c$$$        COMMON /SET1/ U,C,CD,CM,I,J
c$$$        UNIV=U(I)-U(J)
c$$$        IF(UNIV.LT.0.) UNIV=UNIV+1.
c$$$        U(I)=UNIV
c$$$        I=I-1
c$$$        IF(I.EQ.0) I=97
c$$$        J=J-1
c$$$        IF(J.EQ.0) J=97
c$$$        C=C-CD
c$$$        IF(C.LT.0.) C=C+CM
c$$$        UNIV=UNIV-C
c$$$        IF(UNIV.LT.0.) UNIV=UNIV+1
c$$$        RETURN
c$$$        END
c$$$
c$$$C----------------------------------------------------------------------
c$$$C SAVING THE RANDOM NUMBER GENERATOR CONFIGURATION
c$$$
c$$$        SUBROUTINE UWRITE(IUNIT)
c$$$        INTEGER IUNIT
c$$$        REAL U(97)
c$$$        COMMON /SET1/ U,C,CD,CM,I,J
c$$$        OPEN(IUNIT,FILE='TESTRN',FORM='UNFORMATTED')
c$$$        WRITE(IUNIT) U,C,CD,CM,I,J
c$$$        CLOSE(IUNIT)
c$$$        RETURN
c$$$        END
c$$$
c$$$
c$$$C----------------------------------------------------------------------
c$$$C READING IN THE RANDOM NUMBER GENERATOR CONFIGURATION
c$$$
c$$$        SUBROUTINE UREAD(IUNIT)
c$$$        INTEGER IUNIT
c$$$        REAL U(97)
c$$$        COMMON /SET1/ U,C,CD,CM,I,J
c$$$        OPEN(IUNIT,FILE='TESTRN',FORM='UNFORMATTED')
c$$$        READ(IUNIT) U,C,CD,CM,I,J
c$$$        CLOSE(IUNIT)
c$$$        RETURN
c$$$        END
c$$$
c$$$C----------------------------------------------------------------------
c$$$C INITIALIZING THE RANDOM NUMBER GENERATOR
c$$$C TO INITIALIZE CALL RSTART(12,34,56,78)
c$$$
c$$$        SUBROUTINE RSTART(I,J,K,L)
c$$$        REAL U(97)
c$$$        COMMON /SET1/ U,C,CD,CM,ISTART,JSTART
c$$$        IF ((I.LT.0).OR.(I.GT.178 )) STOP 'FIRST SEED .LT.0 OR .GT.178'
c$$$        IF ((J.LT.0).OR.(J.GT.178 )) STOP 'SECOND SEED .LT.0 OR .GT.178'
c$$$        IF ((K.LT.0).OR.(K.GT.178 )) STOP 'THIRD SEED .LT.0 OR .GT.178'
c$$$        IF ((L.LT.0).OR.(L.GT.168 )) STOP 'FOURTH SEED .LT.0 OR .GT.168'
c$$$        IF ( (I.EQ.1).AND.(J.EQ.1).AND.(K.EQ.1) ) STOP
c$$$     &     'FIRST, SECOND AND THIRD SEEDS ARE ALL EQUAL TO 1'
c$$$        ISTART=97
c$$$        JSTART=33
c$$$        IDUM = I
c$$$        JDUM = J
c$$$        KDUM = K
c$$$        LDUM = L
c$$$        DO 2 II=1,97
c$$$        S=0.
c$$$        T=.5
c$$$        DO 3 JJ=1,24
c$$$          M=MOD(MOD(IDUM*JDUM,179)*KDUM,179)
c$$$          IDUM=JDUM
c$$$          JDUM=KDUM
c$$$          KDUM=M
c$$$          LDUM=MOD(53*LDUM+1,169)
c$$$          IF(MOD(LDUM*M,64).GE.32) S=S+T
c$$$3         T=.5*T
c$$$2         U(II)=S
c$$$        C=362436./16777216.
c$$$        CD=7654321./16777216.
c$$$        CM=16777213./16777216.
c$$$        RETURN
c$$$        END



C ======================================================================

      REAL*8 FUNCTION DACOSH(X)
      IMPLICIT REAL*8 (A-H,M-Z)
      DACOSH = LOG(X + SQRT(X**2 -1.D0))
      RETURN
      END

C ======================================================================

















