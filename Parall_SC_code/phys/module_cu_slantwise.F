MODULE module_cu_slantwise

   USE module_wrf_error

!--------------------------------------------------------------------
! Lookup table variables (kept from KF scheme):
      INTEGER, PARAMETER, PRIVATE :: KFNT=300,KFNP=220
      REAL, DIMENSION(KFNT,KFNP),PRIVATE, SAVE :: TTAB,QSTAB
      REAL, DIMENSION(KFNP),PRIVATE, SAVE :: THE0K
      REAL, DIMENSION(200),PRIVATE, SAVE :: ALU
      REAL, PRIVATE, SAVE :: RDPR,RDTHK,PLUTOP
! End of Lookup table variables:

CONTAINS

  SUBROUTINE slantwise_adjust(                                &
               ids,ide, jds,jde, kds,kde                      &
              ,ims,ime, jms,jme, kms,kme                      &
              ,its,ite, jts,jte, kts,kte                      &
              ,DT,DX,SLCUDT                                   &
              ,rho, RAINSCV                                   &
              ,U,V,TH0,T0,W,dz8w,Pcps,pi,f                    & 
              !Pcps should the pressure on the model half level (theta points)
              !the U,V, TH, and T are also all on the theta grid points.
              ,W0AVG,XLV0,XLV1,XLS0,XLS1,CP,R,G,EP1           & !R=r_d
              ,EP2,SVP1,SVP2,SVP3,SVPT0                       &
              ,R_V, XLV                                       &
              ,STEPCU                                         &
              ,QV0, QC0, QR0                                  &
              ,CAPE, SCAPE1D                                  &
              ,ang_SW                                         &
              ,adjust_frac                                    & 
              ,RUSLCUTEN, RVSLCUTEN                           &
              ,RTHSLCUTEN, RQVSLCUTEN                         &
              ,RQVCUTEN,RTHCUTEN,thetaes, SW_calls            )
!

   IMPLICIT NONE
   INTEGER,      INTENT(IN   ) ::                            &
                                  ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte

   INTEGER,      INTENT(IN   ) :: STEPCU

   REAL,         INTENT(IN   ) :: XLV0,XLV1,XLS0,XLS1
   REAL,         INTENT(IN   ) :: CP,R,G,EP1,EP2,R_V,XLV
   REAL,         INTENT(IN   ) :: SVP1,SVP2,SVP3,SVPT0

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         , &
          INTENT(IN   ) ::                                   &
                                                          U, &
                                                          V, &
                                                          W, &
                                                        TH0, &
                                                         T0, &
                                                        QV0, &
                                                        QC0, &
                                                        QR0, &
                                                       dz8w, &
                                                       Pcps, &
                                                        rho, &
                                                         pi
!
   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: f 
!
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         , &
          INTENT(INOUT) ::                                   &
                                                      W0AVG   

   REAL,  INTENT(IN   ) :: DT, DX
   REAL,  INTENT(IN   ) :: SLCUDT
!
   REAL, DIMENSION( ims:ime , jms:jme ),                     &
          INTENT(  OUT) ::                          RAINSCV

   REAL, DIMENSION( ims:ime , jms:jme ),                     &
          INTENT(  OUT) ::                         SW_calls
!
   REAL,  DIMENSION( ims:ime , jms:jme )                   , &
          INTENT(  OUT) ::                                   &
                                                       CAPE, &
                                                    SCAPE1D, & 
                                                adjust_frac 

   REAL, DIMENSION( ims:ime , jms:jme ),                     &
          INTENT(OUT) ::                            ang_SW  

!
! Optional arguments
!

   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),           &
         OPTIONAL,                                           &
         INTENT(IN   ) ::                                    &
                                                   RTHCUTEN, &
                                                   RQVCUTEN 

   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),           &
         OPTIONAL,                                           &
         INTENT(INOUT) ::                                    &
                                                  RUSLCUTEN, &
                                                  RVSLCUTEN, &
                                                 RTHSLCUTEN, &
                                                 RQVSLCUTEN 
!
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),          &
          OPTIONAL,                                          &
          INTENT(INOUT) ::                                   &
                                                    thetaes


! LOCAL VARS
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) ::        &
                                                          T, &
                                                         TH, &
                                                         QV, &
                                                         Tv, &
                                                         RH

   LOGICAL :: twoD_nox_setup, twoD_noy_setup,pathadjustdone
   LOGICAL :: flag_qr, flag_qi, flag_qs, sw_trigger,skipflag

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) ::        &
                                                  UgTEND_IJ, &
                                                   TTEND_IJ, &
                                                   QTEND_IJ

   REAL, DIMENSION( kts:kte ) ::                             &
                                                        U1D, &
                                                        V1D, &
                                                        T1D, &
                                                       DZ1D, &
                                                       QV1D, &
                                                        P1D, &
                                                      RHO1D, &
                                                    W0AVG1D, & 
                                                       Ug1D, &
                                                       Vg1D, &
                                                      vor1D, &
                                                       RUg1 

   REAL, DIMENSION( kts:kte )::                              &
                                                       DQDT, &
                                                       DTDT, &
                                                   DUgDT_RH, &
                                                         Z0, &
                                                      RHavg, &
                                                      QVavg !, &
                                                       !Wavg

   REAL,  DIMENSION( 2 , kts:kte ) ::                       &
                                                       SUg, &
                                                        SP, &
                                                  Sthetaes
                                              
   REAL,  DIMENSION( kts:kte )  ::                          &
                                              SHIFT_IN_AD0, &  
                                               SHIFT_IN_AD, &  
                                                   SLOPE_M, &  
                                                   SLOPE_T, & 
                                                     DUgDT

   REAL, DIMENSION( ims:ime , jms:jme ) ::                   &
                                         RAINincrease

   INTEGER :: ii,jj,kk,overshoot_k,KLCLn1,m, markKbest
   REAL :: SHIFT_IN_M, SHIFT_IN_T, TIMEC, Z1
   REAL :: thetaesLCL, UgMIX, cape0,scape0,avvor,num, adjust_frac1D

   REAL    :: TST,PRS,RHOE,W0,SCR1,DXSQ,tmp,diff,wei1,wei2,wei3,wei4
   integer :: ibegh,iendh,jbegh,jendh
   integer :: istart,iend,jstart,jend,mup,mlow,m1,m2
   INTEGER :: i,j,k,NTST, sign_Lx, sign_Ly, SS,            &
              Iincre, IS, IE, Jincre, JS, JE,              &
              num_Lx, num_Ly, XX, YY, XS, XE, YS, YE,      &
              KLFCsw,KLNBsw,n, kmaxSCAPE1D,                &
              Ip1, Im1, Jp1, Jm1, Kp1, Km1, KPBL, KRHsat
   REAL    :: W0AVGfctr, W0fctr, W0den
   REAL    :: MISSING = -1000.
   REAL    :: PIE = 3.141592654
   REAL    :: mdis = 0.01
   REAL    :: f1D
   REAL    :: ALIQ,BLIQ,CLIQ,DLIQ,ES,QES,Q0
   REAL    :: x_prime, y_prime,dtdx_p, dtdy_p, llx, Ideci, & 
              lly,  Jdeci, dis_IS, dis_IE, dis_JS, dis_JE, &
              Ug_isjs, Ug_isje, Ug_iejs, Ug_ieje, dUgds,   &
              Aff_Lx, Aff_Ly, disa, disb,disc,disd, SLCL,  &
              dis1, dis2, dis3, dis4, disall, Sc_I, Sc_J,  &
              Idis, Jdis ,inc,DZ 
   REAL, PARAMETER :: kapa = 0.286
   REAL, PARAMETER :: Rd = 287.
   REAL, PARAMETER :: Rv = 461.6
   REAL, PARAMETER :: epsl = 0.622
   REAL, PARAMETER :: Cl =  4190
   REAL, PARAMETER :: Lv = 2.5e6
   REAL    ::RAIN, thetaesout, RHout, RHxsys, RHxeys, RHxsye, RHxeye, & 
                                      QVxsys, QVxeys, QVxsye, QVxeye, &
                                       Txsys,  Txeys,  Txsye,  Txeye
                         
   REAL    ::ZLCL,ang_M,scapefrac,GWLCL,VgLCL ,alpha,pathadjust,alphabest
   REAL :: layernum_1, layernum_2, layernum_3, avvor_1, avvor_2, avvor_3,    &
           pthespy_1, pthespy_2, pthespy_3, pthespz_1, pthespz_2, pthespz_3, &
           pMpz_1, pMpz_2, pMpz_3

!==== FOR 2D IDEALIZED TEST CASES, CHANGE BELOW ====
   twoD_nox_setup = .FALSE.
   twoD_noy_setup = .FALSE.
!===================================================
    
   ALIQ = SVP1*1000.
   BLIQ = SVP2
   CLIQ = SVP2*SVPT0
   DLIQ = SVP3

! KEPT FROM KF SCHEME-------------->
   NTST = STEPCU
   TST  = float(NTST*2)
!
   DXSQ=DX*DX

   W0AVGfctr = (TST-1.)
   W0fctr = 1.
   W0den = TST
   !endif
!

  DO J = jts,jte
      DO K=kts,kte
         DO I= its,ite

        IF (K.eq.kts) THEN
         RAINSCV(I,J)=0.
         SW_calls(I,J)=0.     ! If for the low-level air tube originated from
                              ! (I,J) does not show potential for slantwise
                              ! convection, the SW_calls=0. (the scheme is not
                              ! activated at this grid point)
        ENDIF
!
            UgTEND_IJ(I,K,J)=0.
            TTEND_IJ(I,K,J)=0.
            QTEND_IJ(I,K,J)=0.
!
            RUSLCUTEN(I,K,J)=0.
            RVSLCUTEN(I,K,J)=0.
            RTHSLCUTEN(I,K,J)=0.
            RQVSLCUTEN(I,K,J)=0.
!
            Kp1=min(K+1,kte)
            W0=0.5*(w(I,K,J)+w(I,Kp1,J))
!
            W0AVG(I,K,J) = ( W0AVG(I,K,J) * W0AVGfctr + W0 * W0fctr ) / W0den
!
         ENDDO
      ENDDO
   ENDDO
! KEPT FROM KF SCHEME<--------------


! To internally include the CPS tedencies for this given time step 
! to avoid double couting effects :

  DO J = max(jts-2,jds),min(jte+2,jde)
      DO K=kts,kte
         
         DO I= max(its-2,ids),min(ite+2,ide)
!
            TH(I,K,J)=TH0(I,K,J)+RTHCUTEN(I,K,J)*DT
            QV(I,K,J)=QV0(I,K,J)+RQVCUTEN(I,K,J)*DT
            T(I,K,J)=TH(I,K,J)*(Pcps(I,K,J)/1.E5)**kapa

            call RHcalc(QV(I,K,J),T(I,K,J),Pcps(I,K,J),RHout, thetaesout)
            RH(I,K,J)=RHout
            thetaes(I,K,J)=thetaesout
            Tv(I,K,J)=T(I,K,J)*(1.+0.61*QV(I,K,J))
!
         ENDDO
      ENDDO
   ENDDO
!

! For each "tile"-------------------------------------------->

  DO J = jts, jte
      DO I = its, ite

! Initializations

            CAPE(I,J)=0.
            SCAPE1D(I,J)=0.
            sw_trigger=.FALSE. 
            adjust_frac(I,J)=0.
         

            x_prime = 0.
            y_prime = 0.
            dtdx_p = 0.
            dtdy_p = 0.

! Calculating the vertically integrated temperature gradients:
            DO k=kts,kte
!
               Ip1=min(I+1,ide)
               Im1=max(I-1,ids)
               Jp1=min(J+1,jde)
               Jm1=max(J-1,jds)
               Kp1=min(K+1,kde)
               Km1=max(K-1,kds)
               Idis=(Ip1-Im1)*DX
               Jdis=(Jp1-Jm1)*DX
               dtdx_p = (Tv(Ip1,K,J)-Tv(Im1,K,J))/(Idis)-(Pcps(Ip1,K,J)-Pcps(Im1,K,J))/(Idis) &
                          *(Tv(I,Kp1,J)-Tv(I,Km1,J))/(Pcps(I,Kp1,J)-Pcps(I,Km1,J))
               dtdy_p = (Tv(I,K,Jp1)-Tv(I,K,Jm1))/(Jdis)-(Pcps(I,K,Jp1)-Pcps(I,K,Jm1))/(Jdis) &
                          *(Tv(I,Kp1,J)-Tv(I,Km1,J))/(Pcps(I,Kp1,J)-Pcps(I,Km1,J))
               !the new "y" axis (for the 2D plan of CSI, the horizonntal axis
               !of which is now termed as S) is defined as the 
               !negative column-averaged temperature gradient 
               !950-hPa to 350-hPa averages
               IF (Pcps(I,K,J).LE.95000 .AND. Pcps(I,K,J).GE.35000) THEN           
               x_prime = x_prime - dtdx_p
               y_prime = y_prime - dtdy_p
               ENDIF
!
            ENDDO
           
            IF ( twoD_nox_setup ) x_prime=0.
            IF ( twoD_noy_setup ) y_prime=0.
            ang_SW(I,J)=atan2(y_prime,x_prime)
            !The two-argument atan2 function computes the arctangent of y / x
            !given y and x, but with a range of (−π, π]. In other words,
            !atan2(y, x) is the angle between the positive x-axis of a plane
            !and the point (x, y) on it, with positive sign for
            !counter-clockwise angles (upper half-plane, y > 0), and negative
            !sign for clockwise angles (lower half-plane, y < 0). 

! ---------------------------
            sign_Lx=int(sign(1.,x_prime))  
            sign_Ly=int(sign(1.,y_prime))             
!
! Later, We will apply the DUg/DT (momentum adjustment) on the slanted thetaes surface               
! Here, SS=0 corresponds to the location on the Cartesian coordinate (I,J) 
!       SLCL corresponds to the horizontal shift of LCL from SS=0 along S 
!
! First, interpolate thetaes and ug fields in (I,K,J) to (S,K)

                DO SS=1, 2  ! first order forward finite differencing 
!
                  IF (SS.EQ.1) THEN !if SS=1, the original location at (I,J)
                   DO K = kts, kte
                   Sthetaes(SS,K) = thetaes(I,K,J) 
                   SUg(SS,K) = V(I,K,J)*sin(ang_SW(I,J)-0.5*PIE)+U(I,K,J)*cos(ang_SW(I,J)-0.5*PIE)
                   SP(SS,K)  = Pcps(I,K,J) 
                   ENDDO
            !
                  ELSE     
            !
                   llx=float(SS-1)*DX*cos(ang_SW(I,J))
                   Iincre=int(llx/DX)
                   Ideci=MOD(llx,DX)
                        IF (cos(ang_SW(I,J)).lt.0.001) THEN
                           IS=I+Iincre
                           IE=I+Iincre
                        ELSE
                           IS=I+Iincre
                           IE=I+Iincre+sign_Lx  
                        ENDIF
                   lly=float(SS-1)*DX*sin(ang_SW(I,J))
                   Jincre=int(lly/DX)
                   Jdeci=MOD(lly,DX)
                        IF (sin(ang_SW(I,J)).lt.0.001) THEN
                           JS=J+Jincre
                           JE=J+Jincre
                        ELSE
                           JS=J+Jincre
                           JE=J+Jincre+sign_Ly  
                        ENDIF
                   dis_IS = abs(llx-(IS-I)*DX)
                   dis_IE = abs(llx-(IE-I)*DX)
                   dis_JS = abs(lly-(JS-J)*DX)
                   dis_JE = abs(lly-(JE-J)*DX)

                   IF (IS.eq.IE) THEN
                   dis_IS=0.
                   dis_IE=0.
                   ELSEIF (JS.eq.JE) THEN
                   dis_JS=0.
                   dis_JE=0.
                   ENDIF

                        
                   IF (IS.ge.ids .AND. IS.le.ide .AND. &
                       JS.ge.jds .AND. JS.le.jde .AND. &
                       IE.ge.ids .AND. IE.le.ide .AND. &
                       JE.ge.jds .AND. JE.le.jde ) THEN
                         
                    DO K = kts, kte
                  !
                    Ug_isjs=V(IS,K,JS)*sin(ang_SW(I,J)-0.5*PIE)+U(IS,K,JS)*cos(ang_SW(I,J)-0.5*PIE)
                    Ug_ieje=V(IE,K,JE)*sin(ang_SW(I,J)-0.5*PIE)+U(IE,K,JE)*cos(ang_SW(I,J)-0.5*PIE)
                    Ug_isje=V(IS,K,JE)*sin(ang_SW(I,J)-0.5*PIE)+U(IS,K,JE)*cos(ang_SW(I,J)-0.5*PIE)
                    Ug_iejs=V(IE,K,JS)*sin(ang_SW(I,J)-0.5*PIE)+U(IE,K,JS)*cos(ang_SW(I,J)-0.5*PIE)
                  !
                       IF ( (dis_IS .ne. 0 .or. dis_IE .ne. 0 ) .and. &
                            (dis_JS .ne. 0 .or. dis_JE .ne. 0 ) ) THEN
                       !     if (I.EQ.2 .AND. J.EQ.10 .AND. K.EQ.20) print*,'origianl weights'
                  !
                    Sthetaes(SS,K) = (thetaes(IS,K,JS)*dis_IE*dis_JE & 
                                   +  thetaes(IE,K,JE)*dis_IS*dis_JS &
                                   +  thetaes(IS,K,JE)*dis_IE*dis_JS &
                                   +  thetaes(IE,K,JS)*dis_IS*dis_JE)&
                                      /(dis_IS+dis_IE)/(dis_JS+dis_JE)
                   !
                    SUg(SS,K)      = (Ug_isjs*dis_IE*dis_JE & 
                                   +  Ug_ieje*dis_IS*dis_JS &
                                   +  Ug_isje*dis_IE*dis_JS &
                                   +  Ug_iejs*dis_IS*dis_JE)&
                                      /(dis_IS+dis_IE)/(dis_JS+dis_JE)

                    SP(SS,K)       = (Pcps(IS,K,JS)*dis_IE*dis_JE & 
                                   +  Pcps(IE,K,JE)*dis_IS*dis_JS &
                                   +  Pcps(IS,K,JE)*dis_IE*dis_JS &
                                   +  Pcps(IE,K,JS)*dis_IS*dis_JE)&
                                      /(dis_IS+dis_IE)/(dis_JS+dis_JE)
                       ELSE
                          IF ( dis_IS+dis_IE+dis_JS+dis_JE .eq. 0) THEN
                             GOTO 100  
                          ELSEIF ( dis_IS .eq. 0 .and. dis_IE .eq. 0) THEN !IS=IE
                       !     if (I.EQ.2 .AND. J.EQ.10 .AND. K.EQ.20) print*,'loop 1'
                             Sthetaes(SS,K) = (thetaes(IS,K,JS)*dis_JE+ &
                                               thetaes(IS,K,JE)*dis_JS)/(dis_JE+dis_JS)
                             SUg(SS,K)      = (Ug_isjs*dis_JE+ &
                                               Ug_isje*dis_JS)/(dis_JE+dis_JS)
                             SP(SS,K)       = (Pcps(IS,K,JS)*dis_JE+ &
                                               Pcps(IS,K,JE)*dis_JS)/(dis_JE+dis_JS)
                          ELSEIF ( dis_JS .eq. 0 .and. dis_JE .eq. 0) THEN !JS=JE
                       !     if (I.EQ.2 .AND. J.EQ.10 .AND. K.EQ.20) print*,'loop 2'
                             Sthetaes(SS,K) = (thetaes(IS,K,JS)*dis_IE+ &
                                               thetaes(IE,K,JS)*dis_IS)/(dis_IE+dis_IS)
                             SUg(SS,K)      = (Ug_isjs*dis_IE+ &
                                               Ug_iejs*dis_IS)/(dis_IE+dis_IS)
                             SP(SS,K)       = (Pcps(IS,K,JS)*dis_IE+ &
                                               Pcps(IE,K,JS)*dis_IS)/(dis_IE+dis_IS)
                          ENDIF
                       ENDIF    
                  !
                    ENDDO
                   ELSE 
                     GOTO 100
                   ENDIF
            !
                  ENDIF 
                ENDDO ! S-loop
!
!
! Assign variables from 3D grid to 1D local column

            DO K=kts,kte
               U1D(K) =U(I,K,J)
               V1D(K) =V(I,K,J)
               T1D(K) =T(I,K,J)
               RHO1D(K) =rho(I,K,J)
               QV1D(K)=QV(I,K,J)
               P1D(K) =Pcps(I,K,J)
               W0AVG1D(K) =W0AVG(I,K,J)
               DZ1D(K)=dz8w(I,K,J)
               Ug1D(K)= V(I,K,J)*sin(ang_SW(I,J)-0.5*PIE)+U(I,K,J)*cos(ang_SW(I,J)-0.5*PIE)
               Vg1D(K)= V(I,K,J)*sin(ang_SW(I,J))+U(I,K,J)*cos(ang_SW(I,J))

            ! Calculate geostrophic vorticity:
            ! original formula is -partial ug/ partial S on pressure levels 
            ! where ug is the geostrophic wind component
            ! converted to eta levels: 
            ! par ug/ par s|p = (par ug/ par s|eta)-(par ug/ par p)*(par p/ par s|eta)
               if (K.eq. kte) then 
               dUgds = 0.
               else 
               dUgds = (SUg(2,K)-SUg(1,K))/DX-(SUg(1,K+1)-SUg(1,K))/(SP(1,K+1)-SP(1,K))* &
                       (SP(2,K)-SP(1,K))/DX
               endif
               vor1D(K) =  -1.*dUgds !vorticity on pressure levels. 
               RHavg(K)=0.
               QVavg(K)=0.
            ENDDO
        
!
            f1D=f(I,J)
              
! ENTERING A LOCAL 1D MODEL:
! Estimate the CAPE and SCAPE(using the 1D approximation),
! identify the trigger function based on where the low-level "tube" is lifted
! from, calculate the first-guess reference state as a fuction of k, and etc. 
!
            CALL Slantwise1D(I,J,                   &
                 f1D,vor1D, Ug1D, Vg1D,             &
                 U1D,V1D,T1D,QV1D,P1D,DZ1D,W0AVG1D, &
                 DT,DX,DXSQ,RHO1D,                  &
                 XLV0,XLV1,XLS0,XLS1,CP,R,G,        &
                 EP2,SVP1,SVP2,SVP3,SVPT0,          &
                 R_V, XLV,                          &
                 ids,ide, jds,jde, kds,kde,         &
                 ims,ime, jms,jme, kms,kme,         &
                 its,ite, jts,jte, kts,kte,         &
                 !output
                 cape0,scape0,scapefrac, RUg1,         &
                 DUgDT,thetaesLCL,ZLCL,KLFCsw,KLNBsw,  & 
                 KPBL,overshoot_k, KLCLn1,kmaxSCAPE1D, &
                 UgMIX, avvor, adjust_frac1D,          & 
                 sw_trigger,Z0,ang_M,GWLCL, VgLCL, TIMEC,SLCUDT ) 
                
!
! Applying back to the 3D model grid (modification/ fine-tuning is required):

      IF (sw_trigger) THEN

              CAPE(I,J)=cape0
              SCAPE1D(I,J)=scape0
              adjust_frac(I,J)=adjust_frac1D
!              
            !!-->Calculate the averaged thetaes* and Mg slopes
                layernum_1 = 0.
                pthespy_1  = 0.
                pthespz_1  = 0.
                pMpz_1     = 0.
                avvor_1    = 0.
                layernum_2 = 0.
                pthespy_2  = 0.
                pthespz_2  = 0.
                pMpz_2     = 0.
                avvor_2    = 0.
               ! layernum_3 = 0.
               ! pthespy_3  = 0.
               ! pthespz_3  = 0.
               ! pMpz_3     = 0.
               ! avvor_3    = 0.
             DO K = 2, overshoot_k-1
                SLOPE_M(K) = 0.
                SLOPE_T(K) = 0.
              IF (Z0(K).ge.Z0(KPBL) .and. Z0(K).le.Z0(kmaxSCAPE1D)) THEN  
                  layernum_1 = layernum_1+1.                          
                  pthespy_1  = pthespy_1+(Sthetaes(2,K)-Sthetaes(1,K))/DX
                  pthespz_1  = pthespz_1+(Sthetaes(1,K+1)-Sthetaes(1,K-1))/(Z0(K+1)-Z0(K-1))
                  pMpz_1     = pMpZ_1+(Ug1D(K+1)-Ug1D(K-1))/(Z0(K+1)-Z0(K-1))
                  avvor_1    = avvor_1+vor1D(K)
              ENDIF
              IF (Z0(K).ge.Z0(kmaxSCAPE1D) .and. Z0(K).lt.Z0(KLNBsw)) THEN  
                  layernum_2 = layernum_2+1.                          
                  pthespy_2  = pthespy_2+(Sthetaes(2,K)-Sthetaes(1,K))/DX
                  pthespz_2  = pthespz_2+(Sthetaes(1,K+1)-Sthetaes(1,K-1))/(Z0(K+1)-Z0(K-1))
                  pMpz_2     = pMpZ_2+(Ug1D(K+1)-Ug1D(K-1))/(Z0(K+1)-Z0(K-1))
                  avvor_2    = avvor_2+vor1D(K)
              ENDIF
              !IF (Z0(K).ge.Z0(KLFCsw) .and. Z0(K).le. Z0(KLNBsw)) THEN  
              !    layernum_3 = layernum_3+1.                          
              !    pthespy_3  = pthespy_3+(Sthetaes(2,K)-Sthetaes(1,K))/DX
              !    pthespz_3  = pthespz_3+(Sthetaes(1,K+1)-Sthetaes(1,K-1))/(Z0(K+1)-Z0(K-1))
              !    pMpz_3     = pMpZ_3+(Ug1D(K+1)-Ug1D(K-1))/(Z0(K+1)-Z0(K-1))
              !    avvor_3    = avvor_3+vor1D(K)
              !ENDIF
             ENDDO
!
                if (layernum_1.gt.0.) then
                pthespy_1  = pthespy_1/layernum_1
                pthespz_1  = pthespz_1/layernum_1
                pMpz_1     = pMpz_1/layernum_1
                avvor_1    = avvor_1/layernum_1+f1D
                !---
                  IF ( abs(f1D/(f1D+avvor_1)).GT.3.0 ) THEN
                     IF (avvor_1.ge.0.0)then
                      avvor_1=amax1(abs(f1D),abs(avvor_1))
                     ELSE !avor negative
                      avvor_1=-1.*amax1(abs(f1D),abs(avvor_1))
                     ENDIF
                  ELSE
                      avvor_1=avvor_1+f1D
                  ENDIF
                !---
                endif

                if (layernum_2.gt.0.) then
                pthespy_2  = pthespy_2/layernum_2
                pthespz_2  = pthespz_2/layernum_2
                pMpz_2     = pMpz_2/layernum_2
                avvor_2    = avvor_2/layernum_2+f1D
                !---
                  IF ( abs(f1D/(f1D+avvor_2)).GT.3.0 ) THEN
                     IF (avvor_2.ge.0.0)then
                      avvor_2=amax1(abs(f1D),abs(avvor_2))
                     ELSE !avor negative
                      avvor_2=-1.*amax1(abs(f1D),abs(avvor_2))
                     ENDIF
                  ELSE
                      avvor_2=avvor_2+f1D
                  ENDIF
                !---
                endif
                !if (layernum_3.gt.0.) then
                !pthespy_3  = pthespy_3/layernum_3
                !pthespz_3  = pthespz_3/layernum_3
                !pMpz_3     = pMpz_3/layernum_3
                !avvor_3    = avvor_3/layernum_3
                !---
                !  IF ( abs(f1D/(f1D+avvor_1)).GT.3.0 ) THEN
                !     IF (avvor_1.ge.0.0)then
                !      avvor_1=amax1(abs(f1D),abs(avvor_1))
                !     ELSE !avor negative
                !      avvor_1=-1.*amax1(abs(f1D),abs(avvor_1))
                !     ENDIF
                !  ELSE
                !      avvor_1=avvor_1+f1D
                !  ENDIF
                !---
                !endif
            !!Calculate the averaged thetaes* and Mg slopes<---
!            
!!--> Adjustment path estimation (first guess)--- 
        pathadjustdone=.FALSE.
        alpha=0.5
        alphabest=0.5
        pathadjust=0.

300     DO K = 1, overshoot_k

           alpha=AMIN1(alpha,1.)
           alpha=AMAX1(alpha,0.)
           SHIFT_IN_AD0(K)=0.
           inc=0.

             IF (Z0(K).le.Z0(KLFCsw)) THEN

                   SHIFT_IN_AD0(K)=0. ! assume the low-level mixture (usl)
                                      ! remains at the origin of grid points
             ELSE 

                IF ( Z0(K).gt.Z0(KLFCsw) .and. Z0(K).lt.Z0(kmaxSCAPE1D)+0.1*(Z0(KLNBsw)-Z0(kmaxSCAPE1D))) THEN  
                  SLOPE_M(K) = pMpz_1/avvor_1
                  SLOPE_T(K) = -1.*pthespz_1/pthespy_1
                  DZ=Z0(K)-Z0(KLFCsw)
                   IF (SLOPE_T(K).ge.SLOPE_M(K)) THEN
                     inc = (1.0*SLOPE_M(K))*DZ
                   ELSEIF (pthespz_1.le. 5.e-4 ) THEN
                     inc = ((1.3-alpha)*SLOPE_T(K)+(alpha-0.3)*SLOPE_M(K))*DZ 
                   ELSE
                     inc = ((1.-alpha)*SLOPE_T(K)+alpha*SLOPE_M(K))*DZ 
                   ENDIF 
                ELSE
                  SLOPE_M(K) = pMpz_2/avvor_2
                  SLOPE_T(K) = -1.*pthespz_2/pthespy_2
                  DZ=Z0(K)-Z0(KLFCsw)
                   IF (SLOPE_T(K).ge.SLOPE_M(K)) THEN
                     inc = (1.0*SLOPE_M(K))*DZ 
                   ELSEIF (pthespz_2.le. 5.e-4) THEN
                     inc = ((1.3-alpha)*SLOPE_T(K)+(alpha-0.3)*SLOPE_M(K))*DZ 
                   ELSE
                     inc = ((1.-alpha)*SLOPE_T(K)+alpha*SLOPE_M(K))*DZ 
                   ENDIF 
                ENDIF
                !
                IF (inc.lt.0.) THEN
                     SHIFT_IN_AD0(K) = 0.
                ELSE
                     SHIFT_IN_AD0(K) = inc
                ENDIF
                 
             ENDIF 

             IF (K.gt.1 .and. abs(SHIFT_IN_AD0(K)-SHIFT_IN_AD(K-1)).gt.DX) THEN
                SHIFT_IN_AD(K) = 0.5*(SHIFT_IN_AD0(K)+SHIFT_IN_AD(K-1))
             ELSE
                SHIFT_IN_AD(K) = SHIFT_IN_AD0(K)
             ENDIF
!!--> Adjustment path estimation (first guess) --- 
!
! Last, for each K, we need to extrapolate the tendencies back to the four cloestest grid
! points to SHIFT_IN_S(K) in the (I,J) coordinates
!
             IF (cos(ang_SW(I,J)).lt.0.001) THEN
                 llx=0.
                 XS=I
                 XE=I
             ELSE 
                 llx=SHIFT_IN_AD(K)*cos(ang_SW(I,J)) !in meter; the distance
                                                     !of this point to the 
                                                     !origianl starting point(I,J); SS=0
                 ! E.g., if SHIFT_IN_S(K)=5.3 along S, meaning that it is 5.3*DX
                 ! meters away from the grid (I,J). 
                 ! IF and_SW(I,J) is 90-180 degree, llx will be negative (say
                 ! -4.6*DX meters, indicating the horizontal distance from (I) along x.
                 ! Then, we know that the cloest two I along x to this grid S, are I-4 and I-5.
                 ! int(llx/DX)=-4., sign_Lx=-1.
                 ! 
                 XS=I+int(llx/DX) !one of the two closest grid points in the I direction
                 XE=I+int(llx/DX)+sign_Lx !one of the two closest grid points in the I direction
             ENDIF

             IF (sin(ang_SW(I,J)).lt.0.001) THEN
                 lly=0.
                 YS=J
                 YE=J
             ELSE 
                 lly=SHIFT_IN_AD(K)*sin(ang_SW(I,J)) 
                 YS=J+int(lly/DX) !one of the two closest grid points in the I direction
                 YE=J+int(lly/DX)+sign_Ly !one of the two closest grid points in the I direction
             ENDIF
                 !
             ! The real (non-rounded to integers) (I,J) of this center (SHIFT_IN_S(K)) is actually...
             Sc_I=float(I)+llx/DX  !DX is real
             Sc_J=float(J)+lly/DX
             disa=0.
             disb=0.
             disc=0.
             disd=0.
             IF (XS.ge.ids .AND. XS.le.ide .AND. YS.ge.jds .AND. YS.le.jde .AND. &
                 XE.ge.ids .AND. XE.le.ide .AND. YE.ge.jds .AND. YE.le.jde ) THEN
                  disa=abs(float(XS)-Sc_I)
                  disb=1.-disa
                  disc=abs(float(YS)-Sc_J)
                  disd=1.-disc
                  QVxsys= QV0(XS,K,YS)+RQVCUTEN(XS,K,YS)*DT
                  QVxeys= QV0(XE,K,YS)+RQVCUTEN(XE,K,YS)*DT
                  QVxsye= QV0(XS,K,YE)+RQVCUTEN(XS,K,YE)*DT
                  QVxeye= QV0(XE,K,YE)+RQVCUTEN(XE,K,YE)*DT
                  Txsys= (TH0(XS,K,YS)+RTHCUTEN(XS,K,YS)*DT)*(Pcps(XS,K,YS)/1.E5)**kapa
                  Txeys= (TH0(XE,K,YS)+RTHCUTEN(XE,K,YS)*DT)*(Pcps(XE,K,YS)/1.E5)**kapa
                  Txsye= (TH0(XS,K,YE)+RTHCUTEN(XS,K,YE)*DT)*(Pcps(XS,K,YE)/1.E5)**kapa
                  Txeye= (TH0(XE,K,YE)+RTHCUTEN(XE,K,YE)*DT)*(Pcps(XE,K,YE)/1.E5)**kapa
                  !print*,'QVxsys=',QVxsys,'Txsys=',Txsys,'QVxeye=',QVxeye,'Txeye=',Txeye

                  call RHcalc(QVxsys,Txsys,Pcps(XS,K,YS),RHxsys, thetaesout)
                  call RHcalc(QVxeys,Txeys,Pcps(XE,K,YS),RHxeys, thetaesout)
                  call RHcalc(QVxsye,Txsye,Pcps(XS,K,YE),RHxsye, thetaesout)
                  call RHcalc(QVxeye,Txeye,Pcps(XE,K,YE),RHxeye, thetaesout)
                  !print*,'RHxsys=',RHxsys,'RHxeys=',RHxeys,'RHxsye=',RHxsye,'RHxeye=',RHxeye

                  IF ((XS.ne.XE).AND.(YS.ne.YE)) THEN
                    
                    RHavg(K)= RHxsys*disb*disd+RHxeys*disa*disd+ &
                         &    RHxsye*disb*disc+RHxeye*disa*disc

                    QVavg(K)= QVxsys*disb*disd+QVxeys*disa*disd+ &
                         &    QVxsye*disb*disc+QVxeye*disa*disc
                    !Wavg(K) = W0AVG(XS,K,YS)*disb*disd+W0AVG(XE,K,YS)*disa*disd+ &
                    !      &    W0AVG(XS,K,YE)*disb*disc+W0AVG(XE,K,YE)*disa*disc
                  ELSEIF (XS.eq.XE)THEN
                    RHavg(K)= RHxsys*disd+RHxsye*disc
                    QVavg(K)= QVxsys*disd+QVxsye*disc
                    !Wavg(K) = W0AVG(XS,K,YS)*disd+W0AVG(XS,K,YE)*disc
                  ELSEIF (YS.eq.YE)THEN
                    RHavg(K)= RHxsys*disb+RHxeys*disa
                    QVavg(K)= QVxsys*disb+QVxeys*disa
                    !Wavg(K) = W0AVG(XS,K,YS)*disb+W0AVG(XE,K,YS)*disa
                  ENDIF

                  ! Checking if fine-tuning of the adjustment path is necessary
                  ! ..............
                  IF (.not.pathadjustdone .and. Z0(K).gt.Z0(KLFCsw) .and. &
                      Z0(K).le.Z0(KLNBsw) .and. SHIFT_IN_AD(K).ne.0.) THEN  
                  !==== Fine-tuning of the adjustment path is required!!!!! ====
                  IF (K.gt.1 .and. RHavg(K-1).ge.0.9 .and. RHavg(K).lt.0.9) THEN
                    
                    pathadjust=pathadjust+1.

                    if (pathadjust.le.1.1) THEN

                        markKbest=K
                        alphabest=alpha
                        alpha=alpha+0.1
                        goto 300

                    else !adjust more than one time
                    
                      IF (K.gt. markKbest) THEN 
                           markKbest=K
                           alphabest=alpha
                      ENDIF 
                      alpha=alpha+0.1
                              
                      IF (alpha.lt.1.1) THEN 
                         goto 300
                      ELSEIF (alpha.ge.1.1) THEN
                         alpha=alphabest
                         pathadjustdone=.TRUE.
                         goto 300
                      ENDIF

                    endif !(pathadjust.le.1.1)
                  ENDIF
                  ENDIF
                  ! ..............
                  !
             ENDIF !WITHIN DOMAIN
        ENDDO !K Loop 300  
        
        ! ===== Call the Subroutine for parameterizing the latent effects =====
        ! =====        (heating, mositure loss and precipitation)         ===== 
        
             call DTDQV_adjust(DTDT,DQDT,RAIN,scape0,GWLCL,scapefrac,         &
                               ang_M, ZLCL,KLCLn1,KLFCsw,kmaxSCAPE1D,KLNBsw,  &
                               RHavg,QVavg,T1D,P1D,Z0,SLCUDT,                 &
                               DT,ALIQ,BLIQ,CLIQ,DLIQ,XLV,CP,R_V,R,G,EP2,     &
                               kts,kte,ims,ime,jms,jme,I,J, KRHsat,skipflag,  &
                               DUgDT_RH, RUg1, Ug1D, TIMEC, overshoot_k,Z1) 
         
             IF (KRHsat.lt.KLFCsw .or. skipflag) THEN
                 goto 100
             ENDIF

         ! Finally, distribute theta, qv, U, V tendencies to the surrouding
         ! (closest) four model grid points (at a given K) along the adjustment
         ! path. The weight for each grid point is inversely proportional to
         ! their distance to the path.
    
           DO K = 1, overshoot_k
                        
                IF (cos(ang_SW(I,J)).lt.0.001) THEN
                   llx=0.
                   XS=I
                   XE=I
                ELSE 
                   llx=SHIFT_IN_AD(K)*cos(ang_SW(I,J)) !in meter; the distance
                   XS=I+int(llx/DX) !one of the two closest grid points in the I direction
                   XE=I+int(llx/DX)+sign_Lx !one of the two closest grid points in the I direction
                ENDIF
                   !
                IF (sin(ang_SW(I,J)).lt.0.001) THEN
                   lly=0.
                   YS=J
                   YE=J
                ELSE 
                   lly=SHIFT_IN_AD(K)*sin(ang_SW(I,J)) 
                   YS=J+int(lly/DX) !one of the two closest grid points in the I direction
                   YE=J+int(lly/DX)+sign_Ly !one of the two closest grid points in the I direction
                ENDIF
                !the real (non-rounded to integers) (I,J) of this center (SHIFT_IN_S(K)) is actually...
                Sc_I=float(I)+llx/DX  !DX is real
                Sc_J=float(J)+lly/DX
                 
                   wei1=0.
                   wei2=0.
                   wei3=0.
                   wei4=0.
                   disall=0.
                   dis1=0.
                   dis2=0.
                   dis3=0.
                   dis4=0

                IF (XS.ge.ids .AND. XS.le.ide .AND. YS.ge.jds .AND. YS.le.jde .AND. &
                    XE.ge.ids .AND. XE.le.ide .AND. YE.ge.jds .AND. YE.le.jde ) THEN

                    IF ((XS.ne.XE).AND.(YS.ne.YE)) THEN
                     dis1=sqrt((float(XS)-Sc_I)**2.+(float(YS)-Sc_J)**2.)
                     dis2=sqrt((float(XE)-Sc_I)**2.+(float(YS)-Sc_J)**2.)
                     dis3=sqrt((float(XS)-Sc_I)**2.+(float(YE)-Sc_J)**2.)
                     dis4=sqrt((float(XE)-Sc_I)**2.+(float(YE)-Sc_J)**2.)
                       if (dis1.gt.mdis .and. dis2.gt.mdis .and. dis3.gt.mdis .and. dis4.gt.mdis) then
                           disall=float(nint(1./dis1*1.e3))    &
                                 +float(nint(1./dis2*1.e3))    &
                                 +float(nint(1./dis3*1.e3))    &
                                 +float(nint(1./dis4*1.e3))
                           wei1=float(nint(1./dis1*1.e3/disall*100.))/100.
                           wei2=float(nint(1./dis2*1.e3/disall*100.))/100.
                           wei3=float(nint(1./dis3*1.e3/disall*100.))/100.
                           !wei4=float(nint(1./dis4*1000.))/disall
                           wei4=1.-wei1-wei2-wei3
                       else
                         if (dis1.le.mdis) then
                             wei1=1.
                             wei2=0.
                             wei3=0.
                             wei4=0.
                         elseif (dis2.le.mdis) then
                             wei1=0.
                             wei2=1.
                             wei3=0.
                             wei4=0.
                         elseif (dis3.le.mdis) then
                             wei1=0.
                             wei2=0.
                             wei3=1.
                             wei4=0.
                         elseif (dis4.le.mdis) then
                             wei1=0.
                             wei2=0.
                             wei3=0.
                             wei4=1.
                         endif
                       endif
                    ELSEIF (XS.eq.XE) THEN
                      dis1=sqrt((float(YS)-Sc_J)**2.)
                      dis3=sqrt((float(YE)-Sc_J)**2.)
                      disall=float(nint(1./dis1*1.e3))+float(nint(1./dis3*1.e3))
                      if (dis1.ne.0.) THEN
                      wei1=float(nint(1./dis1*1.e3/disall*100.))/100.
                      else
                      wei1=1.
                      endif
                      wei3=1.-wei1
                      wei2=wei1
                      wei4=wei3
                    ELSEIF (YS.eq.YE) THEN
                      dis1=sqrt((float(XS)-Sc_I)**2.)
                      dis2=sqrt((float(XE)-Sc_I)**2.)
                      disall=float(nint(1./dis1*1.e3))+float(nint(1./dis2*1.e3))
                      if (dis1.ne.0.) THEN
                      wei1=float(nint(1./dis1*1.e3/disall*100.))/100.
                      else
                      wei1=1.
                      endif
                      wei2=1.-wei1
                      wei3=wei1
                      wei4=wei2
                    ENDIF
     
                    ! Note that for the standard output, we use DUgDT_RH rather
                    ! than DUgDT. The latter is the value derived following Ma
                    ! (2000), which implicitly assumes that the momentum
                    ! adjustment is applied from near the surface to LSNB.
                    ! On the other hand, DUgDT_RH considers the current
                    ! developing extent in height of slantwise conveciton 
                    ! (based on the grid-resolved RH). Our idealized test
                    ! results suggest the former more realistically reflect the
                    ! neutralization process by the CSI release.

                    UgTEND_IJ(XS,K,YS)=DUgDT_RH(K)*wei1
                    QTEND_IJ(XS,K,YS)=DQDT(K)*wei1
                    TTEND_IJ(XS,K,YS)=DTDT(K)*wei1
                    IF(Z0(K+1).gt.Z1 .AND. Z0(K).le.Z1 ) RAINincrease(XS,YS)=RAIN*wei1
     
                    UgTEND_IJ(XE,K,YS)=DUgDT_RH(K)*wei2
                    QTEND_IJ(XE,K,YS)=DQDT(K)*wei2
                    TTEND_IJ(XE,K,YS)=DTDT(K)*wei2
                    IF(Z0(K+1).gt.Z1 .AND. Z0(K).le.Z1 ) RAINincrease(XE,YS)=RAIN*wei2

                    UgTEND_IJ(XS,K,YE)=DUgDT_RH(K)*wei3
                    QTEND_IJ(XS,K,YE)=DQDT(K)*wei3
                    TTEND_IJ(XS,K,YE)=DTDT(K)*wei3
                    IF(Z0(K+1).gt.Z1 .AND. Z0(K).le.Z1 ) RAINincrease(XS,YE)=RAIN*wei3

                    UgTEND_IJ(XE,K,YE)=DUgDT_RH(K)*wei4 
                    QTEND_IJ(XE,K,YE)=DQDT(K)*wei4
                    TTEND_IJ(XE,K,YE)=DTDT(K)*wei4 
                    IF(Z0(K+1).gt.Z1 .AND. Z0(K).le.Z1 ) RAINincrease(XE,YE)=RAIN*wei4

                    DO XX=XS, XE, sign_Lx
                    DO YY=YS, YE, sign_Ly

                       RUSLCUTEN(XX,K,YY)=RUSLCUTEN(XX,K,YY)+ &
                                          UgTEND_IJ(XX,K,YY)*cos(ang_SW(I,J)-0.5*PIE)   
                       RVSLCUTEN(XX,K,YY)=RVSLCUTEN(XX,K,YY)+ &
                                          UgTEND_IJ(XX,K,YY)*sin(ang_SW(I,J)-0.5*PIE)   
                       RTHSLCUTEN(XX,K,YY)=RTHSLCUTEN(XX,K,YY)+TTEND_IJ(XX,K,YY)
                       RQVSLCUTEN(XX,K,YY)=RQVSLCUTEN(XX,K,YY)+QTEND_IJ(XX,K,YY) 
                       IF(Z0(K+1).gt.Z1 .AND. Z0(K).le.Z1 ) RAINSCV(XX,YY)=RAINSCV(XX,YY)+RAINincrease(XX,YY)
                       IF(K.eq.1) SW_calls(XX,YY)=1.

                    ENDDO
                    ENDDO

                ENDIF !WITHIN DOMAIN
                 
           ENDDO !K-loop

      ENDIF !IF:sw_trigger


100 continue


      ENDDO       ! i-loop
  ENDDO           ! j-loop 
!<---------------------------------------- For each "tile"


  END SUBROUTINE slantwise_adjust

! ***************************************************************************

  SUBROUTINE Slantwise1D  (I, J, f, vor, Ug0, Vg0,              &
                           U0,V0,T0,QV0,P0,DZQ,W0AVG1D,         &
                           DT,DX,DXSQ,rhoe,                     &
                           XLV0,XLV1,XLS0,XLS1,CP,R,G,          &
                           EP2,SVP1,SVP2,SVP3,SVPT0,            &
                           R_V, XLV,                            &
                           ids,ide, jds,jde, kds,kde,           &
                           ims,ime, jms,jme, kms,kme,           &
                           its,ite, jts,jte, kts,kte,           &
                           CAPE,SCAPE1D, SCAPEfrac, RUg1,       &
                           DUgDT,thetaesLCL, ZLCL, KLFCsw,KLNBsw,   &
                           KPBL, overshoot_k, KLCLn1, kmaxSCAPE1D,  &
                           UgMIX, avvor, adjust_frac,sw_trigger,Z0, &
                           ang_M, GWLCL, VgLCL,TIMEC,SLCUDT)
!
! This subroutine mainly estimates CAPE, SCAPE(using the 1D formula), 
! and the reference profile of Ug to elimate SCAPE.
! Note that the above calculation is only done IF the low-level 
! lifted air tube is found "buoyant" (the trigger function is satisfied).
!
      IMPLICIT NONE
      INTEGER, INTENT(IN   ) :: ids,ide, jds,jde, kds,kde, &
                                ims,ime, jms,jme, kms,kme, &
                                its,ite, jts,jte, kts,kte, &
                                I,J

      REAL, DIMENSION( kts:kte ),                          &
            INTENT(IN   ) ::                           U0, &
                                                       V0, &
                                                       T0, &
                                                      QV0, &
                                                       P0, &
                                                     rhoe, &
                                                      DZQ, &
                                                  W0AVG1D, &
                                                      vor, &
                                                      Ug0, &
                                                      Vg0
!
      REAL , INTENT(IN   ) :: f, SLCUDT
      REAL , INTENT(  OUT) :: avvor, adjust_frac,GWLCL,VgLCL
      REAL,  INTENT(IN   ) :: DT,DX,DXSQ
!

      REAL,  INTENT(IN   ) :: XLV0,XLV1,XLS0,XLS1,CP,R,G,XLV
      REAL,  INTENT(IN   ) :: EP2,SVP1,SVP2,SVP3,SVPT0,R_V

!

      REAL, DIMENSION( kts:kte ), INTENT(  OUT) ::         &
                                                    DUgDT, &
                                                       Z0, &
                                                     RUg1

      REAL,  INTENT(INOUT) ::   CAPE, SCAPE1D, thetaesLCL, &  
                                                    UgMIX

      INTEGER,  INTENT(   OUT) :: KLFCsw, KLNBsw, KPBL,    & 
                                  KLCLn1, kmaxSCAPE1D,     &
                                  overshoot_k
      REAL,  INTENT(   OUT)   :: ZLCL, ang_M, SCAPEfrac,TIMEC
      LOGICAL,  INTENT(INOUT) :: sw_trigger
!
!...DEFINE LOCAL VARIABLES...
!
      REAL, DIMENSION( kts:kte ) ::                        &
            Q0,TV0,TU,TVU,QU,QES,WU,                       &
            DZA,DPA,THETEU,QLIQ,QICE                   

      REAL, DIMENSION( kts:kte ) ::                        &
            TVQU,DP,RH,WSPD,check_scape,check_cape

      REAL, DIMENSION( kts:kte ) ::                        &
            CLDHGT,CLDHGTsw

! LOCAL VARS

      REAL    :: P00,T00,RLF,PIE,UgMIX2,overshoot_z
      REAL    :: GDRY,ROCP,ALIQ,BLIQ,CLIQ,DLIQ             
      REAL    :: P300,P150,DPTHMX,THMIX,QMIX,ZMIX,PMIX,    &
                 TMIX,EMIX,TLOG,TDPT,TLCL,TVLCL,           &
                 PLCL,ES,DLP,TENV,QENV,TVEN,               &
                 WKL,WABS,WSIGNE,DTLCL,GDT,WLCL,           &
                 WTW,RHOLCL,WKLCL,QNEWIC,                  &
                 ENTERM,DZZ,DPP,QNEWLQ,VCONV,              &
                 QS,TADVEC,A1,AINC,DPMIN,                  &
                 DILBE, vorLCL, UgLCL,TV0M, TV0M2, DILBE2, &
                 SCAPEres, Ugres, RUg0, SCAPE1D_ad,        &
                 delTvt, delp, sumdelp,rhs,vor1,           &
                 delTvt2, RUgLCL, wslw, midcloudZ
                 
   REAL    ::    ASTRT,TP,VALUE,AINTRP,QFRZ,maxSCAPE1D,Knum
!
   INTEGER :: INDLU,NU,KLFC,KLNB,NKn1 !, midcloudK
   REAL    :: CHMIN,PM15,CHMAX,DTRH,DPPP,P0diff,Z1,Z2
   REAL    :: DTTOT,wmax,slope,Lscale,weightsum,weightsum2

      INTEGER :: m,m1,n,KX,K,KL,K150,KTrop,KLCL
!
      INTEGER :: NCHECK
      INTEGER, DIMENSION (kts:kte) :: KCHECK
      REAL, DIMENSION (kts:kte) ::  RUg, SDILBE,        &
                                    SDILBE2, lapse,     &
                                    TPA  !, RQ,RT 

      INTEGER :: L5,KMIX,LOW,LC,LLFC,                   &
                 NLAYRS,NK,LCL,                         &
                 NK1,Nlayer_klow,                       &
                 Nlayer_kupp,iter,condmaxK,             &
                 kmidSCAPE1D 
      LOGICAL, DIMENSION (kts:kte) :: readjust
      REAL :: u00,qslcl,rhlcl,dqssdt    !jfb
      CHARACTER*1024 message
!
      DATA P00,T00/1.E5,273.16/
      DATA RLF/3.339E5/
      DATA PIE/3.141592654/
!-----------------------------------------------------------
      GDRY=-G/CP
      ROCP=R/CP
      KL=kte
      KX=kte
!
      ALIQ = SVP1*1000.
      BLIQ = SVP2
      CLIQ = SVP2*SVPT0
      DLIQ = SVP3
!
      DPMIN = 5.E3
!...
      P300=P0(1)-30000.
      P150=P0(1)-15000.
      sw_trigger=.FALSE.
      adjust_frac=-1000. !missing value
!
      DO NK=2,KL-1
         IF (P0(NK).ge.1.5e4 .and. P0(NK+1).lt.1.5e4) THEN
             K150=NK
         ENDIF
      ENDDO
!
      DO K=1,KX
!
!  Saturation vapor pressure (ES) is calculated following Buck (1981)
!...IF Q0 IS ABOVE SATURATION VALUE, REDUCE IT TO SATURATION LEVEL...
!
         ES=ALIQ*EXP((BLIQ*T0(K)-CLIQ)/(T0(K)-DLIQ))
         QES(K)=0.622*ES/(P0(K)-ES)
         Q0(K)=AMIN1(QES(K),QV0(K))
         Q0(K)=AMAX1(0.000001,Q0(K))
         RUg1(K)=Ug0(K) !set as the same as original Ug in the beginning 
         RUg(K)=Ug0(K) !set as the same as original Ug in the beginning 
         RH(K) = Q0(K)/QES(K)
         TV0(K)=T0(K)*(1.+0.608*Q0(K))
!   DP IS THE PRESSURE INTERVAL BETWEEN FULL SIGMA LEVELS...
         DP(K)=rhoe(k)*g*DZQ(k)
         CLDHGT(K) = 0.
         CLDHGTsw(K) = 0.
         DUgDT(K) = 0. 
         check_cape(K) = 0. 
         check_scape(K) = 0. 
         IF(P0(K).GE.0.5*P0(1))L5=K
         IF(P0(K).GE.P150)LLFC=K
         SDILBE(K)=0.
         SDILBE2(K)=0.
      ENDDO

!
!...DZQ IS DZ BETWEEN SIGMA SURFACES, DZA IS DZ BETWEEN MODEL HALF LEVEL
        Z0(1)=.5*DZQ(1)
!cdir novector
        DO K=2,KL
          Z0(K)=Z0(K-1)+.5*(DZQ(K)+DZQ(K-1))
          DZA(K-1)=Z0(K)-Z0(K-1)
          DPA(K-1)=P0(K)-P0(K-1)
          lapse(K) = 10. !random initial value 
        ENDDO   
        DZA(KL)=0.
        DPA(KL)=0.
!
! Find the K level for the thermal Tropopause (KTrop)-->
      KTrop = K150
      DO K = 2 ,KL-1
         lapse(K) = -1.*(TV0(K+1)-TV0(K-1))/(Z0(K+1)-Z0(K-1))
      ENDDO
      DO K = KL-1, 2, -1
            IF (lapse(K).GE.2.e-3) THEN
               KTrop=K-1
               EXIT
            ENDIF
      ENDDO
!<-------------------- 
!
      avvor=0.
      Knum=0.
      DO K = 1, KTrop
         IF (Z0(K).lt.8.e3) THEN !only take the average below 8 km height 
         avvor=avvor+vor(K) 
         Knum=Knum+1.
         ENDIF
      ENDDO 
         avvor=avvor/Knum
!   
      IF ( abs(f/(f+avvor)).GT.3.0 ) THEN
        IF (avvor.ge.0.0)then
           avvor=amax1(abs(f),abs(avvor))
        ELSE !avor negative
           avvor=-1.*amax1(abs(f),abs(avvor))
        ENDIF
      ELSE
        avvor=avvor+f
      ENDIF
 
! Below is the "low-level" unstable parcel (or tube to be more specific here) 
! searching method. Taken from the KF scheme.
! ************************************************************************
!  To save time, specify a pressure interval to move up in sequential
!  check of different ~50 mb deep groups of adjacent model layers in
!  the process of identifying updraft source layer (USL).  Note that 
!  this search is terminated as soon as a buoyant parcel is found and 
!  this parcel can produce a cloud greater than specifed minimum depth
!  (CHMIN)...For now, set interval at 15 mb...
!
       NCHECK = 1
       KCHECK(NCHECK)=1
       PM15 = P0(1)-15.E2
       DO K=2,LLFC
         IF(P0(K).LT.PM15)THEN
           NCHECK = NCHECK+1
           KCHECK(NCHECK) = K
           PM15 = PM15-15.E2
         ENDIF
       ENDDO
!
       NU=0
usl:   DO


           NU = NU+1
           IF (NU.GT.NCHECK) THEN
              RETURN    !if exceeding the LLFC (300 hPa from the surface), 
                        !no USL is found. 
           ENDIF
           KMIX = KCHECK(NU)
           LOW=KMIX
!
           LC = LOW
!
!...ASSUME THAT IN ORDER TO SUPPORT A DEEP UPDRAFT YOU NEED A LAYER OF
!...UNSTABLE AIR AT LEAST 50 mb DEEP...TO APPROXIMATE THIS, ISOLATE A
!...GROUP OF ADJACENT INDIVIDUAL MODEL LAYERS, WITH THE BASE AT LEVEL
!...LC, SUCH THAT THE COMBINED DEPTH OF THESE LAYERS IS AT LEAST 50 mb..
!   
           NLAYRS=0
           DPTHMX=0.
           NK=LC-1
           IF ( NK+1 .LT. KTS ) THEN
             WRITE(message,*)'WOULD GO OFF BOTTOM: Slantwise I,J,NK',I,J,NK
             CALL wrf_message (TRIM(message)) 
           ELSE
             DO 
               NK=NK+1   
               IF ( NK .GT. KTE ) THEN
                 WRITE(message,*)'WOULD GO OFF TOP: Slantwise I,J,DPTHMX,DPMIN',I,J,DPTHMX,DPMIN
                 CALL wrf_message (TRIM(message))
                 EXIT
               ENDIF
               DPTHMX=DPTHMX+DP(NK)
               NLAYRS=NLAYRS+1
               IF(DPTHMX.GT.DPMIN)THEN
                 EXIT 
               ENDIF
             END DO    
           ENDIF
           IF(DPTHMX.LT.DPMIN)THEN 
             RETURN
           ENDIF
           KPBL=LC+NLAYRS-1   
!
! ************************************************************************
!...for computational simplicity without much loss in accuracy,
!...mix temperature instead of theta for evaluating convective
!...initiation (triggering) potential...
!
           TMIX=0.
           QMIX=0.
           ZMIX=0.
           PMIX=0.
           UgMIX=0.
!
!...FIND THE THERMODYNAMIC CHARACTERISTICS OF THE LAYER BY
!...MASS-WEIGHTING THE CHARACTERISTICS OF THE INDIVIDUAL MODEL
!...LAYERS...
!
!cdir novector
           DO NK=LC,KPBL
             TMIX=TMIX+DP(NK)*T0(NK)
             QMIX=QMIX+DP(NK)*Q0(NK)
             ZMIX=ZMIX+DP(NK)*Z0(NK)
             PMIX=PMIX+DP(NK)*P0(NK)
             UgMIX=UgMIX+DP(NK)*Ug0(NK)
           ENDDO   
          TMIX=TMIX/DPTHMX
          QMIX=QMIX/DPTHMX
          ZMIX=ZMIX/DPTHMX
          PMIX=PMIX/DPTHMX
          UgMIX=UgMIX/DPTHMX
          EMIX=QMIX*PMIX/(0.622+QMIX)
!
!...FIND THE TEMPERATURE OF THE MIXTURE AT ITS LCL...
!
!        TLOG=ALOG(EMIX/ALIQ)
! ...calculate dewpoint using lookup table...
!
          astrt=1.e-3
          ainc=0.075
          a1=emix/aliq
          tp=(a1-astrt)/ainc
          indlu=int(tp)+1
          value=(indlu-1)*ainc+astrt
          aintrp=(a1-value)/ainc
          tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
          TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
          TLCL=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(TMIX-T00))*(TMIX-TDPT)
          TLCL=AMIN1(TLCL,TMIX)
          TVLCL=TLCL*(1.+0.608*QMIX)
          ZLCL = ZMIX+(TLCL-TMIX)/GDRY

       DO NK = LC, KL
         KLCL = NK
         IF ( ZLCL.LE.Z0(NK) )  EXIT
       END DO

          K=KLCL-1
          KLCLn1=K
          ! KLFC: the K level for the level of free convection
          ! KLNB: the K level for the level of neutral buoyancy
          KLFCsw=1 !before calculating SCAPE, we initialize KLFCsw=1 
          KLNBsw=1 !before calculating SCAPE, we initialize KLNBsw=1 
          KLFC=1 !before calculating CAPE, we initialize KLFC=1 
          KLNB=1 !before calculating CAPE, we initialize KLNB=1 

! calculate DLP using Z instead of log(P)
          DLP=(ZLCL-Z0(K))/(Z0(KLCL)-Z0(K))
          PLCL=P0(K)+(P0(KLCL)-P0(K))*DLP
!     
!...ESTIMATE ENVIRONMENTAL TEMPERATURE AND MIXING RATIO AT THE LCL...
 
          TENV=T0(K)+(T0(KLCL)-T0(K))*DLP
          QENV=Q0(K)+(Q0(KLCL)-Q0(K))*DLP
          TVEN=TENV*(1.+0.608*QENV)
          vorLCL=avvor
          UgLCL=Ug0(K)+(Ug0(KLCL)-Ug0(K))*DLP
          VgLCL=Vg0(K)+(Vg0(KLCL)-Vg0(K))*DLP
!     
! BELOW is the method used in the KF scheme to identify "unstable" 
! (buoyant) low-level air parcels. Keep for now. May be useful to 
! check if it works better than our trigger function
! --------------------------------------------------->     
!...CHECK TO SEE IF CLOUD IS BUOYANT USING FRITSCH-CHAPPELL TRIGGER
!...FUNCTION DESCRIBED IN KAIN AND FRITSCH (1992)...W0 IS AN
!...APROXIMATE VALUE FOR THE RUNNING-MEAN GRID-SCALE VERTICAL
!...VELOCITY, WHICH GIVES SMOOTHER FIELDS OF CONVECTIVE INITIATION
!...THAN THE INSTANTANEOUS VALUE...FORMULA RELATING TEMPERATURE
!...PERTURBATION TO VERTICAL VELOCITY HAS BEEN USED WITH THE MOST
!...SUCCESS AT GRID LENGTHS NEAR 25 km.  FOR DIFFERENT GRID-LENGTHS,
!...ADJUST VERTICAL VELOCITY TO EQUIVALENT VALUE FOR 25 KM GRID
!...LENGTH, ASSUMING LINEAR DEPENDENCE OF W ON GRID LENGTH...
!     
          IF(ZLCL.LT.2.E3)THEN        ! Kain (2004) Eq. 2
            WKLCL=0.02*ZLCL/2.E3
          ELSE
            WKLCL=0.02                ! units of m/s
          ENDIF
          WKL=(W0AVG1D(K)+(W0AVG1D(KLCL)-W0AVG1D(K))*DLP)*DX/25.E3-WKLCL
          IF(WKL.LT.0.0001)THEN
            DTLCL=0.
          ELSE 
            DTLCL=4.64*WKL**0.33      ! Kain (2004) Eq. 1
          ENDIF
!--------
! Higher RH provides an additional delta T (DTRH).

        dtrh = 0.

!...for ETA model, give parcel an extra temperature perturbation based
!...the threshold RH for condensation (U00)...
! as described in Narita and Ohmori (2007, 12th Mesoscale Conf.!
!...for now, just assume U00=0.75...
!...!!!!!! for MM5, SET DTRH = 0. !!!!!!!!
          U00 = 0.75
          IF(U00.lt.1.)THEN
            QSLCL=QES(K)+(QES(KLCL)-QES(K))*DLP
            RHLCL = QENV/QSLCL
            DQSSDT = QMIX*(CLIQ-BLIQ*DLIQ)/((TLCL-DLIQ)*(TLCL-DLIQ))
            IF(RHLCL.ge.0.75 .and. RHLCL.le.0.95)then
              DTRH = 0.25*(RHLCL-0.75)*QMIX/DQSSDT
            ELSEIF(RHLCL.GT.0.95)THEN
              DTRH = (1./RHLCL-1.)*QMIX/DQSSDT
            ELSE
              DTRH = 0.
            ENDIF
          ENDIF   
!<---------------------------------------------

! The current trigger fuction for CSI release: 
! postivie vertical and horizontal initial triggering forces
! AND high relative humidity
! *testing results suggest a non-zero threshold is better for
!  vertical direction to avoid false trigger due to numerical 
!  noises. We set GWLCL > 0.006 m/s for DX=40km.
!  (note that GWLCL is scaled by DX) 
       
      GWLCL=(W0AVG1D(K)+(W0AVG1D(KLCL)-W0AVG1D(K))*DLP)*DX/40.E3
           
trigger1:  IF (GWLCL.GT.0.006 .AND. VgLCL.GT.0. .AND. RHLCL.GT.0.85) THEN   
!             ! Parcel/tube is buoyant
!     
!...COMPUTE EQUIVALENT POTENTIAL TEMPERATURE
!...(THETEU) AND VERTICAL VELOCITY OF THE RISING PARCEL AT THE LCL...
!     
            CALL ENVIRTHT(PMIX,TMIX,QMIX,THETEU(K),ALIQ,BLIQ,CLIQ,DLIQ)
!
!           because the air parcel is saturated, 
!           saturation theta e (thetaes) is conserved for above layers.
             thetaesLCL=THETEU(K)
!
! WLCL is diagnosed from the KF scheme. Not really used in the SC scheme
! We keep this part just in case it may be useful for future testing. 
! WLCL: the (empirically?) "modified" initial vertical velocity of the parcel at its LCL
!
            DTTOT = DTLCL+DTRH
            IF(DTTOT.GT.1.E-4)THEN
              GDT=2.*G*DTTOT*500./TVEN     ! Kain (2004) Eq. 3  (sort of)
              WLCL=1.+0.5*SQRT(GDT)
              WLCL = AMIN1(WLCL,3.)
            ELSE
              WLCL=1.
            ENDIF
            WTW=WLCL*WLCL
!
            TVLCL=TLCL*(1.+0.608*QMIX)
            RHOLCL=PLCL/(R*TVLCL)
!        
            LCL=KLCL
!     
            WU(K)=WLCL
!     
            TU(K)=TLCL
            TVU(K)=TVLCL
            TVQU(K)=TVLCL
            QU(K)=QMIX
            QLIQ(K)=0.
            QICE(K)=0.
!     
!     
!...ENTER THE LOOP FOR UPDRAFT CALCULATIONS...(for the lifted air tube) 
!...CALCULATE UPDRAFT TEMP,MIXING RATIO AT EACH MODEL LEVEL...
!     
!     
            DILBE = 0.
            CAPE = 0. 
            SCAPE1D = 0.

updraft:    DO NK=K,KTrop
              SDILBE(NK)=0.
              NK1=NK+1
              TU(NK1)=T0(NK1)
              THETEU(NK1)=THETEU(NK)
              QU(NK1)=QU(NK)
              QLIQ(NK1)=QLIQ(NK)
              QICE(NK1)=QICE(NK)
             
              call sltpmix2(p0(nk1),theteu(nk1),tu(nk1),qu(nk1),qliq(nk1),        &
                     qice(nk1),qnewlq,qnewic,XLV1,XLV0,ALIQ,BLIQ,CLIQ,DLIQ)
             
              TVQU(NK1)=TU(NK1)*(1.+0.608*QU(NK1)-QLIQ(NK1)-QICE(NK1))
!                        
              call LIFTPARCEL(p0(nk1),p0(nk),TLCL,PLCL,QMIX,TMIX,TPA(nk1),TPA(nk),TVQU(nk1))
!
              TVU(NK1)=TU(NK1)*(1.+0.608*QU(NK1))
!
!
              IF(NK.EQ.K)THEN
                DZZ=Z0(NK1)-ZLCL
                DPP=P0(NK1)-PLCL
              ELSE
                DZZ=DZA(NK)
                DPP=DPA(NK)
              ENDIF
!================================================================
!
!...IF VERT VELOCITY IS LESS THAN ZERO, EXIT THE UPDRAFT LOOP AND,
!...IF CLOUD IS TALL ENOUGH, FINALIZE UPDRAFT CALCULATIONS...
!
              IF(WTW.LT.1.E-3)THEN
                EXIT
              ELSE
                WU(NK1)=SQRT(WTW)
              ENDIF

! Start calculating CAPE and SCAPE: 
! 
! CAPE = integral of (G*{(thetav_p-thetav_e)/thetav_e})*dz 
!      = (G*{(thetav_p/thetav_e)-1})*dz
! Some formula also expresses it using Tv (BOTHE THETAv and Tv can be used to
! approximate the buoyancy in the vertical momentum equation (with differnt
! estimation errors).  
! Here, TVQU is the Tv for the lifted parcel; TV0 for the environmental profiles
!
! For SCAPE, the environmental profiles are the ones along M surfaces (TV0M)
! Estimated using the first-order Tylar expansion and thermal wind approx.
! DPP=DPA(K-1)=P0(K)-P0(K-1)
!
              IF(NK.EQ.K)THEN
                TV0M   = 0.5*(TVEN+TV0(NK1))-0.5*f/R/avvor*0.5*(PLCL+P0(NK1))* &
                           ((Ug0(NK1)-UgMIX)**2.-(UgLCL-UgMIX)**2.)/(-1*(P0(NK1)-PLCL))
                DILBE  = R*((TVLCL+TVQU(NK1))-(TVEN+TV0(NK1)))*-1.*DPP/(PLCL+P0(NK1))
                SDILBE(NK) = R*(0.5*(TVLCL+TVQU(NK1))-TV0M)*-1.*DPP/(0.5*(PLCL+P0(NK1)))
              ELSE
                TV0M   = 0.5*(TV0(NK)+TV0(NK1))-0.5*f/R/avvor*0.5*(P0(NK)+P0(NK1))* &
                           ((Ug0(NK1)-UgMIX)**2.-(Ug0(NK)-UgMIX)**2.)/(-1*(P0(NK1)-P0(NK)))
                DILBE  = R*((TVQU(NK)+TVQU(NK1))-(TV0(NK)+TV0(NK1)))*-1.*DPP/(P0(NK)+P0(NK1)) 
                SDILBE(NK) = R*(0.5*(TVQU(NK)+TVQU(NK1))-TV0M)*-1.*DPP/(0.5*(P0(NK)+P0(NK1))) 
              ENDIF
!
              IF(DILBE.GT.0.)THEN 

                CAPE = CAPE+DILBE
                IF (KLFC.eq.1) KLFC=NK

              ENDIF
!
              IF(SDILBE(NK).GT.0.)THEN 
                SCAPE1D = SCAPE1D+SDILBE(NK)
                IF (KLFCsw.eq.1) THEN 
                KLFCsw=NK  !if a KLFCsw is definded (not equal to zero) for this usl, no need to
                           !replace it (would be inaccurate as there might be layers with negative SDILBE
                           !in the mid levels) 
                ENDIF
              ENDIF
                check_cape(NK)= CAPE
                check_scape(NK)= SCAPE1D
!
            END DO updraft
!
! Some adjustment for SCAPE: IF the during the calculation, there's a thick
! layers (>300 hPa) in which SDILBE become negative, then even if on top of it, 
! there's another layer of positive SDILBE, SCAPE only represents until the
! bottom of where SDILBE becomes negative ------------------>

        IF (SCAPE1D.gt.0) THEN
          DO NK=KLFCsw+1, KTrop
             IF (SDILBE(NK).lt.0) THEN
                 Nlayer_klow=NK
                 Nlayer_kupp=NK !initialization
                 do m = NK+1, KTrop
                    if (SDILBE(m).lt.0) then
                    else
                    Nlayer_kupp=m
                    exit
                    endif
                 enddo
                 P0diff=P0(Nlayer_klow)-P0(Nlayer_kupp) 
                 IF (P0diff.gt.2.e4) THEN !200 hPa
                    KLNBsw=Nlayer_klow ! = Nlayer_klow-1+1
                    SCAPE1D = 0.
                    do m = KLFCsw, KLNBsw
                       IF (SDILBE(m).GT.0) THEN
                       SCAPE1D = SCAPE1D + SDILBE(m)
                       ENDIF
               !        print*,'I,J=',I,J,'******readjustment*******'
                       check_scape(m)= SCAPE1D
                    enddo
                 ENDIF
             ENDIF
          ENDDO
        ENDIF
! <--------------------------------------------
! Find characteristic levels:
 
           DO NK=K,KTrop
             IF (SCAPE1D.gt.0. .and. SCAPE1D.eq. check_scape(NK)) THEN
               KLNBsw=NK+1
             EXIT
             ENDIF
           ENDDO

           DO NK=K,KTrop
             IF (CAPE.gt.0. .and. CAPE.eq. check_cape(NK)) THEN
               KLNB=NK+1
             EXIT
             ENDIF
           ENDDO
!
           IF (SCAPE1D.gt.0.) THEN
           maxSCAPE1D=0.
           kmidSCAPE1D=KLFCsw+int((KLNBsw-KLFCsw)/3.*2.)
           DO m = KLFCsw, kmidSCAPE1D
             IF (SDILBE(m).gt.maxSCAPE1D) THEN
                 maxSCAPE1D=SDILBE(m)
                 kmaxSCAPE1D=m
             ENDIF
           ENDDO
           ENDIF
 
           IF (SCAPE1D.gt.0. )THEN
             IF (KLFCsw .eq. K) THEN 
               CLDHGTsw(LC)=Z0(KLNBsw)-ZLCL
             ELSE
               CLDHGTsw(LC)=Z0(KLNBsw)-Z0(KLFCsw)
             ENDIF
           ENDIF
!
           IF (CAPE.gt.0. )THEN
             IF (KLFC .eq. K) THEN 
               CLDHGT(LC)=Z0(KLNB)-ZLCL
             ELSE
               CLDHGT(LC)=Z0(KLNB)-Z0(KLFC)
             ENDIF
           ENDIF
!
! For now, we set a fixed threshold for the minimum cloud height:

           CHMIN = 2.E3
!     
!...If cloud top height is less than the specified minimum for deep 
!...convection, go back up to check next level...
!     

            SCAPEres=SCAPE1D-CAPE
            if (SCAPE1D.gt.0.) then
            SCAPEfrac=SCAPEres/SCAPE1D 
            else 
            SCAPEfrac=0.
            endif 

!        
! For slantwiase convection, the trigger fuction2 depends on SCAPEres AND
! SCAPEfrac. If the criterion are not met, the current usl will be
! elimiated and another usl one layer above will be checked. 
!
            IF(   (CLDHGTsw(LC).GT.CHMIN) .and.  &
                  ((SCAPEres.GT.100 .and. SCAPEfrac.GT.0.1).or. &
                   (SCAPEres.GT.50  .and. SCAPEfrac.GT.0.5))       ) THEN  
!
              sw_trigger=.TRUE.

              EXIT usl
!
            ENDIF
!
          ENDIF trigger1
        END DO usl

! Other important estimations 
!
   Lscale= abs((Ug0(KLNBsw)-UgMIX)/avvor)
   ang_M=atan(Z0(KLNBsw)/Lscale)
   wmax=sqrt(2*SCAPE1D)*cos(0.5*PIE-ang_M)
   overshoot_z=wmax*3600*3
   DO NK=KLNBsw+1, KTrop
      IF ( (Z0(NK)-Z0(KLNBsw)).le.overshoot_z .and. &
           (Z0(NK+1)-Z0(KLNBsw)).gt.overshoot_z      ) THEN
         overshoot_k=NK
      ELSE
         overshoot_k=KLNBsw
      ENDIF
   ENDDO
   overshoot_k=min(overshoot_k,KTrop)
!
! Starting calculating the reference profile for Ug:
!
    delTvt = 0.
    delp = 0.
    rhs = 0.
    ! Here, we assume RUg0 equals to RUgMIX
    ! RUg1 indicates the first guess of referenced Ug 
    RUg0 = UgMIX

!---- TARGET ON THE LAYER KLFCsw - KLNDsw FIRST ----
! UgMIX is taken for K= LC-KPBL.

    DO NK = KLFCsw, KLNBsw 
           readjust(NK) = .FALSE.
           NK1  = NK+1
           vor1 = avvor
           IF (NK.eq.K) THEN 
             delp   = (P0(NK1)-PLCL)
             delTvt = (TVLCL+TVQU(NK1)-TV0(NK)-TV0(NK1))
             rhs    = R/f*delTvt/(PLCL+P0(NK1))*2*vor1*delp &
                    + (UgLCL-RUg0)**2.   !here assume UgLCL stays fixed
           ELSE
             delp   = (P0(NK1)-P0(NK))
             delTvt = (TVQU(NK)+TVQU(NK1)-TV0(NK)-TV0(NK1))
             rhs    = R/f*delTvt/(P0(NK)+P0(NK1))*2*vor1*delp &
                    + (RUg1(NK)-RUg0)**2.   !RUg1(NK) was initialized with Ug0(NK)
           ENDIF
!
           IF (rhs .GE. 0 ) THEN
               IF (Ug0(NK1).ge.RUg0) THEN
                 RUg1(NK1) = RUg0+sqrt(rhs)
               ELSE
                 RUg1(NK1) = RUg0-sqrt(rhs)
               ENDIF
           ELSE
               RUg1(NK1)=RUg0
               readjust(NK) = .TRUE.
           ENDIF
!    
     ENDDO

! readjustment (to avoid miscalculation when rhs is negative):   

     DO NK = KLNBsw+1, KLFCsw, -1

     IF ( readjust(NK) ) THEN

        DO m = NK+1, KLFCsw+1, -1
           m1=m-1
           vor1= avvor

            IF (m1.eq.K) THEN 
             delTvt = (TVLCL+TVQU(m)-TV0(m1)-TV0(m))
             delp   = (P0(m)-PLCL)
             rhs    = -1.*R/f*delTvt/(P0(m)+PLCL)*2*vor1*delp &
                    + (RUg1(m)-RUg0)**2.   !RUg1(NK) was initialized with Ug0(NK)
            ELSE
             delTvt = (TVQU(m1)+TVQU(m)-TV0(m1)-TV0(m))
             delp   = (P0(m)-P0(m1))
             rhs    = -1.*R/f*delTvt/(P0(m)+P0(m1))*2*vor1*delp &
                    + (RUg1(m)-RUg0)**2.   !RUg1(NK) was initialized with Ug0(NK)
            ENDIF
          
            IF (m1.eq.K) THEN 
              IF (rhs .GE. 0 ) THEN
                IF (UgLCL.ge.RUg0) THEN
                  RUgLCL=RUg0+sqrt(rhs)
                ELSE
                  RUgLCL=RUg0-sqrt(rhs)
                ENDIF
                readjust(m1)=.FALSE.
              ELSE
                RUgLCL=RUg0
                readjust(m1)=.TRUE.
              ENDIF
              RUg1(m1)=(RUgLCL-RUg1(m)*DLP)/(1.-DLP)             
            ELSE
              IF (rhs .GE. 0 ) THEN
                IF (Ug0(m1).ge.RUg0) THEN
                  RUg1(m1)=RUg0+sqrt(rhs)
                ELSE
                  RUg1(m1)=RUg0-sqrt(rhs)
                ENDIF
                readjust(m1)=.FALSE.
              ELSE
                RUg1(m1)=RUg0
                readjust(m1)=.TRUE.
              ENDIF
            ENDIF
        ENDDO
        EXIT

     ENDIF !readjust

     ENDDO
!!------------------------------------------------------ 
! The above first guess RUg0 is obtained by eliminating 
! SCAPE ocmpletly. But the final Reference profile RUg 
! must conserve momentum! A simple way is to adjust the 
! profile by substracting/adding a constant value for 
! every layers...... 
!
    Ugres = 0. 
    sumdelp = 0.
! 
    DO NK = 1, overshoot_k
          Ugres   = Ugres+(RUg1(NK)-Ug0(NK))*((P0(NK+1)-P0(NK))) 
          sumdelp = sumdelp+(P0(NK+1)-P0(NK))
    ENDDO
!

    DO NK= 1, overshoot_k
          RUg(NK)=RUg1(NK)-Ugres/sumdelp 
    ENDDO
!----------------------
 
! Here, we're double checking the 
! SCAPE using the final Reference profile RUg:

    TV0M2=0.
    SCAPE1D_ad=0.
    DILBE2=0.
    UgMIX2=UgMIX-Ugres/sumdelp
    UgLCL=RUg(K)+(RUg(KLCL)-RUg(K))*DLP
    DO NK=K, KTrop
    NK1=NK+1
              IF(NK.EQ.K)THEN
                TV0M2   = 0.5*(TVEN+TV0(NK1))-0.5*f/R/avvor*0.5*(PLCL+P0(NK1))* &
                           ((RUg(NK1)-UgMIX2)**2.-(UgLCL-UgMIX2)**2.)/(-1*(P0(NK1)-PLCL))
                DILBE2  = R*((TVLCL+TVQU(NK1))-(TVEN+TV0(NK1)))*-1.*DPP/(PLCL+P0(NK1))
                SDILBE2(NK) = R*(0.5*(TVLCL+TVQU(NK1))-TV0M2)*-1.*DPP/(0.5*(PLCL+P0(NK1)))
              ELSE
                TV0M2   = 0.5*(TV0(NK)+TV0(NK1))-0.5*f/R/avvor*0.5*(P0(NK)+P0(NK1))* &
                           ((RUg(NK1)-UgMIX2)**2.-(RUg(NK)-UgMIX2)**2.)/(-1*(P0(NK1)-P0(NK)))
                DILBE2  = R*((TVQU(NK)+TVQU(NK1))-(TV0(NK)+TV0(NK1)))*-1.*DPP/(P0(NK)+P0(NK1)) 
                SDILBE2(NK) = R*(0.5*(TVQU(NK)+TVQU(NK1))-TV0M2)*-1.*DPP/(0.5*(P0(NK)+P0(NK1))) 
              ENDIF
              IF(SDILBE2(NK).GT.0.)THEN 
                SCAPE1D_ad = SCAPE1D_ad+SDILBE2(NK)
              ENDIF
    ENDDO

          if (SCAPE1D.gt.0.) adjust_frac=(SCAPE1D-SCAPE1D_ad)/SCAPE1D
!     
!====================================================================
!     
! Prescribe the adjustment time scale for slantwise convection (TIMEC)
!     
! Ma(2000) suggests that the time scale SBI ranges from 0.94~2.2(f*vor)**(-0.5)
! Reuter and Yau's (1990) observational studies suggest 3 hours
!     
        TIMEC=3600*3.          ! unit:sec ; 1 hours  >= TIMEC <= 5 hours
!     
        DO K=1,KX

          DUgDT(K)=(RUg(K)-Ug0(K))/TIMEC

        ENDDO
!
!-----------------------------------------------------------------------
  END SUBROUTINE  Slantwise1D

! ***********************************************************************

  SUBROUTINE DTDQV_adjust(DTHDT_cond,DQDT_cond,RAINincrease,SCAPE1D,GWLCL,  &
                          SCAPEfrac,ang_M,ZLCL,K,KLFCsw,kmaxSCAPE1D,KLNBsw, &
                          RHavg,QVavg,T0,P0,Z0,SLCUDT,DT,ALIQ,BLIQ,CLIQ,    &
                          DLIQ,XLV,CP,R_V,R,G,EP2,kts,kte,ims,ime,jms,jme,  &
                          I,J, KRHsat,skipflag, DUgDT_RH, RUg1,Ug0, TIMEC,  &
                          overshoot_k,Z1)
  IMPLICIT NONE

  INTEGER, INTENT(IN  )  ::kts, kte, ims, ime, jms, jme,I,J  

  REAL, DIMENSION( kts:kte ), INTENT(  OUT) ::  DTHDT_cond, DQDT_cond, DUgDT_RH 

  REAL, INTENT(  OUT)::       RAINincrease,Z1
   
  REAL, DIMENSION( kts:kte ), INTENT(IN  ) ::  RHavg,T0,P0,QVavg,Z0,RUg1, Ug0

  REAL, INTENT(IN  )  :: SCAPE1D, GWLCL, SCAPEfrac, ang_M, & 
                         ZLCL,DT, SLCUDT, TIMEC 
  REAL, INTENT(IN  )  :: ALIQ, BLIQ, CLIQ, DLIQ, XLV,CP,R_V,R,G,EP2

  INTEGER, INTENT(IN   ) :: K, KLFCsw,kmaxSCAPE1D, KLNBsw, overshoot_k
  INTEGER, INTENT(  OUT) :: KRHsat

  ! local
  REAL, DIMENSION( kts:kte ) ::  weight,qvweight,DTDT_cond

  INTEGER :: m, NK, KX, NKn1, condmaxK, kmax,drybot,drytop,KRHsat_low

  REAL    :: wslwmax,w0Z,Z2,subdp2,wslw,pcpslw,weightsum,pp,    &
             ES,Q1,Q2, cond_rate, weightmax,weightsum2,subdp,   &
             lowcdp,DSHDT_cond,sumT,sumQ,DTDT_cond0,asubdqdt,   &
             subdqdt,QVlostpre,sumT2,Z3,ratio,RUg_RH,Ugres,     & 
             sumdelp, drydepth,Zqv,qvweighttotal,qvweightsum

  REAL    :: PIE = 3.141592654, cfrac=0.5 !cfrac: cloud fraction (tunable)

  REAL, PARAMETER :: kapa = 0.286

  LOGICAL, INTENT(  OUT) :: skipflag
!-------------
   KX=kte
   skipflag=.FALSE.

! Now let's deal with condensate, precipitation, latent heat and moisture loss
! Notes: Q1, Q2 are thermodynamic function from RY (1989, Eqns. 7.23 and 7.24)
! A big assumption here is that slantwise convection is partially resolvable
! due to its mesoscale feature! 

    DO NK = 1, KX
         DQDT_cond(NK) = 0.
         DTDT_cond(NK) = 0.
         weight(NK) = 0.
         DUgDT_RH(NK) = 0.
    ENDDO

! Utilize the grid-resolved RH>=90% to identify the extended height of 
! the (developing) slantwise convection
!------------------------------>
    KRHsat_low=KX
    DO NK=1, KLFCsw
       IF (RHavg(NK).ge.0.9) THEN
          KRHsat_low=NK
          goto 400
       ENDIF
    ENDDO

400 IF (KRHsat_low.ne.KX) THEN
    w0Z=Z0(KRHsat_low)
    ELSE
    w0Z=ZLCL
    ENDIF

    KRHsat= 2
    DO NK = KX, 1, -1
       IF (RHavg(NK).ge.0.9) THEN
           KRHsat=NK
           goto 500
       ENDIF 
    ENDDO
!
    
500 drydepth=0.
    NK=KRHsat-1 
    DO WHILE (NK.ge.2) 
       drybot=NK
       drytop=NK
       IF (RHavg(NK).lt.0.9) THEN
           DO m = drytop-1, 1, -1
             IF (RHavg(m).lt.0.9 ) THEN
                 drybot=m 
                 drydepth=Z0(drytop)-Z0(drybot)       
                 IF (drydepth.gt.1.e3 ) THEN
                 KRHsat=drybot
                 ENDIF 
             ELSE
                 goto 100
             ENDIF
           ENDDO
       ENDIF
100    NK=NK-1 
    ENDDO
!<------------------------------

! Parameterizing the vertical motion associated with the slantwise
! convection (scaled by ratio as SCAPE is not going to be released 
! completley) and a parabolic profile for latent heat distribution
!
    ratio=(Z0(KRHsat)-Z0(KLFCsw))/(Z0(KLNBsw)-Z0(KLFCsw))
    ratio=AMIN1(ratio,1.)
    ratio=AMAX1(ratio,0.)
    kmax=min0(KRHsat,overshoot_k)
!
    if (ratio.lt.0.1 .or. Z0(KRHsat).le.Z0(KLFCsw)) then 
    skipflag=.TRUE.
    return
    endif
    IF (ratio.ge.0.9) THEN
    Z1=(Z0(KLFCsw)+0.5*(Z0(KLNBsw)-Z0(KLFCsw)))
    Z2=(Z0(KLFCsw)+0.5*(Z0(KLNBsw)-Z0(KLFCsw)))
    Z3=Z0(kmax)
    ELSE
    Z1=Z0(KRHsat)*1./2.
    Z2=Z0(KRHsat)*1./2.
    Z3=Z0(KRHsat)
    ENDIF    

    wslwmax = 1./0.22*cfrac*GWLCL*sqrt(2.*SCAPE1D)*ratio*SCAPEfrac*cos(0.5*PIE-ang_M)

    subdp2 = 0.
    wslw   = 0.
    pcpslw = 0.
    weightsum = 0.

    DO NK = 1, min0(kmax,KX-1)

       NKn1=max0(NK-1,1)

       IF (Z0(NK).le.w0Z) THEN
         weight(NK)=0.
       ELSEIF (Z0(NK).gt.w0Z .and. Z0(NK).le.Z1) THEN
         pp=1./(4)*(w0Z-Z1)**2.
         weight(NK)=-1./(4*pp)*(Z0(NK)-Z1)**2.+1.
       ELSEIF ( (Z0(NK).gt.Z1) .and. (Z0(NK).le.Z2)) THEN
         weight(NK)=1.
       ELSEIF ( (Z0(NK).gt.Z2) .and. (Z0(NK).le.Z3) ) THEN
         pp=1./(4)*(Z3-Z2)**2.
         weight(NK)=-1./(4*pp)*(Z0(NK)-Z2)**2.+1.
       ELSE
         weight(NK)=0.
       ENDIF

       weightsum = weightsum+weight(NK)
       !
       if ( Z0(NK).gt.w0Z .and. Z0(NK).le.Z3 .and. RHavg(NK).ge.0.9 ) then
       ES = ALIQ*EXP((BLIQ*T0(NK)-CLIQ)/(T0(NK)-DLIQ)) 
           !saturation vapour pressure with respect to water
       Q1 = G/T0(NK)*(XLV/CP/R_V/T0(NK)-1./R)         !R=r_d !Q1:unit(m^-1)
       Q2 = P0(NK)/EP2/ES+XLV**2./CP/R_V/(T0(NK)**2.) !EP2=epsilon=r_d/r_v
       cond_rate = Q1/Q2*wslwmax !condensation rate
       pcpslw = pcpslw+cond_rate*(P0(NK+1)-P0(NKn1))/(NK+1-NKn1)/G
       endif
       !
    ENDDO
    !
       weightmax  = 0.
       condmaxK   = 1
       weightsum2 = 0.
    !
    DO NK= 1, max0(kmax,KX-1)

       if(weightsum.ne.0) then
       weight(NK)=weight(NK)/weightsum
         IF (weight(NK).gt.weightmax) THEN
          weightmax=weight(NK)
          condmaxK=NK
         ENDIF
         weightsum2 = weightsum2+weight(NK)
       else
         return
       endif
    ENDDO

      subdp  = 0.
      lowcdp = 0.
      DSHDT_cond  = 0.
      qvweightsum = 0.
      Zqv         = Z0(KLFCsw)+500.

      DO NK = 1, max0(kmaxSCAPE1D,KLFCsw)
         NKn1=max0(NK-1,1)
         qvweight(NK)=0.
       !A parabolic profile for qv distribution
        IF (Z0(NK).le.0.5*ZLCL) THEN
          qvweight(NK)=-1./((0.-0.5*ZLCL)**2.)*(Z0(NK)-0.5*ZLCL)**2.+1.
        ELSEIF (Z0(NK).gt.0.5*ZLCL.and. Z0(NK).le.0.5*ZLCL) THEN
          qvweight(NK)=1.
        ELSEIF (Z0(NK).gt.0.5*ZLCL .and. Z0(NK).le.Zqv) THEN
          qvweight(NK)=-1./((Zqv-0.5*ZLCL)**2.)*(Z0(NK)-0.5*ZLCL)**2.+1.
        ENDIF
       
       qvweightsum=qvweightsum+qvweight(NK)
       
       ENDDO
    !  
      sumT = 0.
      sumQ = 0.
      DTDT_cond0    = 0.
      qvweighttotal = 0.

      DO NK = 1, max0(kmaxSCAPE1D,KLFCsw)
         NKn1=max0(NK-1,1)
         if (qvweightsum.ne.0) then
           qvweight(NK)=qvweight(NK)/qvweightsum  
         endif 
           qvweighttotal=qvweighttotal+qvweight(NK)
           asubdqdt=pcpslw/((P0(NK+1)-P0(NKn1))/(NK+1-NKn1)/G)*qvweight(NK) !dqv/dt
           QVlostpre= QVavg(NK)-asubdqdt*SLCUDT 
           ! QVlostpre: the Qv after including the Qv tendency from the SC scheme
           ! Check if this value will be come too small (avoid resulting in
           ! negative values!). The threshold we use here 0.2*(original Qv) is 
           ! taken from Ma (2000).
           IF (QVlostpre .LT. 0.2*QVavg(NK)) THEN
              subdqdt=0.8*QVavg(NK)/SLCUDT
           ELSE
              subdqdt=asubdqdt
           ENDIF
           DQDT_cond(NK) = -1.*subdqdt
           DTDT_cond0 = XLV*subdqdt/CP
         !
         sumT= sumT+DTDT_cond0*(P0(NK+1)-P0(NKn1))/(NK+1-NKn1)/G
         sumQ= sumQ+DQDT_cond(NK)*(P0(NK+1)-P0(NKn1))/(NK+1-NKn1)/G
      
      ENDDO

         sumT2=0.

      DO NK= 1, min0(kmax,KX-1)
         NKn1=max0(NK-1,1)
         ! Assumes a parabolic profile! 
         DTDT_cond(NK)=sumT*weight(NK)/((P0(NK+1)-P0(NKn1))/(NK+1-NKn1)/G)
         sumT2= sumT2+DTDT_cond(NK)*(P0(NK+1)-P0(NKn1))/(NK+1-NKn1)/G
         !
      ENDDO
!
! Convert the temperature tendency to theta tendency (assume the pressure 
! remains the same...)
!
      DO NK=1,KX
         DTHDT_cond(NK)=DTDT_cond(NK)*(1.E5/P0(NK))**kapa
      ENDDO
!
!*****************************************************************
      RAINincrease=DT*sumQ  !unit: kg/(m^2)= mm (one kg of water spread 
                            !over m^2 area has 1-mm thickness)
!*****************************************************************
! The momentum tendency profile derived here takes the RUg0 from the 
! 1D slantwise model as well. The momentum conservation is ensured 
! for the column from the lowest level to an estimated extended height 
! of the slantwise convection (intead of the dynamically-defined 
! overhoot_k/LSNB).
! 
! Reason: 
! the slantwise convecetion has a longer developing time 
! than vertical conveciton, and the convection does not develop and 
! reach to LSNB instantly. It doesn't make sense to apply the Ug
! tendency from the near surface to LSNB at every time step when the
! SC scheme trigger function is satisfied.
!
    Ugres   = 0.
    sumdelp = 0.
!
    IF (ratio.ge.0.9)THEN
      kmax=overshoot_k
    ELSE
      kmax=KRHsat
    ENDIF
    
    DO NK = 1, kmax
       Ugres   = Ugres+(RUg1(NK)-Ug0(NK))*((P0(NK+1)-P0(NK)))
       sumdelp = sumdelp+(P0(NK+1)-P0(NK))
    ENDDO

    DO NK= 1, KX
       RUg_RH=0.

       IF(NK.le.kmax)THEN
          RUg_RH=RUg1(NK)-Ugres/sumdelp
          DUgDT_RH(NK)=(RUg_RH-Ug0(NK))/TIMEC
       ENDIF

    ENDDO
       
   return 

  END SUBROUTINE DTDQV_adjust

! ***********************************************************************

  SUBROUTINE LIFTPARCEL(P,P1,Tc,Pc,QMIX,TMIX,TPA,TPA1,TVPA)

  IMPLICIT NONE
  REAL,         INTENT(IN   )   :: P,P1,Tc,Pc,QMIX,TMIX
  REAL,         INTENT(OUT  )   :: TVPA
  REAL,         INTENT(INOUT)   :: TPA, TPA1
  REAL    ::    gama1, gama2, gama3, gama4, gamasp, tpa_c, xx, ws0
  real, parameter :: kapa = 0.286
  real, parameter :: Rd = 287.
  real, parameter :: g = 9.8
  real, parameter :: epsl = 0.622
!-----------------------------------------------------------------------

!P=P at level k, P1=P at level k-1 (lifting from the sfc)
!TPA=TPA at level k, TPA1=TPA at level k-1 (lifting from the sfc)

   IF (P.ge.Pc) THEN

!!below LCL: dry adiabatic (thetapa: mixed potential temperature)
!!below LCL: dry adiabatic (tpa: mixed temperature)
          !TPA=thetpa*(P/p00)**kapa
          TPA=TMIX
   ELSE
          if (P1.ge.Pc) then
!!the surface .or. this level is just below LCL
          call gamas(Tc,Pc,gama1)
          gamasp=gama1/(-g*Pc/(Rd*Tc))
          TPA=Tc+gamasp*(P-Pc)
          call gamas(TPA,P,gama2)
          gamasp=gama2/(-g*P/(Rd*TPA))+gamasp
          TPA=Tc+0.5*gamasp*(P-Pc)
          else
          call gamas(TPA1,P1,gama3)
          gamasp=gama3/(-g*P1/(Rd*TPA1))
          TPA=TPA1+gamasp*(P-P1)
          call gamas(TPA,P,gama4)
          gamasp=gama4/(-g*P/(Rd*TPA))+gamasp
          TPA=TPA1+0.5*gamasp*(P-P1)
          endif
    ENDIF

!!Start calculating virtual temp.-----------
!! 
!!For lifted air parcel:!!
    IF (P.le.Pc) then    !above LCL: saturated

!!ws (below:using the formula shown in the slantwise scheme from Ma 2000)

       tpa_c=TPA-273.15
       xx=6.1078*exp((17.2693882*tpa_c)/(tpa_c+237.3))
       ws0=(621.97*xx/(P/100.-xx))/1000. !--> so that th unit: g g-1

!!tvpa-->using virtual temperature!!
       TVPA=TPA*(1+ws0/epsl)/(1.+ws0) !ws0:saturated mixing ratio
    ELSE

       TVPA=TPA*(1+QMIX/epsl)/(1.+QMIX)
    ENDIF

  END SUBROUTINE LIFTPARCEL

! ***********************************************************************

  SUBROUTINE gamas(temp, press, gamas_ret)

  implicit none
  real, intent(in   ) :: temp, press
  real, intent(  out) :: gamas_ret
  real, parameter :: hl0=2.5e6 !latent heat at 273.15K
  real, parameter :: cpv=1870. !Specific heat of water vapor at constant p
  real, parameter :: t00=273.15
  real, parameter :: g= 9.8
  real, parameter :: epsl= 0.622
  real, parameter :: Cp= 1005.
  real, parameter :: Rd= 287.
  real c0, hl, arg, es, ws, temp0
!
!  temp0=amax1(temp,40.)
  temp0=temp
  c0=4187.  !Specific heat of liquid water at T>273.15 K (> 0 Celsius degree)
!
  if (temp0.lt.t00) c0=c0+0.08*(t00-temp0)/30.  !if T< 0 Celsius degree
  hl=hl0-(c0-cpv)*(temp0-t00)
!!Below are for calculation the es(T) but I don't know which empirical formula
!is used...
  arg=17.67*(temp0-273.15)/(temp0-273.15+243.5)
  es=6.112*exp(arg)*100.
!---
  ws=epsl*es/(press-es)
  gamas_ret=-g/Cp*(1.+hl*ws/(Rd*temp0))/(1.+hl*hl*epsl*ws/(Cp*Rd*temp0*temp0))

  return
  END SUBROUTINE gamas

! ***********************************************************************

  SUBROUTINE RHcalc(QV,T,Pcps,RH,thetaes)
  implicit none
  real, intent(in   ) :: T, QV, Pcps
  real, intent(  out) :: RH, thetaes
  real, parameter :: epsl= 0.622
  real, parameter :: Cp= 1005.
  real, parameter :: Rd= 287.
  REAL, PARAMETER :: Rv = 461.6
  REAL, PARAMETER :: Cl =  4190
  REAL, PARAMETER :: Lv = 2.5e6
  real ES, WS, pdry
!

  if (T.gt.273.15) then
     !  here we use the same code as used in WRF 
      ES = 1000.*0.6112*exp(17.67*(T-273.15)/(T-29.65))
  else
      ES = 1000.*0.6112*exp(21.8745584*(T-273.15)/(T-7.66))
  endif
      WS = Rd/Rv*ES/(Pcps-ES)
            !RH(I,K,J) = Q0/QES
      RH=QV/WS*(1+WS/epsl)/(1+QV/epsl)
      pdry=Pcps/(1.+QV/epsl)
  if (QV.ne.0. .and. T.ne.0.) then
      thetaes=T*((1.E5/pdry)**(Rd/(Cp+Cl*WS)))* &
              (RH**(-1.*WS*Rv/(Cp+Cl*WS)))*(exp(Lv*WS/(Cp+Cl*WS)/T))
  else
      thetaes=0.
  endif
! 
  END SUBROUTINE RHcalc
! ***********************************************************************
  SUBROUTINE SLTPMIX2(p,thes,tu,qu,qliq,qice,qnewlq,qnewic,XLV1,XLV0, & 
                       ALIQ,BLIQ,CLIQ,DLIQ)
!-----------------------------------------------------------------------
   ! Taken from the KF scheme subroutine 
!-----------------------------------------------------------------------
   IMPLICIT NONE
   REAL,         INTENT(IN   )   :: P,THES,XLV1,XLV0,ALIQ,BLIQ,CLIQ,DLIQ
   REAL,         INTENT(OUT  )   :: QNEWLQ,QNEWIC
   REAL,         INTENT(INOUT)   :: TU,QU,QLIQ,QICE
   REAL    ::    TP,QQ,BTH,TTH,PP,T00,T10,T01,T11,Q00,Q10,Q01,Q11,          &
                 TEMP,QS,QNEW,DQ,QTOT,RLL,CPP
   REAL    ::    ES,PI,THTGS,F0,T0,F1,DT
   INTEGER ::    ITCNT

! we need to find temp and qs first:

      ES=ALIQ*EXP((BLIQ*TU-CLIQ)/(TU-DLIQ))
      QS=0.622*ES/(P-ES)
      PI=(1.E5/P)**(0.2854*(1.-0.28*QS))
      THTGS=TU*PI*EXP((3374.6525/TU-2.5403)*QS*(1.+0.81*QS))
    
      F0=THTGS-THES
      TEMP=TU-0.5*F0
      T0=TU
      ITCNT=0

   90 ES=ALIQ*EXP((BLIQ*TEMP-CLIQ)/(TEMP-DLIQ))
      QS=0.622*ES/(P-ES)
      PI=(1.E5/P)**(0.2854*(1.-0.28*QS))
      THTGS=TEMP*PI*EXP((3374.6525/TEMP-2.5403)*QS*(1.+0.81*QS))

      F1=THTGS-THES
      IF(ABS(F1).LT.0.01)GOTO 50
      ITCNT=ITCNT+1
      IF(ITCNT.GT.10) THEN 
       print*,'ITCNT exceeds 10'
       GOTO 50
      ENDIF
      DT=F1*(TEMP-T0)/(F1-F0)
      T0=TEMP
      F0=F1
      TEMP=TEMP-DT
      GOTO 90
!
!
   50 DQ=QS-QU
      IF(DQ.LE.0.)THEN
        QNEW=QU-QS
        QU=QS
      ELSE
!
!   IF THE PARCEL IS SUBSATURATED, TEMPERATURE AND MIXING RATIO MUST BE
!   ADJUSTED...IF LIQUID WATER IS PRESENT, IT IS ALLOWED TO EVAPORATE
! 
        QNEW=0.
        QTOT=QLIQ+QICE
!
!   IF THERE IS ENOUGH LIQUID OR ICE TO SATURATE THE PARCEL, TEMP STAYS AT ITS
!   WET BULB VALUE, VAPOR MIXING RATIO IS AT SATURATED LEVEL, AND THE MIXING
!   RATIOS OF LIQUID AND ICE ARE ADJUSTED TO MAKE UP THE ORIGINAL SATURATION
!   DEFICIT... OTHERWISE, ANY AVAILABLE LIQ OR ICE VAPORIZES AND APPROPRIATE
!   ADJUSTMENTS TO PARCEL TEMP; VAPOR, LIQUID, AND ICE MIXING RATIOS ARE MADE.
!
!...subsaturated values only occur in calculations involving various mixtures of
!...updraft and environmental air for estimation of entrainment and detrainment.
!...For these purposes, assume that reasonable estimates can be given using 
!...liquid water saturation calculations only - i.e., ignore the effect of the
!...ice phase in this process only...will not affect conservative properties...
!
        IF(QTOT.GE.DQ)THEN
          qliq=qliq-dq*qliq/(qtot+1.e-10)
          qice=qice-dq*qice/(qtot+1.e-10)
          QU=QS
        ELSE
          RLL=XLV0-XLV1*TEMP
          CPP=1004.5*(1.+0.89*QU)
          IF(QTOT.LT.1.E-10)THEN
!
!...IF NO LIQUID WATER OR ICE IS AVAILABLE, TEMPERATURE IS GIVEN BY:
            TEMP=TEMP+RLL*(DQ/(1.+DQ))/CPP
          ELSE
!
!...IF SOME LIQ WATER/ICE IS AVAILABLE, BUT NOT ENOUGH TO ACHIEVE SATURATION,
!   THE TEMPERATURE IS GIVEN BY:
!
            TEMP=TEMP+RLL*((DQ-QTOT)/(1+DQ-QTOT))/CPP
            QU=QU+QTOT
            QTOT=0.
            QLIQ=0.
            QICE=0.
          ENDIF
        ENDIF
      ENDIF
      TU=TEMP
      qnewlq=qnew
      qnewic=0.
!
  END SUBROUTINE SLTPMIX2

! ***********************************************************************

  SUBROUTINE ENVIRTHT(P1,T1,Q1,THT1,ALIQ,BLIQ,CLIQ,DLIQ)                       
!-----------------------------------------------------------------------
   ! Taken from the KF scheme subroutine 
!-----------------------------------------------------------------------
   IMPLICIT NONE
   REAL,         INTENT(IN   )   :: P1,T1,Q1,ALIQ,BLIQ,CLIQ,DLIQ
   REAL,         INTENT(INOUT)   :: THT1
   REAL    ::    EE,TLOG,ASTRT,AINC,A1,TP,VALUE,AINTRP,TDPT,TSAT,THT,      &
                 T00,P00,C1,C2,C3,C4,C5
   INTEGER ::    INDLU
   DATA T00,P00,C1,C2,C3,C4,C5/273.16,1.E5,3374.6525,2.5403,3114.834,   &
           0.278296,1.0723E-3/                                          
!                                                                       
!  CALCULATE ENVIRONMENTAL EQUIVALENT POTENTIAL TEMPERATURE...          
!  (Is this for SATURATION theta e?)--> No, i don't think so. If yes, an input of 
!  Q1(mixing ratio) won't be necessary. 
!  But if the parcel is saturated (above LCL), theta e = theta es
!
!  The saturation Q is a function of T and P only.
!
! NOTE: Calculations for mixed/ice phase no longer used...jsk 8/00
!        For example, KF90 Eq. 10 no longer used
!
      EE=Q1*P1/(0.622+Q1)                                             
      TLOG=ALOG(EE/ALIQ)                                              
! ...calculate LOG term using lookup table...
!
!      astrt=1.e-3
!      ainc=0.075
!      a1=ee/aliq
!      tp=(a1-astrt)/ainc
!      indlu=int(tp)+1
!      value=(indlu-1)*ainc+astrt
!      aintrp=(a1-value)/ainc
!      tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
!
      TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)                               
      TSAT=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(T1-T00))*(T1-TDPT) 
      THT=T1*(P00/P1)**(0.2854*(1.-0.28*Q1))                          
      THT1=THT*EXP((C1/TSAT-C2)*Q1*(1.+0.81*Q1))  !Seems from Bolton (1980)                    
!
  END SUBROUTINE ENVIRTHT                                                              

! ***********************************************************************

  SUBROUTINE slantwise_init(RUSLCUTEN,RVSLCUTEN,                    &
                     RTHSLCUTEN, RQVSLCUTEN,                        &
                     swW0AVG,P_QI,P_QS,                             &
                     SVP1,SVP2,SVP3,SVPT0,                          &
                     P_FIRST_SCALAR,restart,allowed_to_read,        &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte                   )
!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           ::  restart,allowed_to_read
   INTEGER , INTENT(IN)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte
   INTEGER , INTENT(IN)           ::  P_QI,P_QS,P_FIRST_SCALAR

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) ::       &
                                                          RUSLCUTEN, &
                                                          RVSLCUTEN, &
                                                         RTHSLCUTEN, &
                                                         RQVSLCUTEN

   REAL ,   DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: swW0AVG

   INTEGER :: i, j, k, itf, jtf, ktf
   REAL, INTENT(IN)    :: SVP1,SVP2,SVP3,SVPT0

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   IF(.not.restart)THEN

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         RUSLCUTEN(i,k,j)=0.
         RVSLCUTEN(i,k,j)=0.
         RTHSLCUTEN(i,k,j)=0.
         RQVSLCUTEN(i,k,j)=0.
      ENDDO
      ENDDO
      ENDDO


      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         swW0AVG(i,k,j)=0.
      ENDDO
      ENDDO
      ENDDO

   endif
 
   CALL SW_LUTAB(SVP1,SVP2,SVP3,SVPT0)

   END SUBROUTINE slantwise_init

! ***********************************************************************

   SUBROUTINE sw_lutab(SVP1,SVP2,SVP3,SVPT0)
!-----------------------------------------------------------------------
   ! Taken from the KF scheme subroutine 
!-----------------------------------------------------------------------
!
!  This subroutine is a lookup table.
!  Given a series of series of saturation equivalent potential 
!  temperatures, the temperature is calculated.
!
!--------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: KP,IT,ITCNT,I
    REAL :: DTH,TMIN,TOLER,PBOT,DPR,                               &
            TEMP,P,ES,QS,PI,THES,TGUES,THGUES,F0,T1,T0,THGS,F1,DT, &
            ASTRT,AINC,A1,THTGS
    REAL    :: ALIQ,BLIQ,CLIQ,DLIQ
    REAL, INTENT(IN)    :: SVP1,SVP2,SVP3,SVPT0
!
! equivalent potential temperature increment
      data dth/1./
! minimum starting temp 
      !data tmin/150./
      data tmin/100./
! tolerance for accuracy of temperature 
      data toler/0.001/
! top pressure (pascals)
      plutop=5000.0
! bottom pressure (pascals)
      pbot=110000.0

      ALIQ = SVP1*1000.
      BLIQ = SVP2
      CLIQ = SVP2*SVPT0
      DLIQ = SVP3
!
! compute parameters
!
! 1._over_(sat. equiv. theta increment)
      rdthk=1./dth
! pressure increment
!
      DPR=(PBOT-PLUTOP)/REAL(KFNP-1)
!      dpr=(pbot-plutop)/REAL(kfnp-1)
! 1._over_(pressure increment)
      rdpr=1./dpr
! compute the spread of thes
!     thespd=dth*(kfnt-1)
!
! calculate the starting sat. equiv. theta
!
      temp=tmin 
      p=plutop-dpr
      do kp=1,kfnp
        p=p+dpr
        es=aliq*exp((bliq*temp-cliq)/(temp-dliq))
        qs=0.622*es/(p-es)
        pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
        the0k(kp)=temp*pi*exp((3374.6525/temp-2.5403)*qs*        &
               (1.+0.81*qs))
      enddo   
!
! compute temperatures for each sat. equiv. potential temp.
!
      p=plutop-dpr
      do kp=1,kfnp
        thes=the0k(kp)-dth
        p=p+dpr
        do it=1,kfnt
! define sat. equiv. pot. temp.
          thes=thes+dth
! iterate to find temperature
! find initial guess
          if(it.eq.1) then
            tgues=tmin
          else
            tgues=ttab(it-1,kp)
          endif
          es=aliq*exp((bliq*tgues-cliq)/(tgues-dliq))
          qs=0.622*es/(p-es)
          pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
          thgues=tgues*pi*exp((3374.6525/tgues-2.5403)*qs*      &
               (1.+0.81*qs))
          f0=thgues-thes
          t1=tgues-0.5*f0
          t0=tgues
          itcnt=0
! iteration loop
          do itcnt=1,11
            es=aliq*exp((bliq*t1-cliq)/(t1-dliq))
            qs=0.622*es/(p-es)
            pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
            thtgs=t1*pi*exp((3374.6525/t1-2.5403)*qs*(1.+0.81*qs))
            f1=thtgs-thes
            if(abs(f1).lt.toler)then
              exit
            endif
!           itcnt=itcnt+1
            dt=f1*(t1-t0)/(f1-f0)
            t0=t1
            f0=f1
            t1=t1-dt
          enddo 
          ttab(it,kp)=t1 
          qstab(it,kp)=qs
        enddo
      enddo   
!
! lookup table for tlog(emix/aliq)
!
! set up intial values for lookup tables
!
       astrt=1.e-3
       ainc=0.075
!
       a1=astrt-ainc
       do i=1,200
         a1=a1+ainc
         alu(i)=alog(a1)
       enddo   
!
  END SUBROUTINE SW_LUTAB

! ***********************************************************************

END MODULE module_cu_slantwise
