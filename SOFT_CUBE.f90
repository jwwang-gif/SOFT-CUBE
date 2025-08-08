PROGRAM SOFT CUBE

      integer,parameter :: imax = 200, jmax = 200, kmax = 100 
      integer,parameter :: imax1 = imax+1,  jmax1 = jmax+1,  kmax1 = kmax+1 
      integer,parameter :: imax2 = imax1+1, jmax2 = jmax1+1, kmax2 = kmax1+1 
      integer:: mm, dd, hh, min_kk, min_jet
      real,dimension(imax2) :: dx, px
      real,dimension(jmax2) :: dy, py
      real,dimension(kmax2) :: dz, pz
      real,dimension(imax) :: x
      real,dimension(jmax) :: y
      real,dimension(kmax) :: z
      real,dimension(imax2, jmax2) :: uu1, vv1, ww1, dtemp1, rtemp1
      real,dimension(imax2, jmax2) :: uu2, vv2, ww2, dtemp2, rtemp2
      real,dimension(imax2, jmax2) :: uu3, vv3, ww3, dtemp3, rtemp3
      real,dimension(imax2, jmax2) :: uu4, vv4, ww4, dtemp4, rtemp4
      real,dimension(imax2, jmax2) :: uu5, vv5, ww5, dtemp5, rtemp5
      real,dimension(imax2, jmax2) :: uu, vv, ww, dtemp, rtemp, ctemp
      real,dimension(imax2, jmax2) :: fuu, fvv, fww
      real,dimension(imax, jmax, kmax) :: u, v, w, ws, t

      real,dimension(imax2,jmax2):: uu11, vv11, ww11, dt11, rt11
      real,dimension(imax2,jmax2):: uu12, vv12, ww12, dt12, rt12
      real,dimension(imax2,jmax2):: uu13, vv13, ww13, dt13, rt13
      real,dimension(imax2,jmax2):: uu14, vv14, ww14, dt14, rt14
      real,dimension(imax2,jmax2):: uu15, vv15, ww15, dt15, rt15
      real,dimension(imax2,jmax2):: uu16, vv16, ww16, dt16, rt16
      real,dimension(imax2,jmax2):: uu17, vv17, ww17, dt17, rt17
      real,dimension(imax2,jmax2):: uu18, vv18, ww18, dt18, rt18

      real,dimension(imax2,jmax2):: uu21, vv21, ww21, dt21, rt21
      real,dimension(imax2,jmax2):: uu22, vv22, ww22, dt22, rt22
      real,dimension(imax2,jmax2):: uu23, vv23, ww23, dt23, rt23
      real,dimension(imax2,jmax2):: uu24, vv24, ww24, dt24, rt24
      real,dimension(imax2,jmax2):: uu25, vv25, ww25, dt25, rt25
      real,dimension(imax2,jmax2):: uu26, vv26, ww26, dt26, rt26
      real,dimension(imax2,jmax2):: uu27, vv27, ww27, dt27, rt27
      real,dimension(imax2,jmax2):: uu28, vv28, ww28, dt28, rt28

      real,dimension(imax2,jmax2):: uu31, vv31, ww31, dt31, rt31
      real,dimension(imax2,jmax2):: uu32, vv32, ww32, dt32, rt32
      real,dimension(imax2,jmax2):: uu33, vv33, ww33, dt33, rt33
      real,dimension(imax2,jmax2):: uu34, vv34, ww34, dt34, rt34
      real,dimension(imax2,jmax2):: uu35, vv35, ww35, dt35, rt35
      real,dimension(imax2,jmax2):: uu36, vv36, ww36, dt36, rt36
      real,dimension(imax2,jmax2):: uu37, vv37, ww37, dt37, rt37
      real,dimension(imax2,jmax2):: uu38, vv38, ww38, dt38, rt38

      real,dimension(imax2,jmax2):: uu41, vv41, ww41, dt41, rt41
      real,dimension(imax2,jmax2):: uu42, vv42, ww42, dt42, rt42
      real,dimension(imax2,jmax2):: uu43, vv43, ww43, dt43, rt43
      real,dimension(imax2,jmax2):: uu44, vv44, ww44, dt44, rt44
      real,dimension(imax2,jmax2):: uu45, vv45, ww45, dt45, rt45
      real,dimension(imax2,jmax2):: uu46, vv46, ww46, dt46, rt46
      real,dimension(imax2,jmax2):: uu47, vv47, ww47, dt47, rt47
      real,dimension(imax2,jmax2):: uu48, vv48, ww48, dt48, rt48

      real,dimension(imax2,jmax2):: uu51, vv51, ww51, dt51, rt51
      real,dimension(imax2,jmax2):: uu52, vv52, ww52, dt52, rt52
      real,dimension(imax2,jmax2):: uu53, vv53, ww53, dt53, rt53
      real,dimension(imax2,jmax2):: uu54, vv54, ww54, dt54, rt54
      real,dimension(imax2,jmax2):: uu55, vv55, ww55, dt55, rt55
      real,dimension(imax2,jmax2):: uu56, vv56, ww56, dt56, rt56
      real,dimension(imax2,jmax2):: uu57, vv57, ww57, dt57, rt57
      real,dimension(imax2,jmax2):: uu58, vv58, ww58, dt58, rt58

      real,dimension(kmax2) :: ldaps_u, ldaps_v, ldaps_t, dir, dt, dws
      real,dimension(kmax2) :: ldaps_uu, ldaps_vv, ldaps_tt
      real,dimension(kmax2) ::ldaps_ws, ldaps_dir, inws_1, inws_2
      real,dimension(kmax2) ::ldaps_ws2, ldaps_dir2
      real,dimension(kmax2) :: inflow_ws_00, inflow_ws_01, inflow_ws_02, inflow_ws_03
      real,dimension(kmax2) :: inflow_ws_04, inflow_ws_05, inflow_ws_06, inflow_ws_07

      REAL :: dt1, dt2, dt3, dt4, rt1, rt2, rt3, rt4, ltemp, rwind, dwd
      REAL :: rws11, rws12, rws13, rws14, rws15, rws16, rws17, rws18
      REAL :: rt91, rt92, rt93, rt94, rt95, rt96, rt97, rt98
      REAL :: alpha, alpha2, pi, roof, road, green, soil, wall
      REAL :: wd_a, wd_b
      REAL :: dwd_a, dwd_b
      REAL :: ws_a, ws_b
      REAL :: dws_a, dws_b
      REAL :: vws_a, vws_b
      REAL :: temp_a, temp_b
      REAL :: roof_temp_a, roof_temp_b
      REAL :: wall_temp_a, wall_temp_b
      REAL :: road_temp_a, road_temp_b
      REAL :: soil_temp_a, soil_temp_b
      REAL :: green_temp_a, green_temp_b
      real,dimension(32) :: dir2
      INTEGER :: angle, angle2
      INTEGER :: temp1, temp2
      INTEGER :: roof_temp1, roof_temp2
      INTEGER :: wall_temp1, wall_temp2
      INTEGER :: road_temp1, road_temp2
      INTEGER :: soil_temp1, soil_temp2
      INTEGER :: green_temp1, green_temp2
      INTEGER :: cfd_ws, cfd_ws2
      INTEGER :: cfd_dws, cfd_dws2
      INTEGER :: cfd_dwd, cfd_dwd2
      INTEGER :: vcfd_ws1, vcfd_ws2
      INTEGER :: SC_yy, SC_mm, SC_dd, SC_hh, change_k

      CHARACTER(999) :: ifile11,ifile12,ifile13,ifile14,ifile15,ifile16,ifile17,ifile18
      CHARACTER(999) :: ifile21,ifile22,ifile23,ifile24,ifile25,ifile26,ifile27,ifile28
      CHARACTER(999) :: ifile31,ifile32,ifile33,ifile34,ifile35,ifile36,ifile37,ifile38
      CHARACTER(999) :: ifile41,ifile42,ifile43,ifile44,ifile45,ifile46,ifile47,ifile48
      CHARACTER(999) :: ifile51,ifile52,ifile53,ifile54,ifile55,ifile56,ifile57,ifile58
      CHARACTER(999) :: ifile61,ifile62,ifile63,ifile64
      CHARACTER(999) :: ifile81,ifile82,ifile83,ifile84,ifile85,ifile86,ifile87,ifile88
      CHARACTER(999) :: ofile1, ofile2, ofile3
      CHARACTER(999) :: l_in_u, l_in_v, l_in_t, l_in_tt, l_in_f

open(70,file='/Beaufort/Inflow_Lv_00.txt',status='old')
open(71,file='/Beaufort/Inflow_Lv_01.txt',status='old')
open(72,file='/Beaufort/Inflow_Lv_02.txt',status='old')
open(73,file='/Beaufort/Inflow_Lv_03.txt',status='old')
open(74,file='/Beaufort/Inflow_Lv_04.txt',status='old')
open(75,file='/Beaufort/Inflow_Lv_05.txt',status='old')
open(76,file='/Beaufort/Inflow_Lv_06.txt',status='old')
open(77,file='/Beaufort/Inflow_Lv_07.txt',status='old')

!c#########################################################################
!c  model cell size

      do 100 i=1,imax2
      dx(i) = 10.
  100 continue
      dx(1)     = 0.
      dx(imax2) = 0.

      do 120 j=1,jmax2
      dy(j) = 10.
  120 continue
      dy(1)     = 0.
      dy(jmax2) = 0. 

      do 140 k=1,kmax2
      dz(k) = 5.
  140 continue
      dz(1)     = 0.
      dz(kmax2) = 0.

!c position of each cell
      px(1) = 0. 
      px(2) = 0.5*dx(2) 
      do 160 i=3,imax1
  160 px(i) = px(i-1)+(dx(i-1)+dx(i))/2. 
      px(imax2) = px(imax1)+0.5*dx(imax1)
 
      py(1) = 0. 
      py(2) = 0.5*dy(2) 
      do 170 j=3,jmax1
  170 py(j) = py(j-1)+(dy(j-1)+dy(j))/2.
      py(jmax2) = py(jmax1)+0.5*dy(jmax1)

      pz(1) = 0.
      pz(2) = 0.5*dz(2)
      do 180 k=3,kmax1
  180 pz(k) = pz(k-1)+(dz(k-1)+dz(k))/2.
      pz(kmax2) = pz(kmax1)+0.5*dz(kmax1)


!c############################################################################
!c  LDAPS or AWS 
!c############################################################################

     do 7000 k = 1, kmax2
     READ(70,*) inflow_ws_00(k)
     READ(71,*) inflow_ws_01(k)
     READ(72,*) inflow_ws_02(k)
     READ(73,*) inflow_ws_03(k)
     READ(74,*) inflow_ws_04(k)
     READ(75,*) inflow_ws_05(k)
     READ(76,*) inflow_ws_06(k)
     READ(77,*) inflow_ws_07(k)
 7000 continue

  770 format(A,I4,I2.2,I2.2,A,I2.2,A)
  777 format(A,I4,I2.2,I2.2,I2.2,A)

     DO SC_yy = syy, syy
     DO SC_mm = smm, smm
     DO SC_dd = sdd, sdd
     DO SC_hh = shh, shh

     WRITE(l_in_f,777) '/RSA/ratio_',SC_yy,SC_mm,SC_dd,SC_hh,'.txt'
     PRINT*, TRIM(l_in_f)
     OPEN(8, file=l_in_f,status='old')
     READ(8,*) roof, road, green, soil, wall

     IF (roof .ge. 0.97 .and. roof .le. 0.995) THEN
     roof_temp1 = 1
     roof_temp2 = 2
     roof_temp_a = roof - 0.97
     roof_temp_b = 0.995 - roof
     ENDIF

     IF (roof .gt. 0.995 .and. roof .le. 1.02) THEN
     roof_temp1 = 2
     roof_temp2 = 3
     roof_temp_a = roof - 0.995
     roof_temp_b = 1.02 - roof
     ENDIF

     IF (roof .gt. 1.02 .and. roof .le. 1.045) THEN
     roof_temp1 = 3
     roof_temp2 = 4
     roof_temp_a = roof - 1.02
     roof_temp_b = 1.045 - roof
     ENDIF

     IF (roof .gt. 1.045 .and. roof .le. 1.07) THEN
     roof_temp1 = 4
     roof_temp2 = 5
     roof_temp_a = roof - 1.045
     roof_temp_b = 1.07 - roof
     ENDIF

     IF (roof .gt. 1.07) THEN
     roof_temp1 = 5
     roof_temp2 = 5
     ENDIF

     IF (roof .lt. 0.97) THEN
     roof_temp1 = 1
     roof_temp2 = 1
     ENDIF


     IF (wall .ge. 0.97 .and. wall .le. 0.985) THEN
     wall_temp1 = 1
     wall_temp2 = 2
     wall_temp_a = wall - 0.97
     wall_temp_b = 0.985 - wall
     ENDIF

     IF (wall .gt. 0.985 .and. wall .le. 1.0) THEN
     wall_temp1 = 2
     wall_temp2 = 3
     wall_temp_a = wall - 0.985
     wall_temp_b = 1.0 - wall
     ENDIF

     IF (wall .gt. 1.0 .and. wall .le. 1.015) THEN
     wall_temp1 = 3
     wall_temp2 = 4
     wall_temp_a = wall - 1.0
     wall_temp_b = 1.015 - wall
     ENDIF

     IF (wall .gt. 1.015 .and. wall .le. 1.03) THEN
     wall_temp1 = 4
     wall_temp2 = 5
     wall_temp_a = wall - 1.015
     wall_temp_b = 1.03 - wall
     ENDIF

     IF (wall .gt. 1.03) THEN
     wall_temp1 = 5
     wall_temp2 = 5
     ENDIF

     IF (wall .lt. 0.97) THEN
     wall_temp1 = 1
     wall_temp2 = 1
     ENDIF


     IF (road .ge. 0.975 .and. road .le. 0.995) THEN
     road_temp1 = 1
     road_temp2 = 2
     road_temp_a = road - 0.975
     road_temp_b = 0.995 - road
     ENDIF

     IF (road .gt. 0.995 .and. road .le. 1.015) THEN
     road_temp1 = 2
     road_temp2 = 3
     road_temp_a = road - 0.995
     road_temp_b = 1.015 - road
     ENDIF

     IF (road .gt. 1.015 .and. road .le. 1.035) THEN
     road_temp1 = 3
     road_temp2 = 4
     road_temp_a = road - 1.015
     road_temp_b = 1.035 - road
     ENDIF

     IF (road .gt. 1.035 .and. road .le. 1.055) THEN
     road_temp1 = 4
     road_temp2 = 5
     road_temp_a = road - 1.035
     road_temp_b = 1.055 - road
     ENDIF

     IF (road .gt. 1.055) THEN
     road_temp1 = 5
     road_temp2 = 5
     ENDIF

     IF (road .lt. 0.975) THEN
     road_temp1 = 1
     road_temp2 = 1
     ENDIF


     IF (green .ge. 0.95 .and. green .le. 0.97) THEN
     green_temp1 = 1
     green_temp2 = 2
     green_temp_a = green - 0.95
     green_temp_b = 0.97 - green
     ENDIF

     IF (green .gt. 0.97 .and. green .le. 0.99) THEN
     green_temp1 = 2
     green_temp2 = 3
     green_temp_a = green - 0.97
     green_temp_b = 0.99 - green
     ENDIF

     IF (green .gt. 0.99 .and. green .le. 1.01) THEN
     green_temp1 = 3
     green_temp2 = 4
     green_temp_a = green - 0.99
     green_temp_b = 1.01 - green
     ENDIF

     IF (green .gt. 1.01 .and. green .le. 1.03) THEN
     green_temp1 = 4
     green_temp2 = 5
     green_temp_a = green - 1.01
     green_temp_b = 1.03 - green
     ENDIF

     IF (green .gt. 1.03) THEN
     green_temp1 = 5
     green_temp2 = 5
     ENDIF

     IF (green .lt. 0.95) THEN
     green_temp1 = 1
     green_temp2 = 1
     ENDIF


     IF (soil .ge. 0.96 .and. soil .le. 0.98) THEN
     soil_temp1 = 1
     soil_temp2 = 2
     soil_temp_a = soil - 0.96
     soil_temp_b = 0.98 - soil

     ENDIF

     IF (soil .gt. 0.98 .and. soil .le. 1.) THEN
     soil_temp1 = 2
     soil_temp2 = 3
     soil_temp_a = soil - 0.98
     soil_temp_b = 1.0 - soil
     ENDIF

     IF (soil .gt. 1. .and. soil .le. 1.02) THEN
     soil_temp1 = 3
     soil_temp2 = 4
     soil_temp_a = soil - 1.0
     soil_temp_b = 1.02 - soil
     ENDIF

     IF (soil .gt. 1.02 .and. soil .le. 1.04) THEN
     soil_temp1 = 4
     soil_temp2 = 5
     soil_temp_a = soil - 1.02
     soil_temp_b = 1.04 - soil
     ENDIF

     IF (soil .gt. 1.04) THEN
     soil_temp1 = 5
     soil_temp2 = 5
     ENDIF

     IF (soil .lt. 0.96 ) THEN
     soil_temp1 = 1
     soil_temp2 = 1
     ENDIF


     WRITE(l_in_u,777) '/LDAPS/U_',gk_yy,gk_mm,gk_dd,gk_hh,'.txt'
     WRITE(l_in_v,777) '/LDAPS/V_',gk_yy,gk_mm,gk_dd,gk_hh,'.txt'
     WRITE(l_in_t,777) '/LDAPS/T_',gk_yy,gk_mm,gk_dd,gk_hh,'.txt'

     PRINT*, TRIM(l_in_u) 
     PRINT*, TRIM(l_in_v) 
     PRINT*, TRIM(l_in_t) 

     OPEN(5, file=l_in_u,status='old')
     OPEN(6, file=l_in_v,status='old')
     OPEN(7, file=l_in_t,status='old')

     DO k = 1, kmax2
     READ(5,*) ldaps_uu(k)
     READ(6,*) ldaps_vv(k)
     READ(7,*) ldaps_tt(k)
     ENDDO


     pi=atan(1.)*4.

     DO k = 1, kmax2
     ldaps_u(k) = ldaps_uu(k)
     ldaps_v(k) = ldaps_vv(k)
     ldaps_t(k) = ldaps_tt(k)

      IF (ldaps_v(k).GT.0.) THEN
      ldaps_dir(k) = ((180./pi)*atan(ldaps_u(k)/ldaps_v(k)))+180.
      ELSE IF (ldaps_u(k) .LT. 0. .AND. ldaps_v(k).LT.0.) THEN
      ldaps_dir(k) = ((180./pi)*atan(ldaps_u(k)/ldaps_v(k)))+0.
      ELSE IF (ldaps_u(k) .GT. 0. .AND. ldaps_v(k).LT.0.) THEN
      ldaps_dir(k) = ((180./pi)*atan(ldaps_u(k)/ldaps_v(k)))+360.
      ELSE
      ldaps_dir(k) = 0.
      END IF
      ldaps_ws(k) = SQRT(ldaps_u(k)**2.0 + ldaps_v(k)**2.0)
      END DO 


     do k = 1, kmax2
      dws(k) = (ldaps_ws(k+1)-ldaps_ws(k-1))*100
      dt(kk) = (ldaps_tt(kk+1)-ldaps_tt(kk-1))*100
     end do

      dws(kmax2) = dws(kmax1)
      dt(kmax2) = dt(kmax1)

     dwd = abs(ldaps_dir(kmax2) - ldaps_dir(21))

     if(dwd .gt. 180.) then
     dwd = 360. - dwd
     end if

     IF (dwd .gt. 0. .and. dwd .le. 45.) THEN
     cfd_dwd = 1
     cfd_dwd2= 2
     dwd_a = 45. - dwd
     dwd_b =  dwd
     ENDIF

     IF (dwd .gt. 45. .and. dwd .le. 90.) THEN
     cfd_dwd = 2
     cfd_dwd2= 3
     dwd_a = 90. - dwd
     dwd_b =  dwd - 45.
     ENDIF

     IF (dwd .gt. 90. .and. dwd .le. 135.) THEN
     cfd_dwd = 3
     cfd_dwd2= 4
     dwd_a = 135. - dwd
     dwd_b = dwd - 90.
     ENDIF

     IF (dwd .gt. 135. .and. dwd .le. 180.) THEN
     cfd_dwd = 4
     cfd_dwd2= 5
     dwd_a = 180. - dwd
     dwd_b = dwd - 135.
     ENDIF

     IF (dwd .eq. 180.) THEN
     cfd_dwd = 5
     cfd_dwd2= 5  ! weighting x 0.0
     ENDIF

     IF (dwd .eq. 0.) THEN
     cfd_dwd = 1
     cfd_dwd2= 1  ! weighting x 0.0
     ENDIF

     
     DO kk = 1, kmax2          !!!!!!!! Vertical Cycle Loop !!!!!!!!!!

     IF (dws(kk) .gt. -10. .and. dws(kk) .le. 0.) THEN
     cfd_dws = 1
     cfd_dws2= 2
     dws_a = dws(kk) + 10.
     dws_b =  - dws(kk)
     ENDIF

     IF (dws(kk) .gt. 0 .and. dws(kk) .le. 10.) THEN
     cfd_dws = 2
     cfd_dws2= 3
     dws_a = 10. - dws(kk)
     dws_b =  dws(kk)
     ENDIF

     IF (dws(kk) .gt. 10. .and. dws(kk) .le. 20) THEN
     cfd_dws = 3
     cfd_dws2= 4
     dws_a = 20. - dws(kk)
     dws_b = dws(kk) - 10.
     ENDIF

     IF (dws(kk) .gt. 20. .and. dws(kk) .le. 30.) THEN
     cfd_dws = 4
     cfd_dws2= 5
     dws_a = 30. - dws(kk)
     dws_b = dws(kk) - 20.
     ENDIF

     IF (dws(kk) .gt. 30.) THEN
     cfd_dws = 5
     cfd_dws2= 5  ! weighting x 0.0
     ENDIF

     IF (dws(kk) .lt. -10.) THEN
     cfd_dws = 1
     cfd_dws2= 1  ! weighting x 0.0
     ENDIF


     IF (ldaps_ws2(kk) .gt. inflow_ws_00(kk) .and. ldaps_ws2(kk) .le. inflow_ws_01(kk)) THEN
     cfd_ws = 00
     cfd_ws2= 01
     inws_1(kk) = inflow_ws_00(kk)
     inws_2(kk) = inflow_ws_01(kk)
     ENDIF

     IF (ldaps_ws2(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws2(kk) .le. inflow_ws_02(kk)) THEN
     cfd_ws = 01
     cfd_ws2= 02
     inws_1(kk) = inflow_ws_01(kk)
     inws_2(kk) = inflow_ws_02(kk)
     ENDIF

     IF (ldaps_ws2(kk) .gt. inflow_ws_02(kk) .and. ldaps_ws2(kk) .le. inflow_ws_03(kk)) THEN
     cfd_ws = 02
     cfd_ws2= 03
     inws_1(kk) = inflow_ws_02(kk)
     inws_2(kk) = inflow_ws_03(kk)
     ENDIF

     IF (ldaps_ws2(kk) .gt. inflow_ws_03(kk) .and. ldaps_ws2(kk) .le. inflow_ws_04(kk)) THEN
     cfd_ws = 03
     cfd_ws2= 04
     inws_1(kk) = inflow_ws_03(kk)
     inws_2(kk) = inflow_ws_04(kk)
     ENDIF

     IF (ldaps_ws2(kk) .gt. inflow_ws_04(kk) .and. ldaps_ws2(kk) .le. inflow_ws_05(kk)) THEN
     cfd_ws = 04
     cfd_ws2= 05
     inws_1(kk) = inflow_ws_04(kk)
     inws_2(kk) = inflow_ws_05(kk)
     ENDIF

     IF (ldaps_ws2(kk) .gt. inflow_ws_05(kk) .and. ldaps_ws2(kk) .le. inflow_ws_06(kk)) THEN
     cfd_ws = 05
     cfd_ws2= 06
     inws_1(kk) = inflow_ws_05(kk)
     inws_2(kk) = inflow_ws_06(kk)
     ENDIF

     IF (ldaps_ws2(kk) .gt. inflow_ws_06(kk) .and. ldaps_ws2(kk) .le. inflow_ws_07(kk)) THEN
     cfd_ws = 06
     cfd_ws2= 07
     inws_1(kk) = inflow_ws_06(kk)
     inws_2(kk) = inflow_ws_07(kk)
     ENDIF

     IF (ldaps_ws2(kk) .gt. inflow_ws_07(kk)) THEN
     cfd_ws = 07
     cfd_ws2= 07  ! weighting x 0.0
     inws_1(kk) = inflow_ws_07(kk)
     inws_2(kk) = inflow_ws_07(kk)
     ENDIF

     IF (ldaps_ws2(kk) .le. inflow_ws_00(kk)) THEN
     cfd_ws = 00
     cfd_ws2= 00  ! weighting x 0.0
     inws_1(kk) = inflow_ws_00(kk)
     inws_2(kk) = inflow_ws_00(kk)
     ENDIF

!!!!!!!!!!!!!!!!! Vertical momentum DB condition !!!!!!!!!!!!!!
     IF (ldaps_ws(21) .gt. inflow_ws_00(21) .and. ldaps_ws(21) .le. inflow_ws_01(21)) THEN
     vcfd_ws1 = 00
     vcfd_ws2 = 01
     vws_a = ldaps_ws(21) -  inws_1(21)
     vws_b = inws_2(21) - ldaps_ws(21)
     ENDIF

     IF (ldaps_ws(21) .gt. inflow_ws_01(21) .and. ldaps_ws(21) .le. inflow_ws_02(21)) THEN
     vcfd_ws1 = 01
     vcfd_ws2 = 02
     vws_a = ldaps_ws(21) -  inws_1(21)
     vws_b = inws_2(21) - ldaps_ws(21)
     ENDIF

     IF (ldaps_ws(21) .gt. inflow_ws_02(21) .and. ldaps_ws(21) .le. inflow_ws_03(21)) THEN
     vcfd_ws1 = 02
     vcfd_ws2 = 03
     vws_a = ldaps_ws(21) -  inws_1(21)
     vws_b = inws_2(21) - ldaps_ws(21)
     ENDIF

     IF (ldaps_ws(21) .gt. inflow_ws_03(21) .and. ldaps_ws(21) .le. inflow_ws_04(21)) THEN
     vcfd_ws1 = 03
     vcfd_ws2 = 04
     vws_a = ldaps_ws(21) -  inws_1(21)
     vws_b = inws_2(21) - ldaps_ws(21)
     ENDIF

     IF (ldaps_ws(21) .gt. inflow_ws_04(21) .and. ldaps_ws(21) .le. inflow_ws_05(21)) THEN
     vcfd_ws1 = 04
     vcfd_ws2 = 05
     vws_a = ldaps_ws(21) -  inws_1(21)
     vws_b = inws_2(21) - ldaps_ws(21)
     ENDIF

     IF (ldaps_ws(21) .gt. inflow_ws_05(21) .and. ldaps_ws(21) .le. inflow_ws_06(21)) THEN
     vcfd_ws1 = 05
     vcfd_ws2 = 06
     vws_a = ldaps_ws(21) -  inws_1(21)
     vws_b = inws_2(21) - ldaps_ws(21)
     ENDIF

     IF (ldaps_ws(21) .gt. inflow_ws_06(21) .and. ldaps_ws(21) .le. inflow_ws_07(21)) THEN
     vcfd_ws1 = 06
     vcfd_ws2 = 07
     vws_a = ldaps_ws(21) -  inws_1(21)
     vws_b = inws_2(21) - ldaps_ws(21)
     ENDIF

     IF (ldaps_ws(21) .gt. inflow_ws_07(21)) THEN
     vcfd_ws1 = 07
     vcfd_ws2 = 07  ! weighting x 0.0
     vws_a = 1
     vws_b = 1
     ENDIF

     IF (ldaps_ws(21) .le. inflow_ws_00(21)) THEN
     vcfd_ws1 = 00
     vcfd_ws2 = 00  ! weighting x 0.0
     vws_a = 1
     vws_b = 1.
     ENDIF

!!!!!!!!!!!!!!!!!!!! Lapse rate DB !!!!!!!!!!!!!!!!!!!!!!!!!!

     IF (dt(kk) .ge. -20. .and. dt(kk) .le. -10.) THEN
     temp1 = 1
     temp2 = 2
     temp_a = dt(kk) + 20.
     temp_b = -10 - dt(kk)
     ENDIF

     IF (dt(kk) .gt. -10 .and. dt(kk) .le. 0.0) THEN
     temp1 = 2
     temp2 = 3
     temp_a = dt(kk) + 10.
     temp_b =  - dt(kk)
     ENDIF

     IF (dt(kk) .gt. 0.0 .and. dt(kk) .le. 10.) THEN
     temp1 = 3
     temp2 = 4
     temp_a = dt(kk)
     temp_b = 10. - dt(kk)
     ENDIF

     IF (dt(kk) .gt. 10. .and. dt(kk) .le. 20.) THEN
     temp1 = 4
     temp2 = 5
     temp_a = dt(kk) - 10
     temp_b = 20 - dt(kk)
     ENDIF

     IF (dt(kk) .gt. 20) THEN
     temp1 = 5
     temp2 = 5
     ENDIF

     IF (dt(kk) .lt. -20) THEN
     temp1 = 1
     temp2 = 1
     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !- TRANS. W.D  ----------------------
     ! AWS, LDAPS -> 0 -> N, 90 -> E ...
     !        CFD -> 0 -> W, 90 -> S ...
     IF (ldaps_dir2(kk) <= 270) THEN
     dir(kk) = 270 -  ldaps_dir2(kk)
     ENDIF

     IF (ldaps_dir2(kk) >= 270) THEN
     dir(kk) = 360 - (ldaps_dir2(kk) - 270)
     ENDIF
     !-------------------------------------

!     PRINT*, 'FOR CFD dir=',dir(k) ! from AWS or LDAPS

     DO i = 1, 32
     ! 0   -> w  ;  90  -> s
     ! 180 -> e  ;  270 -> n
      dir2(i) = i*11.25 - 11.25  
     ENDDO

     ! find out nearest angle -> alpha
     ! ex) dir -> alpha  
     !    2.1 ->  0.00
     !   11.3 -> 11.25
     DO i = 1,32
      IF( dir(kk) >= dir2(i) ) THEN
        alpha = dir2(i)
        angle  = i
        angle2 = i+1
      ENDIF
     ENDDO

     alpha2= dir(kk) - alpha

     ! weighing value -> wd_a, wd_b

     ws_a = ldaps_ws2(kk) -  inws_1(kk)
     ws_b = inws_2(kk) - ldaps_ws2(kk)

     wd_a = (11.25 -alpha2)/11.25
     wd_b = alpha2/11.25


    WRITE(ifile11, 888) '/Roof/Log_',cfd_ws,'/Class_',roof_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile12, 888) '/Roof/Log_',cfd_ws,'/Class_',roof_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile13, 888) '/Roof/Log_',cfd_ws2,'/Class_',roof_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile14, 888) '/Roof/Log_',cfd_ws2,'/Class_',roof_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile15, 888) '/Roof/Log_',cfd_ws,'/Class_',roof_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile16, 888) '/Roof/Log_',cfd_ws,'/Class_',roof_temp2,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile17, 888) '/Roof/Log_',cfd_ws2,'/Class_',roof_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile18, 888) '/Roof/Log_',cfd_ws2,'/Class_',roof_temp2,'/k=',kk,'_moment_',angle2,'.dat'

    WRITE(ifile21, 888) '/Wall/Log_',cfd_ws,'/Class_',wall_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile22, 888) '/Wall/Log_',cfd_ws,'/Class_',wall_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile23, 888) '/Wall/Log_',cfd_ws2,'/Class_',wall_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile24, 888) '/Wall/Log_',cfd_ws2,'/Class_',wall_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile25, 888) '/Wall/Log_',cfd_ws,'/Class_',wall_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile26, 888) '/Wall/Log_',cfd_ws,'/Class_',wall_temp2,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile27, 888) '/Wall/Log_',cfd_ws2,'/Class_',wall_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile28, 888) '/Wall/Log_',cfd_ws2,'/Class_',wall_temp2,'/k=',kk,'_moment_',angle2,'.dat'

    WRITE(ifile31, 888) '/Road/Log_',cfd_ws,'/Class_',road_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile32, 888) '/Road/Log_',cfd_ws,'/Class_',road_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile33, 888) '/Road/Log_',cfd_ws2,'/Class_',road_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile34, 888) '/Road/Log_',cfd_ws2,'/Class_',road_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile35, 888) '/Road/Log_',cfd_ws,'/Class_',road_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile36, 888) '/Road/Log_',cfd_ws,'/Class_',road_temp2,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile37, 888) '/Road/Log_',cfd_ws2,'/Class_',road_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile38, 888) '/Road/Log_',cfd_ws2,'/Class_',road_temp2,'/k=',kk,'_moment_',angle2,'.dat'

    WRITE(ifile41, 888) '/Green/Log_',cfd_ws,'/Class_',green_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile42, 888) '/Green/Log_',cfd_ws,'/Class_',green_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile43, 888) '/Green/Log_',cfd_ws2,'/Class_',green_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile44, 888) '/Green/Log_',cfd_ws2,'/Class_',green_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile45, 888) '/Green/Log_',cfd_ws,'/Class_',green_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile46, 888) '/Green/Log_',cfd_ws,'/Class_',green_temp2,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile47, 888) '/Green/Log_',cfd_ws2,'/Class_',green_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile48, 888) '/Green/Log_',cfd_ws2,'/Class_',green_temp2,'/k=',kk,'_moment_',angle2,'.dat'

    WRITE(ifile51, 888) '/Soil/Log_',cfd_ws,'/Class_',soil_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile52, 888) '/Soil/Log_',cfd_ws,'/Class_',soil_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile53, 888) '/Soil/Log_',cfd_ws2,'/Class_',soil_temp1,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile54, 888) '/Soil/Log_',cfd_ws2,'/Class_',soil_temp1,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile55, 888) '/Soil/Log_',cfd_ws,'/Class_',soil_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile56, 888) '/Soil/Log_',cfd_ws,'/Class_',soil_temp2,'/k=',kk,'_moment_',angle2,'.dat'
    WRITE(ifile57, 888) '/Soil/Log_',cfd_ws2,'/Class_',soil_temp2,'/k=',kk,'_moment_',angle,'.dat'
    WRITE(ifile58, 888) '/Soil/Log_',cfd_ws2,'/Class_',soil_temp2,'/k=',kk,'_moment_',angle2,'.dat'


    WRITE(ifile61, 887) '/Lapse_rate/temperature/',temp1,'/k=',kk,'_moment_',vcfd_ws1,'_02.dat'
    WRITE(ifile62, 887) '/Lapse_rate/temperature/',temp2,'/k=',kk,'_moment_',vcfd_ws1,'_02.dat'
    WRITE(ifile63, 887) '/Lapse_rate/temperature/',temp1,'/k=',kk,'_moment_',vcfd_ws2,'_02.dat'
    WRITE(ifile64, 887) '/Lapse_rate/temperature/',temp2,'/k=',kk,'_moment_',vcfd_ws2,'_02.dat'

    WRITE(ifile81, 999) '/Lapse_rate/windspeed/',cfd_dwd,'/k=',kk,'_moment_',vcfd_ws1,'_',cfd_dws,'.dat'
    WRITE(ifile82, 999) '/Lapse_rate/windspeed/',cfd_dwd,'/k=',kk,'_moment_',vcfd_ws1,'_',cfd_dws2,'.dat'
    WRITE(ifile83, 999) '/Lapse_rate/windspeed/',cfd_dwd,'/k=',kk,'_moment_',vcfd_ws2,'_',cfd_dws,'.dat'
    WRITE(ifile84, 999) '/Lapse_rate/windspeed/',cfd_dwd,'/k=',kk,'_moment_',vcfd_ws2,'_',cfd_dws2,'.dat'
    WRITE(ifile85, 999) '/Lapse_rate/windspeed/',cfd_dwd2,'/k=',kk,'_moment_',vcfd_ws1,'_',cfd_dws,'.dat'
    WRITE(ifile86, 999) '/Lapse_rate/windspeed/',cfd_dwd2,'/k=',kk,'_moment_',vcfd_ws1,'_',cfd_dws2,'.dat'
    WRITE(ifile87, 999) '/Lapse_rate/windspeed/',cfd_dwd2,'/k=',kk,'_moment_',vcfd_ws2,'_',cfd_dws,'.dat'
    WRITE(ifile88, 999) '/Lapse_rate/windspeed/',cfd_dwd2,'/k=',kk,'_moment_',vcfd_ws2,'_',cfd_dws2,'.dat'

      WRITE(ofile1, 778) '/SOFT_CUBE/SOFT_CUBE_',gk_yy,gk_mm,gk_dd,gk_hh,'.dat'

     open(11,file=ifile11,status='old')
     open(12,file=ifile12,status='old')
     open(13,file=ifile13,status='old')
     open(14,file=ifile14,status='old')
     open(15,file=ifile15,status='old')
     open(16,file=ifile16,status='old')
     open(17,file=ifile17,status='old')
     open(18,file=ifile18,status='old')

     open(21,file=ifile21,status='old')
     open(22,file=ifile22,status='old')
     open(23,file=ifile23,status='old')
     open(24,file=ifile24,status='old')
     open(25,file=ifile25,status='old')
     open(26,file=ifile26,status='old')
     open(27,file=ifile27,status='old')
     open(28,file=ifile28,status='old')

     open(31,file=ifile31,status='old')
     open(32,file=ifile32,status='old')
     open(33,file=ifile33,status='old')
     open(34,file=ifile34,status='old')
     open(35,file=ifile35,status='old')
     open(36,file=ifile36,status='old')
     open(37,file=ifile37,status='old')
     open(38,file=ifile38,status='old')

     open(41,file=ifile41,status='old')
     open(42,file=ifile42,status='old')
     open(43,file=ifile43,status='old')
     open(44,file=ifile44,status='old')
     open(45,file=ifile45,status='old')
     open(46,file=ifile46,status='old')
     open(47,file=ifile47,status='old')
     open(48,file=ifile48,status='old')

     open(51,file=ifile51,status='old')
     open(52,file=ifile52,status='old')
     open(53,file=ifile53,status='old')
     open(54,file=ifile54,status='old')
     open(55,file=ifile55,status='old')
     open(56,file=ifile56,status='old')
     open(57,file=ifile57,status='old')
     open(58,file=ifile58,status='old')

     open(61,file=ifile61,status='old')
     open(62,file=ifile62,status='old')
     open(63,file=ifile63,status='old')
     open(64,file=ifile64,status='old')

     open(81,file=ifile81,status='old')
     open(82,file=ifile82,status='old')
     open(83,file=ifile83,status='old')
     open(84,file=ifile84,status='old')
     open(85,file=ifile85,status='old')
     open(86,file=ifile86,status='old')
     open(87,file=ifile87,status='old')
     open(88,file=ifile88,status='old')

     open(60,file=ofile1,status='unknown')

   887 format(A,I1.1,A,I3.3,A,I2.2,A)
   888 format(A,I2.2,A,I1.1,A,I3.3,A,I2.2,A)
   778 format(A,I4,I2.2,I2.2,I2.2,A)
   999 format(A,i1.1,a,I3.3,A,I2.2,A,I2.2,A)

!c############################################################################
!c open wind data
!c############################################################################

      read(61,*) dt1, rt1
      read(62,*) dt2, rt2
      read(63,*) dt3, rt3
      read(64,*) dt4, rt4
      read(81,*) rws11, rt91
      read(82,*) rws12, rt92
      read(83,*) rws13, rt93
      read(84,*) rws14, rt94
      read(85,*) rws15, rt95
      read(86,*) rws16, rt96
      read(87,*) rws17, rt97
      read(88,*) rws18, rt98

      do 6000 j=1,jmax2
      do 6000 i=1,imax2

      read(11,*) uu11(i,j),vv11(i,j),ww11(i,j),dt11(i,j),rt11(i,j)
      read(12,*) uu12(i,j),vv12(i,j),ww12(i,j),dt12(i,j),rt12(i,j)
      read(13,*) uu13(i,j),vv13(i,j),ww13(i,j),dt13(i,j),rt13(i,j)
      read(14,*) uu14(i,j),vv14(i,j),ww14(i,j),dt14(i,j),rt14(i,j)
      read(15,*) uu15(i,j),vv15(i,j),ww15(i,j),dt15(i,j),rt15(i,j)
      read(16,*) uu16(i,j),vv16(i,j),ww16(i,j),dt16(i,j),rt16(i,j)
      read(17,*) uu17(i,j),vv17(i,j),ww17(i,j),dt17(i,j),rt17(i,j)
      read(18,*) uu18(i,j),vv18(i,j),ww18(i,j),dt18(i,j),rt18(i,j)

      read(21,*) uu21(i,j),vv21(i,j),ww21(i,j),dt21(i,j),rt21(i,j)
      read(22,*) uu22(i,j),vv22(i,j),ww22(i,j),dt22(i,j),rt22(i,j)
      read(23,*) uu23(i,j),vv23(i,j),ww23(i,j),dt23(i,j),rt23(i,j)
      read(24,*) uu24(i,j),vv24(i,j),ww24(i,j),dt24(i,j),rt24(i,j)
      read(25,*) uu25(i,j),vv25(i,j),ww25(i,j),dt25(i,j),rt25(i,j)
      read(26,*) uu26(i,j),vv26(i,j),ww26(i,j),dt26(i,j),rt26(i,j)
      read(27,*) uu27(i,j),vv27(i,j),ww27(i,j),dt27(i,j),rt27(i,j)
      read(28,*) uu28(i,j),vv28(i,j),ww28(i,j),dt28(i,j),rt28(i,j)

      read(31,*) uu31(i,j),vv31(i,j),ww31(i,j),dt31(i,j),rt31(i,j)
      read(32,*) uu32(i,j),vv32(i,j),ww32(i,j),dt32(i,j),rt32(i,j)
      read(33,*) uu33(i,j),vv33(i,j),ww33(i,j),dt33(i,j),rt33(i,j)
      read(34,*) uu34(i,j),vv34(i,j),ww34(i,j),dt34(i,j),rt34(i,j)
      read(35,*) uu35(i,j),vv35(i,j),ww35(i,j),dt35(i,j),rt35(i,j)
      read(36,*) uu36(i,j),vv36(i,j),ww36(i,j),dt36(i,j),rt36(i,j)
      read(37,*) uu37(i,j),vv37(i,j),ww37(i,j),dt37(i,j),rt37(i,j)
      read(38,*) uu38(i,j),vv38(i,j),ww38(i,j),dt38(i,j),rt38(i,j)

      read(41,*) uu41(i,j),vv41(i,j),ww41(i,j),dt41(i,j),rt41(i,j)
      read(42,*) uu42(i,j),vv42(i,j),ww42(i,j),dt42(i,j),rt42(i,j)
      read(43,*) uu43(i,j),vv43(i,j),ww43(i,j),dt43(i,j),rt43(i,j)
      read(44,*) uu44(i,j),vv44(i,j),ww44(i,j),dt44(i,j),rt44(i,j)
      read(45,*) uu45(i,j),vv45(i,j),ww45(i,j),dt45(i,j),rt45(i,j)
      read(46,*) uu46(i,j),vv46(i,j),ww46(i,j),dt46(i,j),rt46(i,j)
      read(47,*) uu47(i,j),vv47(i,j),ww47(i,j),dt47(i,j),rt47(i,j)
      read(48,*) uu48(i,j),vv48(i,j),ww48(i,j),dt48(i,j),rt48(i,j)

      read(51,*) uu51(i,j),vv51(i,j),ww51(i,j),dt51(i,j),rt51(i,j)
      read(52,*) uu52(i,j),vv52(i,j),ww52(i,j),dt52(i,j),rt52(i,j)
      read(53,*) uu53(i,j),vv53(i,j),ww53(i,j),dt53(i,j),rt53(i,j)
      read(54,*) uu54(i,j),vv54(i,j),ww54(i,j),dt54(i,j),rt54(i,j)
      read(55,*) uu55(i,j),vv55(i,j),ww55(i,j),dt55(i,j),rt55(i,j)
      read(56,*) uu56(i,j),vv56(i,j),ww56(i,j),dt56(i,j),rt56(i,j)
      read(57,*) uu57(i,j),vv57(i,j),ww57(i,j),dt57(i,j),rt57(i,j)
      read(58,*) uu58(i,j),vv58(i,j),ww58(i,j),dt58(i,j),rt58(i,j)

 5000 format(5f12.6)

 6000 continue


!c############################################################################
!c  Composite wind & temperature field
!c############################################################################
      DO j=1, jmax2
      DO i=1, imax2

!!!!!!!! 1. Roof !!!!!!!!!!!!!!
      IF (ldaps_ws(kk) < inflow_ws_00(kk)) THEN

      if (roof < 0.97) then
      uu1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu12(i,j)
      vv1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv12(i,j)
      ww1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww12(i,j)
      dtemp1(i,j) =  (wd_a*dt11(i,j) + wd_b*dt12(i,j))
      rtemp1(i,j) =  (wd_a*rt11(i,j) + wd_b*rt12(i,j))

      else if (roof > 1.07) then
      uu1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu12(i,j)
      vv1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv12(i,j)
      ww1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww12(i,j)
      dtemp1(i,j) =  (wd_a*dt11(i,j) + wd_b*dt12(i,j))
      rtemp1(i,j) =  (wd_a*rt11(i,j) + wd_b*rt12(i,j))

      else IF (roof >= 0.97 .and. roof <= 1.07) THEN
      uu1(i,j)=  wd_a*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*uu11(i,j) + wd_a*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*uu15(i,j) + wd_b*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*uu12(i,j) + wd_b*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*uu16(i,j)
      vv1(i,j)=  wd_a*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*vv11(i,j) + wd_a*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*vv15(i,j) + wd_b*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*vv12(i,j) + wd_b*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*vv16(i,j)
      ww1(i,j)=  wd_a*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*ww11(i,j) + wd_a*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*ww15(i,j) + wd_b*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*ww12(i,j) + wd_b*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_00(kk)*ww16(i,j)
      dtemp1(i,j)=  wd_a*(roof_temp_b/0.025)*dt11(i,j) + wd_a*(roof_temp_a/0.025)*dt15(i,j) + wd_b*(roof_temp_b/0.025)*dt12(i,j) + wd_b*(roof_temp_a/0.025)*dt16(i,j)
      rtemp1(i,j)=  wd_a*(roof_temp_b/0.025)*rt11(i,j) + wd_a*(roof_temp_a/0.025)*rt15(i,j) + wd_b*(roof_temp_b/0.025)*rt12(i,j) + wd_b*(roof_temp_a/0.025)*rt16(i,j)

      end if
      ENDIF

      IF (ldaps_ws(kk) > inflow_ws_07(kk)) THEN

      if (roof < 0.97) then
      uu1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu12(i,j)
      vv1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv12(i,j)
      ww1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww12(i,j)
      dtemp1(i,j) =  (wd_a*dt11(i,j) + wd_b*dt12(i,j))
      rtemp1(i,j) =  (wd_a*rt11(i,j) + wd_b*rt12(i,j))

      else if (roof > 1.07) then
      uu1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu12(i,j)
      vv1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv12(i,j)
      ww1(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww11(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww12(i,j)
      dtemp1(i,j) =  (wd_a*dt11(i,j) + wd_b*dt12(i,j))
      rtemp1(i,j) =  (wd_a*rt11(i,j) + wd_b*rt12(i,j))

      else IF (roof >= 0.97 .and. roof <= 1.07) THEN
      uu1(i,j)=  wd_a*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*uu11(i,j) + wd_a*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*uu15(i,j) + wd_b*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*uu12(i,j) + wd_b*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*uu16(i,j)
      vv1(i,j)=  wd_a*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*vv11(i,j) + wd_a*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*vv15(i,j) + wd_b*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*vv12(i,j) + wd_b*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*vv16(i,j)
      ww1(i,j)=  wd_a*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*ww11(i,j) + wd_a*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*ww15(i,j) + wd_b*(roof_temp_b/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*ww12(i,j) + wd_b*(roof_temp_a/0.025)*ldaps_ws(kk)/inflow_ws_07(kk)*ww16(i,j)
      dtemp1(i,j)=  wd_a*(roof_temp_b/0.025)*dt11(i,j) + wd_a*(roof_temp_a/0.025)*dt15(i,j) + wd_b*(roof_temp_b/0.025)*dt12(i,j) + wd_b*(roof_temp_a/0.025)*dt16(i,j)
      rtemp1(i,j)=  wd_a*(roof_temp_b/0.025)*rt11(i,j) + wd_a*(roof_temp_a/0.025)*rt15(i,j) + wd_b*(roof_temp_b/0.025)*rt12(i,j) + wd_b*(roof_temp_a/0.025)*rt16(i,j)
      end if
      ENDIF
      !S1_end ---------------------------------------!

      !S2.-  interpolation -------------------------!
      IF (ldaps_ws(kk) >= inflow_ws_00(kk) .and. ldaps_ws(kk) <= inflow_ws_07(kk)) THEN

      if (roof < 0.97) then
      uu1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu14(i,j)
      vv1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv14(i,j)
      ww1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww14(i,j)
      dtemp1(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt14(i,j)
      rtemp1(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt14(i,j)

      else if (roof > 1.07) then
      uu1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu14(i,j)
      vv1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv14(i,j)
      ww1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww14(i,j)
      dtemp1(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt14(i,j)
      rtemp1(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt11(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt13(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt12(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt14(i,j)

      else IF (roof >= 0.97 .and. roof <= 1.07) THEN
      uu1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*uu11(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*uu12(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*uu13(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*uu14(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*uu15(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*uu16(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*uu17(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*uu18(i,j)
      vv1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*vv11(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*vv12(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*vv13(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*vv14(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*vv15(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*vv16(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*vv17(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*vv18(i,j)
      ww1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*ww11(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*ww12(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*ww13(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*ww14(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*ww15(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*ww16(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*ww17(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*ww18(i,j)
      dtemp1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*dt11(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*dt12(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*dt13(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*dt14(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*dt15(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*dt16(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*dt17(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*dt18(i,j)
      rtemp1(i,j)= wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*rt11(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_b/0.025)*rt12(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*rt13(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_b/0.025)*rt14(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*rt15(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(roof_temp_a/0.025)*rt16(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*rt17(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(roof_temp_a/0.025)*rt18(i,j)

      end if
      ENDIF
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! 2. Wall !!!!!!!!!!!!!!!
      !S1.- extrapolation -------------------------!
      IF (ldaps_ws(kk) < inflow_ws_00(kk)) THEN

      if (wall < 0.97) then
      uu2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu22(i,j)
      vv2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv22(i,j)
      ww2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww22(i,j)
      dtemp2(i,j) =  (wd_a*dt21(i,j) + wd_b*dt22(i,j))
      rtemp2(i,j) =  (wd_a*rt21(i,j) + wd_b*rt22(i,j))

      else if (wall > 1.03) then
      uu2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu22(i,j)
      vv2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv22(i,j)
      ww2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww22(i,j)
      dtemp2(i,j) =  (wd_a*dt21(i,j) + wd_b*dt22(i,j))
      rtemp2(i,j) =  (wd_a*rt21(i,j) + wd_b*rt22(i,j))

      else IF (wall >= 0.97 .and. wall <= 1.03) THEN
      uu2(i,j)=  wd_a*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*uu21(i,j) + wd_a*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*uu25(i,j) + wd_b*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*uu22(i,j) + wd_b*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*uu26(i,j)
      vv2(i,j)=  wd_a*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*vv21(i,j) + wd_a*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*vv25(i,j) + wd_b*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*vv22(i,j) + wd_b*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*vv26(i,j)
      ww2(i,j)=  wd_a*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*ww21(i,j) + wd_a*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*ww25(i,j) + wd_b*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*ww22(i,j) + wd_b*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_00(kk)*ww26(i,j)
      dtemp2(i,j)=  wd_a*(wall_temp_b/0.015)*dt21(i,j) + wd_a*(wall_temp_a/0.015)*dt25(i,j) + wd_b*(wall_temp_b/0.015)*dt22(i,j) + wd_b*(wall_temp_a/0.015)*dt26(i,j)
      rtemp2(i,j)=  wd_a*(wall_temp_b/0.015)*rt21(i,j) + wd_a*(wall_temp_a/0.015)*rt25(i,j) + wd_b*(wall_temp_b/0.015)*rt22(i,j) + wd_b*(wall_temp_a/0.015)*rt26(i,j)

      end if
      ENDIF

      IF (ldaps_ws(kk) > inflow_ws_07(kk)) THEN

      if (wall < 0.97) then
      uu2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu22(i,j)
      vv2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv22(i,j)
      ww2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww22(i,j)
      dtemp2(i,j) =  (wd_a*dt21(i,j) + wd_b*dt22(i,j))
      rtemp2(i,j) =  (wd_a*rt21(i,j) + wd_b*rt22(i,j))

      else if (wall > 1.03) then
      uu2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu22(i,j)
      vv2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv22(i,j)
      ww2(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww21(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww22(i,j)
      dtemp2(i,j) =  (wd_a*dt21(i,j) + wd_b*dt22(i,j))
      rtemp2(i,j) =  (wd_a*rt21(i,j) + wd_b*rt22(i,j))

      else IF (wall >= 0.97 .and. wall <= 1.03) THEN
      uu2(i,j)=  wd_a*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*uu21(i,j) + wd_a*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*uu25(i,j) + wd_b*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*uu22(i,j) + wd_b*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*uu26(i,j)
      vv2(i,j)=  wd_a*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*vv21(i,j) + wd_a*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*vv25(i,j) + wd_b*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*vv22(i,j) + wd_b*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*vv26(i,j)
      ww2(i,j)=  wd_a*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*ww21(i,j) + wd_a*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*ww25(i,j) + wd_b*(wall_temp_b/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*ww22(i,j) + wd_b*(wall_temp_a/0.015)*ldaps_ws(kk)/inflow_ws_07(kk)*ww26(i,j)
      dtemp2(i,j)=  wd_a*(wall_temp_b/0.015)*dt21(i,j) + wd_a*(wall_temp_a/0.015)*dt25(i,j) + wd_b*(wall_temp_b/0.015)*dt22(i,j) + wd_b*(wall_temp_a/0.015)*dt26(i,j)
      rtemp2(i,j)=  wd_a*(wall_temp_b/0.015)*rt21(i,j) + wd_a*(wall_temp_a/0.015)*rt25(i,j) + wd_b*(wall_temp_b/0.015)*rt22(i,j) + wd_b*(wall_temp_a/0.015)*rt26(i,j)
      end if
      ENDIF
      !S1_end ---------------------------------------!

      !S2.-  interpolation -------------------------!
      IF (ldaps_ws(kk) >= inflow_ws_00(kk) .and. ldaps_ws(kk) <= inflow_ws_07(kk)) THEN

      if (wall < 0.97) then
      uu2(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu24(i,j)
      vv2(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv24(i,j)
      ww2(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww24(i,j)
      dtemp2(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt24(i,j)
      rtemp2(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt24(i,j)

      else if (wall > 1.03) then
      uu2(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu24(i,j)
      vv2(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv24(i,j)
      ww2(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww24(i,j)
      dtemp2(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt24(i,j)
      rtemp2(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt21(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt23(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt22(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt24(i,j)

      else IF (wall >= 0.97 .and. wall <= 1.03) THEN
      uu2(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*uu21(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*uu22(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*uu23(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*uu24(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*uu25(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*uu26(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*uu27(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*uu28(i,j)
      vv2(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*vv21(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*vv22(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*vv23(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*vv24(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*vv25(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*vv26(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*vv27(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*vv28(i,j)
      ww2(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*ww21(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*ww22(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*ww23(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*ww24(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*ww25(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*ww26(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*ww27(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*ww28(i,j)
      dtemp2(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*dt21(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*dt22(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*dt23(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*dt24(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*dt25(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*dt26(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*dt27(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*dt28(i,j)
      rtemp2(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*rt21(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_b/0.015)*rt22(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*rt23(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_b/0.015)*rt24(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*rt25(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(wall_temp_a/0.015)*rt26(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*rt27(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(wall_temp_a/0.015)*rt28(i,j)

      end if
      ENDIF

!!!!!!!! 3. Road !!!!!!!!!!!!!!!
      !S1.- extrapolation -------------------------!
      IF (ldaps_ws(kk) < inflow_ws_00(kk)) THEN

      if (road < 0.975) then
      uu3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu32(i,j)
      vv3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv32(i,j)
      ww3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww32(i,j)
      dtemp3(i,j) =  (wd_a*dt31(i,j) + wd_b*dt32(i,j))
      rtemp3(i,j) =  (wd_a*rt31(i,j) + wd_b*rt32(i,j))

      else if (road > 1.055) then
      uu3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu32(i,j)
      vv3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv32(i,j)
      ww3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww32(i,j)
      dtemp3(i,j) =  (wd_a*dt31(i,j) + wd_b*dt32(i,j))
      rtemp3(i,j) =  (wd_a*rt31(i,j) + wd_b*rt32(i,j))

      else IF (road >= 0.975 .and. road <= 1.055) THEN
      uu3(i,j)=  wd_a*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu31(i,j) + wd_a*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu35(i,j) + wd_b*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu32(i,j) + wd_b*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu36(i,j)
      vv3(i,j)=  wd_a*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv31(i,j) + wd_a*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv35(i,j) + wd_b*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv32(i,j) + wd_b*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv36(i,j)
      ww3(i,j)=  wd_a*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww31(i,j) + wd_a*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww35(i,j) + wd_b*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww32(i,j) + wd_b*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww36(i,j)
      dtemp3(i,j)=  wd_a*(road_temp_b/0.02)*dt31(i,j) + wd_a*(road_temp_a/0.02)*dt35(i,j) + wd_b*(road_temp_b/0.02)*dt32(i,j) + wd_b*(road_temp_a/0.02)*dt36(i,j)
      rtemp3(i,j)=  wd_a*(road_temp_b/0.02)*rt31(i,j) + wd_a*(road_temp_a/0.02)*rt35(i,j) + wd_b*(road_temp_b/0.02)*rt32(i,j) + wd_b*(road_temp_a/0.02)*rt36(i,j)

      end if
      ENDIF

      IF (ldaps_ws(kk) > inflow_ws_07(kk)) THEN

      if (road < 0.975) then
      uu3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu32(i,j)
      vv3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv32(i,j)
      ww3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww32(i,j)
      dtemp3(i,j) =  (wd_a*dt31(i,j) + wd_b*dt32(i,j))
      rtemp3(i,j) =  (wd_a*rt31(i,j) + wd_b*rt32(i,j))

      else if (road > 1.055) then
      uu3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu32(i,j)
      vv3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv32(i,j)
      ww3(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww31(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww32(i,j)
      dtemp3(i,j) =  (wd_a*dt31(i,j) + wd_b*dt32(i,j))
      rtemp3(i,j) =  (wd_a*rt31(i,j) + wd_b*rt32(i,j))

      else IF (road >= 0.975 .and. road <= 1.055) THEN
      uu3(i,j)=  wd_a*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu31(i,j) + wd_a*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu35(i,j) + wd_b*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu32(i,j) + wd_b*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu36(i,j)
      vv3(i,j)=  wd_a*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv31(i,j) + wd_a*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv35(i,j) + wd_b*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv32(i,j) + wd_b*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv36(i,j)
      ww3(i,j)=  wd_a*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww31(i,j) + wd_a*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww35(i,j) + wd_b*(road_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww32(i,j) + wd_b*(road_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww36(i,j)
      dtemp3(i,j)=  wd_a*(road_temp_b/0.02)*dt31(i,j) + wd_a*(road_temp_a/0.02)*dt35(i,j) + wd_b*(road_temp_b/0.02)*dt32(i,j) + wd_b*(road_temp_a/0.02)*dt36(i,j)
      rtemp3(i,j)=  wd_a*(road_temp_b/0.02)*rt31(i,j) + wd_a*(road_temp_a/0.02)*rt35(i,j) + wd_b*(road_temp_b/0.02)*rt32(i,j) + wd_b*(road_temp_a/0.02)*rt36(i,j)
      end if
      ENDIF
      !S1_end ---------------------------------------!

      !S2.-  interpolation -------------------------!
      IF (ldaps_ws(kk) >= inflow_ws_00(kk) .and. ldaps_ws(kk) <= inflow_ws_07(kk)) THEN

      if (road < 0.975) then
      uu3(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu34(i,j)
      vv3(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv34(i,j)
      ww3(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww34(i,j)
      dtemp3(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt34(i,j)
      rtemp3(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt34(i,j)

      else if (road > 1.055) then
      uu3(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu34(i,j)
      vv3(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv34(i,j)
      ww3(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww34(i,j)
      dtemp3(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt34(i,j)
      rtemp3(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt31(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt33(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt32(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt34(i,j)

      else IF (road >= 0.975 .and. road <= 1.055) THEN
      uu3(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*uu31(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*uu32(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*uu33(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*uu34(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*uu35(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*uu36(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*uu37(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*uu38(i,j)
      vv3(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*vv31(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*vv32(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*vv33(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*vv34(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*vv35(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*vv36(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*vv37(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*vv38(i,j)
      ww3(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*ww31(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*ww32(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*ww33(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*ww34(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*ww35(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*ww36(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*ww37(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*ww38(i,j)
      dtemp3(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*dt31(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*dt32(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*dt33(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*dt34(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*dt35(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*dt36(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*dt37(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*dt38(i,j)
      rtemp3(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*rt31(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_b/0.02)*rt32(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*rt33(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_b/0.02)*rt34(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*rt35(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(road_temp_a/0.02)*rt36(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*rt37(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(road_temp_a/0.02)*rt38(i,j)

      end if
      ENDIF

!!!!!!!! 4. Green !!!!!!!!!!!!!!!
      !S1.- extrapolation -------------------------!
      IF (ldaps_ws(kk) < inflow_ws_00(kk)) THEN

      if (green < 0.95) then
      uu4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu42(i,j)
      vv4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv42(i,j)
      ww4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww42(i,j)
      dtemp4(i,j) =  (wd_a*dt41(i,j) + wd_b*dt42(i,j))
      rtemp4(i,j) =  (wd_a*rt41(i,j) + wd_b*rt42(i,j))

      else if (green > 1.03) then
      uu4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu42(i,j)
      vv4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv42(i,j)
      ww4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww42(i,j)
      dtemp4(i,j) =  (wd_a*dt41(i,j) + wd_b*dt42(i,j))
      rtemp4(i,j) =  (wd_a*rt41(i,j) + wd_b*rt42(i,j))

      else IF (green >= 0.95 .and. green <= 1.03) THEN
      uu4(i,j)=  wd_a*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu41(i,j) + wd_a*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu45(i,j) + wd_b*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu42(i,j) + wd_b*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu46(i,j)
      vv4(i,j)=  wd_a*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv41(i,j) + wd_a*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv45(i,j) + wd_b*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv42(i,j) + wd_b*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv46(i,j)
      ww4(i,j)=  wd_a*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww41(i,j) + wd_a*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww45(i,j) + wd_b*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww42(i,j) + wd_b*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww46(i,j)
      dtemp4(i,j)=  wd_a*(green_temp_b/0.02)*dt41(i,j) + wd_a*(green_temp_a/0.02)*dt45(i,j) + wd_b*(green_temp_b/0.02)*dt42(i,j) + wd_b*(green_temp_a/0.02)*dt46(i,j)
      rtemp4(i,j)=  wd_a*(green_temp_b/0.02)*rt41(i,j) + wd_a*(green_temp_a/0.02)*rt45(i,j) + wd_b*(green_temp_b/0.02)*rt42(i,j) + wd_b*(green_temp_a/0.02)*rt46(i,j)

      end if
      ENDIF

      IF (ldaps_ws(kk) > inflow_ws_07(kk)) THEN

      if (green < 0.95) then
      uu4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu42(i,j)
      vv4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv42(i,j)
      ww4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww42(i,j)
      dtemp4(i,j) =  (wd_a*dt41(i,j) + wd_b*dt42(i,j))
      rtemp4(i,j) =  (wd_a*rt41(i,j) + wd_b*rt42(i,j))

      else if (green > 1.03) then
      uu4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu42(i,j)
      vv4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv42(i,j)
      ww4(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww41(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww42(i,j)
      dtemp4(i,j) =  (wd_a*dt41(i,j) + wd_b*dt42(i,j))
      rtemp4(i,j) =  (wd_a*rt41(i,j) + wd_b*rt42(i,j))

      else IF (green >= 0.95 .and. green <= 1.03) THEN
      uu4(i,j)=  wd_a*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu41(i,j) + wd_a*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu45(i,j) + wd_b*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu42(i,j) + wd_b*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu46(i,j)
      vv4(i,j)=  wd_a*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv41(i,j) + wd_a*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv45(i,j) + wd_b*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv42(i,j) + wd_b*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv46(i,j)
      ww4(i,j)=  wd_a*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww41(i,j) + wd_a*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww45(i,j) + wd_b*(green_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww42(i,j) + wd_b*(green_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww46(i,j)
      dtemp4(i,j)=  wd_a*(green_temp_b/0.02)*dt41(i,j) + wd_a*(green_temp_a/0.02)*dt45(i,j) + wd_b*(green_temp_b/0.02)*dt42(i,j) + wd_b*(green_temp_a/0.02)*dt46(i,j) 
      rtemp4(i,j)=  wd_a*(green_temp_b/0.02)*rt41(i,j) + wd_a*(green_temp_a/0.02)*rt45(i,j) + wd_b*(green_temp_b/0.02)*rt42(i,j) + wd_b*(green_temp_a/0.02)*rt46(i,j)
      end if
      ENDIF
      !S1_end ---------------------------------------!

      !S2.-  interpolation -------------------------!
      IF (ldaps_ws(kk) >= inflow_ws_00(kk) .and. ldaps_ws(kk) <= inflow_ws_07(kk)) THEN

      if (green < 0.95) then
      uu4(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu44(i,j)
      vv4(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv44(i,j)
      ww4(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww44(i,j)
      dtemp4(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt44(i,j)
      rtemp4(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt44(i,j)

      else if (green > 1.03) then
      uu4(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu44(i,j)
      vv4(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv44(i,j)
      ww4(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww44(i,j)
      dtemp4(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt44(i,j)
      rtemp4(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt41(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt43(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt42(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt44(i,j)

      else IF (green >= 0.95 .and. green <= 1.03) THEN
      uu4(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*uu41(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*uu42(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*uu43(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*uu44(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*uu45(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*uu46(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*uu47(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*uu48(i,j)
      vv4(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*vv41(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*vv42(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*vv43(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*vv44(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*vv45(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*vv46(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*vv47(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*uuvv(i,j)
      ww4(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*ww41(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*ww42(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*ww43(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*ww44(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*ww45(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*ww46(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*ww47(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*ww48(i,j)
      dtemp4(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*dt41(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*dt42(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*dt43(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*dt44(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*dt45(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*dt46(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*dt47(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*dt48(i,j)
      rtemp4(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*rt41(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_b/0.02)*rt42(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*rt43(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_b/0.02)*rt44(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*rt45(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(green_temp_a/0.02)*rt46(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*rt47(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(green_temp_a/0.02)*rt48(i,j)

      end if
      ENDIF

!!!!!!!! 5. Soil !!!!!!!!!!!!!!!
      !S1.- extrapolation -------------------------!
      IF (ldaps_ws(kk) < inflow_ws_00(kk)) THEN

      if (soil < 0.96) then
      uu5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu52(i,j)
      vv5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv52(i,j)
      ww5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww52(i,j)
      dtemp5(i,j) =  (wd_a*dt51(i,j) + wd_b*dt52(i,j))
      rtemp5(i,j) =  (wd_a*rt51(i,j) + wd_b*rt52(i,j))

      else if (soil > 1.04) then
      uu5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*uu51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*uu52(i,j)
      vv5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*vv51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*vv52(i,j)
      ww5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_00(kk)*ww51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_00(kk)*ww52(i,j)
      dtemp5(i,j) =  (wd_a*dt51(i,j) + wd_b*dt52(i,j))
      rtemp5(i,j) =  (wd_a*rt51(i,j) + wd_b*rt52(i,j))

      else IF (soil >= 0.96 .and. soil <= 1.04) THEN
      uu5(i,j)=  wd_a*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu51(i,j) + wd_a*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu55(i,j) + wd_b*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu52(i,j) + wd_b*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*uu56(i,j)
      vv5(i,j)=  wd_a*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv51(i,j) + wd_a*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv55(i,j) + wd_b*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv52(i,j) + wd_b*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*vv56(i,j)
      ww5(i,j)=  wd_a*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww51(i,j) + wd_a*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww55(i,j) + wd_b*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww52(i,j) + wd_b*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_00(kk)*ww56(i,j)
      dtemp5(i,j)=  wd_a*(soil_temp_b/0.02)*dt51(i,j) + wd_a*(soil_temp_a/0.02)*dt55(i,j) + wd_b*(soil_temp_b/0.02)*dt52(i,j) + wd_b*(soil_temp_a/0.02)*dt56(i,j)
      rtemp5(i,j)=  wd_a*(soil_temp_b/0.02)*rt51(i,j) + wd_a*(soil_temp_a/0.02)*rt55(i,j) + wd_b*(soil_temp_b/0.02)*rt52(i,j) + wd_b*(soil_temp_a/0.02)*rt56(i,j)

      end if
      ENDIF

      IF (ldaps_ws(kk) > inflow_ws_07(kk)) THEN

      if (soil < 0.96) then
      uu5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu52(i,j)
      vv5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv52(i,j)
      ww5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww52(i,j)
      dtemp5(i,j) =  (wd_a*dt51(i,j) + wd_b*dt52(i,j))
      rtemp5(i,j) =  (wd_a*rt51(i,j) + wd_b*rt52(i,j))

      else if (soil > 1.04) then
      uu5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*uu51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*uu52(i,j)
      vv5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*vv51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*vv52(i,j)
      ww5(i,j)= wd_a*ldaps_ws(kk)/inflow_ws_07(kk)*ww51(i,j) + wd_b*ldaps_ws(kk)/inflow_ws_07(kk)*ww52(i,j)
      dtemp5(i,j) =  (wd_a*dt51(i,j) + wd_b*dt52(i,j))
      rtemp5(i,j) =  (wd_a*rt51(i,j) + wd_b*rt52(i,j))

      else IF (soil >= 0.96 .and. soil <= 1.04) THEN
      uu5(i,j)=  wd_a*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu51(i,j) + wd_a*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu55(i,j) + wd_b*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu52(i,j) + wd_b*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*uu56(i,j)
      vv5(i,j)=  wd_a*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv51(i,j) + wd_a*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv55(i,j) + wd_b*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv52(i,j) + wd_b*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*vv56(i,j)
      ww5(i,j)=  wd_a*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww51(i,j) + wd_a*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww55(i,j) + wd_b*(soil_temp_b/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww52(i,j) + wd_b*(soil_temp_a/0.02)*ldaps_ws(kk)/inflow_ws_07(kk)*ww56(i,j)
      dtemp5(i,j)=  wd_a*(soil_temp_b/0.02)*dt51(i,j) + wd_a*(soil_temp_a/0.02)*dt55(i,j) + wd_b*(soil_temp_b/0.02)*dt52(i,j) + wd_b*(soil_temp_a/0.02)*dt56(i,j)
      rtemp5(i,j)=  wd_a*(soil_temp_b/0.02)*rt51(i,j) + wd_a*(soil_temp_a/0.02)*rt55(i,j) + wd_b*(soil_temp_b/0.02)*rt52(i,j) + wd_b*(soil_temp_a/0.02)*rt56(i,j)
      end if
      ENDIF
      !S1_end ---------------------------------------!

      !S2.-  interpolation -------------------------!
      IF (ldaps_ws(kk) >= inflow_ws_00(kk) .and. ldaps_ws(kk) <= inflow_ws_07(kk)) THEN

      if (soil < 0.96) then
      uu5(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu54(i,j)
      vv5(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv54(i,j)
      ww5(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww54(i,j)
      dtemp5(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt54(i,j)
      rtemp5(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt54(i,j)

      else if (soil > 1.04) then
      uu5(i,j)= wd_a*(ws_b/(ws_a+ws_b))*uu51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*uu53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*uu52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*uu54(i,j)
      vv5(i,j)= wd_a*(ws_b/(ws_a+ws_b))*vv51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*vv53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*vv52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*vv54(i,j)
      ww5(i,j)= wd_a*(ws_b/(ws_a+ws_b))*ww51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*ww53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*ww52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*ww54(i,j)
      dtemp5(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*dt51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*dt53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*dt52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*dt54(i,j)
      rtemp5(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*rt51(i,j) + wd_a*(ws_a/(ws_a+ws_b))*rt53(i,j) + wd_b*(ws_b/(ws_a+ws_b))*rt52(i,j) + wd_b*(ws_a/(ws_a+ws_b))*rt54(i,j)

      else IF (soil >= 0.96 .and. soil <= 1.04) THEN
      uu5(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*uu51(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*uu52(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*uu53(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*uu54(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*uu55(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*uu56(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*uu57(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*uu58(i,j)
      vv5(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*vv51(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*vv52(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*vv53(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*vv54(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*vv55(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*vv56(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*vv57(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*vv58(i,j)
      ww5(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*ww51(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*ww52(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*ww53(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*ww54(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*ww55(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*ww56(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*ww57(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*ww58(i,j)
      dtemp5(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*dt51(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*dt52(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*dt53(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*dt54(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*dt55(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*dt56(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*dt57(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*dt58(i,j)
      rtemp5(i,j)=  wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*rt51(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_b/0.02)*rt52(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*rt53(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_b/0.02)*rt54(i,j) + wd_a*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*rt55(i,j) + wd_b*(ws_b/(ws_a+ws_b))*(soil_temp_a/0.02)*rt56(i,j) + wd_a*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*rt57(i,j) + wd_b*(ws_a/(ws_a+ws_b))*(soil_temp_a/0.02)*rt58(i,j)

      end if
      ENDIF

!c#################### Lapse rate ##################################
      !S1.- extrapolation -------------------------!
      IF (dt(kk) < -20.) then
       if (ldaps_ws(kk) .lt. inflow_ws_01(kk)) Then
      ltemp =  rt1
      else if (ldaps_ws(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws(kk) .le. inflow_ws_07(kk)) Then
      ltemp = (ws_b*rt1 + ws_a*rt3)/(ws_a+ws_b)
      else
      ltemp =  rt1
       end if ! inflow

      ELSE IF (dt(kk) > 20) then
       if (ldaps_ws(kk) .lt. inflow_ws_01(kk)) Then
      ltemp =  rt1
      else if (ldaps_ws(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws(kk) .le. inflow_ws_07(kk)) Then
      ltemp = (ws_b*rt1 + ws_a*rt3)/(ws_a+ws_b)
      else
      ltemp =  rt1
       end if ! inflow

      ELSE IF (dt(kk) >= -20. .and. dt(kk) <= 20.) THEN
       if (ldaps_ws(kk) .lt. inflow_ws_01(kk)) Then
      ltemp =  (temp_b*rt1 + temp_a*rt2)/(temp_a + temp_b)
       else if (ldaps_ws(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws(kk) .le. inflow_ws_07(kk)) Then
      ltemp =  (ws_b*(temp_b*rt1 + temp_a*rt2)/(temp_a + temp_b) + ws_a*(temp_b*rt3 + temp_a*rt4)/(temp_a + temp_b))/(ws_a + ws_b)
       else
      ltemp =  (temp_b*rt1 + temp_a*rt2)/(temp_a + temp_b)
      end if !inflow
      
      END IF  !dt

      !windspeed change ratio---------------------!
      IF(dwd .eq. 0 .or. dwd .eq. 180.) then
      IF (dws(kk) < -10.) then

       if (ldaps_ws(kk) .lt. inflow_ws_01(kk)) Then
      rwind =  rws11
      else if (ldaps_ws(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws(kk) .le. inflow_ws_07(kk)) Then
      rwind = (ws_b*rws11 + ws_a*rws13)/(ws_a+ws_b)
      else
      rwind =  rws11
       end if ! inflow

      else if (dws(kk) > 30) then
       if (ldaps_ws(kk) .lt. inflow_ws_01(kk)) Then
      rwind =  rws11
      else if (ldaps_ws(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws(kk) .le. inflow_ws_07(kk)) Then
      rwind = (ws_b*rws11 + ws_a*rws13)/(ws_a+ws_b)
      else
      rwind =  rws11
       end if ! inflow

      else if (dws(kk) >= -10. .and. dws(kk) <= 30.) THEN
       if (ldaps_ws(kk) .lt. inflow_ws_01(kk)) Then
      rwind =  (dws_b*rws11 + dws_a*rws12)/(dws_a + dws_b)
       else if (ldaps_ws(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws(kk) .le. inflow_ws_07(kk)) Then
      rwind =  (ws_b*(dws_b*rws11 + dws_a*rws12)/(dws_a + dws_b) + ws_a*(dws_b*rws13 + dws_a*rws14)/(dws_a + dws_b))/(ws_a + ws_b)
       else
      rwind =  (dws_b*rws11 + dws_a*rws12)/(dws_a + dws_b)
      end if !inflow

      ENDIF  !dws

      ELSE IF (dwd .gt. 0. .and. dwd .lt. 180.) then
      IF (dws(kk) < -10.) then

      if (ldaps_ws(kk) .lt. inflow_ws_01(kk)) Then
      rwind = (rws11*dwd_b + rws15*dwd_a)/(dwd_a + dwd_b)
      else if (ldaps_ws(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws(kk) .le. inflow_ws_07(kk)) Then
      rwind = (ws_b*(dwd_b*rws11 + dwd_a*rws13)/(dwd_a + dwd_b) + ws_a*(dwd_b*rws15 + dwd_a*rws17)/(dwd_a + dwd_b))/(ws_a + ws_b)
      else
      rwind = (rws11*dwd_b + rws15*dwd_a)/(dwd_a + dwd_b)
       end if ! inflow

      else if (dws(kk) > 30) then
       if (ldaps_ws(kk) .lt. inflow_ws_01(kk)) Then
      rwind = (rws11*dwd_b + rws15*dwd_a)/(dwd_a + dwd_b)
      else if (ldaps_ws(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws(kk) .le. inflow_ws_07(kk)) Then
      rwind = (ws_b*(dwd_b*rws11 + dwd_a*rws13)/(dwd_a + dwd_b) + ws_a*(dwd_b*rws15 + dwd_a*rws17)/(dwd_a + dwd_b))/(ws_a + ws_b)
      else
      rwind = (rws11*dwd_b + rws15*dwd_a)/(dwd_a + dwd_b)
       end if ! inflow

      else if (dws(kk) >= -10. .and. dws(kk) <= 30.) THEN
       if (ldaps_ws(kk) .lt. inflow_ws_01(kk)) Then
      rwind =  (dwd_b*(dws_b*rws11 + dws_a*rws12)/(dws_a + dws_b) + dwd_a*(dws_b*rws15 + dws_a*rws16)/(dws_a + dws_b))/(dwd_a + dwd_b)
       else if (ldaps_ws(kk) .ge. inflow_ws_01(kk) .and. ldaps_ws(kk) .le. inflow_ws_07(kk)) Then
      rwind = (dwd_b*((ws_b*(dws_b*rws11 + dws_a*rws12)/(dws_a + dws_b)) + (ws_a*(dws_b*rws13 + dws_a*rws14)/(dws_a + dws_b)))/(ws_a + ws_b) + dwd_a*((ws_b*(dws_b*rws15 + dws_a*rws16)/(dws_a + dws_b)) + (ws_a*(dws_b*rws17 + dws_a*rws18)/(dws_a + dws_b)))/(ws_a+ws_b))/(dwd_a+dwd_b)
       else
      rwind =  (dwd_b*(dws_b*rws11 + dws_a*rws12)/(dws_a + dws_b) + dwd_a*(dws_b*rws15 + dws_a*rws16)/(dws_a + dws_b))/(dwd_a + dwd_b)
      end if !inflow

      ENDIF  !dws
      END IF !dwd

!c############################################################################

      uu(i,j) = (uu1(i,j) + uu2(i,j) + uu3(i,j) + uu4(i,j) + uu5(i,j))/5.
      vv(i,j) = (vv1(i,j) + vv2(i,j) + vv3(i,j) + vv4(i,j) + vv5(i,j))/5.
      ww(i,j) = (ww1(i,j) + ww2(i,j) + ww3(i,j) + ww4(i,j) + ww5(i,j))/5.
      dtemp(i,j) = ((ldaps_t(kk) + 273.15) * (rtemp1(i,j) + rtemp2(i,j) + rtemp3(i,j) + rtemp4(i,j) + rtemp5(i,j) - 4. )) - 273.15
      ctemp(i,j) = (dtemp(i,j) + 273.15)*ltemp - 273.15    

      if (kk .gt. 20) then !! BL wind 
      fuu(i,j) = uu(i,j)*rwind
      fvv(i,j) = vv(i,j)*rwind
      fww(i,j) = ww(i,j)*rwind
      else  !!RL wind
      fuu(i,j) = uu(i,j)
      fvv(i,j) = vv(i,j)
      fww(i,j) = ww(i,j)
      end if
     
      ENDDO
      ENDDO

      do 350 j= 1,jmax2
      do 350 i= 1,imax2

       write(60, 365) fuu(i,j), fvv(i,j), fww(i,j), dtemp(i,j), ctemp(i,j) 

  350 continue


      CLOSE(5)
      CLOSE(6)
      CLOSE(7)
 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
     
      ENDDO  !! Vertical grid

  365 format(5(1x,f12.6))
  366 format(8(1x,f12.6))
!c############################################################################
      stop
      END PROGRAM





