subroutine proton_temp( height, lat, lon, year, mmdd, hour, res)
implicit none
 real*8 height, lat, lon, hour, dhour, res
 integer*8 year, mmdd, jmag, i

!     month and day: mmdd
!     height in km
!     lat in degree

  real, dimension(20, 500) :: outf
  real, dimension (50) :: oar
  logical, dimension (30) :: jf

integer icalls, nmono, iyearo, idaynro, igino
real rzino, ut0
  COMMON/const2/icalls,nmono,iyearo,idaynro,rzino,igino,ut0

  icalls=0
  nmono=-1
  iyearo=-1
  idaynro=-1
  rzino=-1
  igino=-1
  ut0=-1

  do i=1,50
      oar(i)=-1.0
  enddo

do i=1,30
    jf(i)=.true.
enddo

  jf(1)=.false.                 ! no calculation Ne
  jf(3)=.false.                  ! no ion composition
  jf(5)=.false.               ! URSI foF2 model
  jf(6)=.false.               ! Newest ion composition model
  jf(12)=.false.              ! no konsol output
  jf(23)=.false.              ! TTS Te model is standard
  jf(29)=.false.              ! New Topside options
  jf(30)=.false.              ! NeQuick topside

jmag=1
dhour=hour+25.

call IRI_SUB(jf, jmag, lat, lon, year, mmdd, dhour, height, height, 1., outf, oar)

res = int(outf(3,1)+.5)

end


