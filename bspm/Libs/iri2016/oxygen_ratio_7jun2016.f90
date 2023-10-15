subroutine oxygen_ratio_7jun2016( height, lat, lon, year, mmdd, hour, res)
implicit none
 real*8 height, lat, lon, hour, dhour, res
 integer*8 year, mmdd, jmag, i

!     month and day: mmdd
!     height in km
!     lat in degree

  real :: scid
  real, dimension(20, 1000) :: outf
  real, dimension (100, 1000) :: oar
  logical, dimension (50) :: jf

   call read_ig_rz
   call readapf107
        
   do i=1,100
    oar(i,1)=-1.0
   enddo

  do i=1,50
    jf(i)=.true.
  enddo

  jf(1)=.false.                  ! no Ne
  jf(2)=.false.                  ! no temperature
  jf(3)=.true.                  ! ion composition
  jf(4)=.false.        ! t=B0table f=other models (f)
  jf(5)=.false.               ! URSI foF2 model
  jf(6)=.false.               ! RBV10+TTS03
  jf(12)=.false.              ! no konsol output
  jf(23)=.false.              ! TTS Te model is standard
  jf(28)=.false.	  ! f=spread-F not computed (f)
  jf(29)=.false.              ! New Topside options
  jf(30)=.false.              ! NeQuick topside
  jf(33)=.false. 	  ! f=auroral boundary off (f)
  jf(34)=.false. 	  ! t=messages on 
  jf(35)=.false. 	  ! f=auroral E-storm model off
!          jf(36)=.false. 	  ! t=hmF2 w/out foF2_storm f=with
!          jf(37)=.false. 	  ! t=topside w/out foF2_storm f=with
!          jf(38)=.false. 	  ! t=WRITEs off in IRIFLIP f=on 
  jf(39)=.false. 	  ! new hmF2 models 
!          jf(40)=.false. 	  ! t=AMTB-model, f=Shubin-COSMIC model 
!          jf(41)=.false. 	  ! COV=f(IG12) (IRI before Oct 2015) 
!          jf(42)=.false. 	  ! Te w/o PF10.7 dependance 

jmag=1
dhour=hour+25.

call IRI_SUB(jf, jmag, lat, lon, year, mmdd, dhour, height, height, 1., outf, oar)
      
        scid=1.0E-8
        if(jf(22)) scid=10.
!        jio=INT(OUTF(5,li)*scid+.5)

res=int(outf(5,1)*scid+.5)

end


