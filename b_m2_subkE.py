#C SUBROUTINES AND FUNCTION LIBRARY B-M2-SUB.FOR
#C--------------------------------------------------------------------------
#C SUBROUTINE IN THIS FILE:
#C------------------------
#C   BGRAD
#C FUNCTION IN THIS FILE:
#C------------------------
#C   RGG
#C   BRG
#C   RMAGP
#C--------------------------------------------------------------------------

import numpy as np

CM1=6.
CM2=-24.
CM3=18.
CM4=1728.
CM5=31000.


def RGG(B,CS):
#C      MCILWAIN     UCSD
#C       USES B=6-24*CS+18*CS*CS/(1+(12/R)**3)+31000/R**3
#C      WHERE B=MAGNETIC FIELD IN GAMMAS AND CS=COS(LOCAL TIME)
#C       I.E. MAGNETIC FIELD MODEL M2
#C      TO COMPUTE RGG = THE GEOGRAPHIC RADIAL DISTANCE IN THE MAGNETIC
#C      EQUATORIAL PLANE IN UNITS OF EARTH RADII

#        COMMON /BCOEF/CM1,CM2,CM3,CM4,CM5
#c        DATA  CM1/6./,CM2/-24./,CM3/18./,CM4/1728./,CM5/31000./

    A=CM3*CS*CS
    G=B-CM1-CM2*CS
    H=G-A
    if (H > 0.0) : 
        AA=(CM5/CM4 +G)**2 -4.*CM5/CM4 *A  #GO TO 2
        if (AA < 0.0):
            RGGe=999.0  #GO TO 1
            R3=RGGe**3
        else: 
            R3=0.5*CM4 *(CM5/CM4 -G +np.sqrt(AA))/H
            RGGe=R3**(1./3.)
    
    else:
        RGGe=999.0  #GO TO 1
        R3=RGGe**3
        
    return RGGe


def BRG(R,CS):
#C          MAGNETIC FIELD MODEL M2
#C       R=RADIAL DISTANCE AT MAGNETIC EQUATOR IN EARTH RADII
#C       CS=COS( LOCAL TIME )

#      COMMON /BCOEF/CM1,CM2,CM3,CM4,CM5
    R3=R*R*R
    BRG = CM1+CM2*CS+CM3*CS*CS/(1.0+CM4/R3)+CM5/R3

    return BRG


def RMAGP(CS,FKP):

    RMP1=.028
    RMP2=-.35
    RMP3=-.58
    RMP4=17.87
    RMP5=-233.7
    RMP6=1.057
    RMP7=0.0526

    GS=-CS                  # Noon is zero in orig. coord.
    SN=np.sqrt(1.-GS*GS)
    AA=0.5*(RMP3*SN+RMP4*GS)/(SN*SN+RMP1*SN*GS+RMP2*GS*GS)
    BB=AA*AA-RMP5/(SN*SN+RMP1*SN*GS+RMP2*GS*GS)

    if(BB < 0.0):
        RMAGPe=999.
        return
    else:
        RMAGPe=(-AA+np.sqrt(BB))*(RMP6-RMP7*FKP/(1.+.1*FKP))

    return RMAGPe  # Check!!!


def BGRAD(R,BDR,BDP,CS,SS):
#C        GRADIENTS OF MAG. FIELD MODEL M2
#C         BDR=DB/DR IN GAMMA/EARTH RADII
#C         BDP=DB/DP IN GAMMA/RADIAN

#      COMMON/BCOEF/CM1,CM2,CM3,CM4,CM5
    R2 = R*R
    R3 = R2*R
    R4 = R2*R2
    BDR = -3.*CM5/R4+3.*CM3*CM4*CS*CS/(R2+CM4/R)**2
    BDP = -CM2*SS-CM3*2.*CS*SS/(1.+CM4/R3)*0.26179938
    
    return BDR,BDP  # Check!!!!
