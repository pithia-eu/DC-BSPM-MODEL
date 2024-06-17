def ed_input(file_input):
    
    import numpy as np 
    from b_m2_subkE import RMAGP, RGG
    
#C   program to create an input file for pls_kz program
#c   read in namelist UT0 P0 R0 D0 B0 Dext0 NEP IDZ
#c
#c   create and print new lines at UTmax - DUT ...
#c   until UT0
#c   Write last line -2e22 ...
#c
#c   uto    initial UT
#c   R0    initial R
#c   P0   initial LT
#c   D   intial density...
#c   NEP   first id-number of plasma element .. to be incremented
#c      NEP = 1 implies that the first line will be repeated
#c      NEP =0    >>> not repeated in Reprise.dat
#c   IDZ = 0  means that is with interchange taken into account
#c       = 1  means no interchange  ... MHD
#c   UTmax  max UT to be written
#c   DUT     increment of UT
#c   DELTMA  = increment for next UTmax for following runs

    ut0 = 0.
    Magnetopause = 0
    BMINZ = 8.
    P0min = 18.
    P0max = 30.
    mm1= 4
    P0 = 23.
    R0 = 4.
    R1 = 2.75
    DR1 =  -0.5
    D = 200.0
    DExt = 600.0
    NEP = 1
    IDZ = 0
    
    utmax = 24.000
    dut = 0.125
    DELTMA = 577
    
    INF = -2e21
    
  
    NEPZ = NEP
    RR = R0
    PP = P0
    FKP = 1.

#    inp=[]

#    file_input = 'input.dat'

    with open(file_input, 'w') as output_file:

        while RR >= R1 :  
            t = utmax + dut  # line 71
            NEP=NEPZ
            while  (t >= ut0 + 0.001):
                t = t - dut
#                inp.append([t]+[PP]+[RR]+[D]+[DELTMA]+[DExt]+[NEP]+[IDZ]+[FKP])
                output_file.write(f"{t:16.8e}{PP:16.8e}{RR:16.8e}{D:16.8e}{DELTMA:16.8e}{DExt:16.8e}{NEP:5d}{IDZ:5d}{FKP:6.2f}\n")
                NEP = NEP +1
        
#            inp.append([INF]*6+[NEPZ]+[IDZ]+[FKP])
            output_file.write(6*(f"{INF:16.8e}")+f"{NEPZ:5d}{IDZ:5d}{FKP:6.2f}\n")
            RR=RR+DR1
    
        if (Magnetopause == 0): # GOTO 89
#            inp.append([-INF]*6+[1]+[IDZ]+[FKP])  # Line 89
            output_file.write(6*(f"{INF:16.8e}")+f"{NEPZ:5d}{IDZ:5d}{FKP:6.2f}\n")
            print ('edith',t,PP,RR, NEP,mm1)
#            return inp
            return
    
        PP=P0min
    
        while (PP <= P0max): 
            t = ut0 - dut
            NEP = NEPZ
            CS = np.cos(PP*0.2617993878)
            RR1 = RMAGP(CS,1)            # check functions
            RR2 = RGG(BMINZ,CS)     # check functions
            RR = np.min(RR1,RR2)
    
            while  (t < utmax-0.001):
                t = t + dut
#                inp.append([t]+[PP]+[RR]+[D]+[DELTMA]+[DExt]+[NEP]+[IDZ]+[FKP])
                output_file.write(f"{t:16.8e}{PP:16.8e}{RR:16.8e}{D:16.8e}{DELTMA:16.8e}{DExt:16.8e}{NEP:5d}{IDZ:5d}{FKP:6.2f}\n")
                NEP = NEP +1            
    #        inp.append([INF]*6+[1]+[IDZ]+[FKP])               
            output_file.write(6*(f"{INF:16.8e}")+f"{NEPZ:5d}{IDZ:5d}{FKP:6.2f}\n")
    
            PP = PP + 0.25*mm1

#    return inp
    return

#import pandas as pd
#inp=ed_input()
# year=2015; month=3; day=17
# str_date = "-".join([str(year),str(month), str(day)])
# ed_input('input_'+str_date+'.dat')
#year=2015; month=3; day=17
#str_date = "-".join([str(year),str(month), str(day)])
#df=pd.DataFrame(inp)
#df.to_csv('input_'+str_date+'.dat', sep=' ',float_format='%.4f',header=False,index=False)  
