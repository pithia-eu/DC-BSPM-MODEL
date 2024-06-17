def kp(IP,TU,iedi):
#c   INPUT:
#c   TU   - minimal time
#c   J   - day number
#c
#c   OUTPUT:
#c   AKP   - Kp
#c   APP   - smoothed Kp
   
#      DO 44 I=1,496
#         READ(12,*)M(I),H(I),IP(I),TP(I),TF(I)
#         if (IP(I).eq.0.) then
#             IP(I) = 0.01
#         endif
#         IF (M(I).LE.0) GOTO 5
#  44  CONTINUE
      
    i=int(TU/iedi) # int or not ????
    if (IP[i]==0.): IP[i] = 0.01
    AP=IP[i]
    APP=AP
    TP=1.0
    TF=2.0
    
    if i!=0:
        TU3  = int(TU/3)*3
        TU15 = TU3+TP
        T3U3 = TU3+3
        TUF  = TU3+TF
    
        IM1=i-1
        IP1=i+1
        AP=IP[i]
        
#        print(IP[i],IP[IM1])
    
        if(TU >= TU15 and TU <= T3U3):    # GO TO 2
            if TU <= TUF: AP=IP[i] 
            else: 
#                if IP[i] >= IP[IP1]: # I removed = to obtain the mast smooth !!!
                if IP[i] > IP[IP1]: 
                    AP=IP[i]-(TU-TUF)*(IP[i]-IP[IP1])/(3-TP)
        else:               
            if IP[i] >= IP[IM1]: AP=float(IP[IM1])+(TU-TU3)*(float(IP[i])-float(IP[IM1]))
    
    AKP=AP
    return AKP,APP
