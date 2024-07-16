def fkp_main(dire,fin,fout):
    import pandas as pd
    from kpp import kp
#c      main prog. test the kp
#c
#c   Calculates  Kp index
#c   from a series of real Kp-values
#c   (stored in column 3 of Kp.dat)
#c   The values of tp and tf can be changed
#c   from 0 to 3 hours (tu < 3) to allow 
#c   smoother increase/decrease of fKp
#c
#c   INPUT:
#c   tui    - initial universal time
#c   tuf    - final universal time
#c   m   - day number in a month
#c   H   - hour in day
#c   Ip   - Kp value*10
#c   Tp   - upward time lag       
#c   Tf   - downward time lag       
#c
#c   OUTPUT:
#c   FKP   - smooth Kp   
#c   Ap   - Ap index
#c   tud   - universal time in hours

#    str_date = "-".join([str(year),str(month), str(day)])
#    fn="kp" + str(day) + "-" + str(month) + "-" + str(year) + ".dat"

    INF = -2e+21
    dtu = 0.125
    jf = 10
    
    IFLAG = -1
    tu = -dtu
    tufin = 24*jf  # 92
    
    H=[];IP=[]
    kppe = pd.read_csv(dire+fin,header=None)
        
    for i in range(kppe.shape[0]):
        a=' '.join(list(kppe[0])[i].split(' ')).split()
        H.append(a[1]) ; IP.append(a[2])
    
    IP=[float(i) for i in IP] ; H=[int(i) for i in H]
    
    if (int(H[1])-int(H[0])) == 3 : iedi=3
    else: iedi=1
    
    fkpp=[]
    while tu < tufin:
        try:
            tu=tu+dtu
            tud = 1 + tu/24.
            
    #            print(tu,tud)
            
            fkp,ap = kp(IP,tu,iedi)
        
            if(int(tu)%3)==0:
                if IFLAG==-1:
                    fkpp.append([tud]+[fkp]+[ap])
                    IFLAG=1
                    APOLD = ap
                else:
                    fkpp.append([tud]+[fkp]+[APOLD])
                    fkpp.append([tud]+[fkp]+[ap])
                    APOLD = ap
            else:
                fkpp.append([tud]+[fkp]+[ap])
                APOLD = ap               
        except:pass
        
    print(tud,fkp,ap)    
    fkpp.append([INF]*3)
    df=pd.DataFrame(fkpp)
    df.to_csv(dire+fout, sep=' ',float_format='%.4f',header=False,index=False)  
    
    return 