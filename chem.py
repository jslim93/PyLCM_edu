import numpy as np
from matplotlib import pyplot as plt
from parameters import *

#Dissociation
#    1 SO2(g)            7 S(IV)(aq)
#    2 O3(g)             8 O3(aq)
#    3 H2O2(g)           9 H2O2(aq)
#    4 NH3(g)           10 N(III)(aq)
#    5 HNO3(g)          11 N(V)(aq)
#    6 CO2(g)           12 CO2 (aq)
#                       13 S(VI)(aq)

def eqcon(temp,neql):
    ek = np.zeros(neql)
    delt=(1./temp)-(1./298.)
    
    ek[0]=1.23*np.exp(3150.*delt)       #SO2 (g) <==> SO2 aq
    ek[1]=1.3e-2*np.exp(1960.*delt)        # SO2 aq  <==> HSO3- + H+ !
    ek[2]=6.6e-8*np.exp(1500.*delt)     #HSO3-   <==>  SO32- + H+    !
    ek[3]=1000.                         #S(VI)
    ek[4]=1.02e-2*np.exp(2720.*delt)    #S(VI)
    ek[5]=7.45e4*np.exp(7300.*delt)     #H2O2 (g) <==> H2O2 (aq)
    ek[6]=2.1e5                         #HNO3 (g) <==> HNO3 (aq)
    ek[7]=15.4*np.exp(8700.*delt)       #HNO3     <==> NO3- + H+  !
    ek[8]=1.13e-2*np.exp(2540.*delt)    #O3 (g) <==> O3 (aq)
    ek[9]=62.*np.exp(4110.*delt)        #NH3(g) <==> NH3 (aq)
    ek[10]=1.7e-5*np.exp(-450.*delt)    #NH3 <==> NH4+ + OH-      !
    ek[11]=1.0e-14*np.exp(-6710.*delt)  #H20
    ek[12]=3.4e-2*np.exp(2440.*delt)    #CO2 (g) <==> CO2 (aq)
    ek[13]=4.3e-7*np.exp(-1.e3*delt)    #CO2 (aq) <==> HCO3- + H+  !
    ek[14]=4.68e-11*np.exp(-1.76e3*delt) #HCO3-    <==> CO32- + H+ !

    return ek

def henry(ek,h,ngas): 
# h is the concentration of [H+]. heff is the effective Henry's Law constant 
# ----------------------------------------------------------------------
# Right now 5 Henrys Law coefficients are hardwired in. Input more
# equaitons into this subroutine and change dimension in dropcl
# ----------------------------------------------------------------------
#
# multiply gas partial pressure to get aq. concentraiton. 

    heff = np.zeros(ngas)
    
    #Heff_SO2
    heff[0]=ek[0]*(1+(ek[1]/ h )+(ek[1]*ek[2]/( h **2)))
    #O3 (g) <==> O3 (aq)
    heff[1]=ek[8]
    #H2O2 (g) <==> H2O2 (aq)
    heff[2]=ek[5]
    #NH3
    heff[3]=ek[9]*(1+(ek[10]* h /ek[11]))
    #HNO3 
    heff[4]=ek[6]*((ek[7]/ h ))
    #Heff_CO2
    heff[5]=ek[12]*(1+(ek[13]/ h )+(ek[13]*ek[14]/( h **2)))

    return heff

def accoeff(nspeci,rad):

    xkmt = np.zeros(nspeci)
    eta = np.zeros(nspeci)
    xmfp  = 0.065E-4 #mean free path
    diff  = 0.1

    #--Mass accomodation coeff. SO2, o3, H2O2, HNO3, NH3, CO2
    aw =  [3.5e-2,5.3e-4,1.8e-2,5.0e-2,5.0e-2,5.0e-2]
    
    for i in range(nspeci):
        xkn  = xmfp/rad
        xkni = 1./xkn
        t1=(1.33+(0.71*xkni))/(1.+xkni)   # changed to 0.71 from 0.071 Yiping
        t2=4*(1.-aw[i])/(3.*aw[i])
        t3=(t1+t2)*xkn

        eta[i] = 1./(1.+t3)
    #-- compute mass-transfer coefficients. In this problem, all are equal.
        xkmt[i]=3.* eta[i] * diff / (rad**2)
        
    return eta, xkmt 

def rsource(ya,temp, so2,h,ek,vol):
    d_temp =1./temp-1./298.15
    
    hso3 = ek[0]*so2 /(1.+(h / ek[1]) + (ek[2] / h))
    
    so3  = ek[2] * hso3 / h
    
    so2  = so2 - hso3 - so3

    # o3_aq * 
    r1=ya[7]*((2.4e4*so2)+ (3.7e5*np.exp(-5530.*d_temp)*hso3) + (1.5e9*np.exp(-5280.*d_temp)*so3))  #!!ya(ngas+nbins+l) :: H2O
    
    r2=7.45e7*np.exp(-4430.*d_temp)* h * ya[8]*hso3/(1.+(13.*h)) #!!  ya(ngas+2*nbins+l) :: H2O2 

    rsrc = np.zeros(13)
    
    rsrc[0]=0.
    rsrc[1]=0.
    rsrc[2]=0.
    rsrc[3]=0.
    rsrc[4]=0.
    rsrc[5]=0.       # now CO2
    rsrc[6]=-r1-r2   # SO2(aq)
    rsrc[7]=-r1      # O3 (aq)
    rsrc[8]=-r2      # H2O2 (aq)
    rsrc[9]=0.      # NH3(aq)
    rsrc[10]=0.      # HNO3(aq)
    rsrc[11]=0.      # CO2 (aq)
    rsrc[12]=r1+r2   # S(VI)

    return(rsrc)
from scipy.optimize import fsolve

def elecnuet(ya,ek):
    elecnuet = lambda x: x + x * ya[9] * ek[10]/(ek[10] + ek[10]*x) - ek[7]*ya[10] / (x + ek[7]) - ya[6] * ek[1] * (x + 2*ek[2]) / (x**2 + x*ek[1] + ek[1] * ek[2]) - ek[11]/x- ya[12] * (x + 2* ek[4]) / (x + ek[4])- ya[11] * ek[12] * (x + 2*ek[13]) / (x**2 + x*ek[12] + ek[12]*ek[13])

    a,b = fsolve(elecnuet, [1E-9, 1E1])
    #print("PH = ",-np.log10(b))
    return(a)