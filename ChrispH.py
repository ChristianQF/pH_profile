import numpy
from scipy.integrate import odeint
def tritation(pH, pKa,ConcBuffer):
    AHi=ConcBuffer/(1+10**(pH-pKa))
    Ai=ConcBuffer-AHi
    #AHi=0.09998201
    #Ai=0.0000179855
    
    
    Hi=10.**(-pH)
    OHi=10.**(-(14-pH))
    trit=[]
#Declarando variables
    k1f=10.**(-pKa)
    k1r=1
    kwf=10.**(-14)
    kwr=1
    def model(P, t):
        return [-k1f*P[0]+k1r*P[1]*P[2],
               k1f*P[0]-k1r*P[1]*P[2],
               k1f*P[0]-k1r*P[1]*P[2]+kwf-kwr*P[2]*P[3],
               kwf-kwr*P[2]*P[3]]
    suma=14
    vol=0.01 
    volacumulado=0.0
    for x in range(30):
        ts = numpy.linspace(0, 5000000000, 1000000) 
        P0 = [AHi,Ai,Hi,OHi+vol] 
        Ps = odeint(model, P0, ts) 
        AH = Ps[:,0]
        A = Ps[:,1]
        H = Ps[:,2]
        OH = Ps[:,3]
        AHi=Ps[99999,0]
        Ai=Ps[99999,1]
        Hi=Ps[99999,2]
        OHi=Ps[99999,3]
        pHf=numpy.log10(1/Hi)
        pOHf=numpy.log10(1/OHi)
        volacumulado=volacumulado+vol
        trit.append([x+1,volacumulado,pHf,pOHf+pHf])
    return trit