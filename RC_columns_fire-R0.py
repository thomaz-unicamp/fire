# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 15:07:36 2021

@author: Thomaz Buttignol - University of Campinas

DEFLECTIONS OF RC COLUMNS SUBJECTED TO FIRE - Eurocode 2 and LITS (bilinear) constitutive laws. 

For more information about the design procedures, please consult the paper "SIMPLIFIED DESIGN PROCEDURES FOR THE STRUCTURAL ANALYSIS OF REINFORCED CONCRETE COLUMNS IN FIRE"

About the algorithm:
1. Input data should be inserted manually (materials geometrical and mechanical properties, fire temperature, time interval, load, safety facors, etc.).
2. The type of the concrete class of EC-2 (NSC siliceous, HSC class 1, 3) should be defined manually (see below: EC-2 Constitutive law )
3. The program computes the total displacements of RC columns. Enter the initial eccentricity as equal to 0 (ei =0).
4. The failure (type and time) mode is predicted considering an equivalent axial load. The initial eccentricity (ei) should be specified.
5. Parameters:
gc and gs = concrete and steel safety factors
gf = load magnifying factor
fcT = residual coprssive strength
EcT = residual elastic modulus 
fyT = residual steel yielding stress
fuT = residual steel ultimate stress
EsT = residual steel elastic modulus  
LITS model, EC-2 model:  
sigc, sigc2 = concrete stress
sigs, sigs2 = steel stress
epsc, epsc2 = concrete strain
epss, epss2 = steel strain
Pcr, Pcr2 = buckling load
sc2, ec2 = concrete stress-strain curve from Eurocode 2
Ecs2T = residual secant modulus
RF = residual load carrying capacity
ei = initial eccentricity
c0 = parameter dependent on the first order moment distribution (EC-2 nominal stiffness method)
delta = column deflections
TR = time of fire rewsistance
"""

import numpy as np
import matplotlib.pyplot as plt

######
#INPUT
######
N = 2500*10**3 #Load in N
ei = 3 #initial eccentricity in mm
c0 = 9*np.sqrt(3) #parameter dependent on the first order moment distribution (EC-2 nominal stiffness method)
Lcr=2090 # effective length in mm
Temp=[20,226.094,378.833,494.608,588.149,665.76,731.565,788.36] #ºC - cross-section homogenized temperature
Ts=[20,243,478,630,738,819,883,934,977] #ºC - steel bar temperature
time=[0,30,60,90,120,150,180,210] #time intervals in min
#partial safety factors
gc=1.4 #concrete safety factor
gs=1.15 #steel safety factor
#concrete
h = 500 #height of the cross-section in mm
Ic = 2283333333 #concrete cross-section moment of inetia in mm4
Ac=157587 #cross-sectional area in mm2
fc =36 # concrete compressive strength in MPa
Eci=0.85*5600*np.sqrt(fc) #Concrete tangent modulus in MPa
ag_b=3.11 #ratio of coarse aggregates to binder
binder=331 # binder in kg
#steel
Is=0 #second moment of area of steel bars about the centroid of the concrete cross-section
fyk=420 # steel yielding stress in MPa
fuk=620 #steel ultimate strength in MPa
Es=210000 #steel modulus of elasticity in MPa
As= 12*np.pi*16**2/4 #rebars cross-section area in mm2
epsu=0.2  #ultimate strain
t1=0; t2=0

############
#CALCULATION
############
fc=fc/gc
fyk=fyk/gs
fuk=fuk/gs
fcT=np.zeros(len(Temp)); fyT=np.zeros(len(Temp)); fuT=np.zeros(len(Temp)); EcT=np.zeros(len(Temp)); EsT=np.zeros(len(Temp))
sigc=np.zeros(len(Temp)); sigs=np.zeros(len(Temp)); epsc=np.zeros(len(Temp)); epss=np.zeros(len(Temp)); displ=np.zeros(len(Temp))
sigc2=np.zeros(len(Temp)); sigs2=np.zeros(len(Temp)); epsc2=np.zeros(len(Temp)); epss2=np.zeros(len(Temp)); displ2=np.zeros(len(Temp)); Ec2T=np.zeros(len(Temp)); Ecs2T=np.zeros(len(Temp))
FTS=np.zeros(len(Temp)); Mt=np.zeros(len(Temp)); ni=np.zeros(len(Temp)); et=np.zeros(len(Temp)); ni=np.zeros(len(Temp)); Neq=np.zeros(len(Temp)); ni=np.zeros(len(Temp))
Pcr=np.zeros(len(Temp)); Pcr2=np.zeros(len(Temp))
EI=np.zeros(len(Temp)); EI2=np.zeros(len(Temp))
ec2=np.zeros(100); sc2=np.zeros(100)
RF=np.zeros(len(Temp)); RF2=np.zeros(len(Temp))
EcT[0]=Ec2T[0]=Eci
M1=N*ei #first order moment
delta=np.zeros(len(Temp))
for i in range(len(Temp)):   
    #thermal strain (MC 2010)
    FTS[i]=10**-5*(Temp[i]-20)

#######################
#EC-2 Constitutive law 
#######################
 
    #NSC - siliceous aggregates
    kc=(1.02-0.00019959*Temp[i]-0.0000012179*Temp[i]**2) #concrete strength decay
    if Temp[i]<100: fcT[i]=fc
    else: fcT[i]=fc*kc

    #HSC - class 1
#    kc=(0.793969374522261+0.00153576547877942*Temp[i]-5.37438525581988E-06*Temp[i]**2+3.08586809970606E-09*Temp[i]**3) #concrete strength decay
#    if Temp[i]<100: fcT[i]=fc    
#    else: fcT[i]=fc*kc

    #HSC - class 3
#    kc=(1.4229383892-0.0144388493444795*Temp[i]+0.000110915415137104*Temp[i]**2-4.05059740276654E-07*Temp[i]**3+7.21789961664425E-10*Temp[i]**4-6.16905572722255E-13*Temp[i]**5+2.02863920275789E-16*Temp[i]**6) #concrete strength decay
#    if Temp[i]<100: fcT[i]=fc    
#    else: fcT[i]=fc*kc

   #steel
    if Ts[i]<400:
        fyT[i]=fyk
        fuT[i]=fuk
    else:
        ks=(-3.95071428512473+0.0362770562733001*Ts[i]-0.0000901666666580406*Ts[i]**2+8.818181817333E-08*Ts[i]**3 - 0.0000000000303030303*Ts[i]**4)  #steel strength decay
        fyT[i]=fyk*ks
        fuT[i]=fuk*ks  
    if Ts[i]<100:
        EsT[i]=Es
    elif Ts[i]<800:
        js=(1.51933323794317 - 0.00849959516430642*Ts[i] + 0.0000428965535188524*Ts[i]**2 - 1.0167364122549E-07*Ts[i]**3 + 1.0246501663E-10*Ts[i]**4 - 3.666666E-14*Ts[i]**5) #steel modulus decay
        EsT[i]=Es*js   
    elif Ts[i]<=1000:
        js=0.07 
        EsT[i]=Es*js
  
    #Stress-Strain diagram
    if Temp[i]>500:
        ec1=0.025 #concrete strain at maximum stress of EC-2
    else:
        ec1=0.00178023182954035+0.000042684706560099*Temp[i]+-3.78706531982297E-07*Temp[i]**2+2.32233578752537E-09*Temp[i]**3+-6.52471021333584E-12*Temp[i]**4+6.84980680892169E-15*Temp[i]**5 
    for j in range (100):
        if j==0:
            ec2[j]=0 #concrete strain of EC-2
        else:
            ec2[j]=ec2[j-1]+ec1/100
        sc2[j]=3*ec2[j]*fcT[i]/(ec1*(2+(ec2[j]/ec1)**3)) #concrete stress of EC-2
    for k in range (1,100):
        if sc2[k]<0.3*fcT[i] and Temp[i]>20:
            Ec2T[i]=sc2[k]/ec2[k] #tangent modulus of elasticity of EC-2 

########### 
#LITS model
###########        

#####################################Axial Displacements#########################################
        
    if Temp[i]==20:
        EcT[i]=Eci
    else:
        EcT[i]=-1/((np.log(1+ag_b**3)*(-5.26*10**-5*Temp[i]-9.73*10**-7*Temp[i]**2+3.23*10**-9*Temp[i]**3-4.42*10**-12*Temp[i]**4)+np.log(binder/100)*-np.exp(Temp[i]**0.31)*0.156**5)/1000)      

    #Euler's critical load
    EI[i]=EcT[i]*Ic+EsT[i]*Is
    Pcr[i]=np.pi**2*EI[i]/Lcr**2

    #effective axial load (due to initial eccentricity)   
    ni[i]= 1+(np.pi**2/c0)/(Pcr[i]/N-1)   
    if ni[i]<0: ni[i]=0
    else:
        Mt[i]= M1*ni[i] #total bending moment (first and second order effects)
    et[i]= Mt[i]/N
    Neq[i]= N/np.exp(-2.4*et[i]/h)  
    
    #stresses and strains         
    if sigs[i]<fyT[i]: #check steel bars linear-elastic limit
        sigc[i]=Neq[i]/(Ac+EsT[i]/EcT[i]*As)
        sigs[i]=Neq[i]/(Ac*EcT[i]/EsT[i]+As)        
        if sigc[i]>0.85*fcT[i]: #check concrete crushing
            sigc[i]=0.85*fcT[i]
            sigs[i]=(Neq[i]-sigc[i]*Ac)/As
        RF[i]=0.85*fcT[i]*Ac+fyT[i]*As #reaction force (residual load carrying capacity)
        epss[i]=sigs[i]/EsT[i]
        if sigc[i]==0.85*fcT[i]: epsc[i]=epss[i] #concrete yielding (localized crushing with large displacements)
        else: epsc[i]=sigc[i]/EcT[i]          

    if (sigs[i]>=fyT[i]) and (fyk==fuk):  #check steel bars yielding (elastoplastic model)
        sigc[i]=(Neq[i]-fyT[i]*As)/Ac
        sigs[i]=fyT[i]
        if sigc[i]>0.85*fcT[i]: #check concrete crushing
            sigc[i]=0.85*fcT[i]
            sigs[i]=min(fyT[i],(Neq[i]-sigc[i]*Ac)/As)
        RF[i]=0.85*fcT[i]*Ac+fyT[i]*As #reaction force 
        if sigc[i]==0.85*fcT[i]: #failure - concrete crushing
            epsc[i]=epss[i]=10
        else: 
            epsc[i]=sigc[i]/EcT[i]
            epss[i]=epsc[i]
        
    if (sigs[i]>fyT[i]) and (fyk<fuk): #check steel bars yielding (strain hardening model)
        Ess=(fuk-fyk)/(epsu-fyk/Es) #E's (modulus of elasticity of the hardening branch)
        sigc[i]=(Neq[i]-fyT[i]*As*(1+Ess*js/EsT[i]))/(Ac+Ess*js/EcT[i]*As) 
        sigs[i]=fyT[i]+(sigc[i]/EcT[i]-fyT[i]/EsT[i])*Ess*js        
        if sigc[i]>0.85*fcT[i]: #check concrete crushing
            sigc[i]=0.85*fcT[i]   
            sigs[i]=min(fuT[i],(Neq[i]-sigc[i]*Ac)/As)
        RF[i]=0.85*fcT[i]*Ac+fyT[i]*As #reaction force
        epss[i]=fyT[i]/EsT[i]+(sigs[i]-fyT[i])/(Ess*js)
        if sigc[i]==0.85*fcT[i]: epsc[i]=epss[i] #concrete yielding (localized crushing with large displacements)
        else: epsc[i]=sigc[i]/EcT[i]   
        
    #axial displacement
    displ[i]=(-epsc[i]+FTS[i])*Lcr #subtracting initial displacement

    #failure mode   
    if epsc[i]>(ec1+0.002) and Pcr[i]/Neq[i]<0.99:
        if t1==0:
            t1=time[i]
            sol_LITS=["crushing+buckling", time[i]]
    if RF[i]/Neq[i]<0.99 or epsc[i]>(ec1+0.002):
        if t1==0:
            t1=time[i]
            sol_LITS=["crushing", time[i]]      
    if Pcr[i]/Neq[i]<0.99 and t1==0:
        t1=time[i]
        sol_LITS=["buckling", time[i]] 

############
#Eurocode 2  
############
        
    #stresses and strains
    Ecs2T[i]= Ec2T[i] #1st assumption (tangent modulus equal to the secant modulus)
    if sigs2[i]<fyT[i]: #check steel bars linear-elastic limit         
        sigc2[i]=N/(Ac+EsT[i]/Ecs2T[i]*As) 
        sigs2[i]=N/(Ac*Ecs2T[i]/EsT[i]+As)
        epss2[i]=sigs2[i]/EsT[i]
        if sigc2[i]>fcT[i]: #check concrete crushing
            sigc2[i]=fcT[i]
            sigs2[i]=(N-sigc2[i]*Ac)/As
        if sigc2[i]==fcT[i]: 
            Ecs2T[i]=fcT[i]/ec1 
            epsc2[i]=epss2[i]
        if sigc2[i]<fcT[i]:           
            for k in range (1,100): #finding the secant modulus
                if sc2[k]<sigc2[i] and Temp[i]>20: 
                    Ecs2T[i]=sc2[k]/ec2[k] #secant modulus of elasticity of EC-2
            error=Ec2T[i]/Ecs2T[i] #comparing tangent and secant modulus
            while error>1.01 or error<0.99:
                sigc2[i]=N/(Ac+EsT[i]/Ecs2T[i]*As)
                for k in range (1,100):
                    if sc2[k]<sigc2[i] and Temp[i]>20:
                        ET=sc2[k]/ec2[k] #secant modulus of elasticity of EC-2                
                    error=ET/Ecs2T[i]
                    Ecs2T[i]=ET   
            sigc2[i]=N/(Ac+EsT[i]/Ecs2T[i]*As) 
            sigs2[i]=N/(Ac*Ecs2T[i]/EsT[i]+As)       
            epsc2[i]=sigc2[i]/Ecs2T[i]  
        RF2[i]=fcT[i]*Ac+fyT[i]*As #reaction force       
        if RF2[i]/N<0.99: #check load carrying capacity
            epss2[i]=sigs2[i]/EsT[i]

    if (sigs2[i]>fyT[i]) and (fyk==fuk): #check steel bars yielding (elastoplastic model)
        sigc2[i]=(N-fyT[i]*As)/Ac
        sigs2[i]=fyT[i]        
        if sigc2[i]>fcT[i]: #check concrete crushing
            sigc2[i]=fcT[i]
            sigs2[i]=min(fyT[i],(N-sigc2[i]*Ac)/As)
        if sigc2[i]==fcT[i]: 
            Ecs2T[i]=fcT[i]/ec1 
        else:         
            for k in range (1,100): #finding the secant modulus
                if sc2[k]<sigc2[i] and Temp[i]>20: 
                    Ecs2T[i]=sc2[k]/ec2[k] #secant modulus of elasticity of EC-2
            error=Ec2T[i]/Ecs2T[i] #comparing tangent and secant modulus
            while error>1.01 or error<0.99:
                sigc2[i]=(N-fyT[i]*As)/Ac
                for k in range (1,100):
                    if sc2[k]<sigc2[i] and Temp[i]>20:
                        ET=sc2[k]/ec2[k] #secant modulus of elasticity of EC-2                
                    error=ET/Ecs2T[i]
                    Ecs2T[i]=ET        
        RF2[i]=fcT[i]*Ac+fyT[i]*As #reaction force
        if sigc2[i]==fcT[i]: #failure - concrete crushing
            epsc2[i]=epss2[i]=10
        else: 
            epsc2[i]=sigc2[i]/Ecs2T[i]     
            epss2[i]=epsc2[i]

    if (sigs2[i]>fyT[i]) and (fyk<fuk): #check steel bars yielding (strain hardening model)
        Ess=(fuk-fyk)/(epsu-fyk/Es) #E's (modulus of elasticity of the hardening branch)
        sigc2[i]=(N-fyT[i]*As*(1+Ess*js/EsT[i]))/(Ac+Ess*js/EsT[i]*As)
        sigs2[i]=fyT[i]+(sigc2[i]/Ecs2T[i]-fyT[i]/EsT[i])*Ess*js               
        if sigc2[i]>1.001*fcT[i]: #check concrete crushing
            sigc2[i]=fcT[i]
            sigs2[i]=min(fuT[i],(N-sigc2[i]*Ac)/As)
        if sigc2[i]==fcT[i]: 
            Ecs2T[i]=fcT[i]/ec1 
        else:         
            for k in range (1,100): #finding the secant modulus
                if sc2[k]<sigc2[i] and Temp[i]>20: 
                    Ecs2T[i]=sc2[k]/ec2[k] #secant modulus of elasticity of EC-2
            error=Ec2T[i]/Ecs2T[i] #comparing tangent and secant modulus
            while error>1.01 or error<0.99:
                sigc2[i]=(N-fyT[i]*As+fyT[i]*Ess*js/EsT[i]*As)/(Ac+Ess*js/EsT[i]*As)
                for k in range (1,100):
                    if sc2[k]<sigc2[i] and Temp[i]>20:
                        ET=sc2[k]/ec2[k] #secant modulus of elasticity of EC-2                
                    error=ET/Ecs2T[i]
                    Ecs2T[i]=ET        
            sigc2[i]=(N-fyT[i]*As*(1+Ess*js/EsT[i]))/(Ac+Ess*js/EsT[i]*As)
            sigs2[i]=fyT[i]+(sigc2[i]/Ecs2T[i]-fyT[i]/EsT[i])*Ess*js 
        RF2[i]=fcT[i]*Ac+fyT[i]*As #reaction force  
        epss2[i]=fyT[i]/EsT[i]+(sigs2[i]-fyT[i])/(Ess*js)
        if sigc2[i]==fcT[i]: epsc2[i]=epss2[i] #concrete yielding (localized crushing with large displacements)
        else: epsc2[i]=sigc2[i]/Ecs2T[i]  
    
    #axial displacement
    displ2[i]=(-epsc2[i]+FTS[i])*Lcr #subtracting initial displacement

        
########
#Output 
########
 
#time of fire resistance (TR)
j=np.searchsorted(time, t1)
if j==0: j=i
if RF[j]/Neq[j]<1:
    TR=time[j-1]+30/(RF[j-1]/Neq[j-1]-RF[j]/Neq[j])*(RF[j-1]/Neq[j-1]-1)
if RF[j]/Neq[j]>1 and Pcr[j]/Neq[j]<1:
    TR=time[j-1]+30/(Pcr[j-1]/Neq[j-1]-Pcr[j]/Neq[j])*(Pcr[j-1]/Neq[j-1]-1)
if RF[j]/Neq[j]>1 and Pcr[j]/Neq[j]>1:
    TR=0
print(TR)
print("failure:","time=",t1,"Nr/Neq",RF[j]/Neq[j],"Pcr/Neq=",Pcr[j]/Neq[j])
print(EI[j],0.85*fcT[j],fyT[j],EsT[j],EcT[j])
     
#plot   
#print (epsc/epss)
#print(epsc2/epss2)
plt.ylim((5, -20))
plt.plot(time,displ)
plt.plot(time,displ2)

#export data
d=[displ, displ2]
d=np.transpose(d)
np.savetxt("displ.txt", d, fmt='%.4f', delimiter='        ')

