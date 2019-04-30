import math
import matplotlib.pyplot as plt
import numpy as np

class Kappale:
    #Kolmion muotoinen kappale.
    #Keksitään alkuarvot:
    #xa(cm) = 0, ya(cm) = 0, kulma(alussa) = 0
    #maa tasossa -2
    #xaCM = 0
    #yaCM = 0
    ka = 0
    kl = ka #kulmalopussa
    #   kolmion kärki1 = k1, kärki2 = k2, kärki3 = k3
    #   x(k1) = 0,  y(k1) = 2
    #   x(k2) = -1, y(k2) = -1
    #   x(k3) = 1,  y(k3) = -1

    #massa 1 kg
    m = 1
    #hitausmomentti J = m*r^2 (r=0.75)
    J=m*(0.75)**2

    #vxa(cm) = 5 m/s eli kyseessä siis massakeskipisteen x komponentin suuntainen alkunopeus
    #vya(cm) = 10 m/s
    #wa = 0,2 r/s eli alkupyörimisnopeus

    #alustetaan muutama muuttuja, jotka selitetään koodissa
    w = 0
    rk1 = 0
    rk2 = 0
    rk3 = 0
    WxRK1 = 0
    WxRK2 = 0
    WxRK3 = 0
    vak1 = 0
    vak2 = 0
    vak3 = 0
    I = 0
    RK1xN = 0
    RK2xN = 0
    RK3xN = 0
    maaViivaPieninX = 0 #apumuuttuja, millä piirretään maaviiva
    maaViivaSuurinX = 0

    #konstruktorilla alkuarvot:
    def __init__(self, xaCM, yaCM, vxa, vya, wa):
        self.xaCM = xaCM #x alussa massakeskipiste
        self.xlCM = xaCM
        self.xkoord = [xaCM]
        self.xkoordkulma1 = [xaCM] #kolmion kulma 1
        self.xkoordkulma2 = [xaCM-1]
        self.xkoordkulma3 = [xaCM+1]
        self.yaCM = yaCM
        self.ylCM = yaCM
        self.ykoord = [yaCM]
        self.ykoordkulma1 = [yaCM+2]
        self.ykoordkulma2 = [yaCM-1]
        self.ykoordkulma3 = [yaCM-1]
        #vektori (rk1k2) kulmapisteestä 1 kulmapisteeseen 2
        self.rk1k2 = np.array([(xaCM-1)-(xaCM),(yaCM-1)-(yaCM+2),0])
        self.rk2k3 = np.array([(xaCM+1)-(xaCM-1),(yaCM-1)-(yaCM-1),0])
        self.rk3k1 = np.array([(xaCM)-(xaCM+1),(yaCM+2)-(yaCM-1),0])
        self.vxa = vxa #x suuntainen alkunopeus
        self.vxl = vxa
        self.vya = vya
        self.vyl = vya
        self.wa = wa #pyörimisnopeus
        self.wl = wa


#muita arvoja:
g = 9.81
e = 0.7 #kimmoisuuskerroin
tormaysKappaleeseen1 = False
tormaysKappaleeseen2 = False
tormaysMaahan = False
n = np.array([0,1,0]) #maan törmäysnormaali
nKappaleet = np.array([0,0,0]) #kappaleiden törmäysnormaalin alustus
kYksikkoVektori = np.array([0,0,1]) #yksikkövektori k (eli z suuntainen)
dt = 0.15 #delta t

def tormaysArvot(k, i):
    if(i!=1): #ensimmäisellä iteraatiolla ei käydä läpi
        k.w = np.array([0,0,k.wa])
        #etäisyys rp kärkiin k1,k2,k3
        k.rk1 = np.array([k.xkoordkulma1[-1]-k.xlCM,k.ykoordkulma1[-1]-k.ylCM,0])
        k.rk2 = np.array([k.xkoordkulma2[-1]-k.xlCM,k.ykoordkulma2[-1]-k.ylCM,0])
        k.rk3 = np.array([k.xkoordkulma3[-1]-k.xlCM,k.ykoordkulma3[-1]-k.ylCM,0])
        #ristitulo
        k.WxRK1 = np.cross(k.w,k.rk1)
        k.WxRK2 = np.cross(k.w,k.rk2)
        k.WxRK3 = np.cross(k.w,k.rk3)
        #pisteen nopeus:
        k.vak1 = np.array([k.vxl+k.WxRK1[0],k.vyl+k.WxRK1[1],0+k.WxRK1[2]])
        k.vak2 = np.array([k.vxl+k.WxRK2[0],k.vyl+k.WxRK2[1],0+k.WxRK2[2]])
        k.vak3 = np.array([k.vxl+k.WxRK3[0],k.vyl+k.WxRK3[1],0+k.WxRK3[2]]) 

def tormaysTarkastelu(k, i): #törmäystarkastelu maahan
    if (i!=1):
        #kulmapiste 1 osuu maahan (maa tasossa -2)
        if (k.ykoordkulma1[-1] <= -2) and (np.dot(k.vak1,n) < 0):
            #rp x n
            k.RK1xN = np.cross(k.rk1,n)
            #impulssi
            k.I = -(1+e)*(np.dot(k.vak1,n))/((1/k.m)+((np.linalg.norm(k.RK1xN))**2)/(k.J))
            k.vxl = k.vxl + (k.I/k.m)*(n[0])
            k.vyl = k.vyl + (k.I/k.m)*(n[1])
            k.wl = k.wa + k.I/k.J*(k.RK1xN[2])
            k.xlCM = k.xkoord[-1] + k.vxl*(dt) 
            k.ylCM = k.ykoord[-1] + k.vyl*(dt)
        #kulmapiste 2 osuu maahan
        elif (k.ykoordkulma2[-1] <= -2) and (np.dot(k.vak2,n) < 0):
            #rp x n
            k.RK2xN = np.cross(k.rk2,n)
            #impulssi
            k.I = -(1+e)*(np.dot(k.vak2,n))/((1/k.m)+((np.linalg.norm(k.RK2xN))**2)/(k.J))
            k.vxl = k.vxl + (k.I/k.m)*(n[0])
            k.vyl = k.vyl + (k.I/k.m)*(n[1])
            k.wl = k.wa + k.I/k.J*(k.RK2xN[2])
            k.xlCM = k.xkoord[-1] + k.vxl*(dt) 
            k.ylCM = k.ykoord[-1] + k.vyl*(dt)
        #kulmapiste 3 osuu maahan
        elif (k.ykoordkulma3[-1] <= -2) and (np.dot(k.vak3,n) < 0):
            #rp x n
            k.RK3xN = np.cross(k.rk3,n)
            #impulssi
            k.I = -(1+e)*(np.dot(k.vak3,n))/((1/k.m)+((np.linalg.norm(k.RK3xN))**2)/(k.J))
            k.vxl = k.vxl + (k.I/k.m)*(n[0])
            k.vyl = k.vyl + (k.I/k.m)*(n[1])
            k.wl = k.wa + k.I/k.J*(k.RK3xN[2])
            k.xlCM = k.xkoord[-1] + k.vxl*(dt) 
            k.ylCM = k.ykoord[-1] + k.vyl*(dt)
        #ei osuta maahan
        else:
            k.vxl = k.vxa #loppunopeus = alkunopeus
            k.vyl = k.vya - g*(dt)
            k.xlCM = k.xkoord[-1] + (k.vxa+k.vxl)/2*(dt) 
            k.ylCM = k.ykoord[-1] + (k.vya+k.vyl)/2*(dt)
            k.wl = k.wa #kulmanopeus sama ennen törmäystä

def liikutaKappale(k, i):
    if(i == 1): #ensimmäisellä iteraatiolla:
        k.xkoord[0] = k.xlCM
        k.ykoord[0] = k.ylCM
        #maaviivan piirtämiseen tarvittavat X muuttujat:
        k.maaViivaPieninX = k.xlCM
        k.maaViivaSuurinX = k.xlCM
    else: #muilla iteraatioilla
        k.xkoord.append(k.xlCM)
        k.ykoord.append(k.ylCM)
        #maaviivan piirtämiseen tarvittavat X muuttujat:
        if (k.xlCM < k.maaViivaPieninX):
            k.maaViivaPieninX = k.xlCM
        if (k.xlCM > k.maaViivaSuurinX):
            k.maaViivaSuurinX = k.xlCM
    
    k.kl = k.kl + k.wl*dt *2*(math.pi) #kulma lopussa
    
    #kulmapisteiden alkuarvot    
    kulma1x = 0
    kulma1y = 2
    kulma2x = -1
    kulma2y = -1
    kulma3x = 1
    kulma3y = -1

    #massakeskipisteen avulla kolmion kulmapisteet
    #ensiksi käännetään ja sitten siirretään
    xsiirto = k.xlCM - k.xaCM
    ysiirto = k.ylCM - k.yaCM
    #kulmapisteen 1 kääntö:
    xk1l = kulma1x*math.cos(k.kl)-kulma1y*math.sin(k.kl)
    yk1l = kulma1x*math.sin(k.kl)+kulma1y*math.cos(k.kl)
    #kulmapisteen 1 siirto:
    xk1l = xk1l + xsiirto + k.xaCM
    yk1l = yk1l + ysiirto + k.yaCM
    k.xkoordkulma1.append(xk1l)
    k.ykoordkulma1.append(yk1l)
    #pisteet 2 ja 3
    xk2l = kulma2x*math.cos(k.kl)-kulma2y*math.sin(k.kl)
    yk2l = kulma2x*math.sin(k.kl)+kulma2y*math.cos(k.kl)
    xk2l = xk2l + xsiirto + k.xaCM
    yk2l = yk2l + ysiirto + k.yaCM
    k.xkoordkulma2.append(xk2l)
    k.ykoordkulma2.append(yk2l)
    xk3l = kulma3x*math.cos(k.kl)-kulma3y*math.sin(k.kl)
    yk3l = kulma3x*math.sin(k.kl)+kulma3y*math.cos(k.kl)
    xk3l = xk3l + xsiirto + k.xaCM
    yk3l = yk3l + ysiirto + k.yaCM
    k.xkoordkulma3.append(xk3l)
    k.ykoordkulma3.append(yk3l)

    k.vxa = k.vxl
    k.vya = k.vyl
    k.wa = k.wl
        
#Pisteen (x,y) etäisyys suorasta, joka kulkee pisteestä (x1,y1) -> (x2,y2)
def pisteenEtaisyysSuorasta(x1,y1,x2,y2,x,y):
    d = math.fabs((x2-x1)*(y-y1)-(x-x1)*(y2-y1))/math.sqrt((x2-x1)**2+(y2-y1)**2)
    return d

def tormaysTarkasteluKappaleet(k,k2): #kahden kappaleen törmäys tarkastelu
    #Tarkastetaan onko mikään kappaleen 1 kulmapisteistä kappaleen 2 sisällä
    #Vektorissa (rk2kulma1p1), p1 (eli piste1) on kappaleen 1 tarkasteltava kulmapiste. 
    #Kyseessä on siis vektori kappaleen 2 kulmapisteestä 1 tarkasteltavaan kappaleen 1 kulmapisteeseen.
    rk2kulma1p1 = np.array([k.xkoordkulma1[-1]-k2.xkoordkulma1[-1],k.ykoordkulma1[-1]-k2.ykoordkulma1[-1],0])
    rk2kulma2p1 = np.array([k.xkoordkulma1[-1]-k2.xkoordkulma2[-1],k.ykoordkulma1[-1]-k2.ykoordkulma2[-1],0])
    rk2kulma3p1 = np.array([k.xkoordkulma1[-1]-k2.xkoordkulma3[-1],k.ykoordkulma1[-1]-k2.ykoordkulma3[-1],0])
    #mahdolliset törmäyspisteet 2 ja 3
    rk2kulma1p2 = np.array([k.xkoordkulma2[-1]-k2.xkoordkulma1[-1],k.ykoordkulma2[-1]-k2.ykoordkulma1[-1],0])
    rk2kulma2p2 = np.array([k.xkoordkulma2[-1]-k2.xkoordkulma2[-1],k.ykoordkulma2[-1]-k2.ykoordkulma2[-1],0])
    rk2kulma3p2 = np.array([k.xkoordkulma2[-1]-k2.xkoordkulma3[-1],k.ykoordkulma2[-1]-k2.ykoordkulma3[-1],0])
    rk2kulma1p3 = np.array([k.xkoordkulma3[-1]-k2.xkoordkulma1[-1],k.ykoordkulma3[-1]-k2.ykoordkulma1[-1],0])
    rk2kulma2p3 = np.array([k.xkoordkulma3[-1]-k2.xkoordkulma2[-1],k.ykoordkulma3[-1]-k2.ykoordkulma2[-1],0])
    rk2kulma3p3 = np.array([k.xkoordkulma3[-1]-k2.xkoordkulma3[-1],k.ykoordkulma3[-1]-k2.ykoordkulma3[-1],0])
    
    rk1k2 = np.array([k2.xkoordkulma2[-1]-k2.xkoordkulma1[-1],k2.ykoordkulma2[-1]-k2.ykoordkulma1[-1],0])
    rk2k3 = np.array([k2.xkoordkulma3[-1]-k2.xkoordkulma2[-1],k2.ykoordkulma3[-1]-k2.ykoordkulma2[-1],0])
    rk3k1 = np.array([k2.xkoordkulma1[-1]-k2.xkoordkulma3[-1],k2.ykoordkulma1[-1]-k2.ykoordkulma3[-1],0])

    #Tarkastetaan ristitulon avulla, onko piste 1 kappaleen sisällä.
    tormaysRisti1 = np.cross(rk1k2,rk2kulma1p1)
    tormaysRisti2 = np.cross(rk2k3,rk2kulma2p1)
    tormaysRisti3 = np.cross(rk3k1,rk2kulma3p1)
    #pisteet 2 ja 3
    tormays2Risti1 = np.cross(rk1k2,rk2kulma1p2)
    tormays2Risti2 = np.cross(rk2k3,rk2kulma2p2)
    tormays2Risti3 = np.cross(rk3k1,rk2kulma3p2)
    tormays3Risti1 = np.cross(rk1k2,rk2kulma1p3)
    tormays3Risti2 = np.cross(rk2k3,rk2kulma2p3)
    tormays3Risti3 = np.cross(rk3k1,rk2kulma3p3)
    
    #Törmäys tapahtunut kappaleiden välillä.
    if (tormaysRisti1[2] >= 0) and (tormaysRisti2[2] >= 0) and (tormaysRisti3[2] >= 0):
        #Tarkastetaan mitä sivua lähinnä törmäyspiste on.
        etaisyys1 = pisteenEtaisyysSuorasta(k2.xkoordkulma1[-1],k2.ykoordkulma1[-1],k2.xkoordkulma2[-1],k2.ykoordkulma2[-1],k.xkoordkulma1[-1],k.ykoordkulma1[-1])
        etaisyys2 = pisteenEtaisyysSuorasta(k2.xkoordkulma2[-1],k2.ykoordkulma2[-1],k2.xkoordkulma3[-1],k2.ykoordkulma3[-1],k.xkoordkulma1[-1],k.ykoordkulma1[-1])
        etaisyys3 = pisteenEtaisyysSuorasta(k2.xkoordkulma3[-1],k2.ykoordkulma3[-1],k2.xkoordkulma1[-1],k2.ykoordkulma1[-1],k.xkoordkulma1[-1],k.ykoordkulma1[-1])
        
        if(etaisyys1 < etaisyys2):
            if(etaisyys1 < etaisyys3):
                #etaisyys 1 lyhin
                nKappaleet = np.cross(rk1k2/math.sqrt(rk1k2[0]**2+rk1k2[1]**2),kYksikkoVektori)
            else:
                #etaisyys 3 lyhin
                nKappaleet = np.cross(rk3k1/math.sqrt(rk3k1[0]**2+rk3k1[1]**2),kYksikkoVektori)
        else:
            #etaisyys2 lyhin
            nKappaleet = np.cross(rk2k3/math.sqrt(rk2k3[0]**2+rk2k3[1]**2),kYksikkoVektori)
        
        #kappaleen 2 massakeskipisteen etäisyys törmäyspisteeseen
        rBP = np.array([k.xkoordkulma1[-1]-k2.xkoord[-1],k.ykoordkulma1[-1]-k2.ykoord[-1],0])
        #kappaleen 1
        rAP = np.array([k.xkoordkulma1[-1]-k.xkoord[-1],k.ykoordkulma1[-1]-k.ykoord[-1],0])

        laskeTormaysKappaleet(k,k2,nKappaleet,rBP,rAP)

        return True; #törmäys tapahtui

    elif (tormays2Risti1[2] >= 0) and (tormays2Risti2[2] >= 0) and (tormays2Risti3[2] >= 0):
        etaisyys1 = pisteenEtaisyysSuorasta(k2.xkoordkulma1[-1],k2.ykoordkulma1[-1],k2.xkoordkulma2[-1],k2.ykoordkulma2[-1],k.xkoordkulma2[-1],k.ykoordkulma2[-1])
        etaisyys2 = pisteenEtaisyysSuorasta(k2.xkoordkulma2[-1],k2.ykoordkulma2[-1],k2.xkoordkulma3[-1],k2.ykoordkulma3[-1],k.xkoordkulma2[-1],k.ykoordkulma2[-1])
        etaisyys3 = pisteenEtaisyysSuorasta(k2.xkoordkulma3[-1],k2.ykoordkulma3[-1],k2.xkoordkulma1[-1],k2.ykoordkulma1[-1],k.xkoordkulma2[-1],k.ykoordkulma2[-1])
        
        if(etaisyys1 < etaisyys2):
            if(etaisyys1 < etaisyys3):
                nKappaleet = np.cross(rk1k2/math.sqrt(rk1k2[0]**2+rk1k2[1]**2),kYksikkoVektori)
            else:
                nKappaleet = np.cross(rk3k1/math.sqrt(rk3k1[0]**2+rk3k1[1]**2),kYksikkoVektori)
        else:
            nKappaleet = np.cross(rk2k3/math.sqrt(rk2k3[0]**2+rk2k3[1]**2),kYksikkoVektori)
        
        rBP = np.array([k.xkoordkulma2[-1]-k2.xkoord[-1],k.ykoordkulma2[-1]-k2.ykoord[-1],0])
        rAP = np.array([k.xkoordkulma2[-1]-k.xkoord[-1],k.ykoordkulma2[-1]-k.ykoord[-1],0])
        laskeTormaysKappaleet(k,k2,nKappaleet,rBP,rAP)
        return True;

    elif (tormays3Risti1[2] >= 0) and (tormays3Risti2[2] >= 0) and (tormays3Risti3[2] >= 0):
        etaisyys1 = pisteenEtaisyysSuorasta(k2.xkoordkulma1[-1],k2.ykoordkulma1[-1],k2.xkoordkulma2[-1],k2.ykoordkulma2[-1],k.xkoordkulma3[-1],k.ykoordkulma3[-1])
        etaisyys2 = pisteenEtaisyysSuorasta(k2.xkoordkulma2[-1],k2.ykoordkulma2[-1],k2.xkoordkulma3[-1],k2.ykoordkulma3[-1],k.xkoordkulma3[-1],k.ykoordkulma3[-1])
        etaisyys3 = pisteenEtaisyysSuorasta(k2.xkoordkulma3[-1],k2.ykoordkulma3[-1],k2.xkoordkulma1[-1],k2.ykoordkulma1[-1],k.xkoordkulma3[-1],k.ykoordkulma3[-1])
        
        if(etaisyys1 < etaisyys2):
            if(etaisyys1 < etaisyys3):
                nKappaleet = np.cross(rk1k2/math.sqrt(rk1k2[0]**2+rk1k2[1]**2),kYksikkoVektori)
            else:
                nKappaleet = np.cross(rk3k1/math.sqrt(rk3k1[0]**2+rk3k1[1]**2),kYksikkoVektori)
        else:
            nKappaleet = np.cross(rk2k3/math.sqrt(rk2k3[0]**2+rk2k3[1]**2),kYksikkoVektori)

        rBP = np.array([k.xkoordkulma3[-1]-k2.xkoord[-1],k.ykoordkulma3[-1]-k2.ykoord[-1],0])
        rAP = np.array([k.xkoordkulma3[-1]-k.xkoord[-1],k.ykoordkulma3[-1]-k.ykoord[-1],0])
        laskeTormaysKappaleet(k,k2,nKappaleet,rBP,rAP)
        return True;
    return False; #Törmäystä ei tapahtunut.

def laskeTormaysKappaleet(k,k2,nKappaleet,rBP,rAP):
    k.w = np.array([0,0,k.wa])
    k2.w = np.array([0,0,k2.wa])
    #rBP ristitulo törmäysnormaalin kanssa
    rBPxN = np.cross(rBP,nKappaleet)
    rAPxN = np.cross(rAP,nKappaleet)
    #kappaleen 1 kulmanopeuden ristitulo etäisyyden rAP kanssa
    wAxrAP = np.cross(k.w,rAP)
    wBxrBP = np.cross(k2.w,rBP)
    #törmäyspisteen nopeudet
    vaAP = np.array([k.vxa + wAxrAP[0],k.vya + wAxrAP[1],0]) 
    vaBP = np.array([k2.vxa + wBxrBP[0],k2.vya + wBxrBP[1],0]) 
    #suhteellinen nopeus
    vaAB = vaAP - vaBP
    #impulssi
    I = -(1+e)*(np.dot(vaAB,nKappaleet))/(1/k.m + 1/k2.m + (np.linalg.norm(rAPxN))**2/(k.J) + (np.linalg.norm(rBPxN))**2/(k2.J))
    #uudet arvot törmäyksen jälkeen
    k.vxl = k.vxl + (I/k.m)*(nKappaleet[0])
    k.vyl = k.vyl + (I/k.m)*(nKappaleet[1])
    k.wl = k.wa + I/k.J*rAPxN[2]
    k.xlCM = k.xkoord[-1] + k.vxl*(dt) 
    k.ylCM = k.ykoord[-1] + k.vyl*(dt)
    k2.vxl = k2.vxl - (I/k2.m)*(nKappaleet[0])
    k2.vyl = k2.vyl - (I/k2.m)*(nKappaleet[1])
    k2.wl = k2.wa - I/k2.J*rBPxN[2]
    k2.xlCM = k2.xkoord[-1] + k2.vxl*(dt) 
    k2.ylCM = k2.ykoord[-1] + k2.vyl*(dt)

def kappaleidenTormaysLooppi(): #tässä funktiossa tarkistetaan kappaleiden välinen törmäys
    tormaysKappaleeseen2 = tormaysTarkasteluKappaleet(k1,k2)
    if(tormaysKappaleeseen2 == False):
        tormaysKappaleeseen1 = tormaysTarkasteluKappaleet(k2,k1)
    if (tormaysKappaleeseen2 == True) or (tormaysKappaleeseen1 == True):
        liikutaKappale(k1,i)
        liikutaKappale(k2,i)
    
    tormaysKappaleeseen1 = False;
    tormaysKappaleeseen2 = False;

#Luodaan kappale olio(t)
k1 = Kappale (0,0,7,-3,0.1)
k2 = Kappale (20,0,-9,-5,-0.1)

#lasketaan liike
for i in range(1,25):
    tormaysArvot(k1,i)
    tormaysTarkastelu(k1,i)
    liikutaKappale(k1,i)
    
    tormaysArvot(k2,i)
    tormaysTarkastelu(k2,i)
    liikutaKappale(k2,i)

    kappaleidenTormaysLooppi()


#piirretään:
plt.figure(1)
plt.axis('equal')
plt.title('Kahden kappaleen törmäys')
#massakeskipiste
plt.plot(k1.xkoord,k1.ykoord,"yo",markersize=4) 
plt.plot(k2.xkoord,k2.ykoord,"ro",markersize=4) 
#maaviiva
plt.plot([k1.maaViivaPieninX-2,k2.maaViivaSuurinX+2],[-2,-2],'brown')

def viivanplottaus(x,x2,y,y2,i,vari):
    x1, x2 = x[i], x2[i]
    y1, y2 = y[i], y2[i]
    plt.plot([x1,x2],[y1,y2],vari)

for i in range(1,len(k1.xkoordkulma1)):
    viivanplottaus(k1.xkoordkulma1,k1.xkoordkulma2,k1.ykoordkulma1,k1.ykoordkulma2,i,'-k')
    viivanplottaus(k1.xkoordkulma2,k1.xkoordkulma3,k1.ykoordkulma2,k1.ykoordkulma3,i,'-k')
    viivanplottaus(k1.xkoordkulma3,k1.xkoordkulma1,k1.ykoordkulma3,k1.ykoordkulma1,i,'-k')
    viivanplottaus(k2.xkoordkulma1,k2.xkoordkulma2,k2.ykoordkulma1,k2.ykoordkulma2,i,'-g')
    viivanplottaus(k2.xkoordkulma2,k2.xkoordkulma3,k2.ykoordkulma2,k2.ykoordkulma3,i,'-g')
    viivanplottaus(k2.xkoordkulma3,k2.xkoordkulma1,k2.ykoordkulma3,k2.ykoordkulma1,i,'-g')


plt.figure(2)
plt.axis('equal')
plt.title('Kaksi kappaletta heitetään maahan, josta ne kimpoavat ilmaan\n ja törmäävät toisiinsa. Pelkästään massakeskipisteet piirrettynä.\nAlkupisteet (0,0) ja (20,0)')
plt.plot(k1.xkoord,k1.ykoord,"yo",markersize=4) 
plt.plot(k2.xkoord,k2.ykoord,"ro",markersize=4) 
plt.plot([k1.maaViivaPieninX-2,k2.maaViivaSuurinX+2],[-2,-2],'brown')

plt.figure(3)
plt.title('Kappaleen törmäys maahan')
dt = 0.2
k3 = Kappale (0,0,5,10,-0.2)
for i in range(1,30):
    tormaysArvot(k3,i)
    tormaysTarkastelu(k3,i)
    liikutaKappale(k3,i)
plt.plot(k3.xkoord,k3.ykoord,"yo",markersize=4) 
plt.plot([k3.maaViivaPieninX-2,k3.maaViivaSuurinX+2],[-2,-2],'brown')
for i in range(1,len(k3.xkoordkulma1)):
    viivanplottaus(k3.xkoordkulma1,k3.xkoordkulma2,k3.ykoordkulma1,k3.ykoordkulma2,i,'-k')
    viivanplottaus(k3.xkoordkulma2,k3.xkoordkulma3,k3.ykoordkulma2,k3.ykoordkulma3,i,'-k')
    viivanplottaus(k3.xkoordkulma3,k3.xkoordkulma1,k3.ykoordkulma3,k3.ykoordkulma1,i,'-k')

#plt.axis('scaled')
plt.axis('equal')
plt.show()