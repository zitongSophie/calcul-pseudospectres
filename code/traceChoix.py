import random
import numpy as np
import numpy.matlib
from matplotlib import pyplot as plt
import tkinter as Tk
import sys
import time
from math import *


class ProgressBar:
    
    #This class allows you to make easily a progress bar.
    def __init__(self, steps, maxbar=100, title='Chargement'):
        if steps <= 0 or maxbar <= 0 or maxbar > 200:
            raise ValueError

        self.steps = steps
        self.maxbar = maxbar
        self.title = title

        self.perc = 0
        self._completed_steps = 0

        self.update(False)

    def update(self, increase=True):
        if increase:
            self._completed_steps += 1

        self.perc = floor(self._completed_steps / self.steps * 100)

        if self._completed_steps > self.steps:
            self._completed_steps = self.steps

        steps_bar = floor(self.perc / 100 * self.maxbar)

        if steps_bar == 0:
            visual_bar = self.maxbar * ' '
        else:
            visual_bar = (steps_bar - 1) * '=' + '>' + (self.maxbar - steps_bar) * ' '

        sys.stdout.write('\r' + self.title + ' [' + visual_bar + '] ' + str(self.perc) + '%\n')
        sys.stdout.flush()




def creer_A(taille):
    A1=np.random.randint(taille, size=(taille,taille))
    A2=np.random.randint(taille, size=(taille,taille))
    A=np.array(A1+A2*1j)
    print ("voici la matrice:\n",A)
    return A  



def creer_M(taille):
    X=np.random.rand(taille,taille)
    Z=np.random.rand(taille,taille)
    Y=100*X-50 +(100*Z-50)*1j
    print("\n Y= ",Y)
    return Y


def creer_D(L):
    D=np.diag(L)
    return D
        


###############################################################################
#L'ALGORITHME DE GRID
###############################################################################

def loca_vp(A):
    taille=len(A)
    t=len(A[0])
    som=0
    #verifions si la matrice est carre et est bien inversible
    if(t!=taille):
        print ("erreur matrice non carre")
        return

    if(np.linalg.det(A)==0):
        print ("erreur matrice non inversible")
        return
    
#calculer de disque de gerschgorin
    for i in range(taille):
        som=0
        for j in range(taille):
                if(j!=i):
                        som+=abs(A[i,j])
#calculer du bord au fur et a mesure du parcours dela matrice
        if(i==0):
            xsup= A[i,i].real + som
            xinf= A[i,i].real - som
            ysup= A[i,i].imag + som
            yinf= A[i,i].imag - som
        else:
            if(A[i,i].real+som > xsup):
                xsup= A[i,i].real + som
            if(A[i,i].real-som < xinf):
                xinf= A[i,i].real - som
            if(A[i,i].imag+som > ysup):
                ysup= A[i,i].imag + som
            if(A[i,i].imag-som<yinf):
                yinf= A[i,i].imag - som

    L=[xsup,xinf,ysup,yinf] 
    return  L




def locaPseudoVp(A,e):
    n=len(A)
    #appel de la fonction lacaliser spectre 
    l=loca_vp(A)
    r=[0,0,0,0]
    #modifier le bord contenant les preseudovaleurpropre
    r[0]=l[0]+e*n
    r[1]=l[1]-e*n
    r[2]=l[2]+e*n
    r[3]=l[3]-e*n
    return r


    
def dessinGrid(A,e, pas):
    taille=A.shape[0]
    L=locaPseudoVp(A,e)   
    x=np.arange(L[1],L[0]+pas,pas)    
    y=np.arange(L[3],L[2]+pas,pas)
    print("taille de x: ", x, "\ttaille de y: ", y)
    #creation de la grid
    lx=x.shape[0]
    ly=y.shape[0]
    #print(ly,lx, "L: liste des bornes\n", L)
    XX, YY=np.meshgrid(x,y)
    vpx=np.linalg.eig(A)[0].real
    vpy=np.linalg.eig(A)[0].imag
    #initialisation d une matrice
    sigmin=np.zeros((ly,lx))
    if __name__ == '__main__':
        i = 0
        bar = ProgressBar(100)
        bar.update()
        time.sleep(0.1)
        i=i+1
    progression=lx//90+1
    print("progression: ",progression)
    for m in range(lx):
        #print("m%progr   ",m%progression)
        if(m!=0 and m%progression<1 and bar.perc<98):
            bar.update()
            time.sleep(0.1)
            i=i+1
        for n in range(ly):
            u,sigma,vt = np.linalg.svd((x[m]+1j*y[n])*np.eye(taille)-A)
            sigmin[n,m]=min(sigma)        
    
    plt.title('pseudospectre')
    plt.ylabel('Y-partie imaginaire') 
    plt.xlabel('X-partie reelle')
    
    plt.plot(vpx,vpy,"ob")
    bar.update()
    time.sleep(1)
    i =i+ 1
    c=plt.contour(XX, YY,sigmin,[e])  
    while bar.perc!=100:
        bar.update()
        time.sleep(1)
        i =i+ 1
 
    plt.show()
    plt.close()





#pour le retrace
def nvdessin4(F,yx,ymin,xmax,xmin,nvpas,nve):      
    taille=F.shape[0]   
    x=np.arange(xmin,xmax+nvpas,nvpas)    
    y=np.arange(ymin,yx+nvpas,nvpas)
    #creation de la grid
    lx=x.shape[0]
    ly=y.shape[0]
    #print(ly,lx, "L: liste des bornes\n", L)
    XX, YY= np.meshgrid(x,y) 
    #initialisation d une matrice
    if __name__ == '__main__':
        i = 0
        bar = ProgressBar(100)
        bar.update()
        time.sleep(1)
        vpx=np.linalg.eig(F)[0].real
        vpy=np.linalg.eig(F)[0].imag
        bar.update()
        
    sigmin=np.zeros((ly,lx))
    prog=lx//100+1
    for m in range(lx): 
        if(m!=0 and  m%prog<1 and bar.perc<98):  
            bar.update()
            time.sleep(1)
            i=i+1  
        for n in range(ly):
            u,sigma,vt = np.linalg.svd((x[m]+1j*y[n])*np.eye(taille)-F)
            sigmin[n,m]=min(sigma)
            
    plt.title('pseudospectre mETHODE griD')
    plt.ylabel('Y-partie imaginaire') 
    plt.xlabel('X-partie reelle')
    #plt.plot(vpx,vpy,"ob")
    while(bar.perc<80 ):
        bar.update()
        time.sleep(0.5)
        i=i+ 1
    a=0
    tv=len(vpx)
    for a in range(tv):
        if(vpx[a]<xmax and  vpx[a]>xmin and vpy[a]<yx and vpy[a]>ymin):
            plt.scatter(vpx[a],vpy[a])
    
    #plt.axis([L[1],L[0],L[3],L[2]])
    c=plt.contour(XX, YY,sigmin,[nve])
    while bar.perc!=100:
        bar.update()
        time.sleep(0.1)
        i =i+ 1
           
    plt.show()
    plt.close()
    





###################################################
   #OPTIMISATION DE L'ALGO DE GRID
###################################################

    
def locaOPTI(A,e):
    taille=len(A)
    disque=[]    #calculer les desques pseudovps , (utiliser les disques de gerschgorin) 
    for i in range(taille):#liste des disques original (tuple(centre,rayon)) de chaque valeur propre:disque
        som=0
        d1=[]
        for j in range(taille):
                if(j!=i):
                        som+=abs(A[i,j])                
        d1=[A[i,i],som+taille*e]
        disque.append(d1)
    d=[] #liste des ensembles des disques en relation adjoint directe     
    for d1 in disque:    # verifier la relation adjoint 2 a 2  
        b={(d1[0],d1[1])}
        for d2 in disque:
            if(d2!=d1):
                if(abs(d1[0]-d2[0]) < d1[1]+d2[1]):
                    b.add((d2[0],d2[1]))
        if(b not in d):
            d.append(b)#liste des ensembles des disques en relation adjoint directe
    
    for i in range(len(d)):    #verifier union des disques en relation indirecte
        for j in range(len(d)) :
            if(j!=i):
                if(d[i].isdisjoint(d[j])==False): #verifier union des ensembles en relation indirecte
                    d[i]=d[i].union(d[j])#d[i] contient tout les boucles qui sont adjoint avec d[i]
                    d[j].clear() #ce qui est deja dans d[i] ,on le supprime
        
    borne=[]
    for d1 in d:#des composants :ensemble des boucles adjoints
        i=0
        if(len(d1)==0):#pour qu il nest pas nul
            continue
        for b in d1:#les disques dans chacun composant
            if(i==0):#trouver la borne
                xsup= b[0].real + b[1]
                xinf= b[0].real - b[1]
                ysup= b[0].imag + b[1]
                yinf= b[0].imag - b[1]
                i=1
            else:
                if(b[0].real + b[1] > xsup):
                    xsup= b[0].real + b[1]
                if(b[0].real - b[1] < xinf):
                    xinf= b[0].real - b[1]
                if(b[0].imag + b[1] > ysup):
                    ysup= b[0].imag + b[1]
                if(b[0].imag - b[1]<yinf):
                    yinf= b[0].imag - b[1]
        borne.append([xsup,xinf,ysup,yinf])#on a calcule les bornes par composant
    for i in borne:#cas pour les disques sont disjoints mais les bornes(rectangle) adjointes
        for j in borne:
            if(i!=j):
                if(abs((i[0]+i[1])/2-(j[0]+j[1])/2)<abs((i[0]-i[1])/2+(j[0]-j[1])/2)):
                    if(abs((i[2]+i[3])/2-(j[2]+j[3])/2)<abs((i[2]-i[3])/2+(j[2]-j[3])/2)):
                        i[0]=max(i[0],j[0])
                        i[1]=min(i[1],j[1])
                        i[2]=max(i[2],j[2])
                        i[3]=min(i[3],j[3])
    newborne=[]
    for i in borne:
        if i not in newborne:
            newborne.append(i)
    print ('newborne',newborne)
    return newborne


def dessin(A,borne,e,pas,taille):#borne=[xsup,xinf,ysup,yinf]
    #borne=locaOPTI(A,e) #########################################
    x=np.linspace(borne[1],borne[0],num=int((borne[0]-borne[1])/pas))    
    y=np.linspace(borne[3],borne[2],num=int((borne[2]-borne[3])/pas))
    #creation de la grid
    lx=x.shape[0]
    ly=y.shape[0]
    XX, YY= np.meshgrid(x,y)
    #initialisation d une matrice
    sigmin=np.zeros((ly,lx))
    progression=lx//101+1
    k=1
    for m in range(lx):
        for n in range(ly):
            u,sigma,vt = np.linalg.svd((x[m]+1j*y[n])*np.eye(taille)-A)
            sigmin[n,m]=min(sigma)
    plt.contour(XX, YY,sigmin,[e])


def dessinOpti(A,e,pas):
    taille=A.shape[0]
    L=locaOPTI(A,e)
    vpx=np.linalg.eig(A)[0].real
    vpy=np.linalg.eig(A)[0].imag
    plt.plot(vpx,vpy,"ob")
    plt.title('pseudospectre')
    plt.ylabel('Y-partie imaginaire') 
    plt.xlabel('X-partie reelle')
    for l in L:
        dessin(A,l,e,pas,taille)
    plt.show()
    plt.close()






#########################################################################
#Prediction Correction
#########################################################################

#trouver le premier point 
def prediCPremierPt(eig,e,taille,A):
    z1new=eig+e
    tol=10**(-3) #tolerance
    d=1
    z1new=eig+e*d
    while((abs(np.linalg.svd(z1new*np.eye(taille)-A,compute_uv=False)[taille-1]-e) >tol*e)):   #correction du premier point
        z1old=z1new
        u,sigma,vt = np.linalg.svd(z1old*np.eye(taille)-A)
        umin=u[:,len(u)-1]
        vmint=vt[len(vt)-1]
        vmin=np.conjugate(vmint)
        sigmi=np.linalg.svd(z1new*np.eye(taille)-A,compute_uv=False)[taille-1]
        z1new=eig-d*((sigma[taille-1]-e)/((-d*np.vdot(vmin,umin)).real))   #utiliser le theoreme de newton question:la difference pour umin et umin*
    return z1new


def prediCorrect(A,e,pas):
    eigv=np.linalg.eig(A)[0]    #la liste des valeurs propres
    for i in eigv:              #tracer les points vp
        plt.scatter(i.real,i.imag,color='blue')
    taille=len(A)
    tt=len(eigv)
    if __name__ == '__main__':
        i = 0
        bar = ProgressBar(100)
        bar.update()
        time.sleep(0.1)
        i=i+1
    progression=tt//101+1
    k=1
    #calculer pour chaque valeur propre
    for eig in eigv:
        if(k%progression<0.2 and bar.perc<98 and k!=1):
            bar.update()
            time.sleep(0.1)
            i=i+1
        k=k+1
        ps=[]
        #step 0:initialisation du premier point z1
        z1new=prediCPremierPt(eig,e,taille,A)
        zk=z1new
        plt.scatter(zk.real,zk.imag,color='r')
	#step 1 :prediction
        tk=0.1
        cp=0
        while(abs(zk-z1new)>tk or cp<10):
            u,sigma,vt = np.linalg.svd((zk)*np.eye(taille)-A)
            umin=u[:,len(u)-1]
            vmint=vt[len(vt)-1]
            vmin=np.conjugate(vmint)
            if(np.vdot(vmin,umin) == 0):
                print ('erreur')
                return
            rk=(1j*np.vdot(vmin,umin))/abs(np.vdot(vmin,umin))
            tk=min(0.1,0.5*abs(zk-eig))
            z1=zk+rk*tk
            u,sigma,vt = np.linalg.svd((z1)*np.eye(taille)-A)
            umin=u[:,len(u)-1]
            vmint=vt[len(vt)-1]
            vmin=np.conjugate(vmint)
            zk=z1-((min(sigma)-e)/np.vdot(umin,vmin))
            plt.scatter(zk.real,zk.imag, color='g')
            cp=cp+1
    bar.update()
    time.sleep(0.1)
    i=i+1   
    plt.show()




###########################################################################
#Pseudospectre par composante
###########################################################################
    
def rho2(M):
    taille=len(M)
    s,v=np.linalg.eig(M)
    return max(abs(s))


def ComputRadius(A,e,pas):    
    taille=A.shape[0]
    L=locaPseudoVp(A,e)   
    x=np.arange(L[1],L[0]+pas,pas)    
    y=np.arange(L[3],L[2]+pas,pas)
    vpx=np.linalg.eig(A)[0].real
    vpy=np.linalg.eig(A)[0].imag
    lx=x.shape[0]
    ly=y.shape[0]
    sigmin=np.zeros((ly,lx))
    XX, YY=np.meshgrid(x,y)
    if __name__ == '__main__':
        i = 0
        bar = ProgressBar(100)
        bar.update()
        time.sleep(0.1)
        i=i+1
    k=0
    progression=lx/100 +1
    for i in range (lx):
        k=k+1
        if(k%progression<0.2 and bar.perc<98 and k!=1):
            bar.update()
            time.sleep(0.1)
            i=i+1
        for j in range(ly):
            AM=abs(np.linalg.inv(A-(x[i]+1j*y[j])*np.eye(taille)))
            sigmin[j][i]=rho2(AM*abs(A))
    plt.title('pseudospectre')
    plt.ylabel('Y-partie imaginaire') 
    plt.xlabel('X-partie reelle')
    plt.plot(vpx,vpy,"ob") 
    c=plt.contour(XX, YY,sigmin,[e])
    bar.update()
    time.sleep(0.1)
    i=i+1   
    plt.show()
    plt.close()






########################################
#INTERFACE
########################################

    
def prendreE():
    e=float(e1.get())
    pas=float(p.get())  
    print(F)
    ch=int(choix.get())
    if(ch==2):
        dessinGrid(F,e,pas)
    elif(ch==1):
        prediCorrect(F,e,pas)
    else:
        ComputRadius(F,e,pas)

    
def quitter():
    plt.close()


def connaitre():
    tt=int(e2.get())
    return tt



def nvdessin():
    xmi,xma,ymi,yma=plt.axis()
    npas=float(nvpas.get())
    ne=float(nve.get())
    
    plt.close()
    nvdessin4(F,yma,ymi,xma,xmi,npas,ne)

"""

"""

def tracerOpti():
    e=float(e1.get())
    pas=float(p.get())
    #F=creer_D([1-5*1j, 15-7*1j, 30-3*1j, 10+4*1j])
    #F=creer_D([2-5*1j, 1-7*1j, 2-3*1j, 4*1j,1,-2+1*1j, 30+20*1j,29+20*1j,15+10*1j])
    F=creer_D([2-5*1j, 1-7*1j, 2-3*1j, 4*1j,1,-2+1*1j, 100+20*1j,209+20*1j,15+10*1j,208+19*1j])
    dessinOpti(F,e,pas)
    
        
#####################################################################
    #INTERFACE
#####################################################################

fenetre = Tk.Tk()
fenetre.title("pseudospctre")


iTailleM=Tk.Label(fenetre,text="chosir la taille de la matrice:");
iTailleM.pack()
e2=Tk.Entry()
e2.pack()

b1 = Tk.Button(text='VALIDER', command=connaitre)
b1.pack()

bouton=Tk.Button(fenetre, text="CONTINUER", command=fenetre.quit)
bouton.pack()

fenetre.mainloop()

tailleM=connaitre()
print(tailleM)
#F=creer_M(tailleM)
F=creer_D([1-5*1j, 15-7*1j, 30-3*1j, 10+4*1j])
print("voici F:", F)

#################################################################


l=Tk.Label(fenetre,text="\nchoisir un epsilon >10 PUIS ENSUITE RETRACER AVEC e<<1");
l.pack()
e1 = Tk.Entry()
e1.pack()


l2=Tk.Label(fenetre,text="choisir le pas \n(le pas doit etre <=epsilon");
l2.pack()
p= Tk.Entry()
p.pack()



l3=Tk.Label(fenetre,text="\nchoisir 1:Prediction Correction,  2: algo GRID,    3:Pseudo par composante ");
l3.pack()
choix= Tk.Entry()
choix.pack()


b1 = Tk.Button(text='Valider', command=prendreE)
b1.pack()

z8=Tk.Label(fenetre,text="______________________________________________________________\nTracé iptimisé");
z8.pack()
"""
choix_opti=Tk.Entry()
choix_opti.pack()
"""
opti = Tk.Button(text='tracé optimisé', command=tracerOpti)
opti.pack()

z7=Tk.Label(fenetre,text="\n le nouveau epsilon");
z7.pack()

nve=Tk.Entry()
nve.pack()


z6=Tk.Label(fenetre,text="\nle nouveau pas (<= epsilon )");
z6.pack()

nvpas=Tk.Entry()
nvpas.pack()

z = Tk.Button(text='retracer', command=nvdessin)
z.pack()


b1 = Tk.Button(text='le pseudospectre de depart', command=prendreE)
b1.pack()

bouton=Tk.Button(fenetre, text="Fermer", command=fenetre.quit)
bouton.pack()

fenetre.mainloop()




#########################
#B=creer_A(5)
B=np.array([[3.+4.*1j, 3.+2.*1j, 1.+3.*1j, 0.+4.*1j, 0.+2.*1j], [3.+1.*1j, 0.+1.*1j, 3.+3.*1j, 2.+1.*1j, 3.+0.*1j],[1.+2.*1j, 1.+1.*1j, 1.+3.*1j, 0.+4.*1j, 0.+1.*1j],[2.+1.*1j, 4.+0.*1j, 4.+3.*1j, 3.+4.*1j, 3.+2.*1j],[2.+0.*1j ,3.+4.*1j ,3.+4.*1j ,4.+1.*1j ,2.+4.*1j]]) 



print(locaOPTIv2(B))
#print(locaOPTIv2(M))
#D=creer_diag([1-5.*1j, 150-300.*1j, 300-400.*1j, 1000+40.*1j],4)
#print(locaOPTIv2(D))

#KK=np.diag([1-5.*1j, 150-300.*1j, 300-400.*1j, 1000+40.*1j])
#print(locaOPTIv2(KK))

