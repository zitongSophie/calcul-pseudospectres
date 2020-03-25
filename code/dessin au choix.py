import random
import numpy as np
import numpy.matlib
from matplotlib import pyplot as plt
import tkinter as Tk
def creer_A(taille):
    A1=np.random.random(size=(taille,taille))*10-5
    A2=np.random.random(size=(taille,taille))*10-7
    for i in range(taille):
        for j in range(taille):
            if(random.random()<0.7):
                A1[i][j]=0
                A2[j][i]=0
    A=np.array(A1+A2*1j)
    return A
def loca_vp(A,taille):
    for i in range(taille):
        som=0
        for j in range(taille):
                if(j!=i):
                        som+=abs(A[i,j])  #calculer du bord au fur et a mesure du parcours dela matrice
        if(i==0):
            xsup= A[0,0].real + som
            xinf= A[0,0].real - som
            ysup= A[0,0].imag + som
            yinf= A[0,0].imag - som
        else:
            if(A[i,i].real+som > xsup):
                xsup= A[i,i].real + som
            if(A[i,i].real-som < xinf):
                xinf= A[i,i].real - som
            if(A[i,i].imag+som > ysup):
                ysup= A[i,i].imag + som
            if(A[i,i].imag-som<yinf):
                yinf= A[i,i].imag - som
    return [xsup,xinf,ysup,yinf]  

def locaPseudoVp(A,e,taille):
    l=loca_vp(A,taille) #appel de la fonction lacaliser spectre 
    r=[0,0,0,0] #modifier le bord contenant les preseudovaleurpropre
    r[0]=l[0]+e*taille
    r[1]=l[1]-e*taille
    r[2]=l[2]+e*taille
    r[3]=l[3]-e*taille
    return r
def locaOPTI(A,e,taille):
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
    for d1 in disque:    # verifier la relation adjoint pour chaque d1 avec d'autres 
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
def dessin(A,borne,e,pas,taille): #borne=[xsup,xinf,ysup,yinf]
    x=np.linspace(borne[1],borne[0],num=int((borne[0]-borne[1])/pas))    
    y=np.linspace(borne[3],borne[2],num=int((borne[2]-borne[3])/pas)) #creation de la grid
    lx=x.shape[0]
    ly=y.shape[0]
    XX, YY= np.meshgrid(x,y) #initialisation d une matrice
    sigmin=np.zeros((ly,lx))
    for m in range(lx):
        for n in range(ly):
            u,sigma,vt = np.linalg.svd((x[m]+1j*y[n])*np.eye(taille)-A)
            sigmin[n,m]=sigma[taille-1]
    print (1)
    plt.contour(XX, YY,sigmin,[e])
def dessinOpti(A,e,pas,mode):#mode=true:opti,mode=false:pas opti
    taille=A.shape[0]
    t=len(A[0])#verifions si la matrice est carre et est bien inversible
    if(t!=taille):
        print ("erreur matrice non carre")
        return
    if(np.linalg.det(A)==0):
        print ("erreur matrice non inversible")
        return 
    vpx=np.linalg.eig(A)[0].real
    vpy=np.linalg.eig(A)[0].imag
    plt.plot(vpx,vpy,"ob")
    plt.title('pseudospectre')
    plt.ylabel('Y-partie imaginaire') 
    plt.xlabel('X-partie reelle')
    if(mode):
        L=locaOPTI(A,e,taille)
        for l in L:
            dessin(A,l,e,pas,taille)
    else:
        dessin(A,locaPseudoVp(A,e,taille),e,pas,taille) #dessin(A,borne,e,pas,taille),borne=[xsup,xinf,ysup,yinf]


#f=np.array([[2-1*1j,1+0*1j,1+0*1j,1+0*1j,1+0*1j],[2+3*1j,10+10*1j,2+0*1j,2+0*1j,2+0*1j],[3+0*1j,0+3*1j,20+20*1j,1+3*1j,2+2*1j],[0+0*1j,0+0*1j,0+0*1j,30+0*100j,0+0*1j],[0+0*1j,0+0*1j,0+0*1j,0+0*1j,0+40*1j]])
#f=creer_A(10)
#dessin pour tracer le courbe dans une rectangle mais il fait pas de show
#dessinopti cest  qu on peut chosir le mode opti ou pas
#j ai essaye ,on peut faire show dehors donc
#je sais pas si c est mieux de mettre show dans la fonction ou on le faire nous meme apres 
f=np.array([[1.00001+1j,0,0,30],[0,1.0002+2j,0,0],[0,0,0.9999+3j,0],[0,0,0,0.999+4j]])
dessinOpti(f,0.1,0.1,True)
plt.show()
'''
def prendreE():
    epsilon=float(e.get())
    pas=float(p.get())
    taille=int(t.get())
    A=creer_A(taille)
    dessinOpti(A,epsilon,pas,True)
def nvdessin():
    xmi,xma,ymi,yma=plt.axis()
    npas=float(nvpas.get())
    ne=float(nve.get())  
    plt.close()

    taille=A.shape[0]
    t=len(A[0])#verifions si la matrice est carre et est bien inversible
    if(t!=taille):
        print ("erreur matrice non carre")
        return
    if(np.linalg.det(A)==0):
        print ("erreur matrice non inversible")
        return 
    vpx=np.linalg.eig(A)[0].real
    vpy=np.linalg.eig(A)[0].imag
    plt.plot(vpx,vpy,"ob")
    plt.title('pseudospectre')
    plt.ylabel('Y-partie imaginaire') 
    plt.xlabel('X-partie reelle')
    dessin(b1.get(),[xma,xmi,yma,ymi],ne,npas,)

fenetre = Tk.Tk()
fenetre.title("pseudospctre")


l=Tk.Label(fenetre,text="chosir la taille de la matrice:\nVeuillez commencer par une matrice de taille<10 pour commencer");
l.pack()
t=Tk.Entry()
t.pack()

l2=Tk.Label(fenetre,text="choisir un epsilon\n veuillez choisir un epsilon>=0.1 pour commencer\nAttention si epsilon<0.1 le temps d'attente sera un peu long~5min)\n");
l2.pack()
e=Tk.Entry()
e.pack()
l3=Tk.Label(fenetre,text="choisir le pas \n(le pas doit etre <=epsilon \n");
l3.pack()
p=Tk.Entry()
p.pack()
b=Tk.Button(text='Valider', command=prendreE)
b.pack()
bouton=Tk.Button(fenetre, text="\n vous avez zoomer \n ",command=fenetre.quit)
bouton.pack()
z5=Tk.Label(fenetre,text="\n");
z5.pack()
fenetre.mainloop()
z7=Tk.Label(fenetre,text="\n le nouveau epsilon \n");
z7.pack()
nve=Tk.Entry()
nve.pack()
z6=Tk.Label(fenetre,text="\nle nouveau pas (<= epsilon )\n");
z6.pack()
nvpas=Tk.Entry()
nvpas.pack()
z = Tk.Button(text='retracer', command=nvdessin)
z.pack()
b1 = Tk.Button(text='le pseudospectre de depart\n', command=prendreE)
b1.pack()
bouton=Tk.Button(fenetre, text="Fermer", command=fenetre.quit)
bouton.pack()

fenetre.mainloop()
'''



