import random
import numpy as np
import numpy.matlib
from matplotlib import pyplot as plt
from math import *
def creer_A(taille):
        A1=np.random.random(size=(taille,taille))*10-5
        A2=np.random.random(size=(taille,taille))*50-12

        A=np.array(A1+A2*1j)
        return A
#cet fonction (au dessous) permettre de trouver  y  telque (x+iy)-A
def hamiltonian(x,A,taille,e):
        y=[]
        dsum=np.zeros((taille*2,taille*2),dtype=complex)
        a=np.array(x*np.eye(taille)-np.transpose(np.conjugate(A)))
        b=e*np.eye(taille)
        c=A-x*np.eye(taille)
        dsum[:a.shape[0],:a.shape[1]]=a
        dsum[a.shape[0]:,:a.shape[1]]=b
        dsum[:a.shape[0],a.shape[1]:]=-b
        dsum[a.shape[0]:,a.shape[1]:]=c
        eigs=np.linalg.eigvals(dsum)
        for eig in eigs:        
                if(abs(eig.real)<0.001):
                                y.append(eig.imag)
        return y
def abscissa(A,e):
#etape 1: chercher le valeur propre plus au droite(righteig)
        taille=len(A)
        eigv=np.linalg.eigvals(A)
        for eig in eigv:
                if(eig.real==max(eigv.real)):
                        righteig=eig
#chercher le pseudovaleurpropre plus au droite pour righteig
        d=1
        z1new=righteig+e*d
        print ('righteig',righteig)
        tol=0.1
        while((abs(np.linalg.svd(z1new*np.eye(taille)-A,compute_uv=False)[taille-1]-e) >tol*e) or z1new-righteig<0 ):   #correction for first point
                z1old=z1new
                u,sigma,vt = np.linalg.svd(z1old*np.eye(taille)-A)
                umin=u[:,len(u)-1]
                vmint=vt[len(vt)-1]  
                vmin=np.conjugate(vmint)
                sigmi=np.linalg.svd(z1new*np.eye(taille)-A,compute_uv=False)[taille-1]
                z1new=righteig-d*((sigma[taille-1]-e)/((-np.conjugate(d)*np.vdot(vmin,umin)).real))   #utiliser le theoreme de newton question:la difference pour umin et umin* 
        z1=z1new
        print ('z1',z1)
        zk=z1
        k=0
        while(k==0 or (max(X)-min(X))>0.0001):#convergence
                                                #quand le plus grand et le plus petit sont proche
                Y=[]
                k=k+1
                print ('k=',k,hamiltonian(zk.real,A,taille,e))
                eigs=hamiltonian(zk.real,A,taille,e)
                for eig in eigs:
                        if((np.linalg.svd((zk.real+eig*1j)*np.eye(taille)-A,compute_uv=False)[taille-1])-e<0.00001):
                                Y.append(eig)
                
                Y.sort()
        #trouver midpoint
                midpoint=[]
                for i in range(1,len(Y),2):
                        midpoint.append((Y[i-1]+Y[i])/2)
                X=[]
                print('midpoint',midpoint)
                for y in midpoint:
                        eigs=hamiltonian(-y,A*1j,taille,e)
                        print('eigs',eigs)
                        X.append(max(eigs)+1j*y)
                zk=X[0]
                print('x',X)
                for i in range(1,len(X)):
                        if(zk.real<X[i].real):
                                        zk=X[i]
                
                print('zk',zk)
        print ('abscissa',zk)
        return zk.real

#b=np.array([[2-1*1j,1+0*1j,1+0*1j,1+0*1j,1+0*1j],[2+3*1j,10+10*1j,2+0*1j,2+0*1j,2+0*1j],[3+0*1j,0+3*1j,20+20*1j,1+3*1j,2+2*1j],[0+0*1j,0+0*1j,0+0*1j,30+0*100j,0+0*1j],[0+0*1j,0+0*1j,0+0*1j,0+0*1j,0+40*1j]])
b=creer_A(10)
#b=np.array([[1.00001+1j,0,0,30],[0,1.0002+2j,0,0],[0,0,0.9999+3j,0],[0,0,0,0.999+4j]])
abscissa(b,0.1)

#print(hamiltonian(1.1,b,4,0.1))
'''
	while(i<len(Y)-1):
		if(Y[i]==Y[i+1]):
			interval.add(Y[i])
			i=i+2
		else:
			if(Y[i] not in interval ):
				interval.add(Y[i],Y[i+1])
				i=i+1
			i=i+1



                midinterval=[]
                i=0
                while(i<len(Y)-1):
                        if(Y.count(Y[i])>1):
                                midinterval.append(Y[i])
                                i=i+Y.count(Y[i])
                        else:
                                midinterval.append((Y[i]+Y[i+1])/2)
                                i=i+2
                X=[]
                for y in midinterval:
                        eigs=hamiltonian(-y,A*1j,taille,e)      
                        X.append(max(eigs)+1j*y)
                zk=X[0] 
                for i in range(1,len(X)):
                        if(zk.real<X[i].real):
                                        zk=X[i]


'''


