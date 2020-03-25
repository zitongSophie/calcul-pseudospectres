import random
import numpy as np
import numpy.matlib
from matplotlib import pyplot as plt
from math import *
'''def creer_B(taille):
	b=np.random.random(taille)*500
	B=
'''
def creer_A(taille):
	A1=np.random.random(size=(taille,taille))*10-5
	A2=np.random.random(size=(taille,taille))*50-12

	A=np.array(A1+A2*1j)
	print (A)
	return A

def sigmin(z,A):
	taille=len(A)
	sigma = np.linalg.svd(z*np.eye(taille)-A,compute_uv=False) #svd pour calcuer singular value(sigma) et left and rignt singular value vector u,v;
	return sigma[taille-1]

def prediCorrect(A,e):  #suivi le documentaire page 446 et 447
	eigv=np.linalg.eigvals(A)
	z=[]
	taille=len(A)
	tol=0.0001
	#calculer pour chaque valeur propre
	for eig in eigv:

		plt.scatter(eig.real,eig.imag,color='g')
               #step 0:initialisation du premier point z1
		ps=[] #liste des points pour tracer pseudospectre
		d=1
		z1new=eig+e*d
		while((abs(np.linalg.svd(z1new*np.eye(taille)-A,compute_uv=False)[taille-1]-e) >tol*e) ):   #correction for first point
			z1old=z1new
			u,sigma,vt = np.linalg.svd(z1old*np.eye(taille)-A)
			umin=u[:,len(u)-1]
			vmint=vt[len(vt)-1]  
			vmin=np.conjugate(vmint)
			sigmi=np.linalg.svd(z1new*np.eye(taille)-A,compute_uv=False)[taille-1]
			z1new=eig-d*((sigma[taille-1]-e)/((-np.conjugate(d)*np.vdot(vmin,umin)).real))   #utiliser le theoreme de newton question:la difference pour umin et umin*	
		z1=z1new
		ps.append(z1new)
	#for k=2,3...
		zk=z1
		print (z1)
		plt.scatter(z1.real,z1.imag,color='r')
	#step 1 :prediction
		tk=0.1
		while(abs(zk-z1new)>tk or len(ps)<10):
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
			ps.append(zk)
			plt.scatter(zk.real,zk.imag)
		z.append(ps)
	plt.show() 	
	return z

A=np.array([[2-1*1j,1+0*1j,1+0*1j,1+0*1j,1+0*1j],[2+3*1j,10+10*1j,2+0*1j,2+0*1j,2+0*1j],[3+0*1j,0+3*1j,20+20*1j,1+3*1j,2+2*1j],[0+0*1j,0+0*1j,0+0*1j,30+0*100j,0+0*1j],[0+0*1j,0+0*1j,0+0*1j,0+0*1j,0+40*1j]])
#A=creer_A(5)
b=np.array([[1.00001+1j,0,0,30],[0,1.0002+2j,0,0],[0,0,0.9999+3j,0],[0,0,0,0.999+4j]])
prediCorrect(b,0.1)



