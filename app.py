#!/usr/bin/env python
# -*- coding: utf-8 -*- 
from __future__ import unicode_literals
from __future__ import division#garantit qu'il n'y aura pas de problème avec les convention de division euclidienne/entière et division standard entre les version 2 et 3 de python
import os
from math import *#contient sqrt,exp,log
from random import *#contient uniform 
import matplotlib as mpl#tracers généraux
from mpl_toolkits.mplot3d import Axes3D#tracer 3D
import numpy as np#fonction facilitant la manipulation de listes de listes
import matplotlib.pyplot as plt#tracer 2D
from sys import version# nous dit quelle version de python est utilisée
if version[0]=='3':from tkinter import *#importation de la bibliothèques de construction d'interface graphique : 'tkinter' sur python3 et 'Tkinter' sur le 2
else:from Tkinter import *
from platform import system#Nous dit dans quel os le script est interprété

global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant#toutes les variables global color, plmax,es utilisables par toutes les fonctions du programme
kb=1.38064852*10**(-23)#déclaration des constantes
mproton=1.672622*10**(-27)
plmax=[[],-1,-1,-1]#listemax/T/M/Vini
color=[0,['blue','green','cyan','purple','red','black','orange','magenta','aliceblue','antiquewhite','aqua','aquamarine','azure','beige','brown','deeppink','seashell','sienna','silver','skyblue','slateblue','slategray','snow','springgreen','steelblue','tan','tea','thistle','tomato','turquoise','violet','wheat','white','whitesmoke''yellow''yellowgreen'],0]#rang couleur courante,liste couleurspour tracer 2D,indicateur impéchant de légender deux courbes de même couleur différemment.
def vdeT(T,M):#retourne à une température et une masse données : [vitesse moyenne des particules de ce type de gaz,vitesse la plus probable des particules dans ce gaz
    return([sqrt(8/np.pi*kb*T/M),sqrt(2*kb*T/M)])
def Tdev(v,M):#retourne à une viteese et une massse donnée :la Température du gaz dont les particules vont en moyenne à cette vitesse
    return(np.pi*M*v*v/kb/8)
def listegauss(s):#création de la grosse liste pour tirage gaussien
    if sigma==0 or sigma=='random':return [0]#cas du sigma=0 : rebond miroir, sigma=random : on le fixe à 0 arbitrairement mais le cas random est directement traité dans les méthodes 'randsaut' des classe particule et particuleext
    gauss=list()#on déclare la liste gauss
    K=1.01/exp(-(s*3)**2/(2*s**2))#on fixe la valeur max à représenter dans la liste à 3*sigma
    for i in range(-5000,5000):#on aura un pas de 6*sigma/10000
        I=i*s*3/5000#valeur qui va être écrite dans la liste
        gauss+=int(exp(-(I)**2/(2*s**2))*K)*[I]#pondération de cette valeur par un coef proportionnel à la densité de probabilité gaussienne
    
    return gauss
def listemax(T,M):#création de la grosse liste pour tirage dans une maxwellienne ou dans une seule vitesse fixe prédéfinie
    
    if plmax[0]==[] or plmax[1]!=T or plmax[2]!=M or plmax[3]!=Vini:
        plmax[1],plmax[2],plmax[3]=T,M,Vini
        if Vini==0:#si on choisi de définir la particule par une température
            v=vdeT(T,M)
            plmax[0]=list()
            K=1.01/(((vdeT(T,M)[1]*3)**2)*exp(-M*(vdeT(T,M)[1]*3)**2/(2*kb*T)))
            for i in range(0,10000):
                I=i*vdeT(T,M)[1]*3/10000
                plmax[0]+=int((I**2)*exp(-M*(I)**2/(2*kb*T))*K)*[I]
            
            return(plmax[0])
        return([Vini])#si on choisi de définir la particule par une vitesse.
    else:return(plmax[0])
def listelpm(lpmm):#liste pour piocher une distance avant le prochain contact.
    if lpmm==0:return([0])#Pas de contact
    l=list()
    for i in range(0,10000):
        I=i*(-lpmm*log(0.01))/10000
        l+=int(exp(-I/lpmm)*101)*[I]
        
    return l
def tirageliste(l):#réalise un tirage equiprobable dans une liste
    if l==[]:return(0)
    return l[randint(0,len(l)-1)]
def coefperte(v):#fonction qui retourne le coef de perte en vitesse en fonction de la vitesse actuelle de la particule. Ici on retourne directement la valeur entrée par l'utilisateur mais cette fct est créée pour être adaptée, exemple return(0.0011/v) ou un truc du genre, c'est libre et prêt à être modifié.
    return(perteQ[0])
def th(x):#tangente hyperbolique (inutilisé mais disponible ! (peut apparaitre dans certaines probas)
    return((exp(x)-exp(-x))/(exp(x)+exp(-x)))

def P(v,M):#retourne la probabilité de se coller en fonction d'une vitesse et d'une masse.

    if sticklist[0]==10:#dans ce cas on utilise la loi : https://tel.archives-ouvertes.fr/tel-00438534/document
        
        return(1+(Tdev(v,M)/53/M*mproton))*exp(-(Tdev(v,M)/53/M*mproton))
    else :#dans ce cas on utilise la proba entrée par l'utilisateur
       
        return sticklist[0]



def gaussiennetest():
    lpmm=10**(-6)
    liste=[[i*10*(-lpmm*log(0.05))/10000 for i in range(0,1000)],[0 for i in range(0,1000)]]
    
    l=listelpm(lpmm)
    for i in range(0,100000):
        x=tirageliste(l)
        liste[1][int(-x/10/lpmm/log(0.05)*10000)]+=1
    #liste[1][500]*=0.5 #le int arrondi les 0.0X et les -0.0X à 0 donc il se mange un x2 pas désiré :)
    
    plt.plot(liste[0],liste[1])
    plt.show(block=True)

def ezw(l):#prends une liste l=[a,b,c,d...] et renvoi la chaine de caractère : "a b c d..." (pas d'espaces sur les côtés
    ret=""
    for a in l[0:len(l)-1]:
        ret+=str(a)+" "
    ret+=str(l[len(l)-1])
    return(ret)


def rand11():#renvoi -1 ou 1 equiprobablement
    x=randint(0,1)
    if x==0:x=-1
    return x
def scal(u,v):#u et v des listes de taille 3 : renvoi un float de valeur u produit scalaire v
    return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2])

def multscal(x,u):#retourne la liste de taille 3 correspondant à la multiplication des éléments de la liste de taille 3 u par le flottant x
    return([x*u[0],x*u[1],x*u[2]])
def somvect(u,v):#retourne la liste de taille 3  correspondant à la somme 2 à 2 des éléments des listes de taille 3 u et v
    return([v[0]+u[0],v[1]+u[1],v[2]+u[2]])
def norme(u):#retourne le flottant correspondant à la norme de la liste de taille 3 u
    return sqrt(scal(u,u))
    
def normale(pos):#retourne la normale à la surface (liste unitaire de taille 3)
    if natureQ==1:#cas du cylindre
        if pos[2]<=0 and norme(pos)==r: return(multscal(1/sqrt(2),somvect([0,0,1],multscal(1/norme([-pos[0],-pos[1],0]),[-pos[0],-pos[1],0]))))#cas où on est sur un côté
        elif pos[2]<=0 :return([0,0,1])#cas où on est en bas
        else :return(multscal(1/norme([-pos[0],-pos[1],0]),[-pos[0],-pos[1],0]))
    else :#cas de la sphère
        return(multscal(-1/norme(pos),pos))
def randnormale(pos,vv):#renvoi la normale modifiée
    if sigma==0:return normale(pos)#cas du rebond miroir
    normnorm=normale(pos)
    normnorm[0]=normnorm[0]+tirageliste(gaussliste)
    normnorm[1]=normnorm[1]+tirageliste(gaussliste)
    normnorm[2]=normnorm[2]+tirageliste(gaussliste)
    return (multscal(1/norme(normnorm),normnorm))#cas d'un sigma =/=0
            
    
def tirage(p):#renvoi le booléen True avec une probabilité p
    return( randint(0,1000)/1000 <=p)

def txtipt(nom):#fonction qui va lire les datas ! elle importe les datas du fichier nom.txt en excluant l'en-tête et les sauvegarde dans une liste qu'elle retourne.
    l=list()
    f=open(os.path.join(cheminpy,"Res/"+nom+".txt"),'r')
    ff=f.read()
    f.close()
    pos=0
    while ff[pos]+ff[pos+1]+ff[pos+2]+ff[pos+3]!="ou?\n":pos+=1
    pos+=4
    
    while pos<len(ff):
        texte=str()
        ligne=list()
        encours=1
        while pos<len(ff) and ff[pos]!="\n":
            
            b=ff[pos]
        
            if b==" " and encours==1 and ff[pos+1]!='\n':
                encours=0
                if float(texte)==int(float(texte)):ligne.append(int(float(texte)))
                else :ligne.append(float(texte))
                texte=str()
            if b!=" " and b!="\n" :
                encours=1
                texte=texte+b
                
            pos+=1 
        if float(texte)==int(float(texte)):ligne.append(int(float(texte)))
        else :ligne.append(float(texte))
        l.append(ligne)
        pos+=1
    return (l)
def txtiptet(nom):#importe l'entête du fichier de datas dans une liste et la retourne
    f=open(os.path.join(cheminpy, "Res/"+nom+".txt"),'r')
    ff=f.read()
    f.close()
    pos=0
    
    texte=str()
    ligne=list()
    encours=1
    while ff[pos]!="\n":
        
        b=ff[pos]
        
        if b==" " and encours==1:
            encours=0
            ligne.append(texte)
            texte=str()
        if b!=" " :
            encours=1
            texte=texte+b
            
        pos+=1 
    
    ligne.append(texte)
    
    return (ligne)  

def moyenneur(nom):#traite les datas en fonction de leur type (cf guide) et retourne les données moyennées (à afficher par la fonction grapheur sous forme de données moyennées simples ou de graphe 2D à l'utilisateur, ou encore de graphe 3D de trajectoires)
    l=txtipt(nom)

    nbcol=len(l[0])
    nblignes=len(l)

    if nbcol==3:
        comptes=[0,0,0]
        for a in l:
            comptes[0]+=a[0]
            comptes[1]+=a[1]
            comptes[2]+=a[2]
        return(multscal(1/nblignes,comptes))
    elif nbcol==4  :
        res=list()
        comptes=[0,0,0]
        para=l[0][0]
        nbsimulations=0
        for i in range(0,nblignes):
            para=l[i][0]
            a=l[i][1:]
            comptes[0]+=a[0]
            comptes[1]+=a[1]
            comptes[2]+=a[2]
            nbsimulations+=1
            if i==nblignes-1 or para!=l[i+1][0]:
                res.append([para]+multscal(1/nbsimulations,comptes))
                comptes=[0,0,0]
                nbsimulations=0
        return res

    elif int((nbcol-3)/3)==float((nbcol-3)/3) :
    
        res=list()
        comptes=[0,0,0]
        x,y,z=list(),list(),list()
        xx,yy,zz=list(),list(),list()
        nbsimulations=0
        for i in range(0,nblignes):
            
            a=l[i][0:]
            nbcol=len(a)
            comptes[0]+=a[nbcol-3]
            comptes[1]+=a[nbcol-2]
            comptes[2]+=a[nbcol-1]
            
            for k in range(0,int((nbcol-3)/3)):
                x.append(a[k])
                y.append(a[k+int((nbcol-3)/3)])
                z.append(a[k+2*int((nbcol-3)/3)])
            xx.append(x)
            x=list()
            yy.append(y)
            y=list()
            zz.append(z)
            z=list()
            
            nbsimulations+=1
            if i==nblignes-1:
                res.append(multscal(1/nbsimulations,comptes))
                
        res.append(xx)
        res.append(yy)
        res.append(zz)
        res.append(['t',nbsimulations])
    
        return res
    else:
        
        res=list()
        comptes=[0,0,0]
        para=l[0][nbcol-4]
        x,y,z=list(),list(),list()
        xx,yy,zz=list(),list(),list()
        nbsimulations=0
        nbpara=0
        for i in range(0,nblignes):
            
            a=l[i][0:]
            nbcol=len(a)
            para=a[nbcol-4]
            comptes[0]+=a[nbcol-3]
            comptes[1]+=a[nbcol-2]
            comptes[2]+=a[nbcol-1]
            
            for k in range(0,int((nbcol-3)/3)):
                x.append(a[k])
                y.append(a[k+int((nbcol-3)/3)])
                z.append(a[k+2*int((nbcol-3)/3)])
            xx.append(x)
            x=list()
            yy.append(y)
            y=list()
            zz.append(z)
            z=list()
            
            nbsimulations+=1
            
            if i==nblignes-1 or para!=l[i+1][len(l[i+1][0:])-4]:
                res.append([para]+multscal(1/nbsimulations,comptes))
                comptes=[0,0,0]
                nbsimulations=0
                nbpara+=1
            
        res.append(xx)
        res.append(yy)
        res.append(zz)
        res.append(['tt',nbpara])
        return(res)
def edit():#fonction qui permet d'éditer les parties 'cosmétiques' d'un graphe 2D
    if EE[12].get()!='':para=EE[12].get()
    else:para="T"
    if Ee1.get()!='':nomm=Ee1.get()
    else:nomm="temporaire"
    if Ee0.get()!='':legabs=Ee0.get()
    else:legabs=''
    if EE[6].get()!='':titre=EE[6].get()
    else:titre=nomm
    plt.title(titre)
    if legabs=='':plt.xlabel(para)
    else:plt.xlabel(legabs)
    
def grapheur():#S'occupe d'afficher les résultats selon des paramètres fournis ! cf guide
    
   
    if EE[12].get()!='':para=EE[12].get()
    else:para="T"
    if Ee1.get()!='':nomm=Ee1.get()
    else:nomm="temporaire"
    if Ee2.get()!='':
        if Ee2.get()[0]=='h':
            rhauteur=1
            iii=int(Ee2.get()[1])
        else:iii=int(Ee2.get())
    else:
        rhauteur=1
        iii=3
    if Ee3.get()!='':legg=Ee3.get()
    else:legg='rien'
    if Ee9.get() in ['','1']:
        color[0]+=1
        color[2]=0
    elif color[2]==1:legg='rien'
    if Ee0.get()!='':legabs=Ee0.get()
    else:legabs=''
    if EE[6].get()!='':titre=EE[6].get()
    else:titre=nomm
    
    nom=nomm
    i=iii
    leg=legg
    l=moyenneur(nom)
    lipt=txtipt(nom)


    if (type(l[0])!=list and l==lipt[0]) or (type(l[0])==list and l[len(l)-1][1]==1 and iii!=0):
        if type(l[0])==list:l=l[0]
        Lres["text"]="Temps de résidence dans la cavité :\n{0} secondes\nNombre de contacts :\n{1} contacts\n {2}".format(l[0],l[1],["La particule n'est pas sortie",'La particule est sortie'][int(l[2])])
        
    elif type(l[0])!=list or (l[len(l)-1][0]=='t' and iii!=0):
        if type(l[0])==list:l=l[0]
        Lres["text"]="Temps de résidence dans la cavité moyen :\n{0} secondes\nNombre de contacts moyen :\n{1} contacts\nPart des particules sorties :\n{2}%".format(l[0],l[1],l[2]*100)
        
    elif iii==0 and (l[len(l)-1][0]=='tt' or l[len(l)-1][0]=='t'):
        linfos=txtiptet(nom)
        pos=0
        while linfos[pos]!='natureQ':
            pos+=1
        nat=int(linfos[pos+1])
        while linfos[pos]!='mode2':
            pos+=1
        mmode2=int(linfos[pos+1])
        while linfos[pos]!='r':
            pos+=1
        rr=float(linfos[pos+1])
        while linfos[pos]!='RouL':
            pos+=1
        RR=float(linfos[pos+1])
       


    
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        if l[len(l)-1][0]=='tt':nbpart=len(l[len(l)-2])
        if l[len(l)-1][0]=='t':
            for i in range(0,l[len(l)-1][1]):ax.plot(l[1][i],l[2][i],l[3][i])
        else:
        
            for i in range(0,nbpart):
                ax.plot(l[l[len(l)-1][1]][i],l[l[len(l)-1][1]+1][i],l[l[len(l)-1][1]+2][i])
            
            
        if nat==0 :        
            u = np.linspace(0, 2*np.pi, 50)
            if mmode2==1:v = np.linspace(np.pi-asin(rr/RR), np.pi, 50)
            else:v = np.linspace(0+asin(rr/RR), np.pi, 50)
            ctheta=np.linspace(0,2*np.pi,50)
            xc=r*(np.cos(ctheta))
            yc=r*(np.sin(ctheta))
            if mode2==1:
                ax.plot([0],[0],-RR+2*rr)
                zc=-sqrt(RR*RR-rr*rr)
            else:
                ax.plot([0],[0],R)
                zc=sqrt(RR*RR-rr*rr)
            ax.plot(xc,yc,zc,color='green')
            x = RR * np.outer(np.cos(u), np.sin(v))
            y = RR * np.outer(np.sin(u), np.sin(v))
            z = RR * np.cos(v) #coordonnée sphère en cartésien : z=RcosO,x=rsinOsinP,y=RsinOcosP
            ax.plot_wireframe(x, y, z, rstride=7, cstride=7,color='grey')
            
            
        else:
            u = np.linspace(0, 2*np.pi, 50)
            v = np.linspace(0,RR, 50)
            ctheta=np.linspace(0,2*np.pi,50)
            xc=rr*(np.cos(ctheta))
            yc=rr*(np.sin(ctheta))
            zc=RR
            ax.plot(xc,yc,zc,color='green')
            x = rr *np.cos(u)
            y = rr *np.sin(u)
            if 2*rr>RR:
                ax.plot([0],[0],2*rr)
            else:
                ax.plot([RR/2],[0],0)
                ax.plot([0],[RR/2],0)
                ax.plot([-RR/2],[0],0)
                ax.plot([0],[-RR/2],0)
            deb=0
            if mmode2==1:
                deb=1
                ax.plot(xc,yc,0,color='green')
            
            for a in range(deb,10):
                pas=RR/10
                ax.plot(x,y,a*pas,color='grey')
    
        plt.show(block=True)
    elif (i in [1,2,3]) and l[len(l)-1][0]!='t':
        linfos=txtiptet(nom)
        pos=0
        
        while linfos[pos]!='mode2':
            pos+=1
        mmode2=int(linfos[pos+1])
        
        while linfos[pos]!='RouL':
            pos+=1
        RR=float(linfos[pos+1])
        while linfos[pos]!='para':
            pos+=1
        rhauteur=(linfos[pos+1]=='r')
        x=list()
        y=list()
        if (l[len(l)-1][0]=='tt'):
            for a in l[0:len(l)-4]:
                if rhauteur==1:
                    if mmode2==0:x.append(RR+sqrt(RR*RR-a[0]**2))
                    else:x.append(RR-sqrt(RR*RR-a[0]**2))
                else:x.append(a[0])
                if i==3:y.append(a[i]*100)
                else:y.append(a[i])
        else:
            for a in l:
                if rhauteur==1:
                    if mmode2==0:x.append(RR+sqrt(RR*RR-a[0]**2))
                    else:x.append(RR-sqrt(RR*RR-a[0]**2))
                else:x.append(a[0])
                if i==3:y.append(a[i]*100)
                else:y.append(a[i])
        if leg!='rien':
            color[2]=1
            plt.plot(x,y,label=leg,color=color[1][color[0]])
            plt.legend().draggable()
        else:
            plt.plot(x,y,color=color[1][color[0]])
        plt.title(titre)
        if legabs=='':plt.xlabel(para)
        else:plt.xlabel(legabs)
        plt.ylabel(['Temps passé dans la cavité','Nombre de contacts','Pourcentage de particules sorties'][i-1])  
        plt.gca().get_xaxis().get_major_formatter().set_powerlimits((0, 0))
        plt.show(block=True)
        
        
def grapheurpara(nomm):#s'occupe d'afficher les trajectoires à la suite de la simulation dans le cas d'une simulation de type 3) cf guide
    l=moyenneur(nomm)

    linfos=txtiptet(nomm)
    pos=0
    while linfos[pos]!='natureQ':
        pos+=1
    nat=int(linfos[pos+1])
    while linfos[pos]!='mode2':
        pos+=1
    mmode2=int(linfos[pos+1])
    while linfos[pos]!='r':
        pos+=1
    rr=float(linfos[pos+1])
    while linfos[pos]!='RouL':
        pos+=1
    RR=float(linfos[pos+1])
    while linfos[pos]!='RouL':
        pos+=1
    RR=float(linfos[pos+1])



    fig = plt.figure()
    ax = fig.gca(projection='3d')
    if l[len(l)-1][0]=='tt':nbpart=len(l[len(l)-2])
    if l[len(l)-1][0]=='t':
        for i in range(0,l[len(l)-1][1]):ax.plot(l[1][i],l[2][i],l[3][i])
    else:
        
        for i in range(0,nbpart):
            ax.plot(l[l[len(l)-1][1]][i],l[l[len(l)-1][1]+1][i],l[l[len(l)-1][1]+2][i])
        
        
    if nat==0 :        
        u = np.linspace(0, 2*np.pi, 50)
        if mmode2==1:v = np.linspace(np.pi-asin(rr/RR), np.pi, 50)
        else:v = np.linspace(0+asin(rr/RR), np.pi, 50)
        ctheta=np.linspace(0,2*np.pi,50)
        xc=r*(np.cos(ctheta))
        yc=r*(np.sin(ctheta))
        if mode2==1:
            ax.plot([0],[0],-RR+2*rr)
            zc=-sqrt(RR*RR-rr*rr)
        else:
            ax.plot([0],[0],R)
            zc=sqrt(RR*RR-rr*rr)
        ax.plot(xc,yc,zc,color='green')
        x = RR * np.outer(np.cos(u), np.sin(v))
        y = RR * np.outer(np.sin(u), np.sin(v))
        z = RR * np.cos(v) #coordonnée sphère en cartésien : z=RcosO,x=rsinOsinP,y=RsinOcosP
        ax.plot_wireframe(x, y, z, rstride=7, cstride=7,color='grey')
        
        
    else:
        u = np.linspace(0, 2*np.pi, 50)
        v = np.linspace(0,RR, 50)
        ctheta=np.linspace(0,2*np.pi,50)
        xc=rr*(np.cos(ctheta))
        yc=rr*(np.sin(ctheta))
        zc=RR
        ax.plot(xc,yc,zc,color='green')
        x = rr *np.cos(u)
        y = rr *np.sin(u)
        if 2*rr>RR:
            ax.plot([0],[0],2*rr)
        else:
            ax.plot([RR/2],[0],0)
            ax.plot([0],[RR/2],0)
            ax.plot([-RR/2],[0],0)
            ax.plot([0],[-RR/2],0)
        deb=0
        if mmode2==1:
            deb=1
            ax.plot(xc,yc,0,color='green')
        
        for a in range(deb,10):
            pas=RR/10
            ax.plot(x,y,a*pas,color='grey')
    plt.show(block=True)

def ain(): #fonction appellée lorsqu'on lance une simulation (en cliquant sur le bouton 'simulation' bleu cf guide)
    global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
    #on va ici importer tous les paramètres depuis les zones d'entrée de texte disponible dans l'interface et leur attribuer leur valeur par défaut si le champs est vide.
    if Ee8.get()!='':lpm=float(Ee8.get())*10**(-6)
    else:lpm=0
    if Ee7.get()!='':popv=int(Ee7.get())
    else:popv=0
    if Ee5.get()!='':Ts=float(Ee5.get())
    else:Ts=10
    if Ee6.get()!='':
        A=Ee6.get()
        sticklist=list()
        valeurtempo=str()
        
        for a in range(0,len(A)):
            if A[a]!=',' and A[a]!=';' and A[a]!=' ' and A[a]!='/':
                valeurtempo+=A[a]
            else:
                
                sticklist.append(float(valeurtempo))
                valeurtempo=str()
        sticklist.append(float(valeurtempo))        
        
    else:sticklist=[10]
    if Ee4.get()!='':mode2=int(Ee4.get())
    else:mode2=0
    if EE[18].get()!='':sigma=float(EE[18].get())
    else:sigma='random'
    if EE[0].get()!='':natureQ=int(EE[0].get())
    else:natureQ=1
    if EE[1].get()!='':r=float(EE[1].get())*10**(-6)
    else:r=1*10**(-6)
    if EE[2].get()!='':R=float(EE[2].get())*10**(-6)
    else:R=2*10**(-6)
    L=R
    if popv==1:ext=0
    elif EE[3].get()!='':ext=int(EE[3].get())
    else:ext=0
    if EE[7].get()!='':masse=float(EE[7].get())*mproton
    else:masse=mproton
    if EE[4].get()!='':
        if not(EE[4].get()[0] in ['v','m']):
            T=float(EE[4].get())
            Vini=0
        elif EE[4].get()[0]=='v':
            Vini=float(EE[4].get()[1:])
            T=Tdev(Vini,masse)
        else:
            Vini=vdeT(float(EE[4].get()[1:]),masse)[0]
            T=Tdev(Vini,masse)
    else:T,Vini=350,0
    perteQ=[1,0]
    if EE[5].get()!='':perteQ[0]=float(EE[5].get())
    
    
    if EE[8].get()!='':nbsimu=int(EE[8].get())
    else:nbsimu=1
    
    saveposvvQ=[0,0]#tracer et save
    if EE[9].get()!='':saveposvvQ[0]=int(EE[9].get())
    if EE[10].get()!='':saveposvvQ[1]=int(EE[10].get())
    if EE[11].get()!='':legende=EE[11].get()
    else:legende="temporaire"
    if EE[12].get()!='':para=EE[12].get()
    else:para="T"
    
    paravariant=[0, 10, 710, 100,2]#[Oui/Non,allant de a% de la valeur donnée à,b%,pasen %,arrêt sur valeur pour X simulations]
    if EE[13].get()!='':paravariant[0]=int(EE[13].get())
    if EE[14].get()!='':paravariant[1]=float(EE[14].get())
    if EE[15].get()!='':paravariant[2]=float(EE[15].get())
    if EE[16].get()!='':paravariant[3]=float(EE[16].get())
    if EE[17].get()!='':paravariant[4]=int(EE[17].get())
    paravariant[2]+=paravariant[3]/2
    
    
    class particule:#défini l'objet particule intérieur
        global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
        """proton issue d'un rayon cosmique"""
        def __init__(self,temp,mass=mproton):#initialisation de la vie d'une particule en fonction notamment de sa masse et de sa température/vitesse et des autres paramètres entrés. On va donnerdesn arguments position et vitesse cohérents à la particule
            global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
            
            self.m=mass
            self.T=temp
            self.parcouru=0
            self.lpm=tirageliste(llpm)
            self.v=tirageliste(listemax(temp,mass))
            self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
        
            self.vv=multscal(self.v/norme(self.vv),self.vv)
            
            if natureQ==1:
                if popv==1:
                    xy=[uniform(-1,1),uniform(-1,1)]
                    uu=uniform(0,r)
                    while uu==r:
                        uu=uniform(0,r)
                    normexy=sqrt(xy[0]**2+xy[1]**2)
                    self.pos=[xy[0]*uu/normexy,xy[1]*uu/normexy,uniform(0,L)]
                else: 
                    if mode2==0 and tirage(r/(2*L+r))==1:
                        self.pos=[uniform(-r,r),0,0]
                        self.pos[1]=uniform(-sqrt(r**2-self.pos[0]**2),sqrt(r**2-self.pos[0]**2))
                        self.pos[2]=0
                    else:
                        self.pos=[uniform(-r,r),0,0]
                        x=rand11()
                        self.pos[1]=x*sqrt(r**2-self.pos[0]**2)
                        
                        self.pos[2]=uniform(0,L)
                    
                    if scal(normale(self.pos),self.vv)<0:
                        self.vv=multscal(-1,self.vv)
                        if scal(normale(self.pos),self.vv)<=0:self.vv=multscal(self.v,normale(self.pos))
                    elif scal(normale(self.pos),self.vv)==0:
                        self.vv=multscal(self.v,normale(self.pos))
            else :
                
            
                if mode2==1:
                    if popv==1:
                        xy=[uniform(-R,R),uniform(-R,R),uniform(-R,-sqrt(R*R-r*r))]
                        while norme(xy)>=R or xy[2]==-sqrt(R*R-r*r):
                            xy=[uniform(-R,R),uniform(-R,R),uniform(-R,-sqrt(R*R-r*r))]
                        self.pos=xy
                    else:
                        self.pos=[uniform(-R,-sqrt(R*R-r*r)),0,0]
                        self.pos[1]=uniform(-sqrt(R**2-self.pos[0]**2),sqrt(R**2-self.pos[0]**2))
                        x=rand11()
                        self.pos[2]=x*sqrt(R**2-self.pos[0]**2-self.pos[1]**2)
                        self.pos[0],self.pos[2]=self.pos[2],self.pos[0]
                else:
                    if popv==1:
                        xy=[uniform(-R,R),uniform(-R,R),uniform(-R,sqrt(R*R-r*r))]
                        while norme(xy)>=R or xy[2]==sqrt(R*R-r*r):
                            xy=[uniform(-R,R),uniform(-R,R),uniform(-R,sqrt(R*R-r*r))]
                        self.pos=xy
                    else:
                        """on génère une position extérieure au trou en haut de la sphère de façon equiprobable"""
                        self.pos=[uniform(-R,sqrt(R*R-r*r)),0,0]
                        self.pos[1]=uniform(-sqrt(R**2-self.pos[0]**2),sqrt(R**2-self.pos[0]**2))
                        x=rand11()
                        self.pos[2]=x*sqrt(R**2-self.pos[0]**2-self.pos[1]**2)
                        self.pos[0],self.pos[2]=self.pos[2],self.pos[0]
                
                """on va maintenant rendre cohérent notre vecteur vitesse avec nos conditions physique"""
                if popv==0:
                    if scal(normale(self.pos),self.vv)<0:
                        self.vv=multscal(-1,self.vv)
                    elif scal(normale(self.pos),self.vv)==0:
                        self.vv=multscal(self.v,normale(self.pos))
        def randsaut(self):#méthode qui simule une portion de trajectoire rectiligne de la particule (du point d'arrivé précédent jusqu'à la prochaine collision avec les paroi/collision avec une autre particule ou alors jusqu'à sa sortie par un trou
            global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
            if natureQ==1:#cas du cylindre
                t=max((-2*(self.pos[0]*self.vv[0]+self.pos[1]*self.vv[1])-sqrt(4*(self.pos[0]*self.vv[0]+self.pos[1]*self.vv[1])**2-4*(self.vv[0]**2+self.vv[1]**2)*(self.pos[0]**2+self.pos[1]**2-r**2)))/(2*(self.vv[0]**2+self.vv[1]**2)),(-2*(self.pos[0]*self.vv[0]+self.pos[1]*self.vv[1])+sqrt(4*(self.pos[0]*self.vv[0]+self.pos[1]*self.vv[1])**2-4*(self.vv[0]**2+self.vv[1]**2)*(self.pos[0]**2+self.pos[1]**2-r**2)))/(2*(self.vv[0]**2+self.vv[1]**2)))
                
                if self.vv[2]<0:
                    t2=-self.pos[2]/self.vv[2]
                elif self.vv[2]>0:
                    t2=(L-self.pos[2])/self.vv[2]
                else : t2=t
                if t>t2:t=t2
                
                postemp=multscal(-1,somvect(self.pos,multscal(t,self.vv)))
                
                self.parcouru+=norme(somvect(self.pos,postemp))
                
                if self.lpm==0 or self.parcouru<self.lpm:#cas d'une distance parcouru inférieur à la distance nécessaire pour cogner une autre particule (cf lpm)
                
                    self.pos=somvect(self.pos,multscal(t,self.vv))
                    
                    if sigma=='random':#cas du rebond aléatoire parfait
                        self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
                        self.vv=multscal(coefperte(self.v)*self.v/norme(self.vv),self.vv)
                        self.T=self.T*coefperte(self.v)*coefperte(self.v)
                        if scal(normale(self.pos),self.vv)<0:
                            self.vv=multscal(-1,self.vv)
                            if scal(normale(self.pos),self.vv)<=0:self.vv=multscal(self.v,normale(self.pos))
                        elif scal(normale(self.pos),self.vv)==0:
                            self.vv=multscal(self.v,normale(self.pos))
                    else:#cas du rebond miroir ou aléatoire partiel
                        rn=randnormale(self.pos,self.vv)#génération d'une première normale
                    
                        vv=multscal(coefperte(self.v),somvect(self.vv,multscal(-2*scal(self.vv,rn),rn)))
                        self.T=self.T*coefperte(self.v)*coefperte(self.v)
                        while sigma!=0 and scal(normale(self.pos),vv)<=0 and ((mode2==1 and 0<self.pos[2]<L) or (mode2==0 and self.pos[2]<L)):
                           
                            rn=randnormale(self.pos,self.vv)#génération d'une autre normale si la précédente n'est pas cohérente
                        
                            vv=multscal(coefperte(norme(self.vv)),somvect(self.vv,multscal(-2*scal(self.vv,rn),rn)))
                        self.vv=vv
                    self.v=norme(self.vv)
                
                        
                    return t#renvoie le temps de parcours de cette portion droite de trajectoire
                else :#cas où on cagne une particule
                    
                    savepos=self.pos
                    self.pos=somvect(multscal(-1,postemp),multscal(-(self.parcouru-self.lpm)/self.v,self.vv))
                    t=norme(somvect(self.pos,multscal(-1,savepos)))/self.v
                    self.parcouru=0
                    self.lpm=tirageliste(llpm)
                    self.v=tirageliste(listemax(self.T,self.m))
            
                    self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
                    self.vv=multscal(self.v/norme(self.vv),self.vv)
                    self.v=norme(self.vv)
                    return t#renvoie le temps de parcours de cette portion droite de trajectoire avant la collision
            else :#idem pour la sphère
                
                """met à jour la position et la vitesse et la position de la particule après une trajectoir droite jusqu'à rencontrer un nouveau bord"""
                t=max((-2*scal(self.pos,self.vv)-sqrt(4*scal(self.pos,self.vv)**2-4*(self.v**2)*(norme(self.pos)**2-R**2)))/2/(self.v**2),(-2*scal(self.pos,self.vv)+sqrt(4*scal(self.pos,self.vv)**2-4*(self.v**2)*(norme(self.pos)**2-R**2)))/2/(self.v**2))
                t2=t
                if self.vv[2]>0:
                    if mode2==1:t2=(-sqrt(R*R-r*r)-self.pos[2])/self.vv[2]
                    else:t2=(sqrt(R*R-r*r)-self.pos[2])/self.vv[2]
                if  t>t2:
                    t=t2
                postemp=multscal(-1,somvect(self.pos,multscal(t,self.vv)))
                self.parcouru+=norme(somvect(self.pos,postemp))
                if self.lpm==0 or self.parcouru<self.lpm:
                    self.pos=somvect(self.pos,multscal(t,self.vv))
                    
                    if sigma=='random':
                        self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
                        self.vv=multscal(coefperte(self.v)*self.v/norme(self.vv),self.vv)
                        self.T=self.T*coefperte(self.v)*coefperte(self.v)
                        if scal(normale(self.pos),self.vv)<0:
                            self.vv=multscal(-1,self.vv)
                            if scal(normale(self.pos),self.vv)<=0:self.vv=multscal(self.v,normale(self.pos))
                        elif scal(normale(self.pos),self.vv)==0:
                            self.vv=multscal(self.v,normale(self.pos))
                    else:
                        rn=randnormale(self.pos,self.vv)
                    
                        vv=multscal(coefperte(self.v),somvect(self.vv,multscal(-2*scal(self.vv,rn),rn)))
                        self.T=self.T*coefperte(self.v)*coefperte(self.v)
                        while sigma!=0 and scal(normale(self.pos),vv)<=0 and ((mode2==1 and self.pos[2]<-sqrt(R*R-r*r)) or (mode2==0 and self.pos[2]<sqrt(R*R-r*r))):
                          
                            rn=randnormale(self.pos,self.vv)
                        
                            vv=multscal(coefperte(norme(self.vv)),somvect(self.vv,multscal(-2*scal(self.vv,rn),rn)))
                        self.vv=vv
                    self.v=norme(self.vv)
                
                        
                    return t
                else :
                    savepos=self.pos
                    self.pos=somvect(multscal(-1,postemp),multscal(-(self.parcouru-self.lpm)/self.v,self.vv))
                    t=norme(somvect(self.pos,multscal(-1,savepos)))/self.v
                    self.parcouru=0
                    self.lpm=tirageliste(llpm)
                    self.v=tirageliste(listemax(self.T,self.m))
            
                    self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
                    self.vv=multscal(self.v/norme(self.vv),self.vv)
                    return t
                    
            
        def monosimulation(self):#simule la vie entière de l'objet particule en se servant de la méthode de classe randsaut précédente et sors les résultats et datas désirés durée de la vie/nombre de contacts effectués/booléen de sortie(/coordonnées succéssives en vu d'un futur tracer/valeur du paramètre variant (simu type 3 cf guide )):
            global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
            nbcontacts=0 
            sorti=0
            tempstotal=0
            collepas=1
            
            if saveposvvQ[0]==1 or saveposvvQ[1]==1:#on veut tracer les trajectoires et/ou sauvegarder les positions
                x,y,z=[self.pos[0]],[self.pos[1]],[self.pos[2]]
                aecrire=ezw(self.pos)+" "+ezw(self.vv)
            if natureQ==1:#cas du cylindre
                if mode2==1:#cas du double trou
                    while (collepas==1 or self.parcouru==0) and  (tempstotal==0 or 0<self.pos[2]<L):
                        collepas=not(tirage(P(self.v,self.m)))
                        nbcontacts+=1
                        tempstotal+=self.randsaut()   
                        if saveposvvQ[0]==1 or saveposvvQ[1]==1:
                            x+=[self.pos[0]]
                            y+=[self.pos[1]]
                            z+=[self.pos[2]]
                            aecrire+=" "+ezw(self.pos)+" "+ezw(self.vv)
                        
                        
                    sorti=int(not(0<self.pos[2]<L))
                    if sorti==1:nbcontacts-=1
                else:#cas du simple trou
                    while (collepas==1 or self.parcouru==0) and  (tempstotal==0 or self.pos[2]<L):
                        
                        collepas=not(tirage(P(self.v,self.m)))
                        nbcontacts+=1
                        tempstotal+=self.randsaut()   
                        if saveposvvQ[0]==1 or saveposvvQ[1]==1:
                            x+=[self.pos[0]]
                            y+=[self.pos[1]]
                            z+=[self.pos[2]]
                            aecrire+=" "+ezw(self.pos)+" "+ezw(self.vv)
                        
                        
                    sorti=int(self.pos[2]>=L)
                    if sorti==1:nbcontacts-=1
            else :#cas de la sphère
                if mode2==1:#cas du trou de postion z<0
                    while (collepas==1 or self.parcouru==0) and (tempstotal==0 or self.pos[2]<-sqrt(R*R-r*r)):
                        collepas=not(tirage(P(self.v,self.m)))
                        nbcontacts+=1
                        tempstotal+=self.randsaut()
                        if saveposvvQ[0]==1 or saveposvvQ[1]==1:
                            x+=[self.pos[0]]
                            y+=[self.pos[1]]
                            z+=[self.pos[2]]
                            aecrire+=" "+ezw(self.pos)+" "+ezw(self.vv)
                        
                    sorti=int(not(self.pos[2]<-sqrt(R*R-r*r)))
                    if sorti==1:nbcontacts-=1
                else:#cas du trou de position z>0
                    while (collepas==1 or self.parcouru==0) and (tempstotal==0 or self.pos[2]<sqrt(R*R-r*r)):
                        collepas=not(tirage(P(self.v,self.m)))
                        nbcontacts+=1
                        tempstotal+=self.randsaut()
                        if saveposvvQ[0]==1 or saveposvvQ[1]==1:
                            x+=[self.pos[0]]
                            y+=[self.pos[1]]
                            z+=[self.pos[2]]
                            aecrire+=" "+ezw(self.pos)+" "+ezw(self.vv)
                        
                    sorti=int(not(self.pos[2]<sqrt(R*R-r*r)))
                    if sorti==1:nbcontacts-=1
            dicopara={"v":Vini,"l":lpm,"R":R,"r":r,"Ts":Ts,"T":self.T,"m":self.m,"s":sigma}
            #on retourne les datas de la simulation selon son type et ce que désire l'utilisateur:
            if (saveposvvQ[0]==1  or saveposvvQ[1]==1) and paravariant[0]==0:
                return ([[x,y,z],tempstotal,nbcontacts,sorti])
            elif (saveposvvQ[0]==1  or saveposvvQ[1]==1) and paravariant[0]==1:
                return ([[x,y,z],dicopara[para],tempstotal,nbcontacts,sorti])
            elif paravariant[0]==0:return([tempstotal,nbcontacts,sorti])
            else :return ([dicopara[para],tempstotal,nbcontacts,sorti])
    class particuleext:#idem avec une particule d'origine extérieure !
        global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
        """proton qui rentre par le trou"""
        def __init__(self,temp,mass=mproton):
            global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
            self.m=mass
            self.T=temp
            self.parcouru=0
            self.lpm=tirageliste(llpm)
            self.v=tirageliste(listemax(temp,mass))
            self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
            self.vv=multscal(self.v/norme(self.vv),self.vv)
            if natureQ==1:
                if mode2==1:
                    if tirage(0.5)==1:
                        self.pos=[r,r,L]
                    else:
                        self.pos=[r,r,0]
                        
                else:self.pos=[r,r,L]
            else :
                if mode2==1:self.pos=[r,r,-sqrt(R*R-r*r)]
                else:self.pos=[r,r,sqrt(R*R-r*r)]
            while sqrt(self.pos[0]**2+self.pos[1]**2)>=r:
                self.pos[0]=uniform(-r,r)
                self.pos[1]=uniform(-sqrt(r**2-self.pos[0]**2),sqrt(r**2-self.pos[0]**2))
            if natureQ==1 and self.pos[2]==0:
                if scal([0,0,1],self.vv)<0:
                    self.vv=multscal(-1,self.vv)
                elif scal([0,0,1],self.vv)==0:
                    self.vv=multscal(self.v,[0,0,-1])
            else:
                if scal([0,0,-1],self.vv)<0:
                    self.vv=multscal(-1,self.vv)
                elif scal([0,0,-1],self.vv)==0:
                    self.vv=multscal(self.v,[0,0,-1])
            
        def randsaut(self):
            
            global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
            if natureQ==1:
                t=max((-2*(self.pos[0]*self.vv[0]+self.pos[1]*self.vv[1])-sqrt(4*(self.pos[0]*self.vv[0]+self.pos[1]*self.vv[1])**2-4*(self.vv[0]**2+self.vv[1]**2)*(self.pos[0]**2+self.pos[1]**2-r**2)))/(2*(self.vv[0]**2+self.vv[1]**2)),(-2*(self.pos[0]*self.vv[0]+self.pos[1]*self.vv[1])+sqrt(4*(self.pos[0]*self.vv[0]+self.pos[1]*self.vv[1])**2-4*(self.vv[0]**2+self.vv[1]**2)*(self.pos[0]**2+self.pos[1]**2-r**2)))/(2*(self.vv[0]**2+self.vv[1]**2)))
                if self.vv[2]<0:
                    t2=-self.pos[2]/self.vv[2]
                elif self.vv[2]>0:
                    t2=(L-self.pos[2])/self.vv[2]
                else : t2=t
                if t>t2:t=t2
                postemp=multscal(-1,somvect(self.pos,multscal(t,self.vv)))
                
                self.parcouru+=norme(somvect(self.pos,postemp))
                
                if self.lpm==0 or self.parcouru<self.lpm:
                    self.pos=somvect(self.pos,multscal(t,self.vv))
                    
                    if sigma=='random':
                        self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
                        self.vv=multscal(coefperte(self.v)*self.v/norme(self.vv),self.vv)
                        self.T=self.T*coefperte(self.v)*coefperte(self.v)
                        if scal(normale(self.pos),self.vv)<0:
                            self.vv=multscal(-1,self.vv)
                            if scal(normale(self.pos),self.vv)<=0:self.vv=multscal(self.v,normale(self.pos))
                        elif scal(normale(self.pos),self.vv)==0:
                            self.vv=multscal(self.v,normale(self.pos))
                    else:
                        rn=randnormale(self.pos,self.vv)
                    
                        vv=multscal(coefperte(self.v),somvect(self.vv,multscal(-2*scal(self.vv,rn),rn)))
                        self.T=self.T*coefperte(self.v)*coefperte(self.v)
                        while sigma!=0 and scal(normale(self.pos),vv)<=0 and ((mode2==1 and 0<self.pos[2]<L) or (mode2==0 and self.pos[2]<L)):
                            
                            rn=randnormale(self.pos,self.vv)
                        
                            vv=multscal(coefperte(norme(self.vv)),somvect(self.vv,multscal(-2*scal(self.vv,rn),rn)))
                        self.vv=vv
                    self.v=norme(self.vv)
                
                        
                    return t
                else :
                    savepos=self.pos
                    self.pos=somvect(multscal(-1,postemp),multscal(-(self.parcouru-self.lpm)/self.v,self.vv))
                    t=norme(somvect(self.pos,multscal(-1,savepos)))/self.v
                    self.parcouru=0
                    self.lpm=tirageliste(llpm)
                    self.v=tirageliste(listemax(self.T,self.m))
            
                    self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
                    self.vv=multscal(self.v/norme(self.vv),self.vv)
                    return t
            else :
                
                """met à jour la position et la vitesse et la position de la particule après une trajectoir droite jusqu'à rencontrer un nouveau bord"""
                t=max((-2*scal(self.pos,self.vv)-sqrt(4*scal(self.pos,self.vv)**2-4*(self.v**2)*(norme(self.pos)**2-R**2)))/2/(self.v**2),(-2*scal(self.pos,self.vv)+sqrt(4*scal(self.pos,self.vv)**2-4*(self.v**2)*(norme(self.pos)**2-R**2)))/2/(self.v**2))
                t2=t
                if self.vv[2]>0:
                    if mode2==1:t2=(-sqrt(R*R-r*r)-self.pos[2])/self.vv[2]
                    else:t2=(sqrt(R*R-r*r)-self.pos[2])/self.vv[2]
                if  t>t2:
                    t=t2
                postemp=multscal(-1,somvect(self.pos,multscal(t,self.vv)))
                self.parcouru+=norme(somvect(self.pos,postemp))
                if self.lpm==0 or self.parcouru<self.lpm:
                    self.pos=somvect(self.pos,multscal(t,self.vv))
                    
                    if sigma=='random':
                        self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
                        self.vv=multscal(coefperte(self.v)*self.v/norme(self.vv),self.vv)
                        self.T=self.T*coefperte(self.v)*coefperte(self.v)
                        if scal(normale(self.pos),self.vv)<0:
                            self.vv=multscal(-1,self.vv)
                            if scal(normale(self.pos),self.vv)<=0:self.vv=multscal(self.v,normale(self.pos))
                        elif scal(normale(self.pos),self.vv)==0:
                            self.vv=multscal(self.v,normale(self.pos))
                    else:
                        rn=randnormale(self.pos,self.vv)
                    
                        vv=multscal(coefperte(self.v),somvect(self.vv,multscal(-2*scal(self.vv,rn),rn)))
                        self.T=self.T*coefperte(self.v)*coefperte(self.v)
                        while sigma!=0 and scal(normale(self.pos),vv)<=0 and ((mode2==1 and self.pos[2]<-sqrt(R*R-r*r)) or (mode2==0 and self.pos[2]<sqrt(R*R-r*r))):
                            
                            rn=randnormale(self.pos,self.vv)
                        
                            vv=multscal(coefperte(norme(self.vv)),somvect(self.vv,multscal(-2*scal(self.vv,rn),rn)))
                        self.vv=vv
                    self.v=norme(self.vv)
                
                        
                    return t
                else :
                    savepos=self.pos
                    self.pos=somvect(multscal(-1,postemp),multscal(-(self.parcouru-self.lpm)/self.v,self.vv))
                    t=norme(somvect(self.pos,multscal(-1,savepos)))/self.v
                    self.parcouru=0
                    self.lpm=tirageliste(llpm)
                    self.v=tirageliste(listemax(self.T,self.m))
            
                    self.vv=[rand11()*uniform(1,1000),rand11()*uniform(1,1000),rand11()*uniform(1,1000)]
                    self.vv=multscal(self.v/norme(self.vv),self.vv)
                    return t          
            
        def monosimulation(self):
            global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
            nbcontacts=0 
            sorti=0
            tempstotal=0
            collepas=1
            if saveposvvQ[0]==1 or saveposvvQ[1]==1:
                x,y,z=[self.pos[0]],[self.pos[1]],[self.pos[2]]
                aecrire=ezw(self.pos)+" "+ezw(self.vv)
            if natureQ==1:
                if mode2==1:
                    while (collepas==1 or self.parcouru==0) and  (tempstotal==0 or 0<self.pos[2]<L):
                        collepas=not(tirage(P(self.v,self.m)))
                        nbcontacts+=1
                        tempstotal+=self.randsaut()   
                        if saveposvvQ[0]==1 or saveposvvQ[1]==1:
                            x+=[self.pos[0]]
                            y+=[self.pos[1]]
                            z+=[self.pos[2]]
                            aecrire+=" "+ezw(self.pos)+" "+ezw(self.vv)
                        
                        
                    sorti=int(not(0<self.pos[2]<L))
                    if sorti==1:nbcontacts-=1
                else:
                    while (collepas==1 or self.parcouru==0) and  (tempstotal==0 or self.pos[2]<L):
                        collepas=not(tirage(P(self.v,self.m)))
                        nbcontacts+=1
                        tempstotal+=self.randsaut()   
                        if saveposvvQ[0]==1 or saveposvvQ[1]==1:
                            x+=[self.pos[0]]
                            y+=[self.pos[1]]
                            z+=[self.pos[2]]
                            aecrire+=" "+ezw(self.pos)+" "+ezw(self.vv)
                        
                        
                    sorti=int(self.pos[2]>=L)
                    if sorti==1:nbcontacts-=1
            else :
                if mode2==1:
                    while (collepas==1 or self.parcouru==0) and (tempstotal==0 or self.pos[2]<-sqrt(R*R-r*r)):
                        collepas=not(tirage(P(self.v,self.m)))
                        nbcontacts+=1
                        tempstotal+=self.randsaut()
                        if saveposvvQ[0]==1 or saveposvvQ[1]==1:
                            x+=[self.pos[0]]
                            y+=[self.pos[1]]
                            z+=[self.pos[2]]
                            aecrire+=" "+ezw(self.pos)+" "+ezw(self.vv)
                        
                        
                    sorti=int(not(self.pos[2]<-sqrt(R*R-r*r)))
                    if sorti==1:nbcontacts-=1
                else:
                    
                    while (collepas==1 or self.parcouru==0) and (tempstotal==0 or self.pos[2]<sqrt(R*R-r*r)):
                        collepas=not(tirage(P(self.v,self.m)))
                        nbcontacts+=1
                        tempstotal+=self.randsaut()
                        if saveposvvQ[0]==1 or saveposvvQ[1]==1:
                            x+=[self.pos[0]]
                            y+=[self.pos[1]]
                            z+=[self.pos[2]]
                            aecrire+=" "+ezw(self.pos)+" "+ezw(self.vv)
                    sorti=int(not(self.pos[2]<sqrt(R*R-r*r)))
                    if sorti==1:nbcontacts-=1
            dicopara={"v":Vini,"l":lpm,"R":R,"r":r,"Ts":Ts,"T":self.T,"m":self.m,"s":sigma}
            if (saveposvvQ[0]==1  or saveposvvQ[1]==1) and paravariant[0]==0:
                return ([[x,y,z],tempstotal,nbcontacts,sorti])
            elif (saveposvvQ[0]==1  or saveposvvQ[1]==1) and paravariant[0]==1:
                return ([[x,y,z],dicopara[para],tempstotal,nbcontacts,sorti])
            elif paravariant[0]==0:return([tempstotal,nbcontacts,sorti])
            else :return ([dicopara[para],tempstotal,nbcontacts,sorti])
    
    
    def main():#fonction qui va réaliser la simulation en prenant en compte tout ce que désire l'utilisateur et écrire les datas dans un fichier+graphe facultatif 3D(cf guide pour comprendre sa logique)
        global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
        tittre=ezw(["legende",legende,"natureQ",natureQ,"mode2",mode2,"sigma",sigma,"r",r,"RouL",R,"ext",ext,"Ts",Ts,"T",T,"Vini",Vini,"popvolumique",popv,"coefperte",perteQ,"masse",masse,"nbsimu",nbsimu,"para",para,"saveposvvQ",saveposvvQ,"lpm",lpm,"sticklist",sticklist,"paravariant",paravariant])
        
        f=open(os.path.join(cheminpy,"Res/"+legende+".txt"),'w')
        f.write("")
        f.close()
        f=open(os.path.join(cheminpy, "Res/"+legende+".txt"),'a')
        f.write(tittre+"\nTempsVie(s) ; Nb contacts ; SortieParTrou?")
        gaussliste=listegauss(sigma)
        llpm=listelpm(lpm)
        if saveposvvQ==[0,0] and paravariant[0]==0:#simulation de type 1)
            
            for i in range(0,nbsimu):
                if ext==1:p=particuleext(T,masse)
                else : p=particule(T,masse)
                f.write("\n"+ezw(p.monosimulation()))
            
        if (saveposvvQ[0]==1  or saveposvvQ[1]==1 )and paravariant[0]==0:# simulation de type 2)
            
            if saveposvvQ[0]==1:
                fig = plt.figure()
                ax = Axes3D(fig)
            for i in range(0,nbsimu):
                if ext==1:p=particuleext(T,masse)
                else : p=particule(T,masse)
                pm=p.monosimulation()
                
                if  saveposvvQ[1]==1 :
                
                    f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                else:
                    
                    f.write("\n"+ezw(pm[1:]))
                if saveposvvQ[0]==1 :
                    ax.plot(pm[0][0],pm[0][1],pm[0][2])
                    
            if saveposvvQ[0]==1 and natureQ==0 :        
                u = np.linspace(0, 2*np.pi, 50)
                if mode2==1:v = np.linspace(np.pi-asin(r/R), np.pi, 50)
                else:v = np.linspace(0+asin(r/R), np.pi, 50)
                ctheta=np.linspace(0,2*np.pi,50)
                xc=r*(np.cos(ctheta))
                yc=r*(np.sin(ctheta))
                if mode2==1:
                    ax.plot([0],[0],-R+2*r)
                    zc=-sqrt(R*R-r*r)
                else:
                    ax.plot([0],[0],R)
                    zc=sqrt(R*R-r*r)
                ax.plot(xc,yc,zc,color='green')
                x = R * np.outer(np.cos(u), np.sin(v))
                y = R * np.outer(np.sin(u), np.sin(v))
                z = R * np.cos(v) #coordonnée sphère en cartésien : z=RcosO,x=rsinOsinP,y=RsinOcosP
                ax.plot_wireframe(x, y, z, rstride=7, cstride=7,color='grey')
                
                
            elif saveposvvQ[0]==1:
                u = np.linspace(0, 2*np.pi, 50)
                v = np.linspace(0,L, 50)
                ctheta=np.linspace(0,2*np.pi,50)
                xc=r*(np.cos(ctheta))
                yc=r*(np.sin(ctheta))
                zc=L
                ax.plot(xc,yc,zc,color='green')
                x = r *np.cos(u)
                y = r *np.sin(u)
                if 2*r>L:
                    ax.plot([0],[0],2*r)
                else:
                    ax.plot([L/2],[0],0)
                    ax.plot([0],[L/2],0)
                    ax.plot([-L/2],[0],0)
                    ax.plot([0],[-L/2],0)
                deb=0
                if mode2==1:
                    deb=1
                    ax.plot(xc,yc,0,color='green')
                
                for a in range(deb,10):
                    pas=L/10
                    ax.plot(x,y,a*pas,color='grey')
                
            plt.show(block=True)
        if saveposvvQ==[0,0] and paravariant[0]==1:#simulation de type 3) simple
            
            if para=='R':
                R=paravariant[1]*10**(-6)
                L=paravariant[1]*10**(-6)
                while R<=paravariant[2]*10**(-6):
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        f.write("\n"+ezw(p.monosimulation()))
                    R+=paravariant[3]*10**(-6)
                    L+=paravariant[3]*10**(-6)
                    
            if para=='r':
                r=paravariant[1]*10**(-6)
                while r<=paravariant[2]*10**(-6):
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        f.write("\n"+ezw(p.monosimulation()))
                    r+=paravariant[3]*10**(-6)
                    
            if para=='T':
                T=paravariant[1]
                while T<=paravariant[2]:
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        f.write("\n"+ezw(p.monosimulation()))
                    T+=paravariant[3]
            if para=='Ts':
                Ts=paravariant[1]
                while Ts<=paravariant[2]:
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        f.write("\n"+ezw(p.monosimulation()))
                    Ts+=paravariant[3]        
            if para=='s':
                sigma=paravariant[1]
                while sigma<=paravariant[2]:
                    gaussliste=listegauss(sigma)
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        f.write("\n"+ezw(p.monosimulation()))
                    sigma+=paravariant[3]
            if para=='l':
                lpm=paravariant[1]
                while lpm<=paravariant[2]:
                    llpm=list()
                    
                    llpm=listelpm(lpm)
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        f.write("\n"+ezw(p.monosimulation()))
                    lpm+=paravariant[3]
                    
            if para=='m':
                masse=paravariant[1]*mproton
                while masse<=paravariant[2]*mproton:
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        f.write("\n"+ezw(p.monosimulation()))
                    masse+=paravariant[3]*mproton
            if para=='v':
                Vini=paravariant[1]
                while Vini<=paravariant[2]:
                
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    Vini+=paravariant[3]                
                    
        if saveposvvQ==[0,1]  and paravariant[0]==1:#simulation de type 3) et sauvegarde des postions
            
            if para=='R':
                R=paravariant[1]*10**(-6)
                L=paravariant[1]*10**(-6)
                while R<=paravariant[2]*10**(-6):
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    R+=paravariant[3]*10**(-6)
                    L+=paravariant[3]*10**(-6)
                    
            if para=='r':
                r=paravariant[1]*10**(-6)
                while r<=paravariant[2]*10**(-6):
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    r+=paravariant[3]*10**(-6)
                    
            if para=='T':
                T=paravariant[1]
                while T<=paravariant[2]:
                
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    T+=paravariant[3]
            if para=='Ts':
                Ts=paravariant[1]
                while Ts<=paravariant[2]:
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    Ts+=paravariant[3]        
            if para=='s':
                sigma=paravariant[1]
                while sigma<=paravariant[2]:
                    gaussliste=listegauss(sigma)
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    sigma+=paravariant[3]
            if para=='l':
                lpm=paravariant[1]
                while lpm<=paravariant[2]:
                    llpm=list()
                    
                    llpm=listelpm(lpm)
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    lpm+=paravariant[3]
                    
            if para=='m':
                masse=paravariant[1]*mproton
                while masse<=paravariant[2]*mproton:
                
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    masse+=paravariant[3]*mproton   
            if para=='v':
                Vini=paravariant[1]
                while Vini<=paravariant[2]:
                
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    Vini+=paravariant[3]        
                            
        if saveposvvQ[0]==1 and paravariant[0]==1:#simulation de type 3) et tracer 3D
            if para=='R':
                R=paravariant[1]*10**(-6)
                L=paravariant[1]*10**(-6)
                while R<=paravariant[2]*10**(-6):
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    R+=paravariant[3]*10**(-6)
                    L+=paravariant[3]*10**(-6)
                    
            if para=='r':
                r=paravariant[1]*10**(-6)
                while r<=paravariant[2]*10**(-6):
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    r+=paravariant[3]*10**(-6)
                    
            if para=='T':
                T=paravariant[1]
                while T<=paravariant[2]:
                
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    T+=paravariant[3]
            if para=='Ts':
                Ts=paravariant[1]
                while Ts<=paravariant[2]:
                    
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    Ts+=paravariant[3]        
            if para=='s':
                sigma=paravariant[1]
                while sigma<=paravariant[2]:
                    gaussliste=listegauss(sigma)
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    sigma+=paravariant[3]
            if para=='l':
                lpm=paravariant[1]
                while lpm<=paravariant[2]:
                    llpm=list()
                    
                    llpm=listelpm(lpm)
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    lpm+=paravariant[3]
                    
            if para=='m':
                masse=paravariant[1]*mproton
                while masse<=paravariant[2]*mproton:
                
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    masse+=paravariant[3]*mproton 
            if para=='v':
                Vini=paravariant[1]
                while Vini<=paravariant[2]:
                
                    for i in range(0,paravariant[4]):
                        if ext==1:p=particuleext(T,masse)
                        else : p=particule(T,masse)
                        pm=p.monosimulation()
                        f.write("\n"+ezw(pm[0][0])+" "+ezw(pm[0][1])+" "+ezw(pm[0][2])+" "+ezw(pm[1:]))
                    Vini+=paravariant[3]
            
        f.close()#fermeture du fichier de datas
        if saveposvvQ[0]==1 and paravariant[0]==1:grapheurpara(legende)#si l'on a sélectionné une simu type 3) avec tracer
    main()

def start():#récupère les données de la première fenêtre d'acceuil la détruit
    global color, plmax, Vini, cheminpy, legabs, titre, llpm, lpm, popv, Ts,sticklist, mode2, nomm,iii,legg, kb,mproton, gaussliste,sigma,natureQ,r,R,L,ext,T,perteQ,masse,nbsimu,para,saveposvvQ,legende,paravariant
    if FE.get()!='':
        cheminpy=FE.get()
    else:
        cheminpy = './data'
    try:
        import tkinter.filedialog
    except:
        pass
    F.destroy()
    
#Coeur de l'interface graphique
#fenêtre d'acceuil:
F=Tk()
F.configure(bg="#99aaaa")
F.title("IASB - particle sticking on icy grains simulation")
FL=Label(F,bg="#99aaaa",text="Entrez le chemin\nabsolu vers le répertoire data\n(exemple:'C:/data/' or '/home/data/')")
FL.pack()
FE=Entry(F)
FE.pack()
FB=Button(F,text='lancer',command=start,bg="#999999")
FB.pack()
F.mainloop()
#Interface principale:
f=Tk()
Ck=PhotoImage(file=os.path.join(cheminpy,"checked.gif"))
Gp=PhotoImage(file=os.path.join(cheminpy,"graphr.gif"))
Fd=PhotoImage(file=os.path.join(cheminpy,"fond.gif"))
f.configure(bg="#99aaaa",width=1200,height=700)
c = Label(f,width=1200,height=700,image=Fd)

c.place(y=-2,x=-2)
f.title("Stage IASB")
L0=Label(f,justify=LEFT,text="[random] écart type de\nl'aléatoire du rebond (0 pour miroir)")
L0.place(y=200,x=400)
E0=Entry(f)
E0.place(y=200,x=600)
L1=Label(f,justify=LEFT,text="[1]Type de Géométrie\n(0 pour la Sphère;1 pour le Cylindre)")
L1.place(y=10,x=0)
E1=Entry(f)
E1.place(y=10,x=200)
L2=Label(f,justify=LEFT,text="[1]Rayon du trou\nvers l'extérieur\n(en micromètres)")
L2.place(y=50,x=0)
E2=Entry(f)
E2.place(y=55,x=200)
L3=Label(f,justify=LEFT,text="[2]Rayon de la Sphère ou\nHauteur du Cylindre (en micromètres)")
L3.place(y=60,x=400)
E3=Entry(f)
E3.place(y=60,x=600)
L4=Label(f,justify=LEFT,text="[0]Particules\narrivent-elles par le trou?")
L4.place(y=150,x=0)
E4=Entry(f)
E4.place(y=150,x=200)
L5=Label(f,justify=LEFT,text="[350K]Température des particules(K)\n ou vitesse (m/s) ('v'+Vitesse ou\n'm'+Température si l'on veut la vitesse moyenne associée)")
L5.place(y=100,x=0)
E5=Entry(f)
E5.place(y=100,x=200)
L6=Label(f,justify=LEFT,text="[1]Coefficient de perte\nde vitesse après un\nchoc (1 si vitesse conservée)")
L6.place(y=250,x=0)
E6=Entry(f)
E6.place(y=250,x=200)
L7=Label(f,text="[Nom du fichier]Titre Graphe")
L7.place(y=330,x=800)
E7=Entry(f)
E7.place(y=330,x=1000)
L8=Label(f,justify=LEFT,text="[1]Masse des particules\n(en masse de proton)")
L8.place(y=200,x=0)
E8=Entry(f)
E8.place(y=200,x=200)
L9=Label(f,justify=LEFT,text="[1]Nombre de\nparticules à simuler")
L9.place(y=300,x=400)
E9=Entry(f)
E9.place(y=300,x=600)
L10=Label(f,justify=LEFT,text="[0]Tracer la(es) trajectoire(s)?")
L10.place(y=350,x=0)
E10=Entry(f)
E10.place(y=350,x=200)
L11=Label(f,justify=LEFT,text="[0]Faut-il sauvegarder les\ntrajectoires?")
L11.place(y=350,x=400)
E11=Entry(f)
E11.place(y=350,x=600)
L12=Label(f,justify=LEFT,text="[temporaire]Nom à donner à\nla sauvegarde(sans l'extension)")
L12.place(y=400,x=400)
E12=Entry(f)
E12.place(y=400,x=600)
L13=Label(f,justify=LEFT,text="[T]Paramètre à faire varier\nparmi r,R,s,T,m,Ts,l")
L13.place(y=450,x=400)
E13=Entry(f)
E13.place(y=450,x=600)
L14=Label(f,justify=LEFT,text="[0]Faire Varier un paramètre?")
L14.place(y=450,x=0)
E14=Entry(f)
E14.place(y=450,x=200)
L15=Label(f,justify=LEFT,text="[10]de")
L15.place(y=500,x=100)
E15=Entry(f)
E15.place(y=500,x=200)
L16=Label(f,justify=LEFT,text="[710]à")
L16.place(y=500,x=500)
E16=Entry(f)
E16.place(y=500,x=600)
L17=Label(f,justify=LEFT,text="[100]par pas de")
L17.place(y=550,x=70)
E17=Entry(f)
E17.place(y=550,x=200)
L18=Label(f,justify=LEFT,text="[2]Nombre de simulations à\neffectuer à chaque pallier")
L18.place(y=550,x=400)
E18=Entry(f)
E18.place(y=550,x=600)
EE=[E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15,E16,E17,E18,E0]
LL=[L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15,L16,L17,L18,L0]
Ll1=Label(f,text="[temporaire]nom du fichier\nà représenter (sans l'extension)",bg="#99aaaa")
Ee1=Entry(f)
Ll1.place(y=370,x=800)
Ee1.place(y=380,x=1000)
Ll0=Label(f,text="[Paramètre] Légende abscisses",bg="#99aaaa")
Ee0=Entry(f)
Ll0.place(y=530,x=800)
Ee0.place(y=530,x=1000)
Ll2=Label(f,justify=LEFT,text="[h3]Que tracer (entrer 0,1,2ou3) ?\nrien,1,2,3:moyennes des simulations\nGraphe 2D:\n(h)1>temps,(h)2>contacts,(h)3>%sorties\nGraphe 3D : 0>Tracer des trajectoires",bg="#99aaaa")
Ee2=Entry(f)
Ll2.place(y=420,x=800)
Ee2.place(y=450,x=1000)
Ll3=Label(f,text="[Aucune]Légende de la courbe",bg="#99aaaa")
Ee3=Entry(f)
Ll3.place(y=560,x=800)
Ee3.place(y=560,x=1000)
Ll4=Label(f,justify=LEFT,text="[0]Option de Géométrie activée?\n(Double Trou pour le Cylindre\nHauteur trou<0 pour la Sphère)",bg="#99aaaa")
Ee4=Entry(f)
Ll4.place(y=2,x=400)
Ee4.place(y=10,x=600)
Ll5=Label(f,justify=LEFT,text="[10]Température de la\nsurface de la cavité(K)",bg="#99aaaa")
Ee5=Entry(f)
Ll5.place(y=100,x=400)
Ee5.place(y=100,x=600)
Ll6=Label(f,justify=LEFT,text="[Proba de Bergeron pour p(v,M)]\nProbabilité de se coller",bg="#99aaaa")
Ee6=Entry(f)
Ll6.place(y=300,x=0)
Ee6.place(y=300,x=200)
Ll7=Label(f,justify=LEFT,text="[0]Type de zone d'apparition\nuniforme :\n0 =Surfacique, 1=Volumique?",bg="#99aaaa")
Ee7=Entry(f)
Ll7.place(y=150,x=400)
Ee7.place(y=150,x=600)
Ll8=Label(f,justify=LEFT,text="[0]libre parcours\nmoyen, 0 si infini(en micromètre)",bg="#99aaaa")
Ee8=Entry(f)
Ll8.place(y=250,x=400)
Ee8.place(y=250,x=600)
Ll9=Label(f,justify=LEFT,text="[1]changer de couleur ?",bg="#a66464")
Ee9=Entry(f)
Ll9.place(y=500,x=800)
Ee9.place(y=500,x=1000)

Lres=Label(f,bg="#99aaaa")
Lres.place(y=600,x=800)


for a in LL:
    a["bg"]="#99aaaa"
for a in [Ll1,Ll2,Ll3,L7,Ll0]:
    a["bg"]="#a66464"
Lres["bg"]="#c1888a"
B=Button(f,image=Ck,command=ain,bg="#999999")
B.place(y=600,x=200)
B2=Button(f,image=Gp,command=grapheur,bg="#a66464")
B2.place(y=600,x=1100)
Bedit=Button(f,text='éditer\nabscisse\n et titre',command=edit,bg="#a66464")
Bedit.place(y=600,x=1045)

f.mainloop()#loop de l'interface principale
