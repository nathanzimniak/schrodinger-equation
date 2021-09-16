#!/usr/bin/env python3
# coding: utf-8

'''
Programme pour résoudre l'équation de Schrödinger dépendante du temps par la méthodes des différences finies
'''

__author__ = 'Nathan Zimniak'


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import time
start_time = time.time()




#Initialisation des constantes
Lx = 201                 #Longueur spatiale
Ly = 201                 #Largeur spatiale
Nbi = 40001            #Nombre d'itérations temporelles
dx = 1/(Lx-1)           #Pas spatial suivant x
dy = 1/(Ly-1)           #Pas spatial suivant y
dt = 1e-7               #Pas temporel
hb = 1                  #Constante de Planck réduite
m = 1                   #Masse



#Création du potentiel
X, Y = np.meshgrid(np.linspace(0, 1, Lx), np.linspace(0, 1, Ly))
mux = 1/2
muy = 1/2
sigmax = 1/10
sigmay = 1/10
V = -1e4*np.exp(-((X-mux)**2/(2*sigmax**2)+(Y-muy)**2/(2*sigmay**2)))



#Création du tableau des solutions
psi0 = np.zeros((Nbi, Ly, Lx))
psii = np.sin(np.pi*X)*np.sin(np.pi*Y)     #Condition initiale
psi0[0, :, :] = psii



def finite_difference_method(Z):
    ''' Calcule la fonction d'onde pour chaque itération temporelle
        ----------
        :param Z: 3D array, tableau des solutions vide
        :return: Z, tableau des solutions
        ----------
    '''
    for k in range(0, Nbi-1):
        for i in range(1, Ly-1):
            for j in range(1, Lx-1):
                Z[k + 1, i, j] = 1j*(hb*dt)/(2*m) * ((Z[k][i+1][j] + Z[k][i-1][j] - 2*Z[k][i][j])/dx**2  + (Z[k][i][j+1] + Z[k][i][j-1] - 2*Z[k][i][j])/dy**2) + Z[k][i][j]*(1 - 1j * (V[i][j]*dt)/hb)
    return Z



psi = finite_difference_method(psi0.astype(complex))



#Calcul de la densité de probabilité normalisée
psicarre = np.absolute(psi)**2

for i in range(0, Nbi):
    psicarre[i, :, :] = psicarre[i,:, :]/np.sum(psicarre[i, :, :])





#Affiche le résultat

#Plot (Image 2D)

plt.style.use('dark_background')
plt.figure()
X, Y = np.meshgrid(np.arange(0, Lx), np.arange(0, Ly))
imageNbi = 30000
plt.contourf(X, Y, psicarre[imageNbi, :, :], 100, cmap = plt.cm.viridis)
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Densité de probabilité à t = " + str(round(imageNbi*dt,5)) + " s")
plt.savefig('2D_Schrodinger_Equation.png')



#Plot (Animation 2D)

plt.style.use('dark_background')
fig = plt.figure()

def animate(k):
    k=k*100
    plt.clf()
    plt.pcolormesh(psicarre[k, :, :], cmap = plt.cm.viridis)
    plt.colorbar()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Densité de probabilité à t = {k*dt:.5f} s")
    return

anim = animation.FuncAnimation(fig, animate, frames = int(Nbi/100), interval = 50, repeat = True)

#plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'
#Writer = animation.writers['ffmpeg']
#writermp4 = Writer(fps=30, bitrate=1800)
#anim.save("2D_Schrodinger_Equation.mp4", writer=writermp4)
writergif = animation.PillowWriter(fps=30)
writergif.setup(fig, "2D_Schrodinger_Equation.gif")
anim.save("2D_Schrodinger_Equation.gif", writer=writergif)



#Plot (Image 3D)
    
##plt.style.use('dark_background')
##fig = plt.figure()
##ax = plt.axes(projection = '3d')
##X, Y = np.meshgrid(np.arange(0, Lx), np.arange(0, Ly))
##imageNbi = 30000
##surf = ax.plot_surface(X, Y, psicarre[imageNbi, :, :], cmap=plt.cm.viridis)
##ax.set_xlabel("x")
##ax.set_ylabel("y")
##ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
##ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
##ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
##ax.grid(False)
##plt.title("Densité de probabilité à t = " + str(round(imageNbi*dt,5)) + " s")
##plt.savefig('3D_Schrodinger_Equation.png')



#Plot (Animation 3D)

##plt.style.use('dark_background')
##fig = plt.figure()
##ax = plt.axes(projection = '3d')
##X, Y = np.meshgrid(np.arange(0, Lx), np.arange(0, Ly))
##V = 2e-6*V
##
##def Animate3D(k):
##    k=k*100
##    ax.clear()
##    ax.set_zlim3d(0, np.max(psicarre))
##    #ax.set_zlim3d(np.min(V), np.max(psicarre))
##    #ax.plot_surface(X, Y, V, cmap=plt.cm.gray)
##    ax.plot_surface(X, Y, psicarre[k, :, :], cmap=plt.cm.viridis)
##    ax.set_xlabel("x")
##    ax.set_ylabel("y")
##    #fig.set_facecolor('black')
##    #ax.set_facecolor('black')
##    ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
##    ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
##    ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
##    ax.grid(False)
##    plt.title(f"Densité de probabilité à t = {k*dt:.5f} s")
##    ax.view_init(azim=k/100)
##    return
##
##anim3D = animation.FuncAnimation(fig, Animate3D, frames = int(Nbi/100), interval = 50, blit = False, repeat = True)
##
####plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'
####Writer = animation.writers['ffmpeg']
####writermp4 = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
####anim3D.save("3D_Schrodinger_Equation.mp4", writer=writermp4)
##writergif = animation.PillowWriter(fps=30)
##anim3D.save("3D_Schrodinger_Equation.gif", writer=writergif)


print("--- %s seconds ---" % (time.time() - start_time))
plt.show()
