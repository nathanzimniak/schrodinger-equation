#!/usr/bin/env python3
# coding: utf-8

'''
Programme pour résoudre l'équation de Schrödinger dépendante du temps par la méthodes des différences finies
'''

__author__ = 'Nathan Zimniak'



import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import eigsh
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import time
start_time = time.time()




#Initialisation des constantes
N = 100
xmax = 0.5
xmin = -0.5
dx = 0.05
X, Y = np.meshgrid(np.linspace(xmin,xmax,N), np.linspace(xmin,xmax,N))
hb = 1
m = 1
omega = 2*np.pi



#Création du potentiel
def get_potential(x, y):
    return (1/2)*(omega**2)*(x**2 + y**2) #np.exp(-(x-0.3)**2/(2*0.1**2))*np.exp(-(y-0.3)**2/(2*0.1**2))
V = get_potential(X,Y)



#Création du Hamiltonien
diags = np.array([np.ones(N), -2*np.ones(N), np.ones(N)])
D = sparse.spdiags(diags, (-1,0,1), N, N)
T = (-hb**2)/(2*m) * sparse.kronsum(D,D)/dx**2
U = sparse.diags(V.reshape(N**2), (0))
H = T+U



#Calcul des valeurs et vecteurs propres du Hamiltonien
n=11
E, psi = eigsh(H, k=n, which='SM')

def nth_psi(n):
    npsi = psi.T[n].reshape((N,N))
    return npsi

def nth_psicarre(n):
    psicarre = np.absolute(nth_psi(n))**2
    return psicarre





#Affiche le résultat

#Plot (Animation 2D)

##plt.style.use('dark_background')
##fig = plt.figure()
##r=1000
##def animate(k):
##    plt.clf()
##    plt.contourf(X, Y, nth_psicarre(int(k*n/r)), 100, cmap = plt.cm.viridis)
##    plt.colorbar()
##    plt.xlabel("x")
##    plt.ylabel("y")
##    plt.title("Densité de probabilité pour l'état n = " + str(int(k*n/r)))
##    return
##
##anim = animation.FuncAnimation(fig, animate, frames = r, interval = 500, repeat = True)
##
##plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'
##Writer = animation.writers['ffmpeg']
##writermp4 = Writer(fps=30, bitrate=1800)
##anim.save("2D_Time_Independant_Schrodinger_Equation.mp4", writer=writermp4)
##writergif = animation.PillowWriter(fps=30)
##writergif.setup(fig, "2D_Time_Independant_Schrodinger_Equation.gif")
##anim.save("2D_Time_Independant_Schrodinger_Equation.gif", writer=writergif)



#Plot (Animation 3D)

plt.style.use('dark_background')
fig = plt.figure()
ax = plt.axes(projection = '3d')
r=1000
def Animate3D(k):
    ax.clear()
    ax.set_zlim3d(0, np.max(nth_psicarre(int(k*n/r))))
    ax.zaxis.set_rotate_label(False)
    #ax.set_zlim3d(np.min(V), np.max(psicarre))
    #ax.plot_surface(X, Y, V, cmap=plt.cm.gray)
    ax.plot_surface(X, Y, nth_psicarre(int(k*n/r)), cmap=plt.cm.viridis)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("$|\Psi|^2$", rotation=0)
    #fig.set_facecolor('black')
    #ax.set_facecolor('black')
    ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.grid(False)
    plt.title("Densité de probabilité pour l'état n = " + str(int(k*n/r)))
    ax.view_init(azim=k)
    return

anim3D = animation.FuncAnimation(fig, Animate3D, frames = r, interval = 500, blit = False, repeat = True)

plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'
Writer = animation.writers['ffmpeg']
writermp4 = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
anim3D.save("3D_Time_Independant_Schrodinger_Equation.mp4", writer=writermp4)
writergif = animation.PillowWriter(fps=30)
anim3D.save("3D_Time_Independant_Schrodinger_Equation.gif", writer=writergif)


print("%s secondes" % (time.time() - start_time))
plt.show()
