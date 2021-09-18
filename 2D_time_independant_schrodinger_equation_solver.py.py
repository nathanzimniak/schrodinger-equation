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



#Initialisation des constantes
N = 100
xmax = 0.5
xmin = -0.5
dx = 0.05
X, Y = np.meshgrid(np.linspace(xmin,xmax,N), np.linspace(xmin,xmax,N))
hb = 1
m = 1
omega = 2*np.pi

def get_potential(x, y):
    return (1/2)*(omega**2)*(x**2 + y**2) #np.exp(-(x-0.3)**2/(2*0.1**2))*np.exp(-(y-0.3)**2/(2*0.1**2))

V = get_potential(X,Y)

diags = np.array([np.ones(N), -2*np.ones(N), np.ones(N)])
D = sparse.spdiags(diags, (-1,0,1), N, N)
T = (-hb**2)/(2*m) * sparse.kronsum(D,D)/dx**2

U = sparse.diags(V.reshape(N**2), (0))

H = T+U

n=21
eigenvalues, eigenvectors = eigsh(H, k=n, which='SM')

def ev(n):
    evn = eigenvectors.T[n].reshape((N,N))
    return evn

def evc(n):
    evcn = np.absolute(ev(n))**2
    return evcn

##plt.style.use('dark_background')
##plt.figure()
##for i in range(0,2):
##    plt.contourf(X, Y, evc(i)**2, 100, cmap = plt.cm.viridis)
##    plt.colorbar()
##    plt.xlabel("x")
##    plt.ylabel("y")
##    plt.title("Densité de probabilité à n = " + str(i) + "\n d'énergie E = " + str(eigenvalues[i]))
##    plt.show()



#Plot (Animation 2D)

plt.style.use('dark_background')
fig = plt.figure()

def animate(k):
    plt.clf()
    plt.contourf(X, Y, evc(k)**2, 100, cmap = plt.cm.viridis)
    plt.colorbar()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Densité de probabilité à n = " + str(k) + "\n d'énergie E = " + str(eigenvalues[k]))
    return

anim = animation.FuncAnimation(fig, animate, frames = int(n), interval = 500, repeat = True)

plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'
Writer = animation.writers['ffmpeg']
writermp4 = Writer(fps=1, bitrate=1800)
anim.save("2D_Schrodinger_Equation.mp4", writer=writermp4)
writergif = animation.PillowWriter(fps=1)
writergif.setup(fig, "2D_Time_Independant_Schrodinger_Equation.gif")
anim.save("2D_Time_Independant_Schrodinger_Equation.gif", writer=writergif)

plt.show()
