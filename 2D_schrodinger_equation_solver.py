##RÉSOLUTION DE L'ÉQUATION DE LA CHALEUR 2D



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import time

start_time = time.time()


#Initialisation des constantes et des variables

Lx = 31                 #Longueur de l'espace
Ly = 31                 #Largeur de l'espace
Nbi = 100000              #Nombre d'itérations temporelles
dx = 1/(Lx-1)            #Pas spatial suivant x
dy = 1/(Ly-1)            #Pas spatial suivant y
dt = 1e-7                #Pas temporel

X, Y = np.meshgrid(np.linspace(0, 1, Lx), np.linspace(0, 1, Ly))
mux, sigmax = 1/2, 1/10
muy, sigmay = 1/2, 1/10
V = -1e4*np.exp(-((X-mux)**2/(2*sigmax**2)+(Y-muy)**2/(2*sigmay**2)))

psii = np.sin(np.pi*X)*np.sin(np.pi*Y)     #Condition initiale

psi0 = np.zeros((Nbi, Ly, Lx))    #Création du tableau des conditions initiales
psi0[0, :, :] = psii


#Création d'une fonction qui retourne le tableau des solutions en appliquant la méthode des différences finies

def mdf(y):
    for k in range(0, Nbi-1):
        for i in range(1, Ly-1):
            for j in range(1, Lx-1):
                y[k + 1, i, j] = y[k][i][j] + 1j/2 * dt * ((y[k][i+1][j] + y[k][i-1][j] - 2*y[k][i][j])/dx**2  + (y[k][i][j+1] + y[k][i][j-1] - 2*y[k][i][j])/dy**2) - 1j * dt * V[i][j] * y[k][i][j]
    return y



#Création du tableau des solutions

psi = mdf(psi0.astype(complex))
psicarre = np.absolute(psi)**2

for i in range(0, Nbi):
    psicarre[i, :, :] = psicarre[i,:, :]/np.sum(psicarre[i, :, :])


print("%s secondes" % (time.time() - start_time))



###Plot (Image 2D)
##
##plt.figure(1)
##X, Y = np.meshgrid(np.arange(0, Lx), np.arange(0, Ly))
##imageNbi = 6000
##plt.contourf(X, Y, psicarre[imageNbi, :, :], 100, cmap = plt.cm.viridis)
##plt.colorbar()
##plt.xlabel("x")
##plt.ylabel("y")
##plt.title("Densité de probabilité \n de la fonction d'onde à t = " + str(round(imageNbi*dt,5)) + " s")





###Plot (Animation 2D)
##
##fig = plt.figure()
##
##def animate(k):
##    k=k*100
##    plt.clf()
##    plt.pcolormesh(psicarre[k, :, :], cmap = plt.cm.viridis)
##    plt.colorbar()
##    plt.xlabel("x")
##    plt.ylabel("y")
##    plt.title(f"Densité de probabilité \n de la fonction d'onde à t = {k*dt:.5f} s")
##    return
##
##anim = animation.FuncAnimation(fig, animate, frames = 100, interval = 50, repeat = True)
##
##writergif = animation.PillowWriter(fps=30)
##anim.save("2DSchrodingerEquation.gif", writer=writergif)





###Plot (Image 3D)
##
##fig = plt.figure(3)
##ax1 = plt.axes(projection = '3d')
##X, Y = np.meshgrid(np.arange(0, Lx), np.arange(0, Ly))
##imageNbi = 9999
##surf = ax1.plot_surface(X, Y, psicarre[imageNbi, :, :], cmap=plt.cm.viridis)
##fig.colorbar(surf, shrink=0.5)
##ax1.set_xlabel("x")
##ax1.set_xlabel("y")
##plt.title("Densité de probabilité \n de la fonction d'onde à t = " + str(round(imageNbi*dt,5)) + " s")





#Plot (Animation 3D)

fig = plt.figure()
ax = plt.axes(projection = '3d')

X, Y = np.meshgrid(np.arange(0, Lx), np.arange(0, Ly))

def Animate3D(k):
    k=k*100
    ax.clear()
    ax.set_zlim3d(0, np.max(psicarre))
    ax.plot_surface(X,Y,psicarre[k, :, :],cmap=plt.cm.viridis)
    ax.set_xlabel("x")
    ax.set_xlabel("y")
    plt.title(f"Densité de probabilité \n de la fonction d'onde à t = {k*dt:.5f} s")
    return
anim3D = animation.FuncAnimation(fig, Animate3D, frames = int(Nbi/100), interval = 50, blit = False, repeat = True)
writergif = animation.PillowWriter(fps=30)
anim3D.save("3DSchrodingerEquation.gif", writer=writergif)






plt.show()
