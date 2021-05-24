import pandas as pd
import numpy as np
from numpy.linalg import solve
import matplotlib.pylab as plt
import os

#Parameters
S = 1
T1 = 1
T2 = 5

#savefolder
savefolder = r'C:\Users\Admin\Desktop\Kai\Figure 2_Hete-isotropic-transient'
if not os.path.exists(savefolder):
    os.makedirs(savefolder)
os.chdir(savefolder)

#Domain
Nx = 51
Ny = 51
Lx = 1
Ly = 1
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)

# Time step
dt = 0.01

# Coeffient
A1 = (S*dx*dy)/(T1*dt)
A2 = (S*dx*dy)/(T2*dt)
#B2 = (dx*dy)/T2 #source/sink

#Grid
x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)
X, Y = np.meshgrid(x, y)

#Initial matrices 
Nt = 30
N = Nx*Ny #of unknowns
M = np.zeros((N,N)) #N rows, N columns
Bm = np.zeros(N) # N elements 
phi_vec = np.zeros(N) # N elements
Head = []
t_hour = []

#Set up interior nodes
for k in range (Nt):
	for i in range (1, Nx-1):
		for j in range (2, Ny):
			n = i+(j-1)*(Nx) # convert ij grid point to the nth grid point
			M [n,n] = -(4+A1) #main diogonal
			M [n,n-1] = 1 # off diagonal to the left
			M [n,n+1] = 1 # off diagonal to the right
			M [n, n-Nx] = 1 #far off diagonal to the left
			M [n, n+Nx] = 1 #far off diagonal to the right
			Bm [n] = -A1*phi_vec[n]

	for i in range (int(Nx/2), Nx-1):
		for j in range (int(Ny/2), Ny):
			n = i+(j-1)*(Nx) # convert ij grid point to the nth grid point
			M [n,n] = -(4+A2) #main diogonal
			M [n,n-1] = 1 # off diagonal to the left
			M [n,n+1] = 1 # off diagonal to the right
			M [n, n-Nx] = 1 #far off diagonal to the left
			M [n, n+Nx] = 1 #far off diagonal to the right
			Bm [n] = -A2*phi_vec[n]
			
			#Souce/sink
			# if (i==45 and j==45):
			# 	Bm[n] = -A2*phi_vec[n]-0.01/(dx*dy)
			# 	#Bm[n] = -A2*phi_vec[n]-500*B2
			# else:
			# 	Bm[n] = -A2*phi_vec[n]
				
	#Left	
	i = 0
	for j in range (1,Ny+1):
		n = i+(j-1)*Nx # nth row for this ij
		M[n,n] = 1 #main diagonal
		Bm[n] = 100

	#Right
	i = Nx-1
	for j in range (1,Ny+1):
		n = i+(j-1)*Nx # nth row for this ij
		M[n,n] = 1 #main diagonal
		Bm[n] = 10


	#Bot
	j = 1
	for i in range (0,Nx):
		n = i+(j-1)*Nx # nth row for this ij
		M[n,n] = -1 #main diagonal
		M[n,n+Nx] = 1
		Bm[n] = 0
    
	#Top
	j = Ny
	for i in range (0, Nx):
		n=i+(j-1)*Nx # nth row for this ij
		M[n,n] = 1 #main diagonal
		M[n,n-Nx] = -1
		Bm[n] = 0

	#Solve the matrix
	new_phi_vec = solve(M, Bm)
	phi_vec = new_phi_vec #create phi_vec to update value Bm in the interior loop
	phi = np.zeros((Nx,Ny))
	
    #Rearrange valu in domain
	for j in range(1,Ny+1): 
		for i in range(1, Nx + 1):
			n = j + (i-1)*Nx-1
			phi[i-1, j-1] = phi_vec[n]
	
	#Plot
	fig = plt.figure(figsize=[5,4])
	ax = fig.add_subplot(111)
	colorinterpolation = 50
	colourMap = plt.cm.jet
	contour = ax.contourf(X, Y, phi, colorinterpolation, cmap=colourMap)
	ax.set_title('Time at = {}'.format(round(dt*k,3)))
	plt.colorbar(contour, label='Hydraulic head (m)', shrink=0.75)
	plt.tight_layout()
	plt.savefig(os.path.join(savefolder, 'PDEs_Solver_{}.png'.format(str(k+1).zfill(3))), dpi=100, facecolor='white', edgecolor='k', transparent=False)
	plt.close()
# end loop
