import pandas as pd
import numpy as np
from numpy.linalg import solve
import matplotlib.pylab as plt
from matplotlib.font_manager import FontProperties
import pylab as pl
import os

#Parameters
S = 0.00035 #S trong sách dm: trang 220 S=Ssxb=(10^-5ft)x35ft=3.5x10^-4
T = 150 # K = 50m/day; b thickness = 30n => T=Kb = 50x30 = 150m^2/day

#savefolder
savefolder = r'C:\Users\Admin\Desktop\Kai\Figure 1_Homo_transient'
if not os.path.exists(savefolder):
    os.makedirs(savefolder)
os.chdir(savefolder)

#Domain
Nx = 51
Ny = 51
Lx = 1000
Ly = 1000
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)

#delta steps
dt = 0.01

# Coeffient
A = (S*dx*dx)/(T*dt)

#Grid
x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)
X, Y = np.meshgrid(x, y)

#Initial matrices 
Nt = 90 #time steps
N = Nx*Ny # unknowns
M = np.zeros((N,N)) #N rows, N columns
Bm = np.zeros(N) # N elements 
phi_vec = np.zeros(N) # N elements
result = []
#Set up interior nodes
for k in range (Nt):
	for i in range (1, Nx-1):
		for j in range (2, Ny):
			n = i+(j-1)*(Nx) # convert ij grid point to the nth grid point
			M [n,n] = -(4+A) #main diogonal
			M [n,n-1] = 1 # off diagonal to the left
			M [n,n+1] = 1 # off diagonal to the right
			M [n, n-Nx] = 1 #far off diagonal to the left
			M [n, n+Nx] = 1 #far off diagonal to the right
			Bm [n] = -A*phi_vec[n]
			# if (i==25 and j==25):
			# 	Bm[n] = -(A*phi_vec[n]+67500/(dx*dy))
			# 	#Bm[n] = -A*phi_vec[n]+500/(dx*dy)*50
			# else:
			# 	Bm[n] = -A*phi_vec[n]
			
	#Left
	i = 0
	for j in range (1,Ny+1):
		n = i+(j-1)*Nx # nth row for this ij
		M[n,n] = 1 #main diagonal
		Bm[n] = 50

	#Right
	i = Nx-1
	for j in range (1,Ny+1):
		n = i+(j-1)*Nx # nth row for this ij
		M[n,n] = 1 #main diagonal
		Bm[n] = 5

	#Bottom
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
	phi_vec = new_phi_vec # create phi_vec to update the value of Bm in the loop (interior nodes loop)
	phi = np.zeros((Nx,Ny))
	
	
	# Rearrange value in domain
	for j in range(1,Ny+1): 
		for i in range(1, Nx + 1):
			n = j + (i-1)*Nx-1
			phi[i-1, j-1] = phi_vec[n]
			
	result.append(phi[Nx-1,:])

			
	#Plot
	fig = plt.figure(figsize=[5,4])
	ax = fig.add_subplot(111)
	colorinterpolation = 50
	colourMap = plt.cm.jet
	contour = ax.contourf(X, Y, phi, colorinterpolation, cmap=colourMap)
	ax.set_title('Time at = {}'.format(round(dt*k*100,3)))
	plt.colorbar(contour, label='Hydraulic head (m)', shrink=0.75)
	plt.xlabel('Distance (m)')
	plt.ylabel('Distance (m)')
	plt.tight_layout()
	plt.savefig(os.path.join(savefolder, 'PDEs_Solver_{}.png'.format(str(k+1).zfill(3))), dpi=100, facecolor='white', edgecolor='k', transparent=False)
	plt.close()
# end loop
#phi1 = phi[Nx-1,:]
# for p in range (Nt):
# 	#fig, ax = pl.subplots(1, 1)
# 	plt.plot(x,result[p], label = 'Classical method')
# 	plt.xlabel('Distance (km)')
# 	plt.ylabel('Groundwater head (m)')
# 	plt.legend()
# 	plt.show()

	
	#print(result[p])
fontP = FontProperties()
fontP.set_size('xx-small')

plt.plot(x,result[0], label = 't = 1 day')
plt.plot(x,result[9], label = 't = 10 day')
plt.plot(x,result[19], label = 't = 20 day')
plt.plot(x,result[29], label = 't = 30 day')
plt.plot(x,result[39], label = 't = 40 day')
plt.plot(x,result[49], label = 't = 50 day')
plt.plot(x,result[59], label = 't = 60 day')
plt.plot(x,result[69], label = 't = 70 day')
plt.plot(x,result[79], label = 't = 80 day')
plt.plot(x,result[89], label = 't = 90 day')
plt.xlabel('Distance (m)')
plt.ylabel('Groundwater head (m)')
plt.legend()
#plt.legend(title='Legend', bbox_to_anchor=(1.0, 1), loc='upper left')
plt.show()