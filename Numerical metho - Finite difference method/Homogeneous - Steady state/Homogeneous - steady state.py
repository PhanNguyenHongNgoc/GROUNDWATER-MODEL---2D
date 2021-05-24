
import numpy as np
import scipy 
from numpy.linalg import solve
import matplotlib.pylab as plt

#Inputs
Nx = 51
Ny = 51
Lx = 1000
Ly = 1000
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)

#Grid
x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)
X, Y = np.meshgrid(x, y)

#Initial matrices 
N = Nx*Ny #of unknowns
N_ma = np.zeros((N,N))
M = np.zeros((N,N)) #N rows, N columns
B = np.zeros((N,1)) # N rows, 1 column

#Set up interior nodes
for i in range (1, Nx-1):
	for j in range (2,Ny):
		n = i+(j-1)*(Nx) # convert ij grid point to the nth grid point
		M [n,n] = -4 #main diogonal
		M [n,n-1] = 1 # off diagonal to the left
		M [n,n+1] = 1 # off diagonal to the right
		M [n, n-Nx] = 1 #far off diagonal to the left
		M [n, n+Nx] = 1 #far off diagonal to the right
		B [n] = 0 
		if (i== 20 and j== 20):
			B[n] = -(67500/(dx*dy))
		else:
			B[n] = 0
		
#Left (Dirichlet boundary)
i = 0 
for j in range (1,Ny+1):
	n=i+(j-1)*Nx # nth row for this ij
	M[n,n]=1 #main diagonal
	B[n]=50
	
#Right (Dirichlet boundary)
i = Nx-1
for j in range (1,Ny+1):
	n=i+(j-1)*Nx # nth row for this ij
	M[n,n]=1 #main diagonal
	B[n]=5
	
#Bottom (Neuman boundary for no flow boundary)
j = 1
for i in range (0,Nx):
	n=i+(j-1)*Nx # nth row for this ij
	M[n,n]=-1 #main diagonal
	M[n,n+Nx]=1
	B[n]=0

#Top (Neuman boundary for no flow boundary)
j = Ny
for i in range (0, Nx):
	n=i+(j-1)*Nx # nth row for this ij
	M[n,n]=1 #main diagonal
	M[n,n-Nx]=-1
	B[n]= 0

#phi = scipy.sparse.linalg.spsolve(M,B)
phi_vec = solve(M, B)
phi = np.zeros((Nx,Ny))

#rearrange value in the domain
for j in range(1,Ny+1): 
    for i in range(1, Nx + 1): 
    	n = j + (i-1)*Nx-1
    	phi[i-1, j-1] = phi_vec[n]

#Plot
colorinterpolation = 50
colourMap = plt.cm.jet
fig, ax = plt.subplots()


CS_classical = plt.contour(X,Y, phi, colors = 'white')
plt.clabel(CS_classical, inline=0.01, fontsize=10 )

plt.contourf( X, Y, phi, colorinterpolation, cmap=colourMap)
ax.set_title('Classical method')


plt.colorbar(label='Hydraulic head (m)')
plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')
plt.show()

positions = np.vstack([Y.ravel(), X.ravel()])

is_perimeter = (positions[0] == np.min(y)) | (positions[0] == np.max(y)) | \
               (positions[1] == np.min(x)) | (positions[1] == np.max(x))

plt.scatter(*positions[::-1], c=is_perimeter)
plt.show()

#SREAMLINE
u = np.zeros((Nx,Ny))
v = np.zeros((Nx,Ny))

for i in range(1,Nx-1): 
    for j in range(1, Ny-1): 
    	u[i,j] = -(phi[i,j+1]-phi[i,j-1])/(2*dx)
    	v[i,j] = -(phi[i+1,j]-phi[i-1,j])/(2*dy)


    	# u[i,j] = -(phi[i+1,j]-phi[i-1,j])/(2*dx)
    	# v[i,j] = -(phi[i,j+1]-phi[i,j-1])/(2*dy)

plt.figure(figsize=(7,3.85))
# plt.streamplot(X,Y,u,v) 
plt.streamplot(X,Y,u.transpose(),v.transpose()) 
#plt.quiver(X,Y,u.transpose(),v.transpose())
# plt.streamplot(X,Y,u.T,v.T)
plt.show()