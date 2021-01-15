import numpy as np
from prox_ml1 import prox_ml1
import time
from matplotlib import pyplot as plt
from matplotlib import cm

#%% Script for testing the proximal operator of the matrix induced l1 norm


#% Dummy test
# toto = [ 3, 2, 1, 4.1 ; 2, 1, 0 , 3; 1, 0, 0, 2 ];
toto = np.random.randn(100,10)
[U, t, nu]= prox_ml1(toto,10,tol=1e-10)
print(U,t,nu)


out = np.zeros([10,10])
# Parameters
#lamb = 10

for n in range(1,11):
    print(1000*n)
    for m in range(1,11):
        temp = np.zeros(5)
        for i in range(5):
            # Data
            X = np.random.randn(1000*n,10*m)
            # 50% max lambda
            lamb = 0.5*np.sum(np.max(np.abs(X),axis=0))
            # Prox vs sum
            clock_alg=time.time()
            prox_ml1(X,lamb,tol=1e-8)
            t=time.time() - clock_alg
            # output
            temp[i] = t#/(10*m*(np.log2(100*n)**2))
        out[n-1,m-1] = np.mean(temp)

fig = plt.figure()
ax = fig.gca(projection='3d')
Xgrid = np.arange(1000, 11000, 1000)
Ygrid = np.arange(10, 110, 10)
Xgrid, Ygrid = np.meshgrid(Xgrid, Ygrid)
surf = ax.plot_surface(Xgrid,Ygrid,out.T, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#ax.axis([100, 1000, 10, 100])
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('n')
ax.set_ylabel('m')
ax.set_zlabel('t')
#fig.title('Effective computation time for the prox of l11')
plt.show()


#c = colorbar;
#c.Label.String = 'log_{10}(t/nm log_2(n))';
#title('Effective computation time of prox of l11')
