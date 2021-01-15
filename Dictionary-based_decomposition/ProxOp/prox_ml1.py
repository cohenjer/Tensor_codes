import numpy as np
from numpy.matlib import repmat as repmat

def ml1(X):
    '''
    Computes the induced matrix l1 norm of X
    '''
    return np.max(np.sum(np.abs(X), axis=0))

def SoftT(x, lamb):
    '''
    Computes the Soft-Thresholding of input vector x, with coefficient lamb
    '''
    return np.maximum(np.abs(x) - lamb, 0)*np.sign(x)

def prox_ml1(X, lamb, tol=1e-10):
    '''
    Computes the proximal operator of the matrix induced l1 norm

    Inputs:

    X    :  numpy array, input of the proximal operator
    lamb :  float, regularization parameter
    tol  :  float, small tolerance on the value of the maximal columnwise 1 norm

    Outputs:

    U    :  numpy array, the proximal operator of X
    t    :  maximum l1 norm of the columns of U
    nu_t :  optimal dual parameters


    Credits to Jeremy E. Cohen, v0 06/01/2021
    Reference: ``Computing the proximal operator of the l1 induced matrix norm, J.E.Cohen, arxiv:2005.06804v2''.
    '''

    # Remove single column case
    if X.ndim==1:
        out = SoftT(B, lamb)
        return out, np.sum(np.abs(out)), 1

    #%% Precomputations

    # Constants
    n, m  = np.shape(X)
    ps    = np.linspace(1,n,n)
    Xsort = np.sort(np.abs(X), axis=0)
    Xsort = np.flip(Xsort, axis=0)
    Xcumsum=np.cumsum(Xsort, axis=0)

    # Maximal value of t (prox is identity)
    Xsum  = np.sum(Xsort, axis=0)  # Xsort for abs value
    tmax = np.max(Xsum)
    tmin = 0
    t = tmax/2  # initial value in the middle of the admissible interval

    # Find the order of visited columns in X
    sorted_sums = np.flip(np.sort(Xsum))
    order = np.flip(np.argsort(Xsum))

    # Deal with the lamb>lamb_max case in advance to avoid bug
    if lamb>=np.sum(Xsort[0,:]):
        U = np.zeros([n,m])
        t = 0
        nu_t = Xsort[0,:]/lamb
        return U, t, nu_t

    #%% Starting bisection
    while tmax-tmin > tol:
        # Compute the current active set and its size
        I = order[sorted_sums>t]
        l = len(I)

        # Process commulative sums
        Xcumsum_t = (Xcumsum[:, I]-t)/lamb

        # Compute the candidate nu values
        Ps = np.transpose(repmat(ps, l, 1))  # matrix of sequences from 1 to n in with i columns
        nu = Xcumsum_t/Ps

        nu_t = np.zeros(m)
        N = Xsort[:, I] - lamb*nu
        # temp is reverse sorted
        N = np.flip(N, axis=0)  # temp is sorted
        for j in range(l):
            # Find the largest index for which the condition described in the paper is satisfied
            # i.e.  xsort(p) - \lambda\nu(p) >= 0
            idx = np.searchsorted(N[:,j], 0, side='left')
            idx = len(N[:,j]) - idx - 1  # counter effect flip
            nu_t[I[j]] = nu[idx, j]

        # end j loop

        # Checking dual condition 1< or 1> to move t
        if np.sum(nu_t) < 1:
            # t must be decreased
            tmax = t
            t = (t + tmin)/2
        else:
            # t must be increased
            tmin = t
            t = (t + tmax)/2

    # Final step, thresholding vectors that need to be
    U = SoftT(X, lamb*nu_t)

    return U, t, nu_t

if __name__ == '__main__':
    A = np.array([[1, 2, 5],[3, 4, 0], [1, 3, -1], [0.4, 2, 1]])
    out = prox_ml1(A, 11.9)
    print(out)

    B = np.array([[1, 2, 3],[2,3,4]])
    out2 = prox_ml1(B, 2)
    print(out2)
