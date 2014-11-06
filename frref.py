import numpy as np
import numpy.linalg as nlg
import scipy
import scipy.sparse as sp
import scipy.linalg as splg

def frref(A):
    m = A.shape[0]
    n = A.shape[1]
    
    r = nlg.qr(A, mode = 'r')
    indep_rows = (r!=0).max(axis=1)
    indep_rows = np.ravel(indep_rows)

    i_dep = np.nonzero(indep_rows)[0]
    i_indep = np.nonzero(1-indep_rows)[0]

    L_indep = np.min(i_indep.shape)
    L_dep   = np.min(i_dep.shape)

    F = sp.lil_matrix((m,n))
    if L_indep != 0:
        tmp = nlg.lstsq(r[indep_rows,:][:,i_dep],\
                        r[indep_rows,:][:,i_indep])[0]
        F[indep_rows,i_indep] = tmp
    if L_dep != 0:
        F[indep_rows,i_dep] = np.eye(L_dep)
    
    if sp.issparse(A):
        F.to_csr()
    else:
        F = F.todense()
    return F

X1 = np.array([[1,1,0],[1,1,0],[0,0,1]])
X2 = np.array([[1,2,3],[1,2,4],[3,2,1]])
X3 = sp.lil_matrix((10,10))
X3[0,0:2]  = np.array([1,1])
X3[1,0:2]  = np.array([1,1])
X3[2,0:3]  = np.array([0,0,1])
X3[3,1:4]  = np.array([1,2,3])
X3[4,2:5]  = np.array([1,2,3])
X3[5,3:6]  = np.array([1,2,3])
X3[6,4:7]  = np.array([1,2,3])
X3[7,5:8]  = np.array([1,2,3])
X3[8,6:9]  = np.array([1,2,3])
X3[9,7:10] = np.array([1,2,3])

print X3.todense()
print X3.shape


Y1 = frref(X1)
print 'X1 -> Y1'
print X1
print Y1


Y2 = frref(X2)
print 'X2 -> Y2'
print X2
print Y2

Y3 = frref(X3)
print 'X3 -> Y3'
print X3
print Y3
