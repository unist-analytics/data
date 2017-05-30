from math import sqrt
import numpy as np
from numpy.linalg import inv

def norm(x):

    return sqrt(sum([x_i**2 for x_i in x]))

def Q(A):

    m, n = np.shape(A) #m x n matrix
    col = []
    q = []

    for i in range(n):
        col.append(A[:, i:i+1]) # divide matrix into each column

    #print col
    for j in range(len(col)): # making orthogonal vector
        temp = col[j]
        if j == 0:
            q.append(col[j])
            continue

        for k in range(len(q)):
            #print(temp)
            temp = temp - (q[k].T).dot(col[j]) / (q[k].T).dot(q[k]) * q[k] # gram-schmidt
            #print (q[k].T).dot(col[j])
            #print (q[k].T).dot(q[k])
            #print((q[k].T).dot(col[j])[0][0] / (q[k].T).dot(q[k])[0][0])
        q.append(temp)
    #print q

    norm_q = []
    for a in range(len(q)): # Normalization
        temp1 = q[a] / norm(q[a])
        norm_q.append(temp1)
    #print(norm_q)

    value = []
    for f in range(len(norm_q)): #Reshape Q
        for g in range(len(norm_q[f])):
            value.append(norm_q[g][f][0])
    Q = np.array(value).reshape((m,n))



    return Q

def R(Q):

    Qinv = inv(Q)
    R = np.dot(Qinv, A) # A = QR

    return R

def Backsubstitution(Q,R,b):

    Rinv = inv(R)
    x_1 = np.dot(Rinv, Q.T)
    x = np.dot(x_1, b)

    return x


if __name__ == '__main__':

    #Ax = b, find x
    A = np.array([[-3,2,1],
                  [-1,4,1],
                  [1,-3,1]], dtype=np.float64)
    b = np.array([[4],
                  [3],
                  [2]], dtype=np.float64)


    print("QR Factorization\n")
    print("A:")
    print A
    print("Q:")
    print Q(A).round(3)
    print('R:')
    print R(Q(A)).round(3)
    print('x:')
    print Backsubstitution(Q(A), R(Q(A)), b).round(3)

    #print(np.shape(x)[0])
    #print b