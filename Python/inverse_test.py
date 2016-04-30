import numpy as np
from copy import deepcopy
import random
from numpy.linalg import det

def not_in_list(item_list, item):
    return item not in item_list

def in_list(item_list, item):
    return item in item_list

def inverse(Q, N, excl_rows, excl_cols, err):
    C = deepcopy(Q)
    n = N
    j = 0
    while (j < N):

        #Exclude columns that are to be excluded
        if in_list(excl_cols, j):
            j += 1
            continue

        #IF... add a non-zero row to row j
        if abs(C[j][j]) < err:
            k = None
            #Find the non-zero k
            ##pragma omp parallel for
            for new_col in xrange(j+1,N):
                if not_in_list(excl_cols, new_col) and abs(C[new_col][j])>err:
                    k = new_col
                    break

            #Add the row k to row j
            for row in xrange(0, 2*N):
                if not_in_list(excl_rows, row):
                    C[j][row] +=  C[k][row]

        ajj = C[j][j]
        
        #Divide out ajj
        ##pragma omp parallel for
        for row in xrange(0, 2*N):
            if not_in_list(excl_rows, row):
                C[j][row] /= ajj

        #Subtract row j multiplied by appropriate constant from other rows
        for col in xrange(0,N):
            if (col!=j) and not_in_list(excl_cols, col):
                aij = C[col][j]
                for row in xrange(0, 2*N):
                    if not_in_list(excl_rows, row):
                        C[col][row] = C[col][row] - aij*C[j][row]
        j += 1
        #print "finished row: ", C[j]

    return C

def get_inv_from_C(C, N):
    A_inv = [[] for i in xrange(0,N)]
    for i in xrange(0,N):
        for j in xrange(0,N):
            A_inv[i].append(C[i][j+N])

    return A_inv

def generate_matrix(N,a,b):
    return [[float(random.randint(a,b)) for i in range(N)] for j in range(N)]

def my_inv_test(A, N, my_inv, excl, err):
    a = np.array(A)
    I = a.dot(np.array(my_inv))

    error = False
    for i in xrange(0,N):
        if i in excl:
            continue

        for j in xrange(0,N):
            if j in excl:
                continue

            if j == i:
                if abs(I[i][j] - 1) > err:
                    print "ERROR for 1 - (i,j):", i, j
                    print I[i][j]
                    return None
                    
            else:
                if abs(I[i][j] - 0) > err:
                    print "ERROR for 0 - (i,j):", i, j
                    print I[i][j]
                    return None

    if not error:
        print "INV TEST SUCCEEDED"
    return I

def augment(A,N):
    C = [[] for i in xrange(0, len(A))]
    for i in xrange(0, len(A)):
        for j in xrange(0, len(A)*2):
            if j<len(A):
                C[i].append(A[i][j])
            elif j == len(A) + i:
                C[i].append(1.)
            else:
                C[i].append(0.)
    return C

def test(num_tests, N, a, b, err):
    for i in xrange(0, num_tests):
        A = generate_matrix(N,a,b)
        excl = [0,1,N-1]
        npA = np.array(deepcopy(A))
        npA = np.delete(npA, (excl), axis=0)
        npA = np.delete(npA, (excl), axis=1)

        if abs(det(np.array(npA))) < err:
            print "0 DETERMINANT"
            continue

        C = augment(A,N)
        my_inv = get_inv_from_C(inverse(C, len(A), excl, excl, err), N)
        if my_inv_test(A, N, my_inv, excl, err) is None:
            return
        

def det_test():
    A = [[1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0], [1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0], [1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0], [1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0]]
    print det(A)

#det_test()
test(200,40,0,1,10**-6)