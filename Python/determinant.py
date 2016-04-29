import numpy as np
import random
from numpy.linalg import det
from copy import deepcopy

def determinant(A,N):
    det = 1
    for rc in xrange(0,N):
        #swap rows
        if A[rc][rc] == 0:
            #print "swapping rows"
            new_col = None
            for i in xrange(rc+1,N):
                if A[i][rc] != 0:
                    new_col = i
                    break

            if new_col != None:
                det = det * -1.
                for i in xrange(0,N):
                    temp = A[rc][i]
                    A[rc][i] = A[new_col][i]
                    A[new_col][i] = temp

        for row_below in xrange(rc+1,N):
            if A[row_below][rc] != 0:
                val = A[row_below][rc]/A[rc][rc]
                for i in xrange(0,N):
                    A[row_below][i] -= val*A[rc][i]

    for i in xrange(0,N):
        det *= A[i][i]

    return det

def generate_matrix(N,a,b):
    return [[float(random.randint(a,b)) for i in range(N)] for j in range(N)]

def test(num_tests, N, a, b, err):
    for i in xrange(0, num_tests):
        A = generate_matrix(N,a,b)

        my_det = determinant(deepcopy(A),N)
        np_det = det(np.array(deepcopy(A)))
        
        val1 = abs(my_det) < err
        val2 = abs(np_det) < err
        #print my_det
        #print np_det
        if val1 != val2:
            print "TEST FAILED"
            print A
            print my_det
            print np_det
            break
        else:
            print "TEST SUCCEEDED"


def specific_test():
    A = [[1.0, 1.0, 0.0, 1.0], [1.0, 0.0, 0.0, 0.0], [1.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0]]
    np_A = deepcopy(np.array(A))
    print determinant(A, len(A))
    print det(np_A)

#test(1000,40,0,2000,10**-15)         

