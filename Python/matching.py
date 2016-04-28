import numpy as np
import random
from numpy.linalg import det
from determinant import determinant
from inverse_test import inverse, augment, generate_matrix

def test(num_tests, N, m, a, b):
    for test in xrange(0, num_tests):
        A = generate_matrix(N,0,0)
        for counter in xrange(0,m):
            i = random.randint(0, N-1)
            j = random.randint(0, N-1)
            val = float(random.randint(0,2*m))
            if i < j:
                A[i][j] = val
                A[j][i] = -1.*val
            else:
                A[j][i] = -1.*val
                A[i][j] = val

        print determinant(A,N)


test(10,88,87,0,24)





