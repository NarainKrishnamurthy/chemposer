from copy import deepcopy
import numpy as np
import random
from numpy.linalg import det, inv
from determinant import determinant
from inverse_test import inverse, augment, generate_matrix, get_inv_from_C

def not_in_list(item_list, item):
    return item not in item_list

def in_list(item_list, item):
    return item in item_list

def set_weights(A,N,m):
    counter = 0
    while counter < m:
        i = random.randint(0, N-1)
        j = random.randint(0, N-1)

        if A[i][j] == 0:
            val = float(random.randint(0,2*m))
            if i < j:
                A[i][j] = val
                A[j][i] = -1.*val
            elif i > j:
                A[i][j] = -1.*val
                A[j][i] = val
            counter += 1

    return A


def np_matching(np_T, size_T, err):
    E = deepcopy(np_T)
    T = np.array(np_T)
    if det(T) < err:
        print "NUMPY - ZERO DET"
        return None

    M = []
    true_cols = range(0,size_T)
    true_rows = range(0, size_T)

    while len(M) < size_T/2:
        if det(T) == 0:
            print "SINGULAR"
            return None

        N = inv(T)

        true_row = true_rows[0]
        j_col = None
        for col in xrange(0,len(N)):
            if N[0][col] != 0 and E[true_row][true_cols[col]]!= 0:
                j_col = col
                break

        #print true_row, j_col, true_cols[j_col]
        e = (true_row, true_cols[j_col])
        M.append(e)
        print "\n", e
        print T
        T = np.delete(T, (0), axis=0)
        T = np.delete(T, (0), axis=1)
        T = np.delete(T, (j_col-1), axis=0)
        T = np.delete(T, (j_col-1), axis=1)
        print T
        del true_cols[j_col]
        del true_cols[0]
        del true_rows[j_col]
        del true_rows[0]
        
    return M



def check_matching(M, N):
    v_dict = [0]*N
    for i in xrange(0,N):
        v_dict[i] = 0

    for (v1,v2) in M:
        v_dict[v1] = 1
        v_dict[v2] = 1

    v_sum = reduce(lambda x,y: x+y, v_dict)
    if v_sum != N:
        print "MATCHING FAILED"
    else:
        print "MATCHING SUCCEDED"

def test(num_tests, N, m, err):
    for test in xrange(0, num_tests):
        T = set_weights(generate_matrix(N,0,0), N, m)
        M = np_matching(T, N, err)
        print M
        if M is not None:
            check_matching(M,N)

N = 10
test(1, N, N*4, 10**-4)



'''

def matching(T, N, err):

    if abs(determinant(deepcopy(T),N)) < err:
        return None

    if det(np.array(deepcopy(T))) < err:
        print "NUMPY - ZERO DET"

    M = []
    excl_rows = set()
    excl_cols = set()

    T_C = augment(deepcopy(T),N)

    while len(M) < N/2:
        N_C = inverse(deepcopy(T_C), N, excl_rows, excl_cols)
        T_N = get_inv_from_C(N_C, N)

        first_row = -1
        for row in xrange(0, N):
            if not_in_list(excl_rows, row):
                first_row = row
                break

        j_col = -1
        for col in xrange(0,N):
            if not_in_list(excl_cols, col) and T_N[first_row][col] != 0 \
                and T[first_row][col] != 0:
                j_col = col
                break

        e = (first_row, j_col)
        M.append(e)

        excl_rows.add(first_row)
        excl_rows.add(j_col)
        excl_cols.add(first_row)
        excl_cols.add(j_col)

    return M'''





