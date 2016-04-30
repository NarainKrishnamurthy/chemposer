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

def set_weights(A,N):
    counter = 0
    m = random.randint(N,N*(N-1)/2)
    while counter < m:
        i = random.randint(0, N-1)
        j = random.randint(0, N-1)

        if A[i][j] == 0:
            val = float(random.randint(1,N**2))
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
    det_T = deepcopy(np_T)
    if determinant(det_T,size_T,err) < err:
        print "NUMPY - ZERO DET"
        return []

    M = []
    true_cols = range(0,size_T)
    true_rows = range(0, size_T)

    while len(M) < size_T/2:
        N = inv(T)

        true_row = true_rows[0]
        j_col = None
        for col in xrange(0,len(N)):
            if abs(N[0][col]) > err  and E[true_row][true_cols[col]]!= 0:
                j_col = col
                break

        e = (true_row, true_cols[j_col])
        M.append(e)
        T = np.delete(T, (0), axis=0)
        T = np.delete(T, (0), axis=1)
        T = np.delete(T, (j_col-1), axis=0)
        T = np.delete(T, (j_col-1), axis=1)
        del true_cols[j_col]
        del true_cols[0]
        del true_rows[j_col]
        del true_rows[0]

    return M

def matching(np_T, size_T, err):
    E = np_T
    if determinant(deepcopy(np_T),size_T,err) < err:
        print "NUMPY - ZERO DET"
        return []

    excl = set()
    C_T = augment(np_T, size_T)
    M = []

    while len(M) < size_T/2:
        N = get_inv_from_C(inverse(deepcopy(C_T), size_T, excl, excl, err), size_T)
        first_row = None
        for row in xrange(0, size_T):
            if not row in excl:
                first_row = row
                break

        j_col = None
        for col in xrange(0,size_T):
            if not col in excl and abs(N[first_row][col]) > err and E[first_row][col]!= 0:
                j_col = col
                break

        e = (first_row, j_col)
        M.append(e)

        excl.add(first_row)
        excl.add(j_col)

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
        return 1
    else:
        print "MATCHING SUCCEDED"
        return 0

def test(num_tests, N, err):
    fail_counter = 0
    for test in xrange(0, num_tests):
        T = set_weights(generate_matrix(N,0,0), N)
        M = matching(T, N, err)
        #print M
        if M is not None:
            fail_counter += check_matching(M,N)

    print "\n", fail_counter, " our of ", num_tests, " matchings failed"

#hard_coded()
N = 40
test(100, N, 10**-6)

def hard_coded():
    N = 10
    T = [[  0.,   0.,  61.,  32.,  99.,  83.,   0.,  11.,   0.,  30.],
         [  0.,   0.,  77.,   0.,   0.,   0.,   0.,  20.,   0.,   0.],
         [-61., -77.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.],
         [-32.,   0.,   0.,   0.,   0.,  80.,  32.,   0.,  37.,  74.],
         [-99.,   0.,   0.,   0.,   0.,   0.,   0.,  29.,   0.,  26.],
         [-83.,   0.,   0., -80.,   0.,   0.,  35.,   6.,  38.,   0.],
         [  0.,   0.,   0., -32.,   0., -35.,   0.,  10.,   0.,  28.],
         [-11., -20.,   0.,   0., -29.,  -6., -10.,   0.,  96.,  50.],
         [  0.,   0.,   0., -37.,   0., -38.,   0., -96.,   0.,  63.],
         [-30.,   0.,   0., -74., -26.,   0., -28., -50., -63.,   0.]]
    M = np_matching(T, N, 10**-4)
    #print M
    if M is None:
        return
    elif len(M) != 0:
        if check_matching(M,N) is None:
            return


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





