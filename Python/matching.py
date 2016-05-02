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

def matching(np_T, n, err):
    E = np_T
    if abs(determinant(deepcopy(np_T),n,err)) < err:
        print "NUMPY - ZERO DET"
        return []

    excl = set()
    C_T = augment(np_T, n)
    M = []

    while len(M) < n/2:
        N = inverse(deepcopy(C_T), n, excl, excl, err)
        first_row = None
        for row in xrange(0, n):
            if not row in excl:
                first_row = row
                break

        j_col = None
        for col in xrange(n,2*n):
            if not (col-n) in excl and abs(N[first_row][col]) > err and E[first_row][col-n]!= 0:
                j_col = col - n
                break

        e = (first_row, j_col)
        M.append(e)

        excl.add(first_row)
        excl.add(j_col)

    return M

def is_blue(val):
    return val == "blue"

def is_red(val):
    return val == "red"

def color(p):
    x = random.uniform(0,1)
    if p < x:
        return "blue"
    else:
        return "red"

def is_colored(val):
    return is_blue(val) or is_red(val)

def select(V):
    if len(V) == 0:
        return "inf"
    else:
        return V[0]

def greedy_maximal_matching(A, N, p):
    pi = [0]*N
    sigma = [0]*N
    for v in xrange(0,N):
        pi[v] = "blue"

    done = False
    while not done:
        #Assign Vertex Colors
        done = True
        for v in xrange(0,N):
            if is_colored(pi[v]):
                done = False
                pi[v] = color(p)

        #Blue Vertices Propose to Red Vertices
        for v in xrange(0,N):
            if is_blue(pi[v]):
                colored_set = filter(lambda j: is_colored(pi[j]), A[v])
                if len(colored_set) == 0:
                    sigma[v] = "dead"
                else:
                    select_set = filter(lambda j: is_red(pi[j]), A[v])
                    sigma[v] = select(select_set)
            else:
                sigma[v] = "inf"

        #Red vertices respond to blue vertices
        for v in xrange(0,N):
            if is_red(pi[v]):
                colored_set = filter(lambda j: is_colored(pi[j]), A[v])
                if len(colored_set) == 0:
                    sigma[v] = "dead"
                else:
                    select_set = filter(lambda j: is_blue(pi[j]) and sigma[j] == v, A[v])
                    sigma[v] = select(select_set)

        #Match mutual proposals
        for v in xrange(0,N):
            if sigma[v] == "dead":
                pi[v] = "dead"
            elif not sigma[v] == "inf":
                if sigma[sigma[v]] == v:
                    pi[v] = min(v, sigma[v])

        return pi

def matrix_to_adjlist(A, N):
    G = []
    for i in xrange(N):
        G.append([])

    for r in xrange(N):
        for c in xrange(N):
            if A[r][c] != 0:
                G[r].append(c)

        random.shuffle(G[r])
    return G

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

def greedy_las_vegas(G,N,p):
    num_unmatched = 99
    counter = 0
    pi = []
    while num_unmatched != 0 and counter < 3000:
        pi = greedy_maximal_matching(G,N,p)
        unmatched = filter(lambda x: is_blue(x) or is_red(x), pi)
        num_unmatched = len(unmatched)
        counter += 1

    return pi


def test(num_tests, N, err, p):
    fail_counter = 0
    for test in xrange(0, num_tests):
        T = set_weights(generate_matrix(N,0,0), N)
        if abs(det(np.array(T))) < err:
            print "ZERO DET"
        else:
            G = matrix_to_adjlist(T, N)
            pi = greedy_las_vegas(G, N, 0.5)
            print pi


#test(1, 20, 10**-6, 0.53406)
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
    
    G = matrix_to_adjlist(T, N)
    #print G
    pi = greedy_las_vegas(G, N, 0.5)
    print pi

#hard_coded()

'''

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

    return M'''





