#!/usr/bin/env python3
"""
Implementation of the Four Russians Algorithm for CSE 60445.

Author: Quinn Bardwell, William Dennen, and Tommy O'Connor.
Date: 5/12/21
"""
import math
import numpy as np
import random

def FourRussians(A, B):
    """
    Return the result of AxB.
    
    A & B are expected to be square matrixs.
    """
    assert len(A) > 0 and len(A[0]) == len(A)
    assert len(B) > 0 and len(B[0]) == len(A[0]) #Check if A & B are square
    n = len(A) #pre padded
    m = math.floor(math.log(n, 2)) #Could be natural log (?)
    padLength = m - (n % m)
    if (n%m != 0):
        for pad in range(padLength):
            B.append([0] * n)
        [row.extend([0]*padLength) for row in A]
    
    pretty_print(A)
    pretty_print(B)

    C = [[0] * n for _ in range(n)] # check maybe

    for i in range(math.ceil(n/m)): # use i, start at 0??
        RS = [[0] * n for _ in range(2**m+1)] # 2^m x n
        bp = 1
        k = 0
        for j in range(1, 2**m + 1):
            RS[j] = [x | y for x, y in zip(RS[j - 2**k], B[m*i + m - (k + 1)])]
            if bp == 1:
                bp = j + 1
                k = k + 1
            else:
                bp = bp - 1

        Ci = [[0] * n for _ in range(n)] # nxn all 0's
        for j in range(n):
            Ci[j] = RS[listToDecimal(A[j][m*i:m*i + m])]
        #Or bitwise C and Ci
        for x in range(n):
            for y in range(n):
                C[x][y] = C[x][y] | Ci[x][y]
    return C    

def FourRussians_energy(A, B):
    """
    Return the result of AxB.
    
    A & B are expected to be square matrixs.
    """
    assert len(A) > 0 and len(A[0]) == len(A)
    assert len(B) > 0 and len(B[0]) == len(A[0]) #Check if A & B are square
    n = len(A) #pre padded
    #Assuming that cache is sized correctly with input matrix size
    m = math.floor(math.log(n, 2)) 
    padLength = m - (n % m)
    #Ignoring energy this code for energy analysis, for now
    if (n%m != 0):
        for pad in range(padLength):
            B.append([0] * n)
        [row.extend([0]*padLength) for row in A]
    
    L2 = 0
    L1 = 0
    OR = 0

    C = [[0] * n for _ in range(n)] # Initialize all 0s output matrix

    for i in range(math.ceil(n/m)):
        #Set aside 2^m x n storage in L1 for RS table
        RS = [[0] * n for _ in range(2**m+1)] # 2^m x n

        #Maybe take in B part

        bp = 1
        k = 0
        for j in range(1, 2**m + 1):
            
            RS[j] = [x | y for x, y in zip(RS[j - 2**k], B[m*i + m - (k + 1)])]
            L2 += n #Read row of B from L2 cache
            L1 += n #Read a previous row of RS from L1 cache
            L1 += n #Write row for RS[j] to L1 cache
            OR += n #n bit-wise XOR, basically nothing

            if bp == 1:
                bp = j + 1
                k = k + 1
            else:
                bp = bp - 1

        #Set aside nxn matrix in L1
        Ci = [[0] * n for _ in range(n)] # nxn all 0's

        for j in range(n):
            Ci[j] = RS[listToDecimal(A[j][m*i:m*i + m])]
            L1 += 2*n #Read and write an entire row to the L1 cache for Ci[j] matrix
            L2 += m #Read in m row for A[j] in from L2
            L1 += n #Read row of RS from L1 cache

        #Or bitwise C and Ci
        for x in range(n):
            for y in range(n):
                C[x][y] = C[x][y] | Ci[x][y]
    L2 += n*n #Write a matrix to L2 final matrix value
    L1 += n*n #Read a matrix of Ci[x] from L1
    #No OR done in hardware, Ci holds final value
    return round(100*L2 + 10*L1)

def FourRussians_delay(A, B):
    """
    Return the result of AxB.
    
    A & B are expected to be square matrixs.
    """
    assert len(A) > 0 and len(A[0]) == len(A)
    assert len(B) > 0 and len(B[0]) == len(A[0]) #Check if A & B are square
    n = len(A) #pre padded
    #Assuming that cache is sized correctly with input matrix size
    m = math.floor(math.log(n, 2)) 
    padLength = m - (n % m)
    #Ignoring energy this code for energy analysis, for now
    if (n%m != 0):
        for pad in range(padLength):
            B.append([0] * n)
        [row.extend([0]*padLength) for row in A]
    
    time = 0

    C = [[0] * n for _ in range(n)] # Initialize all 0s output matrix

    for i in range(math.ceil(n/m)):
        #Set aside 2^m x n storage in L1 for RS table
        RS = [[0] * n for _ in range(2**m+1)] # 2^m x n

        #Maybe take in B part

        bp = 1
        k = 0
        for j in range(1, 2**m + 1):
            time += 1
            #some unit of time for a read/write of a row of data from L2 and L1
            RS[j] = [x | y for x, y in zip(RS[j - 2**k], B[m*i + m - (k + 1)])]

            if bp == 1:
                bp = j + 1
                k = k + 1
            else:
                bp = bp - 1

        #Set aside nxn matrix in L1
        Ci = [[0] * n for _ in range(n)] # nxn all 0's

        for j in range(n):
            time += 1
            #some unit of time for a read/write of a row of data from L2 and L1
            Ci[j] = RS[listToDecimal(A[j][m*i:m*i + m])]

        #Or bitwise C and Ci
        for x in range(n):
            for y in range(n):
                C[x][y] = C[x][y] | Ci[x][y]
    
    time += n
    #n units time for n read/writes of a row data from L1 to L2
    return time

def pretty_print(A):
    [print(" ".join([str(y) for y in x])) for x in A]
    print()

def listToDecimal(L):
    """
    Binary list input.
    Returns values when joined.
    """
    binary = "".join([str(x) for x in L])
    decimal = int(binary, 2)
    return decimal
    
def matrix_mult(A, B):
    """
    Input A, B are square matrices.
    Return A*B.
    """
    assert len(A) > 0 and len(A[0]) == len(A)
    assert len(B) > 0 and len(B[0]) == len(A[0]) #Check if A & B are square
    n = len(A) #Size of matrix
    C = [[0] * n for _ in range(n)] # Make empty output

    for i in range(n):
        for j in range(n):
            for k in range(n):
                C[i][j] += A[i][k]* B[k][j]

    return C

def matrix_mult_energy(A, B):
    """
    Input A, B are square matrices.
    Return A*B.
    """
    assert len(A) > 0 and len(A[0]) == len(A)
    assert len(B) > 0 and len(B[0]) == len(A[0]) #Check if A & B are square
    n = len(A) #Size of matrix
    C = [[0] * n for _ in range(n)] # Make empty output
    
    L2 = 0
    L1 = 0
    MACs = 0

    for i in range(n):
        L2 += n #Cachse A[i] into L1
        for j in range(n):
            #Set aside space for C[i][j] in L1 cache
            for k in range(n):
                L1 += 1 #Read A[i][k] from cache
                L1 += 2 #Read C[i][j] from cache
                L2 += 1 #Read B[k][j]
                MACs += 1 
                C[i][j] += A[i][k]* B[k][j]
            #Read C[i][j] from L1, write to L2, don't need to use again
            L1 += 1
            L2 += 1
    return round(100*L2 + 10*L1 + 0.1*MACs)

def matrix_mult_delay(A, B):
    """
    Input A, B are square matrices.
    Return A*B.
    """
    assert len(A) > 0 and len(A[0]) == len(A)
    assert len(B) > 0 and len(B[0]) == len(A[0]) #Check if A & B are square
    n = len(A) #Size of matrix
    C = [[0] * n for _ in range(n)] # Make empty output
    
    time = 0

    for i in range(n):
        time += 1
        #Read/write an entire row to L1 from L2
        for j in range(n):
            time += 1
            for k in range(n):
                C[i][j] += A[i][k]* B[k][j]
    return time

def test_algo(n, algo): #Tests for correctness/debuggers
    """
    Creates two random nxn binary matrices and uses function "algo" to get their product.
    Prints the matrices, their product, and whether it is correct or not.
    No return.
    """
    sparcity = 10

    A = np.array([[1 if random.randint(0, 100) <= sparcity else 0 for __ in range(n)] for _ in range(n) ])
    B = np.array([[1 if random.randint(0, 100) <= sparcity else 0 for __ in range(n)] for _ in range(n)])
    #A = np.mod(np.random.permutation(n*n).reshape(n, n), 2)
    #B = np.mod(np.random.permutation(n*n).reshape(n, n), 2)

    #Get rid of values > 1 in a bit of an ugly way
    C_answer = [[min(1, y) for y in x] for x in np.matmul(A, B).tolist()]
    C = algo(A.tolist(), B.tolist())
    print("Matrix A:") 
    pretty_print(A.tolist())
    print("Matrix B:")
    pretty_print(B.tolist())
    print("Matrix C:") 
    pretty_print(C)
    print("Matrix C_answer:") 
    pretty_print(C_answer)
    if (C_answer == C):
        print("Matrix C equals correct matrix.\n")
    else:
        print("Matrix C equals correct matrix.\n")

def test_metric(n, metric):
    sparcity = 10

    A = np.array([[1 if random.randint(0, 100) <= sparcity else 0 for __ in range(n)] for _ in range(n) ])
    B = np.array([[1 if random.randint(0, 100) <= sparcity else 0 for __ in range(n)] for _ in range(n)])
    #A = np.mod(np.random.permutation(n*n).reshape(n, n), 2)
    #B = np.mod(np.random.permutation(n*n).reshape(n, n), 2)

    #Get rid of values > 1 in a bit of an ugly way
    return metric(A.tolist(), B.tolist())

def mass_test_metric(n, metric):
    results = [test_metric(size, metric) for size in n]
    return results

def EDP(results1, results2):
    edp = [x*y for x, y in zip(results1, results2)]
    return edp

def main():
    #test_algo(60, FourRussians)
    #test_algo(10, matrix_mult)
    n = [5, 10, 50, 100]#, 500, 1000]

    energy1 = mass_test_metric(n, matrix_mult_energy)
    print(energy1)
    energy2 = mass_test_metric(n, FourRussians_energy)
    print(energy2)
    
    delay1 = mass_test_metric(n, matrix_mult_delay)
    print(delay1)
    delay2 = mass_test_metric(n, FourRussians_delay)
    print(delay2)

    edp1 = EDP(energy1, delay1)
    print(edp1)
    edp2 = EDP(energy2, delay2)
    print(edp2)

if __name__ == "__main__":
    main()