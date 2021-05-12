#!/usr/bin/env python3
"""
Implementation of the Four Russians Algorithm for CSE 60445.

Author: Quinn Bardwell, William Dennen, and Tommy O'Connor.
Date: 5/12/21
"""
import math

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
        print("RS")
        pretty_print(RS)
        Ci = [[0] * n for _ in range(n)] # nxn all 0's
        for j in range(n):
            Ci[j] = RS[listToDecimal(A[j][m*i:m*i + m])]
        #Or bitwise C and Ci
        for x in range(n):
            for y in range(n):
                C[x][y] = C[x][y] | Ci[x][y]
    return C    

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
    

def main():
    A = [[1, 1, 0, 0, 0]
        , [0, 0, 1, 1, 1]
        , [1, 0, 0, 1, 0]
        , [1, 0, 0, 1, 1]
        , [1, 0, 1, 0, 1]]
    B = [[0, 1, 0, 0, 1]
        , [0, 0, 0, 0, 0]
        , [1, 1, 0, 0, 1]
        , [1, 0, 1, 0, 0]
        , [1, 1, 0, 1, 0]]
    
    C_answer =[[0, 1, 0, 0, 1]
            , [1, 1, 1, 1, 1]
            , [1, 1, 1, 0, 1]
            , [1, 1, 1, 1, 1]
            , [1, 1, 0, 1, 1]]

    C = FourRussians(A, B)
    print("Matrix A, B:") 
    pretty_print(A)
    pretty_print(B)
    print("Matrix C:") 
    pretty_print(C)
    print("Matrix C_answer:") 
    pretty_print(C_answer)

if __name__ == "__main__":
    main()