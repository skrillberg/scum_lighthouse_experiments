#!/bin/python

import math
import os
import random
import re
import sys

#
# Complete the 'largestMatrix' function below.
#
# The function is expected to return an INTEGER.
# The function accepts 2D_INTEGER_ARRAY arr as parameter.
#

def largestMatrix(arr):
    # Write your code here
    return largestMatrixBruteForce(arr)

def largestMatrixSearch(arr,row,col):
    #base case
    if size(arr) == 1 and aize(arr[0]) == 1:
        if arr[0][0] == 1:
            return True 
    else:
        return False
    
    #increase array search     
    largestMatrixSearch(arr)

def largestMatrixBruteForce(arr):
    largestLen= 0
    
    for i in range(0,len(arr)):
        for j in range(0, len(arr)):
            matSize = 0
            #search matrix whose upper left corner is i,j
         
            for side_len in range(1,min(len(arr)-i, len(arr)-j)):
                #boolean that keeps track of all ones
                allOnes = True
                #print(i,j,side_len)
                for row in arr[i:i+side_len]:
                    for column in row[j:j+side_len]:
                        #print(column)
                        if column != 1:
                            allOnes = False
                if allOnes:
                    matSize = side_len
            
            if matSize > largestLen:
                largestLen = matSize
                
    return largestLen        
                
                    


if __name__ == '__main__':
    fptr = open(os.environ['OUTPUT_PATH'], 'w')

    arr_rows = int(raw_input().strip())
    arr_columns = int(raw_input().strip())

    arr = []

    for _ in xrange(arr_rows):
        arr.append(map(int, raw_input().rstrip().split()))

    result = largestMatrix(arr)

    fptr.write(str(result) + '\n')

    fptr.close()
