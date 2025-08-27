import numpy as np
import matplotlib.pyplot as plt

def printMatrix(M,dim):
    for i in range(dim):
        row = ""
        for j in range(dim):
            row += f"{M[i][j]}"
            row += " "
        print(row)
        

def prod(x,depth):
    res = 1
    for _ in range(depth):
        res *= x
        x = x-1
    return res

def poly(x,deg,derr):
    if derr == 0:
        return x**deg
    else:
        coeff = prod(deg,derr)
        return coeff*x**(abs(deg-derr))

# generate hermite basis matrix
def genA(n):
    A = np.zeros((n,n))
    for i in range(0,n):
        for j in range(0,n):
            A[i][j] = poly(i%2,j,i//2)
    return A

def getPoly(index,dim,M):
    p = np.zeros(dim)
    for i in range(dim):
        p[i] = M[i][index]
    return p

def HotPlot(M,dim,N):
    xlist = np.linspace(0,1,N)
    ylist = np.zeros((dim,N))
    for i in range(dim): # what poly
        p = getPoly(i,dim,M)
        for j in range(N): # Which point x        
            ylist[i][j] = sum([p[k]*(xlist[j])**k for k in range(dim)])
    for i in range(dim):
        plt.plot(xlist,ylist[i],label=f"H{i}")
    
    plt.legend()
    plt.title("Quintic Splines")
    plt.show()

def round_matrix(matrix,n):
    rounded_matrix = []
    for row in matrix:
        rounded_row = [round(elem, n) for elem in row]
        rounded_matrix.append(rounded_row)
    return rounded_matrix
    
def main():
    dim = 4
    # A is hermite basis matrix
    A = genA(dim)
    printMatrix(A,dim)

    # invert A to get M where columns in M are hermite basis functions
    M = np.linalg.inv(A)
    M = round_matrix(M,5)
    print("\n")
    printMatrix(M,dim)
    # get polynomial from M where first argument decides which polynomial
    N = 1000
    HotPlot(M,dim,N)
    

if __name__ == '__main__':
    main()