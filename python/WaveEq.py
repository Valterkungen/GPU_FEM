
import numpy as np
import matplotlib.pyplot as plt
from interfaces import *
from fem_hat import FEMHatDiscretization
from solver_cg import ConjugateGradientSolver
from solver_gaussian import GaussianSolver
from scipy.stats import norm

def gauss(x):
    rstar = 0.05
    return np.exp(-((0.5-x)/rstar)**2)
def gaussian(x, mu, sig):
    return (
        1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
    )
def U0(x):
    return np.sin((np.pi*x)/L)
def U1(x):
    return np.zeros_like(x)
def f(x):
    return 0

def simulation(L, N, T, solver = "Hat", boundry = BoundaryCondition.periodic, Plot = True):
    mean = 0
    std_dev = 1
    x_values = np.linspace(-3, 3, N)

    # Gaussian function
    x_vec = gaussian(x_values, mean, std_dev)

    #x_vec = gauss(np.linspace(0, 0.1, N))
    t_vec, h = np.linspace(0, T, 10*N, retstep=True)
    if solver == "Hat":
        FEM = FEMHatDiscretization(1, boundry, x_vec, f)
        A = FEM.A
        M = FEM.M
        F = FEM.F
    print(A)
    M_inv = np.linalg.inv(M)
    eigen = np.linalg.eigvals(M_inv@A)
    spectral = np.abs(eigen).max()
    dt = 1 / spectral

    if solver == "Hat":
        
        cg_solver = GaussianSolver(FEM)

        eps = cg_solver.step(40*N, x_vec,dt)
        
        print(eps)

        if Plot:
            fix, ax = plt.subplots()
            [line1] = ax.plot(x_vec, eps[0, :])
            plt.title("Wave Equation with Hat functions")
            ax.set_ylim([-2,2])
            plt.draw()
            plt.pause(0.5)
            
            for tdx, solution in enumerate(eps):
                if tdx%10 == 0:
                    line1.set_ydata(solution)
                    plt.draw()
                    plt.pause(1e-6)
            plt.show()
    return


if __name__ == "__main__":
    a = 0
    b = 1
    L = b-a
    N = 5
    T = 100
    simulation(L, N, T, solver="Hat", Plot=False)




