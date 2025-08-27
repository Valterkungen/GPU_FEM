from fem_hat import fem_hat
from fem_hermite import fem_hermite
import numpy as np
import matplotlib.pyplot as plt

def gauss(x):
    rstar = 0.05
    return np.exp(-((0.5-x)/rstar)**2)
def U0(x):
    return np.sin((np.pi*x)/L)
def U1(x):
    #return np.cos((np.pi*x)/L)
    return np.zeros_like(x)
def f(x):
    return 0

def simulation(L, N, T, solver = "Hat", boundry = "Periodic", Plot = True):
    
    x_vec = np.linspace(0, L, N)
    t_vec, h = np.linspace(0, T, 1000*N, retstep=True)
    if solver == "Hermite":
        M, A, F = fem_hermite(x_vec, f, boundry)
    if solver == "Hat":
        M, A, F = fem_hat(x_vec, f, boundry)
    M_inv = np.linalg.inv(M)
    eigen = np.linalg.eigvals(M_inv@A)
    spectral = np.abs(eigen).max()
    dt = 1 / spectral

    def euler_forward(eps, epsder, h):
        eps = eps + h*epsder
        epsder = epsder + h * M_inv@(F - A@eps)
        return eps, epsder
    
    def euler_back():
        return

    if solver == "Hermite":
        x_vec = np.linspace(0, L, 2*N)
        eps = np.zeros((len(t_vec), len(x_vec)))
        eps[0, :] = gauss(x_vec)
        epsder = np.zeros((len(t_vec), len(x_vec)))
        epsder[0, :] = U1(x_vec)
        
        for tdx in range(len(t_vec)-1):
            eps[tdx + 1,:], epsder[tdx + 1,:] = euler_forward(eps[tdx, :], epsder[tdx, :], h)


        if Plot:
            fix, ax = plt.subplots()
            [line1] = ax.plot(x_vec[0::2], eps[0, :][0::2])
            plt.title("Wave Equation with Hermite functions")
            ax.set_ylim([-2,2])
            plt.draw()
            plt.pause(0.5)
            
            for tdx, solution in enumerate(eps):
                if tdx%100 == 0:
                    line1.set_ydata(solution[0::2])
                    plt.draw()
                    plt.pause(1e-6)
            plt.show()

    if solver == "Hat":
        eps = np.zeros((len(t_vec) + 1, len(x_vec)))
        eps[0, :] = gauss(x_vec)
        epsder = np.zeros((len(t_vec) + 1, len(x_vec)))
        epsder[0, :] = U1(x_vec)
        


        #for tdx in range(len(t_vec)-1):
        for tdx, _ in enumerate(t_vec):
            eps[tdx + 1,:], epsder[tdx + 1,:] = euler_forward(eps[tdx, :], epsder[tdx, :], dt)


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
    N = 100
    T = 30
    simulation(L, N, T, solver="Hat", Plot=True)




