from interfaces import *
from fem_hat import FEMHatDiscretization

class ConjugateGradientSolver(Solver):
    def __init__(self, discretization: FEMDiscretization):
        """
        Conjugate gradient solver. 

        ## Parameters
         - ``discretization`` A FEM Discretization. 
        """
        self.discretization = discretization

        # Do other setup if requried here

    def step(self, n_steps: int, x: np.ndarray | csr_matrix, T: int, dt: float) -> np.ndarray | csr_matrix:
        """
        Step solution ``n_steps`` steps with timestep ``dt``. 

        ## Parameters
         - ``n_steps`` Number of steps. 
         - ``dt`` Timestep. 
        """
        M, A, F = self.discretization.matrices

        def conjugate_gradient_solve(A, b, x0 = None, tol = 1e-6):
            
            #Solves Ax=b
            #If inital guess oterwise zero

            if x0 is None:
                x_new = 0*b
            else:
                x_new = x0
  
            #residual r_0, r_0^Tr_0, p_0 = r_0: 
            r = b - A@x_new
            rho = np.linalg.norm(r)**2
            p = np.copy(r)
            err = 2*tol

            while err > tol:    
                x = x_new
                w = A@p
                Anorm_p_squared = np.dot(p, w)

                if Anorm_p_squared == 0:
                    break
            
                alpha = rho/Anorm_p_squared
                x_new = x + alpha*p
                r -= alpha * w
                rho_prev = rho
                rho = np.linalg.norm(r) ** 2
                p = r + (rho/rho_prev) * p
                err = np.linalg.norm(x_new - x)/np.linalg.norm(x_new)

            return x_new


       
        def rk4_step(x, v, dt, b, A, M, cg_solver):
            
            def dxdt(x, v):
                return v
            
            def dvdt(x, v):
                Ax =  A @ x
                rhs = b - Ax
                return cg_solver(M, rhs)

            k1_x = dxdt(x, v)
            k1_v = dvdt(x, v)

            k2_x = dxdt(x + 0.5*dt*k1_x, v + 0.5*dt*k1_v)
            k2_v = dvdt(x + 0.5*dt*k1_x, v + 0.5*dt*k1_v)

            k3_x = dxdt(x + 0.5*dt*k2_x, v + 0.5*dt*k2_v)
            k3_v = dvdt(x + 0.5*dt*k2_x, v + 0.5*dt*k2_v)

            k4_x = dxdt(x + dt*k3_x, v + dt*k3_v)
            k4_v = dvdt(x + dt*k3_x, v + dt*k3_v)

            x_next = x + dt*(k1_x + 2*k2_x + 2*k3_x + k4_x) / 6
            v_next = v + dt*(k1_v + 2*k2_v + 2*k3_v + k4_v) / 6
        
            return x_next, v_next

        
        def gauss(x):
            rstar = 0.05
            return np.exp(-((0.2-x)/rstar)**2)

        with Timer('Conjugate gradient solver'):
            #Set initial conditions 
            x_vec = x
            t_vec = np.linspace(0, T, n_steps)

            eps = np.zeros((len(t_vec) + 1, len(x_vec)))
            eps[0, :] = gauss(x_vec)

            epsder = np.zeros((len(t_vec) + 1, len(x_vec)))
            epsder[0, :] = np.zeros_like(x_vec)


            #Yield? 
            for tdx in range(n_steps):
                eps[tdx + 1, :], epsder[tdx + 1, :] = rk4_step(eps[tdx, :], epsder[tdx, :], dt, F, A, M, conjugate_gradient_solve)


        return eps, epsder



