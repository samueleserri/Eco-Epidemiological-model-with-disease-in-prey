from EcoEpiModel import EcoEpiModel
import numpy as np
import matplotlib.pyplot as plt
class DiscreteTimeModel(EcoEpiModel):

    """
    This class implements the discrete-time version of the model studied in:
    Complex dynamical behaviors in a discrete eco-epidemiological model with disease in :
    Hu, Z., Teng, Z., Jia, C. et al. Complex dynamical behaviors in a discrete eco-epidemiological model with disease in prey. Adv Differ Equ 2014, 265 (2014). 
    https://doi.org/10.1186/1687-1847-2014-265

    The model is defined by the following system of difference equations:
    S(t + 1) = S(t)*exp{r*(1-(S(t) + I(t))/K) - beta*I(t)}
    I(t + 1) = I(t)*exp{beta*S(t) - c - b*Y(t)/(m*Y(t) + I(t))}
    Y(t + 1) = Y(t)*exp{(k*b*I(t))/(m*Y(t) + I(t)) - d}
    """


    def __init__(self, param):
        super().__init__(param)
        self.R0 = self.beta * self.K / self.c
    
    def next_step(self, S, I, Y):
        S_next = self.__S_next(S, I, Y)
        I_next = self.__I_next(S, I, Y)
        Y_next = self.__Y_next(S, I, Y)
        return np.array([S_next, I_next, Y_next])
    
    def run(self, init, max_iter: int = 10**4) -> dict:
        S_values = np.zeros(max_iter)
        I_values = np.zeros(max_iter)
        Y_values = np.zeros(max_iter)

        S_values[0] = init["S0"]
        I_values[0] = init["I0"]
        Y_values[0] = init["Y0"]

        for t in range(1, max_iter):
            S_next, I_next, Y_next = self.next_step(S_values[t-1], I_values[t-1], Y_values[t-1])
            S_values[t] = S_next
            I_values[t] = I_next
            Y_values[t] = Y_next

        return {'S': S_values, 'I': I_values, 'Y': Y_values}


    def plot_phase_space(self, init, res = None, max_iter: int = 10**4):
        """
        Plot the phase space trajectory of the system in 3D.
        
        Args:
            S0: Initial susceptible prey population
            I0: Initial infected prey population
            Y0: Initial predator population
            max_iter: Number of iterations to simulate
        """
        if res is not None:
            results = res
        else:
            results = self.run(init, max_iter)
        
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        ax.plot(results['S'], results['I'], results['Y'], 'b-', linewidth=1.5, alpha=0.7)
        ax.scatter(results['S'][0], results['I'][0], results['Y'][0], c='green', s=100, label='Start', marker='o')
        ax.scatter(results['S'][-1], results['I'][-1], results['Y'][-1], c='red', s=100, label='End', marker='s')
        
        ax.set_xlabel('Susceptible Prey (S)')
        ax.set_ylabel('Infected Prey (I)')
        ax.set_zlabel('Predator (Y)')
        ax.set_title('Phase Space Trajectory')
        ax.legend()
        
        plt.tight_layout()
        plt.show()
    
    
    def plot_populations_vs_parameter(self, param_name: str, lower_bound: float, upper_bound: float, 
                                     init: dict, num_points: int = 100 , max_iter: int = 10**4):
        """
        Plot steady-state populations as a function of a parameter.
        
        Args:
            param_name: Name of the parameter to vary (e.g., 'r', 'beta', 'K', 'b', 'd', etc.)
            lower_bound: Lower bound for the parameter range
            upper_bound: Upper bound for the parameter range
            init: Dictionary with initial conditions {'S0': ..., 'I0': ..., 'Y0': ...}
            num_points: Number of points to sample in the parameter range
            max_iter: Number of iterations per simulation (last values taken as steady-state)
        """
        if not hasattr(self, param_name):
            print(f"Parameter '{param_name}' not found in model")
            return
        
        param_values = np.linspace(lower_bound, upper_bound, num_points) # *generate num_points values between lower_bound and upper_bound equally spaced
        S_steady = np.zeros(num_points)
        I_steady = np.zeros(num_points)
        Y_steady = np.zeros(num_points)
        
        original_value = getattr(self, param_name) # save original value to restore later
        
        for i, param_val in enumerate(param_values):
            # print(f"Simulating for {param_name} = {param_val:.4f} ({i+1}/{num_points})")
            setattr(self, param_name, param_val)
            results = self.run(init, max_iter)
            S_steady[i] = results['S'][-1] # last value as equilibrium
            I_steady[i] = results['I'][-1]
            Y_steady[i] = results['Y'][-1]
        
        setattr(self, param_name, original_value) # restore original parameter value
        
        fig, axes = plt.subplots(3, 1, figsize=(10, 10))
        
        axes[0].plot(param_values, S_steady, 'b-', linewidth=2, label='Susceptible Prey (S)')
        axes[0].set_ylabel('Population')
        axes[0].set_title(f'Susceptible Prey vs {param_name}')
        axes[0].grid(True, alpha=0.3)
        axes[0].legend()
        
        axes[1].plot(param_values, I_steady, 'r-', linewidth=2, label='Infected Prey (I)')
        axes[1].set_ylabel('Population')
        axes[1].set_title(f'Infected Prey vs {param_name}')
        axes[1].grid(True, alpha=0.3)
        axes[1].legend()
        
        axes[2].plot(param_values, Y_steady, 'g-', linewidth=2, label='Predator (Y)')
        axes[2].set_ylabel('Population')
        axes[2].set_xlabel(param_name)
        axes[2].set_title(f'Predator vs {param_name}')
        axes[2].grid(True, alpha=0.3)
        axes[2].legend()
        
        plt.tight_layout()
        plt.show()






    def __S_next(self, S, I, Y):
        return S * np.exp(self.r * (1 - (S + I) / self.K) - self.beta * I)
    
    def __I_next(self, S, I, Y):
        return I * np.exp(self.beta * S - self.c - (self.b * Y) / (self.m * Y + I))
    
    def __Y_next(self, S, I, Y):
        return Y * np.exp((self.k * self.b * I) / (self.m * Y + I) - self.d)

    


    

    