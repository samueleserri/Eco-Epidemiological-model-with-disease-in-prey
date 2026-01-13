from EcoEpiModel import EcoEpiModel
import numpy as np
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
    
    def next_step(self, S, I, Y):
        S_next = self.__S_next(S, I, Y)
        I_next = self.__I_next(S, I, Y)
        Y_next = self.__Y_next(S, I, Y)
        return np.array([S_next, I_next, Y_next])
    
    def run(self, S0, I0, Y0, max_iter: int) -> dict:
        S_values = np.zeros(max_iter)
        I_values = np.zeros(max_iter)
        Y_values = np.zeros(max_iter)

        S_values[0] = S0
        I_values[0] = I0
        Y_values[0] = Y0

        for t in range(1, max_iter):
            S_next, I_next, Y_next = self.next_step(S_values[t-1], I_values[t-1], Y_values[t-1])
            S_values[t] = S_next
            I_values[t] = I_next
            Y_values[t] = Y_next

        return {'S': S_values, 'I': I_values, 'Y': Y_values}
    
    
    def __S_next(self, S, I, Y):
        return S * np.exp(self.r * (1 - (S + I) / self.K) - self.beta * I)
    
    def __I_next(self, S, I, Y):
        return I * np.exp(self.beta * S - self.c - (self.b * Y) / (self.m * Y + I))
    
    def __Y_next(self, S, I, Y):
        return Y * np.exp((self.k * self.b * I) / (self.m * Y + I) - self.d)

    


    

    