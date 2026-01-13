import numpy as np

from EcoEpiModel import EcoEpiModel

class ContinuousTimeModel(EcoEpiModel):

    """
    This class implements the model studied in:
    Xiao, Y, Chen, L: A ratio-dependent predator-prey model with disease in the prey. Appl. Math. Comput. 131, 397-414 (2002)
    
    the model is defined by the following system of ordinary differential equations:
        dS/dt = r*S*(1 - (S + I)/K) - beta*S*I
        dI/dt = beta*S*I - c*I - bIY/(mY + I)
        dY/dt = -dY + k*bIY/(mY + I)
    
        where S is the susceptible prey population density, I is the infected prey population density, and Y is the predator population density.
    """

    def __init__(self, param):
        super().__init__(param)

    def dSdt(self, S, I, Y):
        return self.r * S * (1 - (S + I) / self.K) - self.beta * S * I

    def dIdt(self, S, I, Y):
        return self.beta * S * I - self.c * I - (self.b * I * Y) / (self.m * Y + I)
    

    def dYdt(self, S, I, Y):
        return -self.d * Y + (self.k * self.b * I * Y) / (self.m * Y + I)
    
    def derivatives(self, S, I, Y):
        dS = self.dSdt(S, I, Y)
        dI = self.dIdt(S, I, Y)
        dY = self.dYdt(S, I, Y)
        return np.array([dS, dI, dY])