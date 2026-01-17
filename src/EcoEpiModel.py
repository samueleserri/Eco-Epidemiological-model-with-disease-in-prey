class EcoEpiModel:
    """
    This is an abstract base class for an ecological epidemiological model
    with diseases in the prey population.

    Parameters:
    ----------
    r: intrinsic growth rate of the prey population
    K: carrying capacity of the prey population
    beta: transmission rate of the disease
    c: death rate due to the disease
    m: ratio-dependent rate
    b: predation coefficient
    k: conversion efficiency of prey into predator offspring
    d: natural death rate of the predator population
    """

    def __init__(self, param: dict):
        self.r = param['r']
        self.K = param['K']
        self.beta = param['beta']
        self.c = param['c']
        self.m = param['m']
        self.b = param['b']
        self.k = param['k']
        self.d = param['d']

    