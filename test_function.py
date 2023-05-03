#####################
# Test function for NSYSU Swarm intelligence
# Editor: yi cheng yang
#####################

import numpy as np
 

def Ackley(x):
    x = np.array(x)
    temp1 = x * x
    temp2 = np.cos(2.0 * np.pi * x)
    return -20.0 * np.exp(-0.2 * np.sqrt(np.sum(temp1) / len(x))) - np.exp(np.sum(temp2) / len(x)) + 20.0 + np.exp(1)

def Griewank(x):
    x = np.array(x)
    y = np.array([i + 1 for i in range(len(x))])
    temp1 = (x * x) / 4000.0
    temp2 = np.cos(x / np.sqrt(y))
    return np.sum(temp1) - np.prod(temp2) + 1.0

def BentCigar(x):
    x = np.array(x)
    temp1 = x * x
    return temp1[0] + pow(10.0, 6) * np.sum(temp1[1:])


def Michalewicz(x):
    x = np.array(x)
    y = np.array([i + 1 for i in range(len(x))])
    temp1 = np.sin(x) * np.power(np.sin((y * x * x) / np.pi), 20)
    return -np.sum(temp1)

def Rosenbrock(x):
    x = np.array(x)
    temp1 = (x[1:] - (x[:-1] * x[:-1])) * (x[1:] - (x[:-1] * x[:-1]))
    temp2 = (x - 1) * (x - 1)
    return np.sum(100.0 * temp1 + temp2[:-1])

def Schwefel(x):
    x = np.array(x)
    temp1 = x * np.sin(np.sqrt(np.abs(x)))
    return 418.9829 * len(x) - np.sum(temp1)

def Zakharov(x):
    x = np.array(x)
    y = np.array([i + 1 for i in range(len(x))])
    temp1 = x * x
    temp2 = 0.5 * y * x
    return np.sum(temp1) + np.power(np.sum(temp2), 2) + np.power(np.sum(temp2), 4)
    

def HappyCat(x):
    x = np.array(x)
    temp1 = x * x
    sum1 = np.sum(temp1)
    sum2 = np.sum(x)
    return np.power(np.abs(sum1 - len(x)), 0.25) + (sum1 * 0.5 + sum2) / len(x) + 0.5
    

def Rastrigin(x):
    x = np.array(x)
    temp1 = x * x - 10.0 * np.cos(2.0 * np.pi * x)
    return 10.0 * len(x) + np.sum(temp1)
    

def HGBat(x):
    x = np.array(x)
    temp1 = x * x
    sum1 = np.sum(temp1)
    sum2 = np.sum(x)
    return np.power(np.abs(pow(sum1, 2) - pow(sum2, 2)), 0.5) + (sum1 * 0.5 + sum2) / len(x) + 0.5
    


function_table = {1:Ackley, 2:Griewank, 3:BentCigar, 4:Michalewicz, 5:Rosenbrock, 6:Schwefel, 7:Zakharov, 8:HappyCat, 9:Rastrigin, 10:HGBat}
search_range = {1:[-32.768, 32.768], 2:[-600.0, 600.0], 3:[-100.0, 100.0], 4:[0.0, np.pi], 5:[-10.0, 10.0], 6:[-500.0, 500.0], 7:[-10.0, 10.0], 8:[-20.0, 20.0], 9:[-5.12, 5.12], 10:[-15.0, 15.0]}



def cal_function_objective(x, func_num):
    return function_table[func_num](x)

def set_search_range(func_num):
    return search_range[func_num][0], search_range[func_num][1]
