#!/bin/ python
"""
This source code aims to evolve the Oscillon profile from the analytical solution.

This code is written in Python 3
Author: Athul Muralidhar
Date: 06 June 2018
Due credits given to S.Evangelos for his guidence and  utmost patience, and helping
the author with the code
"""
# initial imports
import matplotlib.pyplot as plt
import math

if __name__ == '__main__':
    #initial Conditions
    Dx = 0.1
    Dt = 0.1*Dx
    eps = 0.1
    timesteps=int(10000/Dt)
    xsteps=int(2000/Dx)
    sumin = 0

    TS_1=[]
    TS_2=[]
    TS_3=[]

    TS_1.append(0.0) # first xtep value
    TS_2.append(0.0)

    for i in range(1,xsteps-1):
        value = 2*eps*math.sqrt(2/3)*1/math.cosh(eps*i*Dx)
        TS_1.append(value)
        TS_2.append(value)
    TS_1.append(0.0)
    TS_2.append(0.0)

    # initial energy and radius of localization
    for i in range(xsteps-1):
        sumin = sumin+Dx*(0.5*math.pow((TS_1[i]-TS_2[i])/Dt,2)+0.5*math.pow((TS_2[i+1]-TS_2[i])/Dx,2) +0.5*math.pow(TS_2[i],2)-0.25*math.pow(TS_2[i],4)+(1/(6*eps**2))*math.pow(TS_2[i],6))
    print(sumin)

    # test plotting
    step_x = [i for i in range(xsteps)]
    plt.plot(step_x,TS_1)
    plt.show()
