import numpy as np
import Bai3

def Euler(a, air, top, h, NumOfStep):
    t = 0
    i = 1
    AIR = 0
    TOP = 0
    while i <= NumOfStep:
        t += h
        AIR = air + h * a[0].dot([ air, top, 1 ])
        TOP = top + h * a[1].dot([ air, top, 1 ])
        air = AIR
        top = TOP
        i += 1
    return np.array([AIR, TOP])

def rk4(a, air, top, h, NumOfStep):
    t = 0
    i = 1
    AIR = 0
    TOP = 0
    while i <= NumOfStep:
        t += h
        k1 = a[0].dot([ air, top, 1 ])
        k2 = a[0].dot([ air+k1*h/2, top+k1*h/2, 1 ])
        k3 = a[0].dot([ air+k2*h/2, top+k2*h/2, 1 ])
        k4 = a[0].dot([ air+k3*h,   top+k3*h,   1 ])
        AIR = air + (k1 + 2*k2 + 2*k3 + k4)*h/6
        
        k1 = a[1].dot([ air, top, 1 ])
        k2 = a[1].dot([ air+k1*h/2, top+k1*h/2, 1 ])
        k3 = a[1].dot([ air+k2*h/2, top+k2*h/2, 1 ])
        k4 = a[1].dot([ air+k3*h,   top+k3*h,   1 ])
        TOP = top + (k1 + 2*k2 + 2*k3 + k4)*h/6
        
        air = AIR
        top = TOP
        i += 1
    return np.array([AIR, TOP])
    
# A = Bai3.dx()
A = np.array([[-4.84599885e-02,  4.80036055e-02,  6.62436475e-01],
              [ 2.32594789e-01, -2.32594790e-01,  6.19319974e-08]])

# print('\n================================')

CO2air = 1078.205
CO2top = 1141.206

B = np.array([ CO2air, CO2top, 1 ])

print(A @ B)
print()
# h = 1
# SoStep = 900

# print('\n================================')
# B = Euler(A, CO2air, CO2top, h, SoStep)
# print('Approx. using Euler:\t', end='')
# print(B)   

# print('\n================================')
# B = rk4(A, CO2air, CO2top, h, SoStep)
# print('Approx. using rk4:  \t', end='')
# print(B)


# x = 599
# y = 634
# h = 1
# NumOfStep = 1

# B = Euler(A, x, y, h, NumOfStep)
# print(B)