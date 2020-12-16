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
    # AIR = AIR[0]        #Lay gia tri CO2-air tai thoi diem cuoi
    # TOP = TOP[0]        #Lay gia tri CO2-top tai thoi diem cuoi
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
    # AIR = AIR[0]        #Lay gia tri CO2-air tai thoi diem cuoi
    # TOP = TOP[0]        #Lay gia tri CO2-top tai thoi diem cuoi
    return np.array([AIR, TOP])
    
A = Bai3.dx()
CO2air = 1
CO2top = 1
h = 1
SoStep = 1
print('================================')
B = Euler(A, CO2air, CO2top, h, SoStep)
print('EULER:\t', end='')
print(B)   
print('================================')
B = rk4(A, CO2air, CO2top, h, SoStep)
print('RK4: \t', end='')
print(B)
