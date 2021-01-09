import numpy as np
from matplotlib import pyplot as plt


e = 2.7

cst = 6.4 * 10 ** (-9)

#####

def cap_VP_Air(M_Water, h_Air, R, T_Air):
    return (M_Water * h_Air) / (R * (T_Air + 273.15))

def cap_VP_Top(M_Water, h_Top, R, T_Top):
    return (M_Water * h_Top) / (R * (T_Top + 273.15))
    
def rho_Air(rho_Air0, g, M_Air, h_elevation, R):
    return rho_Air0 * e ** ((g * M_Air * h_elevation) /(293.15 * R))
    
#####
def VP_Can(T_Can):
    return 610.78 * e ** ( 17.2694 * T_Can / (T_Can + 238.3) )

def VP_Air(T_Air, RH_Air):
    return 610.78 * e ** ( 17.2694 * T_Air / (T_Air + 238.3) ) * RH_Air / 100

def VEC_CanAir(rho_Air, c_pAir, LAI, Delta_H, gamma, r_b, r_s):
    return 2 * rho_Air * c_pAir * LAI / (Delta_H * gamma * (r_b + r_s))

def MV_CanAir(VEC_CanAir, VP_Can, VP_Air):
    return VEC_CanAir * (VP_Can - VP_Air)

#####

def f_Pad(U_Pad, phi_Pad, A_Flr):
    return (U_Pad * phi_Pad / A_Flr)

def MV_PadAir(rho_Air, f_Pad, eta_Pad, chi_Pad, chi_Out):
    return f_Pad * (eta_Pad * (chi_Pad - chi_Out) + chi_Out) * rho_Air
    
def MV_AirOut_Pad(f_Pad, M_Water, R, VP_Air, T_Air):
    return f_Pad * (M_Water / R) * ( VP_Air / (T_Air + 273.15)) *rho_Air
    
def MV_AirOut_Pad_Lin(f_Pad, M_Water, R, T_Air):
    return np.array([f_Pad * (M_Water / R) * ( 1 / (T_Air + 273.15)), 0, 0])
    
#####

def MV_FogAir(U_Fog, phi_Fog, A_Flr):
    return (U_Fog * phi_Fog / A_Flr)
    
#####
   
def MV_BlowAir(eta_HeatVap, U_Blow, P_Blow, A_Flr):
    return (eta_HeatVap * U_Blow * P_Blow / A_Flr)

#####

def f_ThScr(U_ThScr, K_ThScr, T_Air, T_Top, g, rho_Mean_Air, rho_Air, rho_Top):
    result = U_ThScr * K_ThScr * abs(T_Air - T_Top)**(2/3) + (1 - U_ThScr) * (g * (1 - U_ThScr) / (2 * rho_Mean_Air) * abs(rho_Air - rho_Top))**(1/2)
    return result

def f_VentRoofSide(C_d, A_Flr, U_Roof, U_Side, A_Roof, A_Side, g, h, h_SideRoof, T_Air, T_Out, T_Mean_Air, C_w, v_Wind):
    if (A_Roof == 0 and A_Side == 0):
        return 0
    UA_Roof = U_Roof * A_Roof
    UA_Side = U_Side * A_Side
    first_ex = UA_Roof**2 * UA_Side**2 * 2 * g * h_SideRoof * (T_Air - T_Out) / (T_Mean_Air * (UA_Roof**2 + UA_Side**2))
    second_ex = ((UA_Roof + UA_Side) / 2)**2 * C_w * v_Wind**2
    result = C_d / A_Flr * (first_ex + second_ex)**(1/2)
    return result
    
def eta_InsScr(zeta_InsScr):
    return zeta_InsScr * (2 - zeta_InsScr)

def f_leakage(v_Wind, c_leakage):
    if v_Wind < 0.25: 
        v_Wind = 0.25
    return v_Wind * c_leakage
    
def f_VentSide(eta_Side, eta_Side_Thr, eta_InsScr, f_VentRoofSide_0, f_VentRoofSide, U_ThScr, f_leakage):
    VentSide = 0
    if eta_Side >= eta_Side_Thr:
        VentSide = f_VentRoofSide_0
    else:
        VentSide = U_ThScr * f_VentRoofSide_0 + (1 - U_ThScr) * f_VentRoofSide * eta_Side
    return eta_InsScr * VentSide + 0.5 * f_leakage

def f_VentForced(eta_InsScr, U_VentForced, phi_VentForced, A_Flr):
    return eta_InsScr * U_VentForced * phi_VentForced / A_Flr
    
def ff_VentRoof(C_d, U_Roof, A_Roof, A_Flr, g, h_Roof, T_Air, T_Out, T_Mean_Air, C_w, v_Wind):
    return C_d * U_Roof * A_Roof / (2 * A_Flr) * (g * h_Roof * (T_Air - T_Out) / (2 * T_Mean_Air) + C_w * v_Wind**2)**(1/2)

def f_VentRoof(eta_Roof, eta_Roof_Thr, eta_Side, eta_InsScr, ff_VentRoof, U_ThScr, f_leakage):
    if eta_Roof >= eta_Roof_Thr:
        VentRoof = ff_VentRoof
    else:
        VentRoof = U_ThScr * ff_VentRoof + (1 - U_ThScr) * f_VentRoofSide * eta_Side
    return eta_InsScr * VentRoof * 0.5 * f_leakage
    
#####

def MV_TopOut(M_Water, R, f_VentRoof, VP_Top, VP_Out, T_Top, T_Out):
    return (M_Water * f_VentRoof / R) * ((VP_Top/(T_Top + 273.15)) - (VP_Out/(T_Out + 273.15)))
    
def MV_TopOut_Lin(M_Water, R, f_VentRoof, VP_Out, T_Top, T_Out):
    return np.array([0,(M_Water * f_VentRoof / R) * (1/(T_Top + 273.15)),-(M_Water * f_VentRoof / R) * (VP_Out/(T_Out + 273.15))])  

def MV_AirTop(M_Water, R, f_ThScr, VP_Air, VP_Top, T_Air, T_Top):
    return (M_Water * f_ThScr / R) * ((VP_Air/(T_Air + 273.15)) - (VP_Top/(T_Top + 273.15)))

def MV_AirTop_Lin(M_Water, R, f_ThScr, T_Air, T_Top):
    return np.array([(M_Water * f_ThScr / R) * (1/(T_Air + 273.15)),-(M_Water * f_ThScr / R) * (1/(T_Top + 273.15)),0])  
    
def MV_AirOut(M_Water, R, f_VentSide, f_VentForced, VP_Air, VP_Out, T_Air, T_Out):
    return (M_Water * (f_VentSide + f_VentForced) / R) * ((VP_Air/(T_Air + 273.15)) - (VP_Out/(T_Out + 273.15)))
    
def MV_AirOut_Lin(M_Water, R, f_VentSide, f_VentForced, VP_Out, T_Air, T_Out):
    return np.array([(M_Water * (f_VentSide + f_VentForced) / R) * (1/(T_Air + 273.15)), 0, -(M_Water * (f_VentSide + f_VentForced) / R) * (VP_Out/(T_Out + 273.15))])
    
#####

def HEC_MechAir(U_MechCool, COP_MechCool, P_MechCool, T_Air, T_MechCool, Delta_H, VP_Air, VP_MechCool):
    upper = U_MechCool * COP_MechCool * P_MechCool
    lower = T_Air - T_MechCool + cst * Delta_H * (VP_Air - VP_MechCool)
    result = upper / lower
    return result

def MV_AirMech(VP_Air, VP_MechCool, HEC_MechAir):
    if VP_Air > VP_MechCool:
        result = cst * HEC_MechAir * (VP_Air - VP_MechCool)
    else:
        result = 0
    return result
    
#####
    
def HEC_AirThScr(U_ThScr, T_Air, T_ThScr):
    return 1.7 * U_ThScr * abs(T_Air - T_ThScr) ** (0.33)
    
def MV_AirThScr(VP_Air, VP_ThScr, HEC_AirThScr):
    if VP_Air > VP_ThScr:
        result = cst * HEC_AirThScr * (VP_Air - VP_ThScr)
    else:
        result = 0
    return result

#####

def HEC_Top_Cov_in(c_HECin, T_Top, T_Cov_in, A_Cov, A_Flr):
    return c_HECin * (T_Top - T_Cov_in) ** (0.33) * A_Cov / A_Flr

def MV_Top_Cov_in(VP_Top, VP_Cov_in, HEC_Top_Cov_in):
    if VP_Top > VP_Cov_in:
        result = cst * HEC_Top_Cov_in * (VP_Top - VP_Cov_in)
    else:
        result = 0
    return result

#####

#dx

def dx():
    #cap
    
    M_Water             = float(input('Enter M_Water             : '))
    h_Air               = float(input('Enter h_Air               : '))
    R                   = float(input('Enter R                   : '))
    T_Air               = float(input('Enter T_Air               : '))
    RH_Air              = float(input('Enter RH_Air              : '))
    h_Top               = float(input('Enter h_Top               : '))
    T_Top               = T_Air + 1
    
    cap_air = cap_VP_Air(M_Water, h_Air, R, T_Air)
    
    cap_top = cap_VP_Top(M_Water, h_Top, R, T_Top)
    
    #Pad
    
    A_Flr               = float(input('Enter A_Flr               : '))
    rho_Air0            = float(input('Enter rho_Air0            : '))
    g                   = 9.81
    M_Air               = float(input('Enter M_Air               : '))
    h_elevation         = float(input('Enter h_elevation         : '))
    U_Pad               = 0.5
    phi_Pad             = float(input('Enter phi_Pad             : '))
    eta_Pad             = 0
    chi_Pad             = float(input('Enter chi_Pad             : '))
    chi_Out             = float(input('Enter chi_Out             : '))
    
    Rho_Air = rho_Air(rho_Air0, g, M_Air, h_elevation, R)
    
    F_Pad   = f_Pad(U_Pad, phi_Pad, A_Flr)
    
    PadAir  =  np.array([0,0, MV_PadAir(Rho_Air, F_Pad, eta_Pad, chi_Pad, chi_Out)])
    
    AirOut_Pad = MV_AirOut_Pad_Lin(F_Pad, M_Water, R, T_Air)
    
    #Can
    
    VP_Air_0            = VP_Air(T_Air, RH_Air)
    VP_Can_             = VP_Can(T_Air + 1)
    r_smin              = float(input('Enter r_smin              : '))
    c_pAir              = float(input('Enter c_p,Air             : '))
    LAI                 = float(input('Enter LAI                 : '))
    Delta_H             = float(input('Enter Delta_H             : '))
    gamma               = float(input('Enter gamma               : '))
    r_b                 = float(input('Enter r_b                 : '))

    r_s                 = r_smin
    VEC                 = VEC_CanAir(Rho_Air, c_pAir, LAI, Delta_H, gamma, r_b, r_s)
    
    CanAir = np.array([-VEC, 0, VEC * VP_Can_])
    
    #FogAir
    
    U_Fog               = 0.5
    phi_Fog             = float(input('Enter phi_Fog             : '))
    
    FogAir = np.array([0, 0, MV_FogAir(U_Fog, phi_Fog, A_Flr)])
    
    #BlowAir
    
    eta_HeatVap         = float(input('Enter eta_HeatVap         : '))
    U_Blow              = 0.5
    P_Blow              = float(input('Enter P_Blow              : '))
    
    
    BlowAir = np.array([0, 0, MV_BlowAir(eta_HeatVap, U_Blow, P_Blow, A_Flr)])
    
    #Air - Top - Out
    
    #TopOut
    
    U_ThScr             = 0.5
    K_ThScr             = float(input('Enter K_ThScr             : '))
    
    rho_Mean_Air        = float(input('Enter rho_Mean_Air        : '))
    rho_Top             = float(input('Enter rho_Top             : '))
    
    C_d                 = float(input('Enter C_d                 : '))
    U_Roof              = 0.5
    U_Side              = 0.5
    A_Roof              = float(input('Enter A_Roof              : '))
    A_Side              = float(input('Enter A_Side              : '))
    h                   = float(input('Enter h                   : '))
    h_SideRoof          = float(input('Enter h_SideRoof          : '))
    T_Out               = float(input('Enter T_Out               : '))
    T_Mean_Air          = float(input('Enter T_Mean_Air          : '))
    C_w                 = float(input('Enter C_w                 : '))
    v_Wind              = float(input('Enter v_Wind              : '))
    
    zeta_InsScr         = float(input('Enter zeta_InsScr         : '))
    eta_InsulScr        = eta_InsScr(zeta_InsScr)

    c_leakage           = float(input('Enter c_leakage           : '))
    f_Leakage           = f_leakage(v_Wind, c_leakage)
    
    h_Roof              = float(input('Enter h_Roof              : '))
    eta_Roof            = float(input('Enter eta_Roof            : '))
    eta_Roof_Thr        = float(input('Enter eta_Roof_Thr        : '))
    VP_Out              = float(input('Enter VP_Out              : '))
    eta_Side            = float(input('Enter eta_Side            : '))
    eta_Side_Thr        = float(input('Enter eta_Side_Thr        : '))
    U_VentForced        = 0.5
    phi_VentForced      = float(input('Enter phi_VentForced      : '))
    
    ff_VentRoof_ = ff_VentRoof(C_d, U_Roof, A_Roof, A_Flr, g, h_Roof, T_Air, T_Out, T_Mean_Air, C_w, v_Wind)
    f_VentRoof_  = f_VentRoof(eta_Roof, eta_Roof_Thr, eta_Side, eta_InsulScr, ff_VentRoof_, U_ThScr, f_Leakage)
    
    #####
    TopOut = MV_TopOut_Lin(M_Water, R, f_VentRoof_, VP_Out, T_Top, T_Out)
    
    
    
    #AirTop
    
    f_ThermalScr = f_ThScr(U_ThScr, K_ThScr, T_Air, T_Top, g, rho_Mean_Air, Rho_Air, rho_Top)
    
    AirTop = MV_AirTop_Lin(M_Water, R, f_ThermalScr, T_Air, T_Top)
    
    #AirOut
    
    f_VentRoofSide_ = f_VentRoofSide(C_d, A_Flr, U_Roof, U_Side, A_Roof, A_Side, g, h, h_SideRoof, T_Air, T_Out, T_Mean_Air, C_w, v_Wind)
    f_VentRoofSide0 = f_VentRoofSide(C_d, A_Flr, U_Roof, U_Side, 0, A_Side, g, h, h_SideRoof, T_Air, T_Out, T_Mean_Air, C_w, v_Wind)

    
    f_VentSide_  = f_VentSide(eta_Side, eta_Side_Thr, eta_InsulScr, f_VentRoofSide0, f_VentRoofSide_, U_ThScr, f_Leakage)
    
    f_VentForced_ = f_VentForced(eta_InsulScr, U_VentForced, phi_VentForced, A_Flr)
    #####
    AirOut = MV_AirOut_Lin(M_Water, R, f_VentSide_, f_VentForced_, VP_Out, T_Air, T_Out)
    
    #Air - Mech - ThScr - TopCov,input
    
    #####
    
    U_MechCool          = float(input('Enter U_MechCool          : '))
    COP_MechCool        = float(input('Enter COP_MechCool        : '))
    P_MechCool          = float(input('Enter P_MechCool          : '))
    T_MechCool          = T_Air + 1
    VP_MechCool         = 0
    
    if VP_Air_0 > VP_MechCool:
        AirMech = np.array([0, 0, 0])
    else:
        HEC_MechAir_    = HEC_MechAir(U_MechCool, COP_MechCool, P_MechCool, T_Air, T_MechCool, Delta_H, VP_Air_0, VP_MechCool)
        # AirMech_        = MV_AirMech(VP_Air_0, VP_MechCool, HEC_MechAir_)
        AirMech         = np.array([cst * HEC_MechAir_, 0, -cst * HEC_MechAir_ * VP_MechCool]) 
    
#####

    T_ThScr             = T_Air + 1
    VP_ThScr            = 0
    
    HEC_AirThScr_   = HEC_AirThScr(U_ThScr, T_Air, T_ThScr)
    AirThScr_       = MV_AirThScr(VP_Air_0, VP_ThScr, HEC_AirThScr_)
    if AirThScr_ == 0:
        AirThScr = np.array([0, 0, 0])
    else:
        AirThScr = np.array([cst * HEC_AirThScr_, 0, -cst * HEC_AirThScr_ * VP_ThScr])

#####
    
    c_HECin             = float(input('Enter c_HECin             : '))
    T_Cov_in            = T_Air + 1
    VP_Cov_in           = 0
    A_Cov               = float(input('Enter A_Cov               : '))
    VP_Top_0            = float(input('VP_Top_0                  : '))
    
    HEC_Top_Cov_in_ = HEC_Top_Cov_in(c_HECin, T_Top, T_Cov_in, A_Cov, A_Flr)
    Top_Cov_in_     = MV_Top_Cov_in(VP_Top_0, VP_Cov_in, HEC_Top_Cov_in_)
    if Top_Cov_in_ == 0:
        Top_Cov_in = np.array([0, 0, 0])
    else:
        Top_Cov_in = np.array([0, cst * HEC_Top_Cov_in_, -cst * HEC_AirThScr_ * VP_ThScr])
    
    #####
    
    result = np.array([(CanAir + PadAir + FogAir + BlowAir - AirThScr - AirTop - AirOut - AirOut_Pad - AirMech)/cap_air , (AirTop - Top_Cov_in - TopOut)/cap_top])
    return result
    
    #####
    
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
    
A = dx()
print('================================')
print(A)   

VPair = 490
VPtop = 490
h = 1
SoStep = 1

print('================================')
B = Euler(A, VPair, VPtop, h, SoStep)
print('EULER:\t', end='')
print(B)   
print('================================')
B = rk4(A, VPair, VPtop, h, SoStep)
print('RK4: \t', end='')
print(B)
    
############ PLOTTING ###############

# require a 'static' A

# eu_result = []
# rk_result = []

# VPair_e = VPair_r = VPtop_e = VPtop_r = 490

# for nstep in range(300):
#     B = Euler(A, VPair_e, VPtop_e, h, nstep)
#     VPair_e, VPtop_e = B
#     eu_result.append(VPtop_e)
#     B = rk4(A, VPair_r, VPtop_r, h, nstep)
#     VPair_r, VPtop_r = B
#     rk_result.append(VPtop_r)

# plt.plot(eu_result, '-', label='Euler', linewidth=2)
# plt.plot(rk_result, 'r--', label='rk4', linewidth=2)
# plt.legend()
# plt.show()



