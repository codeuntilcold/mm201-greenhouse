import numpy as np

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
