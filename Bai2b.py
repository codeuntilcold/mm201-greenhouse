import numpy as np

#####

def MC_BlowAir(eta_HeatCO2, U_Blow, P_Blow, A_Flr):
    return (eta_HeatCO2 * U_Blow * P_Blow / A_Flr)

def MC_ExtAir(U_ExtCO2, phi_ExtCO2, A_Flr):
    return (U_ExtCO2 * phi_ExtCO2 / A_Flr)

def MC_PadAir(U_Pad, phi_Pad, A_Flr, CO2_Out, CO2_Air):
    return (U_Pad * phi_Pad / A_Flr) * (CO2_Out - CO2_Air)

def MC_PadAir_Lin(U_Pad, phi_Pad, A_Flr, CO2_Out):
    return (U_Pad * phi_Pad / A_Flr) * np.array([ -1, 0, CO2_Out ])

#####

def f_ThScr(U_ThScr, K_ThScr, diff_T_AirTop, g, rho_Mean_Air, diff_rho_AirTop):
    result = U_ThScr * K_ThScr * abs(diff_T_AirTop)**(2/3) + \
             (1 - U_ThScr) * (g * (1 - U_ThScr) * abs(diff_rho_AirTop) / (2 * rho_Mean_Air))**(1/2)
    return result

def MC_AirTop(f_ThScr, CO2_Air, CO2_Top):
    return (f_ThScr * (CO2_Air - CO2_Top))

def MC_AirTop_Lin(f_ThScr):
    return f_ThScr * np.array([1, -1, 0])

#####

def f_VentRoofSide(C_d, A_Flr, U_Roof, U_Side, A_Roof, A_Side, g, h_SideRoof, diff_T_AirOut, T_Mean_Air, C_w, v_Wind):
    if (A_Roof == 0 and A_Side == 0):
        return 0
    UA_Roof = U_Roof * A_Roof
    UA_Side = U_Side * A_Side
    first_ex = UA_Roof**2 * UA_Side**2 * 2 * g * h_SideRoof * (diff_T_AirOut) / (T_Mean_Air * (UA_Roof**2 + UA_Side**2))
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

def MC_AirOut(f_VentSide, f_VentForced, CO2_Air, CO2_Out):
    return (f_VentSide + f_VentForced) * (CO2_Air - CO2_Out)

def MC_AirOut_Lin(f_VentSide, f_VentForced, CO2_Out):
    return (f_VentSide + f_VentForced) * np.array([ 1, 0, -CO2_Out ])

#####

def ff_VentRoof(C_d, U_Roof, A_Roof, A_Flr, g, h_Roof, diff_T_AirOut, T_Mean_Air, C_w, v_Wind):
    return C_d * U_Roof * A_Roof / (2 * A_Flr) * (g * h_Roof * (diff_T_AirOut) / (2 * T_Mean_Air) + C_w * v_Wind**2)**(1/2)

def f_VentRoof(eta_Roof, eta_Roof_Thr, eta_Side, eta_InsScr, ff_VentRoof, U_ThScr, f_leakage):
    if eta_Roof >= eta_Roof_Thr:
        VentRoof = ff_VentRoof
    else:
        VentRoof = U_ThScr * ff_VentRoof + (1 - U_ThScr) * f_VentRoofSide * eta_Side
    return eta_InsScr * VentRoof * 0.5 * f_leakage

def MC_TopOut(f_VentRoof, CO2_Top, CO2_Out):
    return f_VentRoof * (CO2_Top - CO2_Out)

def MC_TopOut_Lin(f_VentRoof, CO2_Out):
    return np.array([0, f_VentRoof, -f_VentRoof * CO2_Out])

#####

def Photo_Lin(Res):
    return np.array([1 / (3 * Res), 0, 0])

def MC_AirCan_Lin(M_CH2O, h_CBuf, P, R = 0.0):
    return M_CH2O * h_CBuf * (P - np.array([0, 0, R]))    # R = 0 or R = 1%P

##### Cong thuc 28, 24, 25

# e = 2.7

# def k_T(T, T_0, k_T0, H_a, R, LAI):
#     return LAI * k_T0 * e**( (-H_a/R) * (1/T - 1/T_0) ) 
    
# def f_T(T, T_0, H_d, R, S):
#     return ( 1 + e**(-H_d/R * (1/T_0 - S/H_d))) / ( 1 + e**(-H_d/R * (1/T - S/H_d))) 

# def P_Max(k_T, f_T):
#     return k_T * f_T

# ##### Cong thuc 27, 29

# def L(L0, K, LAI, m = 0.1):
#     return L0 * (1 - (K * e **(-K * LAI)) / (1 - m))

# def P_Max_MM(P_MLT, P_Max, L, L05):
#     return P_MLT * P_Max * L / (L + L05)

