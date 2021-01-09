import numpy as np
import Bai2b 

# dx
def dx():
    cap_CO2_Air     = float(input('Enter cap_CO2_Air     : '))
    cap_CO2_Top     = float(input('Enter cap_CO2_Top     : '))

    eta_HeatCO2     = float(input('Enter eta_HeatCO2     : '))
    U_Blow          = float(input('Enter U_Blow          : '))
    P_Blow          = float(input('Enter P_Blow          : '))
    A_Flr           = float(input('Enter A_Flr           : '))
    #####
    BlowAir         = np.array([0, 0, Bai2b.MC_BlowAir(eta_HeatCO2, U_Blow, P_Blow, A_Flr)])
    print('\n[BlowAir]\t' + str(BlowAir) + '\n')

    U_ExtCO2        = float(input('Enter U_ExtCO2        : '))
    phi_ExtCO2      = float(input('Enter phi_ExtCO2      : '))
    #####
    ExtAir          = np.array([0, 0, Bai2b.MC_ExtAir(U_ExtCO2, phi_ExtCO2, A_Flr)])
    print('\n[ExtAir]\t' + str(ExtAir) + '\n')

    U_Pad           = float(input('Enter U_Pad           : '))
    phi_Pad         = float(input('Enter phi_Pad         : '))
    CO2_Out         = float(input('Enter CO2_Out         : '))
    #####
    PadAir          = Bai2b.MC_PadAir_Lin(U_Pad, phi_Pad, A_Flr, CO2_Out)
    print('\n[PadAir]\t' + str(PadAir) + '\n')


    U_ThScr         = float(input('Enter U_ThScr         : '))
    K_ThScr         = float(input('Enter K_ThScr         : '))
    diff_T_AirTop   = float(input('Enter diff_T_AirTop   : '))
    g               = float(input('Enter g               : '))
    rho_Mean_Air    = float(input('Enter rho_Mean_Air    : '))
    diff_rho_AirTop = float(input('Enter diff_rho_AirTop : '))
    # rho_Top         = float(input('Enter rho_Top      : '))
    f_ThermalScr    = Bai2b.f_ThScr(U_ThScr, K_ThScr, diff_T_AirTop, g, rho_Mean_Air, diff_rho_AirTop)
    #####
    AirTop          = Bai2b.MC_AirTop_Lin(f_ThermalScr)
    print('\n[AirTop]\t' + str(AirTop) + '\n')


    C_d             = float(input('Enter C_d             : '))
    U_Roof          = float(input('Enter U_Roof          : '))
    U_Side          = float(input('Enter U_Side          : '))
    A_Roof          = float(input('Enter A_Roof          : '))
    A_Side          = float(input('Enter A_Side          : '))
    h_SideRoof      = float(input('Enter h_SideRoof      : '))
    diff_T_AirOut   = float(input('Enter diff_T_AirOut   : '))
    T_Mean_Air      = float(input('Enter T_Mean_Air      : '))
    C_w             = float(input('Enter C_w             : '))
    v_Wind          = float(input('Enter v_Wind          : '))
    f_VentRoofSide_ = Bai2b.f_VentRoofSide(C_d, A_Flr, U_Roof, U_Side, A_Roof, A_Side, g, h_SideRoof, diff_T_AirOut, T_Mean_Air, C_w, v_Wind)
    f_VentRoofSide0 = Bai2b.f_VentRoofSide(C_d, A_Flr, U_Roof, U_Side, 0, A_Side, g, h_SideRoof, diff_T_AirOut, T_Mean_Air, C_w, v_Wind)

    zeta_InsScr     = float(input('Enter zeta_InsScr     : '))
    eta_InsulScr    = Bai2b.eta_InsScr(zeta_InsScr)

    c_leakage       = float(input('Enter c_leakage       : '))
    f_Leakage       = Bai2b.f_leakage(v_Wind, c_leakage)

    eta_Side        = float(input('Enter eta_Side        : '))
    eta_Side_Thr    = float(input('Enter eta_Side_Thr    : '))
    U_ThScr         = float(input('Enter U_ThScr         : '))
    f_VentSide_     = Bai2b.f_VentSide(eta_Side, eta_Side_Thr, eta_InsulScr, f_VentRoofSide0, f_VentRoofSide_, U_ThScr, f_Leakage)
    
    U_VentForced    = float(input('Enter U_VentForced    : '))
    phi_VentForced  = float(input('Enter phi_VentForced  : '))
    f_VentForced_   = Bai2b.f_VentForced(eta_InsulScr, U_VentForced, phi_VentForced, A_Flr)
    AirOut          = Bai2b.MC_AirOut_Lin(f_VentSide_, f_VentForced_, CO2_Out)
    print('\n[AirOut]\t' + str(AirOut) + '\n')


    M_CH2O          = float(input('Enter M_CH2O          : '))
    h_CBuf          = float(input('Enter h_CBuf          : '))
    Res             = float(input('Enter Res             : '))
    P               = Bai2b.Photo_Lin(Res)
    AirCan          = Bai2b.MC_AirCan_Lin(M_CH2O, h_CBuf, P)
    print('\n[AirCan]\t' + str(AirCan) + '\n')


    h_Roof          = float(input('Enter h_Roof          : '))
    eta_Roof        = float(input('Enter eta_Roof        : '))
    eta_Roof_Thr    = float(input('Enter eta_Roof_Thr    : '))
    ff_VentRoof_    = Bai2b.ff_VentRoof(C_d, U_Roof, A_Roof, A_Flr, g, h_Roof, diff_T_AirOut, T_Mean_Air, C_w, v_Wind)
    f_VentRoof_     = Bai2b.f_VentRoof(eta_Roof, eta_Roof_Thr, eta_Side, eta_InsulScr, ff_VentRoof_, U_ThScr, f_Leakage)
    #####
    TopOut          = Bai2b.MC_TopOut_Lin(f_VentRoof_, CO2_Out)
    print('[TopOut]\t' + str(TopOut) + '\n')

    result = np.array([(BlowAir + ExtAir + PadAir - AirOut - AirTop - AirCan) / cap_CO2_Air, (AirTop - TopOut) / cap_CO2_Top])    
    return result

