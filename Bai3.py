import numpy as np
import Bai2b 

# dx
def dx():
    eta_HeatCO2 = eval(input('Enter eta_HeatCO2  : '))
    U_Blow      = eval(input('Enter U_Blow       : '))
    P_Blow      = eval(input('Enter P_Blow       : '))
    A_Flr       = eval(input('Enter A_Flr        : '))
    #####
    BlowAir = np.array([0, 0, Bai2b.MC_BlowAir(eta_HeatCO2, U_Blow, P_Blow, A_Flr)])
    print('\n[BlowAir]\t', str(BlowAir))
    print()

    U_ExtCO2    = eval(input('Enter U_ExtCO2     : '))
    phi_ExtCO2  = eval(input('Enter phi_ExtCO2   : '))
    #####
    ExtAir = np.array([0, 0, Bai2b.MC_ExtAir(U_ExtCO2, phi_ExtCO2, A_Flr)])
    print('\n[ExtAir]\t', str(ExtAir))
    print()

    U_Pad       = eval(input('Enter U_Pad        : '))
    phi_Pad     = eval(input('Enter phi_Pad      : '))
    CO2_Out     = eval(input('Enter CO2_Out      : '))
    #####
    PadAir = Bai2b.MC_PadAir_Lin(U_Pad, phi_Pad, A_Flr, CO2_Out)
    print('\n[PadAir]\t', end='')
    print(PadAir)
    print()

    U_ThScr         = eval(input('Enter U_ThScr      : '))
    K_ThScr         = eval(input('Enter K_ThScr      : '))
    T_Air           = eval(input('Enter T_Air        : '))
    T_Top           = eval(input('Enter T_Top        : '))
    g               = eval(input('Enter g            : '))
    rho_Mean_Air    = eval(input('Enter rho_Mean_Air : '))
    rho_Air         = eval(input('Enter rho_Air      : '))
    rho_Top         = eval(input('Enter rho_Top      : '))
    f_ThermalScr = Bai2b.f_ThScr(U_ThScr, K_ThScr, T_Air, T_Top, g, rho_Mean_Air, rho_Air, rho_Top)
    #####
    AirTop = Bai2b.MC_AirTop_Lin(f_ThermalScr)
    print('\n[AirTop]\t', end='')
    print(AirTop)
    print()


    C_d             = eval(input('Enter C_d          : '))
    U_Roof          = eval(input('Enter U_Roof       : '))
    U_Side          = eval(input('Enter U_Side       : '))
    A_Roof          = eval(input('Enter A_Roof       : '))
    A_Side          = eval(input('Enter A_Side       : '))
    h               = eval(input('Enter h            : '))
    h_SideRoof      = eval(input('Enter h_SideRoof   : '))
    T_Out           = eval(input('Enter T_Out        : '))
    T_Mean_Air      = eval(input('Enter T_Mean_Air   : '))
    C_w             = eval(input('Enter C_w          : '))
    v_Wind          = eval(input('Enter v_Wind       : '))
    f_VentRoofSide_ = Bai2b.f_VentRoofSide(C_d, A_Flr, U_Roof, U_Side, A_Roof, A_Side, g, h, h_SideRoof, T_Air, T_Out, T_Mean_Air, C_w, v_Wind)
    f_VentRoofSide0 = Bai2b.f_VentRoofSide(C_d, A_Flr, U_Roof, U_Side, 0, A_Side, g, h, h_SideRoof, T_Air, T_Out, T_Mean_Air, C_w, v_Wind)

    zeta_InsScr     = eval(input('Enter zeta_InsScr  : '))
    eta_InsulScr  = Bai2b.eta_InsScr(zeta_InsScr)

    c_leakage       = eval(input('Enter c_leakage    : '))
    f_Leakage   = Bai2b.f_leakage(v_Wind, c_leakage)

    eta_Side        = eval(input('Enter eta_Side     : '))
    eta_Side_Thr    = eval(input('Enter eta_Side_Thr : '))
    U_ThScr         = eval(input('Enter U_ThScr      : '))
    f_VentSide_  = Bai2b.f_VentSide(eta_Side, eta_Side_Thr, eta_InsulScr, f_VentRoofSide0, f_VentRoofSide_, U_ThScr, f_Leakage)
    
    U_VentForced    = eval(input('Enter U_VentForced : '))
    phi_VentForced  = eval(input('Enter phi_VentForced: '))
    f_VentForced_ = Bai2b.f_VentForced(eta_InsulScr, U_VentForced, phi_VentForced, A_Flr)
    #####
    AirOut = Bai2b.MC_AirOut_Lin(f_VentSide_, f_VentForced_, CO2_Out)
    print('\n[AirOut]\t', end='')
    print(AirOut)
    print()


    M_CH2O          = eval(input('Enter M_CH2O       : '))
    h_CBuf          = eval(input('Enter h_CBuf       : '))
    #####
    P = 0
    AirCan = np.array([0, 0, Bai2b.MC_AirCan(M_CH2O, h_CBuf, P)])
    print('\n[AirCan]\t', str(AirCan))
    print()

    h_Roof          = eval(input('Enter h_Roof       : '))
    eta_Roof        = eval(input('Enter eta_Roof     : '))
    eta_Roof_Thr    = eval(input('Enter eta_Roof_Thr : '))
    ff_VentRoof_ = Bai2b.ff_VentRoof(C_d, U_Roof, A_Roof, A_Flr, g, h_Roof, T_Air, T_Out, T_Mean_Air, C_w, v_Wind)
    f_VentRoof_ = Bai2b.f_VentRoof(eta_Roof, eta_Roof_Thr, eta_Side, eta_InsulScr, ff_VentRoof_, U_ThScr, f_Leakage)
    #####
    TopOut = Bai2b.MC_TopOut_Lin(f_VentRoof_, CO2_Out)
    print('[TopOut]\t', end='')
    print(TopOut)
    print()

    result = np.array([BlowAir + ExtAir + PadAir - AirOut - AirTop - AirCan, AirTop - TopOut])
    return result




