from pytank.aquifer import aquifer_carter_tracy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import least_squares
import cvxpy as cp
import math
from pytank.pvt_interp import interp_pvt_matbal
from pytank.pvt_correlations import Bo_bw, comp_bw_nogas

#%%

df_ta = pd.read_csv("mbal_Dataframe3.csv")
df_ta2 = df_ta[df_ta["Tank"] == "tank_center"]

df_pvt = pd.read_csv("../tests/data_for_tests/full_example_1/pvt.csv")
df_pvt = df_pvt.fillna(method="ffill")
# time_step = time
ppvt_col = "Pressure"
oil_fvf_col = "Bo"
gas_fvf_col = "Bg"
gas_oil_rs_col = "GOR"
df_ta2['gas_prod_cum']= df_ta2['gas_prod_cum']*89
def mbal(p, pi, Np, wp, bo, cf, cw, sw0, boi, N, we,bw):
    #if rsi == rs:
    Eo = (bo - boi)
    Efw = boi * (((cw * sw0) + cf) / (1 - sw0)) * (pi - p)
    F = (Np * bo) + (wp * bw)
    funcion_P = (N * (Eo + Efw)) + (we * bw) - F
    # else:
    # Eo = (bo - boi) + (rsi-rs)*bg
    # Efw = (1+m) * boi * (((cw * sw0) + cf) / (1 - sw0)) * (pi - p)
    # Eg = boi*((bg/bgi)-1)
    # Gp= (gp-(Np*rs))*bg
    # F = (Np * bo) + (wp * bw) + Gp
    # funcion_P = (N * (Eo + Efw + m*Eg)) + (we * bw) - F

    return funcion_P


def aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p, theta, k, water_visc, p_anterior, cum, pi):
    delta_t = 365
    wi = (math.pi / 5.615) * (aq_radius ** 2 - res_radius ** 2) * aq_thickness * aq_por
    f = theta / 360
    wei = ct * wi * pi * f
    rd = aq_radius / res_radius
    j = (0.00708 * k * aq_thickness * f) / (water_visc * (math.log(abs(rd))))
    pa = pi * (1 - (cum / wei))
    pr_avg = (p_anterior + p) / 2
    we = (wei / pi) * \
         (1 - np.exp((-1 * j * pi * delta_t) / wei)) * \
         (pa - pr_avg)
    cum_water_influx = cum + we
    return cum_water_influx


def press(p, Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc,
          p_anterior, cum, pi, sw0, N, boi):
    if p < 0:
        return 1e10

    bo = interp_pvt_matbal(df_pvt, ppvt_col, oil_fvf_col, p)

    bw = Bo_bw(p, t, salinity, unit=1)
    cw = comp_bw_nogas(p, t, salinity, unit=1)
    ct = cw + cf
    we = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p, theta, k, water_visc, p_anterior, cum,
                           pi)
    return mbal(p, pi, Np, wp, bo, cf, cw, sw0, boi, N, we,bw)


cf = 4.5e-6
t = 200
salinity = 30000
aq_radius = 15000
res_radius = 2000
aq_thickness = 14
aq_por = 0.25
theta = 300
k = 25
water_visc = 0.35
cum = 0
pi = 4500
sw0 = 0.25
boi = interp_pvt_matbal(df_pvt, ppvt_col, oil_fvf_col, pi)
print(boi)
N = 70e+6
x0 = 3600
P_calculada = [pi]
for i in range(len(df_ta2['Pressure'])):
    Np = df_ta2['oil_prod_cum'][i]
    wp = df_ta2['water_prod_cum'][i]
    gp = df_ta2['gas_prod_cum'][i]
    p_anterior = P_calculada[i]
    # Calculate current reservoir pressure given all other material balance variables through numeric solving.
    presion = fsolve(
        press, x0, args=(
            Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc,
            p_anterior, cum, pi, sw0, N, boi
        )
    )[0]
    print(f"Calculated Reservoir Pressure: {presion}")
    x0 = presion

    P_calculada.append(presion)
    cw = comp_bw_nogas(presion, t, salinity, unit=1)
    ct = cf + cw
    cum = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, presion, theta, k, water_visc, p_anterior,
                            cum, pi)
    #print(f"Wei:{cum}")

# def op(parametros,df,cf,t,salinity,df_pvt,sw0,cum,x0,k,water_visc,pi):
#     N,aq_radius,res_radius,aq_thickness,theta=parametros
#     P3 = pd.DataFrame({'P_calculada': pi}, index=[0])
#     for i in range(len(df['DATE'])):
#         Np = df['oil_prod_cum'][i]
#         wp = df['water_prod_cum'][i]
#         p_anterior = P3['P_calculada'][i]
#         presion = fsolve(press, x0, args=(
#         Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc, p_anterior,
#         cum, pi, sw0, N))
#         x0 = presion[0]
#         P3 = P3._append({'P_calculada': presion[0]}, ignore_index=True)
#         p_nueva = P3['P_calculada'][i + 1]
#         cw = comp_bw_nogas(p_nueva, t, salinity, unit=1)
#         ct = cf + cw
#         cum = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p_nueva, theta, k, water_visc,
#                                 p_anterior, cum, pi)
#     nueva_fila = pd.DataFrame({'Tank': 'tank_center', 'Pressure': 3700.00},
#                               index=[0])
#     df = pd.concat([nueva_fila, df]).reset_index(drop=True)
#     df.iloc[0] = df.iloc[0].fillna(0)
#     n_total=len(df['Pressure'])
#     error = (((df['Pressure']-P3['P_calculada'])**2).sum())/n_total
#     return error
#
#
# df=df_ta2
#
# bnds=[(1e+6,100e+6),(0,None),(0,None),(0,None),(20,360)]
# result = minimize(op, x0=[N,aq_radius,res_radius,aq_thickness,theta],args=(df,cf,t,salinity,df_pvt,sw0,cum,x0,k,water_visc,pi),bounds=bnds)
# Nop,aq_radiusop,res_radiusop,thicknessop, thetaop = result.x
#
# print("Valor de sTOIIP que minimiza la función:", Nop)
# print("Valor de aq radius que minimiza la función:", aq_radiusop)
# print("Valor de res radiusque minimiza la función:", res_radiusop)
# print("Valor de aq_thickness que minimiza la función:", thicknessop)
# print("Valor de theta que minimiza la función:",  thetaop)
#
# N=Nop
# aq_radius=aq_radiusop
# res_radius=res_radiusop
# theta=thetaop
# aq_thickness=thicknessop
#
# P4 = pd.DataFrame({'P_calculada': 3700.00},index=[0])
# for i in range(len(df_ta2['DATE'])):
#     Np = df_ta2['oil_prod_cum'][i]
#     wp = df_ta2['water_prod_cum'][i]
#     p_anterior = P4['P_calculada'][i]
#     presion = fsolve(press, x0, args=(
#         Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc, p_anterior,
#         cum, pi, sw0, N))
#     x0 = presion[0]
#     P4 = P4._append({'P_calculada': presion[0]}, ignore_index=True)
#     p_nueva = P4['P_calculada'][i + 1]
#     cw = comp_bw_nogas(p_nueva, t, salinity, unit=1)
#     ct = cf + cw
#     cum = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p_nueva, theta, k, water_visc,
#                             p_anterior, cum, pi)
#%%
nueva_fila = pd.DataFrame({'Date': '1987-09-01', 'Pressure': 4000.00, 'oil_fvf': 1.1},
                          index=[0])
df_ta2 = pd.concat([nueva_fila, df_ta2]).reset_index(drop=True)
df_ta2["Date"] = pd.to_datetime(df_ta2["Date"])
df_ta2.iloc[0] = df_ta2.iloc[0].fillna(0)
#%%
fig, ax = plt.subplots(figsize=(15, 10))
ax.scatter(df_ta2['Date'].dt.year, df_ta2['Pressure'], label='Presion Observada')
plt.plot(df_ta2['Date'].dt.year, P_calculada, c='g', label='Presion Calculada')
# plt.plot(df_ta2['Date'], P4, c='r', label='Presion Calculada')
plt.title('Gráfico P vs t ', fontsize=15)
plt.xlabel("Tiempo (Años)", fontsize=15)
plt.ylabel('Presion (psia)', fontsize=15)
#ax.set_ylim([400, 4000])
# plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax.grid(axis='both', color='gray', linestyle='dashed')
plt.legend(fontsize=15)
plt.show()




