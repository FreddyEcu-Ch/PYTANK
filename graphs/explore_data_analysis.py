import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from utilities.utilities import days_in_month, interp_from_dates, interp_dates_row
from scipy.interpolate import interp1d
from material_balance.material_balance import underground_withdrawal, pressure_vol_avg

# Constants for column names
DATE_COL = "START_DATETIME"
WELL_NAME_COL = "ITEM_NAME"
TANK_NAME_COL = "Tank"
OIL_CUM_COL, WATER_CUM_COL, GAS_CUM_COL = "OIL_CUM", "WATER_CUM", "GAS_CUM"
OIL_RATE_COL, WATER_RATE_COL, GAS_RATE_COL = "oil_rate", "water_rate", "gas_rate"
LIQUID_RATE_COL, LIQUID_CUM_COL = "liquid_rate", "liquid_cum"
PRESS_COL, PRESS_TYPE_COL = "PRESSURE_DATUM", "TEST_TYPE"
LIQUID_VOL_COL = "liquid_vol"
OIL_FVF_COL, GOR_COL, GAS_FVF_COL = "Bo", "GOR", "Bg"
CAL_DAY_COL = "cal_day"
UW_COL = "UW"

# Formatter for tick labels
formatter = ticker.EngFormatter()


class ExploreDataAnalysis:
    def __init__(self, production_file, pressure_file, pvt_file):
        # Read CSV files into DataFrames
        self.df_prod = pd.read_csv(production_file)
        self.df_press = pd.read_csv(pressure_file)
        self.df_pvt = pd.read_csv(pvt_file)

        # Process the data automatically upon object creation
        self._process_data()
    def _process_data(self):
        self._cast_date_column()
        self._calculate_rates()
        self._create_dataframe()
        self._interpolate_pvt_info()
        self._interpolate_cumulatives_into_press_df()
        self._calculate_underground_withdrawal()
        self._calculate_pressure_volumetric_avg()

    def _cast_date_column(self):
        # Convert date columns to datetime format
        self.df_prod[DATE_COL] = pd.to_datetime(self.df_prod[DATE_COL])

        # Rename columns in the pressure file
        self.df_press = self.df_press.rename(columns={"WELLBORE": WELL_NAME_COL, "DATE": DATE_COL})
        self.df_press[DATE_COL] = pd.to_datetime(self.df_press[DATE_COL])

    def _create_dataframe(self):
        # Creation of dataframes
        df_tank_cols = ["START_DATETIME", "Tank", "liquid_vol"]
        self.df_tank = (
            self.df_prod[df_tank_cols]
            .groupby(df_tank_cols[:-1]).sum()
            .groupby("Tank").cumsum()
            .reset_index()
        )
        self.df_tank.rename(columns={"liquid_vol": "liquid_cum"}, inplace=True)

    def _calculate_rates(self):
        # Calculate production rates
        self.df_prod[CAL_DAY_COL] = self.df_prod[DATE_COL].map(lambda date: days_in_month(date))

        cols_input = [OIL_CUM_COL, WATER_CUM_COL, GAS_CUM_COL]
        cols_output = [OIL_RATE_COL, WATER_RATE_COL, GAS_RATE_COL]

        df_input = self.df_prod[[WELL_NAME_COL, *cols_input]]
        self.df_prod[cols_output] = (df_input.groupby(WELL_NAME_COL).diff().fillna(df_input)
                                     .div(self.df_prod[CAL_DAY_COL], axis=0))

        self.df_prod[LIQUID_RATE_COL] = self.df_prod[OIL_RATE_COL] + self.df_prod[WATER_RATE_COL]
        self.df_prod[LIQUID_CUM_COL] = self.df_prod[OIL_CUM_COL] + self.df_prod[WATER_CUM_COL]

        self.df_prod[LIQUID_VOL_COL] = self.df_prod[LIQUID_RATE_COL] * self.df_prod[CAL_DAY_COL]
        df_field = self.df_prod.groupby(DATE_COL)[LIQUID_VOL_COL].sum().cumsum().reset_index()
        df_field.rename(columns={LIQUID_VOL_COL: LIQUID_CUM_COL}, inplace=True)

        self.df_press[LIQUID_VOL_COL] = interp_from_dates(self.df_press[DATE_COL],
                                                          df_field[DATE_COL],
                                                          df_field[LIQUID_CUM_COL])


    def _interpolate_pvt_info(self):
        # Interpolate PVT information
        #self.df_pvt = self.df_pvt.drop_duplicates(subset=PRESS_COL)  # Eliminar duplicados
        oil_fvf_interp = interp1d(self.df_pvt["Pressure"], self.df_pvt[OIL_FVF_COL], fill_value="extrapolate")
        gas_oil_rs_interp = interp1d(self.df_pvt["Pressure"], self.df_pvt[GOR_COL], fill_value="extrapolate")
        gas_fvf_interp = interp1d(self.df_pvt["Pressure"], self.df_pvt[GAS_FVF_COL], fill_value="extrapolate")

        # Apply the functions to the pressure DataFrame
        self.df_press[OIL_FVF_COL] = oil_fvf_interp(self.df_press[PRESS_COL])
        self.df_press[GOR_COL] = gas_oil_rs_interp(self.df_press[PRESS_COL])
        self.df_press[GAS_FVF_COL] = gas_fvf_interp(self.df_press[PRESS_COL])

    def _interpolate_cumulatives_into_press_df(self):
        # Interpolate oil, gas, and water cumulatives into the pressure DataFrame
        for col in [OIL_CUM_COL, WATER_CUM_COL, GAS_CUM_COL]:
            self.df_press[col] = self.df_press.apply(
                lambda x: interp_dates_row(x, DATE_COL, self.df_prod, DATE_COL,
                                           col, WELL_NAME_COL, WELL_NAME_COL,
                                           left=0.0), axis=1
            )
            # For wells not available in the production DataFrame, fill NaNs with 0
            self.df_press[col].fillna(0, inplace=True)

    def _calculate_underground_withdrawal(self):
        # Calculate underground withdrawal for each well
        self.df_press[UW_COL] = underground_withdrawal(self.df_press, OIL_CUM_COL, WATER_CUM_COL,
                                                        GAS_CUM_COL, OIL_FVF_COL, 1,
                                                        GAS_FVF_COL, GOR_COL, 0)

    def _calculate_pressure_volumetric_avg(self):
        # Calculate the pressure volumetric average per tank
        avg_freq = "12MS"
        self.df_press_avg = self.df_press.groupby(TANK_NAME_COL).apply(
            lambda g: pressure_vol_avg(g, WELL_NAME_COL, DATE_COL, PRESS_COL, UW_COL,
                                       avg_freq, "end")
        ).reset_index(0)


class RatePerWell(ExploreDataAnalysis):
    def __init__(self, production_file, pressure_file, pvt_file):
        super().__init__(production_file, pressure_file, pvt_file)

        # Automatically plot production rate per well upon object creation
        self.plot_production_rate_per_well()

    def plot_production_rate_per_well(self):
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

        self.df_prod.pivot_table(OIL_RATE_COL, DATE_COL, WELL_NAME_COL).plot(colormap="Greens_r", lw=1, ax=ax1, legend=False)
        self.df_prod.pivot_table(WATER_RATE_COL, DATE_COL, WELL_NAME_COL).plot(colormap="Blues_r", lw=1, ax=ax2, legend=False)

        fig.suptitle("Production Rate per Well")
        ax1.set_ylabel("Oil Rate (STB/D)")
        ax2.set_ylabel("Water Rate (STB/D)")
        ax2.set_xlabel("Date")
        plt.show()


class RatePerTank(ExploreDataAnalysis):
    def __init__(self, production_file, pressure_file, pvt_file):
        super().__init__(production_file, pressure_file, pvt_file)

        # Automatically plot production rate per tank upon object creation
        self.plot_production_rate_per_tank()

    def plot_production_rate_per_tank(self):
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

        df_prod_tank = (self.df_prod.groupby([DATE_COL, TANK_NAME_COL])[["oil_rate", "water_rate"]]
                        .sum().reset_index())
        df_prod_tank.pivot_table("oil_rate", DATE_COL, TANK_NAME_COL).plot(lw=1, ax=ax1)
        df_prod_tank.pivot_table("water_rate", DATE_COL, TANK_NAME_COL).plot(lw=1, ax=ax2, legend=False)

        ax1.legend(fontsize=6)
        fig.suptitle("Production Rate per Tank")
        ax1.set_ylabel("Oil Rate (STB/D)")
        ax2.set_ylabel("Water Rate (STB/D)")
        ax2.set_xlabel("Date")
        plt.show()


class PressureData_Vs_LiquidCum(ExploreDataAnalysis):
    def __init__(self, production_file, pressure_file, pvt_file):
        super().__init__(production_file, pressure_file, pvt_file)

        # Automatically plot pressure vs. liquid cumulative upon object creation
        self.plot_pressure_vs_liquid_cum()

    def plot_pressure_vs_liquid_cum(self):
        fig, (ax1, ax2) = plt.subplots(2, 1)

        self.df_press.pivot_table(PRESS_COL, DATE_COL, "TEST_TYPE").plot(style="o", ax=ax1, ms=2)
        self.df_press.pivot_table(PRESS_COL, LIQUID_VOL_COL, "TEST_TYPE").plot(style="o", ax=ax2, ms=2, legend=False)

        ax1.set_title("Pressure data vs. Date")
        ax1.set_xlabel("Date")
        ax1.set_ylabel("Pressure (psia)")
        ax1.tick_params(axis="x", labelsize=8)
        ax1.legend(fontsize=8)
        ax1.yaxis.set_major_formatter(formatter)

        ax2.set_title("Pressure data vs. Liquid Cumulative")
        ax2.set_xlabel("Liquid Cumulative (STB)")
        ax2.set_ylabel("Pressure (psia)")
        ax2.xaxis.set_major_formatter(formatter)
        ax2.yaxis.set_major_formatter(formatter)

        plt.tight_layout()
        plt.show()


class LiquidCumulativesPerTank(ExploreDataAnalysis):
    def __init__(self, production_file, pressure_file, pvt_file):
        super().__init__(production_file, pressure_file, pvt_file)

        # Automatically plot liquid cumulatives per tank upon object creation
        self.plot_LiquidCumulativePerTank()

    def plot_LiquidCumulativePerTank(self):
        ax1 = (
            self.df_tank.pivot_table(LIQUID_CUM_COL, DATE_COL, TANK_NAME_COL)
            .fillna(method="ffill")
            .plot()
        )

        ax1.set_title("Liquid Cumulatives per Tank")
        ax1.set_ylabel("Liquid Cum (STB/D)")
        ax1.set_xlabel("Date")
        ax1.yaxis.set_major_formatter(formatter)
        plt.show()


class PressureVsDate_PressureVsLC(ExploreDataAnalysis):
    def __init__(self, production_file, pressure_file, pvt_file):
        super().__init__(production_file, pressure_file, pvt_file)

        # Automatically plot pressure vs. date and pressure vs. liquid cumulative upon object creation
        self.plot_PVsDate_PVsLc()

    def plot_PVsDate_PVsLc(self):
        fig_1, (ax1, ax2) = plt.subplots(2, 1)

        self.df_press.pivot_table(PRESS_COL, DATE_COL, TANK_NAME_COL).plot(ax=ax1, style="o")
        self.df_press.pivot_table(PRESS_COL, LIQUID_VOL_COL, TANK_NAME_COL).plot(ax=ax2, style="o", legend=False)

        ax1.set_title("Pressure data vs. Date")
        ax1.set_xlabel("Date")
        ax1.set_ylabel("Pressure (psia)")
        ax1.tick_params(axis="x", labelsize=8)
        ax1.legend(fontsize=8)
        ax1.yaxis.set_major_formatter(formatter)

        ax2.set_title("Pressure data vs. Liquid Cumulative")
        ax2.set_xlabel("Liquid Cumulative (STB)")
        ax2.set_ylabel("Pressure (psia)")
        ax2.xaxis.set_major_formatter(formatter)
        ax2.yaxis.set_major_formatter(formatter)

        plt.tight_layout()
        plt.show()


class AVG_Data(ExploreDataAnalysis):
    def __init__(self, production_file, pressure_file, pvt_file):
        super().__init__(production_file, pressure_file, pvt_file)

        # Automatically plot average data upon object creation
        self.plot_avg_data()

    def plot_avg_data(self):
        df_press_avg_tank = self.df_press_avg.loc[self.df_press_avg[TANK_NAME_COL] == "tank_center",
                                                  [DATE_COL, PRESS_COL]]

        df_press_tank = self.df_press.loc[self.df_press[TANK_NAME_COL] == "tank_center",
                                          [DATE_COL, PRESS_COL]]
        fig_1, ax1 = plt.subplots(1, 1)

        df_press_avg_tank.plot(x=DATE_COL, y=PRESS_COL, ax=ax1, style="bo", label="avg")
        df_press_tank.plot(x=DATE_COL, y=PRESS_COL, ax=ax1, style="ro", label="data")
        plt.show()