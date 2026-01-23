import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import cmath

def load_data():
    nk_df = pd.read_csv('nk_data.csv')
    nk_df['wavelength_nm'] = nk_df['wl'] * 1000

    sun_df = pd.read_csv('sun_data.csv')
    sun_df['wavelength_nm'] = sun_df['Wavelength, microns'] * 1000
    sun_df.rename(columns={'E-490 W/m2/micron': 'irradiance_per_micron'}, inplace=True)
    sun_df['irradiance_per_nm'] = sun_df['irradiance_per_micron'] / 1000

    return nk_df, sun_df

def spectral_properties(wavelength_nm, deg, n_data, k_data):
    theta_rad = np.radians(deg)
    #qi = np.cos(theta_rad)

    n = n_data(wavelength_nm)
    k = k_data(wavelength_nm)

    #eps = (n**2 - k**2) + 1j * (2 * n * k)

    rs_comp = (1 - (n +1j*k)) / (1 + (n +1j*k))
    R_spectral = abs(rs_comp)**2
    """
    입사각이 0이 아닌 일반적인 경우에는
    hi = cmath.sqrt(epsilon - np.sin(theta_i_rad)**2)
    rs_complex = (qi - hi) / (qi + hi)
    Rs = abs(rs_complex)**2
    rp_complex = (epsilon * qi - hi) / (epsilon * qi + hi)
    Rp = abs(rp_complex)**2
    R_spectral = 0.5 * (Rs + Rp)
    """
    A_spectral = 1 - R_spectral

    return R_spectral, A_spectral

def calculate_RA(deg, nk_df, sun_df):
    n_interpolate = interp1d(nk_df['wavelength_nm'], nk_df['n'], fill_value="extrapolate")
    k_interpolate = interp1d(nk_df['wavelength_nm'], nk_df['k'], fill_value="extrapolate")

    wavelengths_nm = sun_df['wavelength_nm'].values
    irradiances_nm = sun_df['irradiance_per_nm'].values

    R_spectral_values = []
    A_spectral_values = []
    for wl in wavelengths_nm:
        R_s, A_s = spectral_properties(wl, deg, n_interpolate, k_interpolate)
        R_spectral_values.append(R_s)
        A_spectral_values.append(A_s)
        
    R_spectral_values = np.array(R_spectral_values)
    A_spectral_values = np.array(A_spectral_values)

    numerator_R = np.trapz(R_spectral_values * irradiances_nm, wavelengths_nm)
    numerator_A = np.trapz(A_spectral_values * irradiances_nm, wavelengths_nm)
    denominator = np.trapz(irradiances_nm, wavelengths_nm)
    
    if denominator == 0:
        return 0, 0
        
    R_total = numerator_R / denominator
    A_total = numerator_A / denominator
    
    return R_total, A_total


if __name__ == "__main__":
    
    nk_data, sun_data = load_data()
    
    R0, A0 = calculate_RA(deg=0.0, nk_df=nk_data, sun_df=sun_data)
    
    print("\n--- Optical Constants ---")
    print(f"Total Reflectance (R₀): {R0:.6f}")
    print(f"Total Absorptance (A₀): {A0:.6f}")
    print(f"Check (R₀ + A₀): {R0 + A0:.6f}\n")

