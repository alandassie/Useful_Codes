"""
    Created on February 2026 by Alan D.K. for comparison.

    This code follows the equation XI.36 from the book
    "A. Messiah, Quantum Mechanics" to calculate the Rutherford CM
    cross-section for a particle of charge Z1 on a target of charge Z2.

    The cross-section is given in mb/sr. The energy in CM is given in MeV.
    The angles CM are given in degrees.
"""

import numpy as np
from scipy.special import gamma

# -------------------------------
# Physical constants
# -------------------------------
e2 = 1.44              # MeV·fm
hbar_c = 197.3269804   # MeV·fm
amu = 931.494          # MeV/c^2
fm2_to_mb = 10.0       # 1 fm^2 = 10 mb

# -------------------------------
# Coulomb point-like cross-section
# -------------------------------

def rutherford_cross_section(theta_cm,
                             E_cm,
                             Z1,
                             Z2,
                             theta_in_degrees=True):
    """
    Calculate Rutherford differential cross section.

    Parameters
    ----------
    theta_cm : float
        Center-of-mass scattering angle.
        Degrees if theta_in_degrees=True, otherwise radians.
    E_cm : float
        Center-of-mass energy in MeV.
    Z1, Z2 : int
        Charges of projectile and target.
    theta_in_degrees : bool, optional
        If True, theta_cm is given in degrees (default: True).

    Returns
    -------
    float
        Differential cross section in millibarn per steradian (mb/sr).
    """

    # Convert angle to radians if needed
    if theta_in_degrees:
        theta = np.deg2rad(theta_cm)
    else:
        theta = theta_cm

    # Rutherford formula in fm^2/sr
    prefactor = (Z1 * Z2 * e2 / (4.0 * E_cm))**2
    dsigma_domega_fm2 = prefactor / (np.sin(theta / 2.0)**4)

    # Convert fm^2 to mb (1 fm^2 = 10 mb)
    dsigma_domega_mb = dsigma_domega_fm2 * fm2_to_mb

    return dsigma_domega_mb


# Example: p particle on 6Li
# theta = 90.03      # degrees
# E_cm = 3.2        # MeV
# Z1 = 1            # proton particle
# Z2 = 3            # 6Li nucleus

# Calculation:
# Charges
Z1 = 2
Z2 = 2
# CM Energy
E_cm = np.linspace(0.001, 3.501, 100)
# Angles
# theta = [54.77, 63.45, 73.95, 90.03, 104.6, 116.6, 125.3, 140.8]
theta = [0.00001]

# sigma = rutherford_cross_section(theta, E_cm, Z1, Z2)

# print(f"dσ/dΩ = {sigma:.6f} mb/sr")

# elastic p on 6Li
# Define angles
# theta = [57.3, 68.3, 79, 90.75, 126.1, 159.12]

# # Define energies
# E_cm = np.linspace(0.001, 3.501, 100)

# for theta_cm in theta:
#     outputfile = open("rutherford_p_elastic_6Li_Theta%s.out"% (theta_cm), "a")
#     for E in E_cm:
#         sigma = rutherford_cross_section(theta_cm, E, 1, 3)
#         outputfile.write("%f %f\n" % (E, sigma))
#     outputfile.close()

for theta_cm in theta:
    outputfile = open("rutherford_Z1_%s_Z2_%s_Theta%s.out"% (Z1,Z2,theta_cm), "a")
    for E in E_cm:
        sigma = rutherford_cross_section(theta_cm, E, 3, 4)
        outputfile.write("%f %e\n" % (E, sigma))
    outputfile.close()