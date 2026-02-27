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

def cm_to_lab_cs(theta_cm,
                 M1,
                 M2,
                 rutherford_cm_cs,
                 theta_in_degrees=True):
    """
    Calculate Rutherford differential cross section
    in Lab frame from data in CM frame.

    Parameters
    ----------
    theta_cm : float
        Center-of-mass scattering angle.
        Degrees if theta_in_degrees=True, otherwise radians.
    rutherford_cm_cs : array-like
        Rutherford CS in the center-of-mass frame in mb/sr.
    Z1, Z2 : int
        Masses of projectile and target.
    theta_in_degrees : bool, optional
        If True, theta_cm is given in degrees (default: True).

    Returns
    -------
    float
        Differential [LAB] cross section in millibarn per steradian (mb/sr).
    """

    # Convert angle to radians if needed
    if theta_in_degrees:
        theta = np.deg2rad(theta_cm)
    else:
        theta = theta_cm

    # Rutherford formula in fm^2/sr
    rho = M1 / M2
    numerator = (1 + (rho)**2 + 2 * rho * np.cos(theta))**(3/2)
    denominator = abs(1 + rho * np.cos(theta))
    
    rutherford_lab_cs = rutherford_cm_cs * numerator / denominator    

    return rutherford_lab_cs

def cm_to_lab_energy(E_cm, M1, M2):
    return E_cm * (1+M1 / M2)


# Example: p particle on 6Li
# theta = 90.03      # degrees
# E_cm = 3.2        # MeV
# Z1 = 1            # proton particle
# Z2 = 3            # 6Li nucleus

# Calculation:
# Charges and Masses
Z1 = int(input("Enter charge of projectile: "))
M1 = int(input("Enter mass A_p of projectile: "))
Z2 = int(input("Enter charge of target: "))
M2 = int(input("Enter mass A of target: "))
# CM Energy
e_cm_min = float(input("Enter minimum CM energy in MeV: "))
e_cm_max = float(input("Enter maximum CM energy in MeV: "))
E_cm = np.linspace(e_cm_min, e_cm_max, 200)
E_lab = cm_to_lab_energy(E_cm, M1, M2)
# CM Angles
angles_aux = input("Enter CM angles in degrees (comma-separated): ")
thetas_in_cm = [float(angle) for angle in angles_aux.split(",")]

for theta_cm in thetas_in_cm:
    outputfile = open("rutherford_Z1_%s_Z2_%s_Theta%s.out"% (Z1,Z2,theta_cm), "w")
    for i, E in enumerate(E_cm):
        energy_lab = E_lab[i]
        sigma = rutherford_cross_section(theta_cm, E, Z1, Z2)
        sigma_lab = cm_to_lab_cs(theta_cm, M1, M2, sigma)
        outputfile.write("%9.5f %9.5f %9.5f %9.5f\n" % (E, energy_lab, sigma, sigma_lab))
    outputfile.close()