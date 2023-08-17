
import numpy as np
from scipy.optimize import newton_krylov as JFNK


def BreakingPropagation(H1, T1, DIR1, h1, ANGbati):
    Bcoef = 0.55
    DIRrel = rel_angle_cartesian(nauticalDir2cartesianDir(DIR1), ANGbati)
    
    h2l0 = H1 / Bcoef
    
    H2 = np.zeros_like(H1)
    DIR2 = np.zeros_like(DIR1)
    h2 = np.zeros_like(H1)
    
    mask_h2l0_ge_h1 = h2l0 >= h1
    H2[mask_h2l0_ge_h1] = H1[mask_h2l0_ge_h1]
    DIR2[mask_h2l0_ge_h1] = DIR1[mask_h2l0_ge_h1]
    h2[mask_h2l0_ge_h1] = h2l0[mask_h2l0_ge_h1]
    
    mask_H1_le_0_1 = H1 <= 0.1
    H2[mask_H1_le_0_1] = H1[mask_H1_le_0_1]
    DIR2[mask_H1_le_0_1] = DIR1[mask_H1_le_0_1]
    h2[mask_H1_le_0_1] = h2l0[mask_H1_le_0_1]
    
    propProf = np.where((np.abs(DIRrel) <= 90) & (H1 > 0.1) & (h2l0 < h1))[0]
    
    if propProf.size > 0:
        # 
        # def myFun(x):
        #     return LinearShoalBreak_Residual(x, H1[propProf], T1[propProf], DIR1[propProf], h1, ANGbati, Bcoef)[0]
        myFun = lambda x: LinearShoalBreak_Residual(x, H1[propProf], T1[propProf], DIR1[propProf], h1, ANGbati, Bcoef)[0]
        
        h2l = JFNK(myFun, h2l0[propProf])
        _, H2l, DIR2l = LinearShoalBreak_Residual(h2l, H1[propProf], T1[propProf], DIR1[propProf], h1, ANGbati, Bcoef)
        
        H2[propProf] = H2l
        DIR2[propProf] = DIR2l
        h2[propProf] = h2l
    
    return H2, DIR2, h2

def abs_angle_cartesian(relD, batiD):
    """
    Absolute angle in cartesian notation, angle between [180,-180],
    0 is in EAST and positive counterclockwise.
    From a relative angle from wave and bathymetry.
    The same as rel_angle_cartesian(relD, -1 * batiD)

    Parameters:
    relD (float or array-like): Relative wave angle between wave and bathymetry, 0 is the bathymetry and positive counterclockwise.
    batiD (float or array-like): Bathymetry angle (normal to the shoreline) in Cartesian notation.

    Returns:
    waveD (float or array-like): Wave angle in Cartesian notation.
    """
    waveD = relD + batiD
    waveD[waveD > 180] = waveD[waveD > 180] - 360
    waveD[waveD < -180] = waveD[waveD < -180] + 360
    return waveD

def ADEAN(D50):
    """
    Dean parameter, D50 in meters.

    Parameters:
    D50 (float): D50 value in meters.

    Returns:
    A (float): Dean parameter.
    """
    A = 0.51 * wMOORE(D50) ** 0.44
    return A

def cartesianDir2nauticalDir(cDir):
    """
    Convert angles from Cartesian convention (0 in East, positive counterclockwise) to Nautical convention (0 in North, positive clockwise).

    Parameters:
    cDir (numpy.array or float): Angle(s) in Cartesian convention.

    Returns:
    nDir (numpy.array or float): Angle(s) in Nautical convention.
    """
    nDir = 90. - cDir
    nDir[nDir < 0] = 360 + nDir[nDir < 0]
    return nDir

def GroupCelerity(L, T, h):
    """
    Calculate the group celerity of waves.

    Parameters:
    L (float): Wave length.
    T (float): Wave period.
    h (float): Depth of wave conditions.

    Returns:
    Cg (float): Group celerity of waves.
    """
    c = L / T
    k = 2 * np.pi / L
    N = 1 + 2 * k * h / np.sinh(2 * k * h)
    Cg = c / 2 * N
    return Cg

def hunt(d, T):
    """
    Calculate wave length using the Hunt method.

    Parameters:
    d (float): Water depth.
    T (float): Wave period.

    Returns:
    L (float): Wave length.
    L0 (float): Deep water wave length.
    """
    # Constantes
    g = 9.81  # [m/s^2]

    # Calculos
    L0 = g * T ** 2 / (2 * np.pi)

    G = 2 * np.pi * (d / L0)

    p = [0.067, 0.0864, 0.4622, 0.6522, 1]

    F = G + 1 / np.polyval(p, G)

    L = T * (g * d / F) ** 0.5

    return L, L0

def LinearShoal(H1, T1, DIR1, h1, h2, ANGbati):
    """
    Wave shoaling and refraction applying linear theory with parallel, rectilinear bathymetry.

    Parameters:
    H1 (float): Initial wave height.
    T1 (float): Wave period.
    DIR1 (float): Initial wave direction in nautical convention.
    h1 (float): Initial depth of wave conditions.
    h2 (float): Final depth of wave conditions.
    ANGbati (float): Bathymetry angle (normal to the shoreline) in Cartesian convention.

    Returns:
    H2 (float): Wave height during breaking. Wave period is assumed invariant due to linear theory.
    DIR2 (float): Wave direction during breaking in nautical convention.
    """
    relDir1 = rel_angle_cartesian(nauticalDir2cartesianDir(DIR1), ANGbati)
    L1, _ = hunt(h1, T1)
    L2, _ = hunt(h2, T1)
    CG1 = GroupCelerity(L1, T1, h1)
    CG2 = GroupCelerity(L2, T1, h2)
    relDir2 = Snell_Law(L1, L2, relDir1)
    KS = np.sqrt(CG1 / CG2)
    KR = np.sqrt(np.cos(relDir1 * np.pi / 180.) / np.cos(relDir2 * np.pi / 180.))
    H2 = H1 * KS * KR
    DIR2 = cartesianDir2nauticalDir(abs_angle_cartesian(relDir2, ANGbati))

    return H2, DIR2

def LinearShoalBreak_Residual(h2l, H1, T1, DIR1, h1, ANGbati, Bcoef):
    """
    Calculate the residual of the linear shoaling and breaking equation.

    Parameters:
    h2l (float): Depth value for calculation.
    H1 (float): Initial wave height.
    T1 (float): Wave period.
    DIR1 (float): Initial wave direction in nautical convention.
    h1 (float): Initial depth of wave conditions.
    ANGbati (float): Bathymetry angle (normal to the shoreline) in Cartesian convention.
    Bcoef (float): Breaking coefficient.

    Returns:
    res (float): Residual value.
    H2l (float): Calculated wave height during breaking.
    DIR2l (float): Calculated wave direction during breaking in nautical convention.
    """
    H2l, DIR2l = LinearShoal(H1, T1, DIR1, h1, h2l, ANGbati)
    H2comp = h2l * Bcoef
    res = H2l - H2comp
    # res = res.T  # Transpose if necessary

    return res, H2l, DIR2l

def nauticalDir2cartesianDir(nDir):
    """
    Convert nautical convention angle to cartesian convention angle.

    Parameters:
    nDir (float or array-like): Nautical convention angle(s) in degrees.

    Returns:
    cDir (float or array-like): Cartesian convention angle(s) in degrees.
    """
    cDir = 90 - nDir
    cDir[cDir < -180] = 360 + cDir[cDir < -180]
    
    return cDir

def rel_angle_cartesian(waveD, batiD):
    """
    Calculate the relative angle (in degrees) between wave direction and bathymetry.
    
    Parameters:
    waveD (float or array-like): Wave angle in Cartesian notation.
    batiD (float or array-like): Bathymetry angle (normal to the shoreline) in Cartesian notation.
    
    Returns:
    relD (float or array-like): Relative wave angle between wave and bathymetry, 0 is the bathymetry and positive counterclockwise.
    """
    relD = waveD - batiD
    relD[relD > 180] = relD[relD > 180] - 360
    relD[relD < -180] = relD[relD < -180] + 360
    
    return relD

def Snell_Law(L1, L2, alpha1):
    """
    Calculate wave refraction using Snell's law.
    
    Parameters:
    L1 (float): Initial wave length.
    L2 (float): Final wave length.
    alpha1 (float): Initial wave direction in Cartesian notation.
    
    Returns:
    alpha (float): Final wave direction in Cartesian notation.
    """
    alpha = np.arcsin(L2 * np.sin(alpha1 * np.pi / 180.) / L1) * 180. / np.pi
    
    return alpha

def wast(hb, D50):
    """
    Calculate the width of the active surf zone.
    
    Parameters:
    hb (float): Depth of closure.
    D50 (float): Mean sediment grain size (meters).
    
    Returns:
    wsf (float): Width of the active surf zone.
    """
    wsf = (hb / ADEAN(D50)) ** (3 / 2)
    
    return wsf

def wMOORE(D50):
    """
    Calculate fall velocity using Moore's equation (1982).
    
    Parameters:
    D50 (numpy.ndarray or float): Mean sediment grain size (meters).
    
    Returns:
    ws (numpy.ndarray): Fall velocity values.
    """

    if D50 <= 0.1 * 1e-3:
        ws = 1.1 * 1e6 * D50 ** 2
    elif 0.1 * 1e-3 < D50 < 1 * 1e-3:
        ws = 273 * D50 ** 1.1
    elif D50 >= 1 * 1e-3:
        ws = 4.36 * D50 ** 0.5
    
    return ws
