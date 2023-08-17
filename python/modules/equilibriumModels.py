import numpy as np
from scipy.signal import convolve
from waveslib import BreakingPropagation, wMOORE, wast
from numba import jit

class EBSM:
    def __init__(self, time, Hs, SS, AT, Tp, d50):
        self.time = time
        self.Hs = Hs
        self.SS = SS
        self.AT = AT
        self.Tp = Tp
        self.d50 = d50
        self.Wast = None
        self.Hb = None
        self.depthb = None
        self.Omega = None
        self.MD = None
        self.Y09 = None
        self.SF = None
    
    def LinearBreak(self, theta, depth, angleBathy):
        self.Hb, _, self.depthb = BreakingPropagation(self.Hs, self.Tp, theta, depth, angleBathy)
        self.depthb[self.Hb == 0] = 0.5 / 0.78
        self.Hb[self.Hb == 0] = 0.5
        self.Omega = self.Hb / (wMOORE(self.d50) * self.Tp)
        self.Wast = wast(self.depthb, self.d50)
    
    def MillerDean(self, dt, Yi, Hberm, flagP):
        self.MD = {'dt': dt, 'Yi': Yi, 'Hberm': Hberm, 'flagP': flagP}
        self.MD['sl'] = self.AT + self.SS
    
    def Yates09(self, dt, S0):
        self.Y09 = {'dt': dt, 'S0': S0, 'E': self.Hb ** 2}
    
    def ShoreFor(self, dt, Sini):
        self.SF = {'dt': dt, 'Sini': Sini}
        self.SF['P'] = 1 / 16 * 1025 * 9.81 * self.Hb ** 2 * (9.81 * self.Hb / 0.78) ** 0.5

def millerDean04(settings, params):
    DY0 = params[0]
    cacr = params[1]
    cero = params[2]
    
    kero = np.zeros(len(settings.Hb))
    kacr = np.zeros(len(settings.Hb))
    
    if settings.MD["flagP"] == 1:
        kero[:] = cero
        kacr[:] = cacr
    elif settings.MD["flagP"] == 2:
        kero[:] = cero * settings.Hb ** 2
        kacr[:] = cacr * settings.Hb ** 2
    elif settings.MD["flagP"] == 3:
        kero[:] = cero * settings.Hb ** 3
        kacr[:] = cacr * settings.Hb ** 3
    elif settings.MD["flagP"] == 4:
        kero[:] = cero * settings.Omega
        kacr[:] = cacr * settings.Omega
    
    Y = np.full(len(settings.Hb), np.nan)
    wl = 0.106 * settings.Hb + settings.MD["sl"]
    Wast = settings.Wast
    yeq = DY0 - Wast * wl / (settings.MD["Hberm"] + settings.depthb)
    Y[0] = settings.MD["Yi"]
    
    for i in range(1, len(settings.MD["sl"])):
        r = yeq[i] - Y[i-1] > 0
        k = kacr[i] * r + kero[i] * (not r)
        A = k * settings.MD["dt"] * 0.5
        Y[i] = (Y[i-1] + A * (yeq[i] + yeq[i-1] - Y[i-1])) / (1 + A)
   
    return Y

def yates09(settings, params):
    a = -params[0]
    b = params[1]
    cacr = -params[2]
    cero = -params[3]
    
    Seq = (settings.Y09["E"] - b) / a
    
    S = np.full(settings.Y09["E"].shape, np.nan)
    S[0] = settings.Y09["S0"]
    
    for i in range(len(S) - 1):
        r = S[i] < Seq[i + 1]
        k = cacr * r + cero * (not r)
        S[i + 1] = ((S[i] - Seq[i + 1]) * 
                    np.exp(-1. * a * k * (settings.Y09["E"][i + 1] ** 0.5) * settings.Y09["dt"])) + Seq[i + 1]
    
    return S

def shorefor(settings, params):
    phi = params[0]
    c = params[1]

    D = 2 * phi

    ii = np.arange(0, (D - 1) * 24 + 1, settings.SF["dt"])

    phivecP = 10 ** (-np.abs(ii) / (phi * 24))

    IDX = len(phivecP)

    phivecP = np.concatenate((np.zeros(IDX), phivecP))

    vent = phivecP / np.sum(phivecP)

    OmegaEQ = convolve(settings.Omega - np.mean(settings.Omega), vent, mode='same') + np.mean(settings.Omega)

    F = (settings.SF["P"] ** 0.5) * (OmegaEQ - settings.Omega) / np.std(OmegaEQ)

    F[:IDX - 1] = 0

    S = np.full(len(settings.Omega), np.nan)
    rero = F < 0
    racr = F >= 0
    S[0] = params[2]

    r = np.abs(np.sum(F[racr]) / np.sum(F[rero]))

    r_rero_F = r * rero[1:] * F[1:]
    racr_F = racr[1:] * F[1:]
    r_rero_F_prev = r * rero[:-1] * F[:-1]
    racr_F_prev = racr[:-1] * F[:-1]
    S[1:] = 0.5 * settings.SF["dt"] * c * np.cumsum(r_rero_F + racr_F + r_rero_F_prev + racr_F_prev) + S[0]
    
    return S

def Objective(model, settings, params, MetObj, ENS):
    if model == "Y09":
        # Y = yates09(settings, np.exp(params))
        Y = yates09_jit(settings.Y09["E"], settings.Y09["S0"], settings.Y09["dt"], params)
    elif model == "MD":
        # Y = millerDean04(settings, np.exp(params))
        Y = millerDean04_jit(settings.Hb, settings.MD["sl"], settings.depthb, settings.Omega, 
                             settings.MD["flagP"], settings.Wast, settings.MD["Hberm"], settings.MD["dt"],
                             settings.MD["Yi"], params)
    elif model == "SF":
        Y = shorefor(settings, np.exp(params))
        # Y = shorefor_jit(settings.Omega, settings.SF["P"], settings.SF["dt"], params)
    
    YYsl = Y[ENS["indexes"]]

    if MetObj == "Pearson":
        metVal = 1 - np.abs(np.sum((YYsl - np.mean(YYsl)) * (ENS["Yobs"] - np.mean(ENS["Yobs"]))) / 
                          (np.sqrt(np.sum((YYsl - np.mean(YYsl))**2) * np.sum((ENS["Yobs"] - np.mean(ENS["Yobs"]))**2))))
    elif MetObj == "RMSE":
        metVal = np.sqrt(np.mean((YYsl - ENS["Yobs"])**2))
    elif MetObj == "MSS":
        metVal = np.sum((YYsl - ENS["Yobs"])**2) / len(YYsl) / (np.var(YYsl) + np.var(ENS["Yobs"]) + (np.mean(YYsl) - np.mean(ENS["Yobs"]))**2)
    elif MetObj == "BSS":
        YYref = np.mean(YYsl)  # You need to define YYref here based on your requirements
        metVal = (np.mean((YYsl - ENS["Yobs"])**2) - np.mean((YYref - ENS["Yobs"])**2)) / np.mean((YYref - ENS["Yobs"])**2)
    
    return metVal


################################################
################################################
################################################
######Using numba JIT to speed up the code######
################################################
################################################
################################################

@jit(nopython = True)
def yates09_jit(E, S0, dt, params):
    a = -params[0]
    b = params[1]
    cacr = -params[2]
    cero = -params[3]
    
    Seq = (E - b) / a
    
    S = np.full(E.shape, np.nan)
    S[0] = S0
    
    for i in range(len(S) - 1):
        if S[i] < Seq[i + 1]:
            k = cacr
        else:
            k = cero

        S[i + 1] = ((S[i] - Seq[i + 1]) * 
                    np.exp(-1. * a * k * (E[i + 1] ** 0.5) * dt)) + Seq[i + 1]
    
    return S

@jit(nopython = True)
def millerDean04_jit(Hb, sl, depthb, Omega, flagP, Wast, Hberm, dt, Yi, params):
    DY0 = params[0]
    cacr = params[1]
    cero = params[2]
    
    kero = np.zeros(len(Hb))
    kacr = np.zeros(len(Hb))
    
    if flagP == 1:
        kero[:] = cero
        kacr[:] = cacr
    elif flagP == 2:
        kero[:] = cero * Hb ** 2
        kacr[:] = cacr * Hb ** 2
    elif flagP == 3:
        kero[:] = cero * Hb ** 3
        kacr[:] = cacr * Hb ** 3
    elif flagP == 4:
        kero[:] = cero * Omega
        kacr[:] = cacr * Omega
    
    Y = np.full(len(Hb), np.nan)
    wl = 0.106 * Hb + sl
    yeq = DY0 - Wast * wl / (Hberm + depthb)
    Y[0] = Yi
    
    for i in range(1, len(sl)):
        if Y[i-1] < yeq[i]:
            k = kacr[i]
        else:
            k = kero[i]
        A = k * dt * 0.5
        Y[i] = (Y[i-1] + A * (yeq[i] + yeq[i-1] - Y[i-1])) / (1 + A)
   
    return Y

@jit(nopython = True)
def shorefor_jit(Omega, P, dt, params):
    phi = params[0]
    c = params[1]

    D = 2 * phi

    ii = np.arange(0, (D - 1) * 24 + 1, dt)

    phivecP = 10 ** (-np.abs(ii) / (phi * 24))

    IDX = len(phivecP)

    phivecP = np.concatenate((np.zeros(IDX), phivecP))

    vent = phivecP / np.sum(phivecP)

    OmegaEQ = np.convolve(Omega - np.mean(Omega), vent, mode='same') + np.mean(Omega)

    F = (P ** 0.5) * (OmegaEQ - Omega) / np.std(OmegaEQ)

    F[:IDX - 1] = 0

    S = np.full(len(Omega), np.nan)
    rero = F < 0
    racr = F >= 0
    S[0] = params[2]

    r = np.abs(np.sum(F[racr]) / np.sum(F[rero]))

    r_rero_F = r * rero[1:] * F[1:]
    racr_F = racr[1:] * F[1:]
    r_rero_F_prev = r * rero[:-1] * F[:-1]
    racr_F_prev = racr[:-1] * F[:-1]
    S[1:] = 0.5 * dt * c * np.cumsum(r_rero_F + racr_F + r_rero_F_prev + racr_F_prev) + S[0]
    
    return S