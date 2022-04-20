import numpy as np
import numpy.random as rn
import scipy.stats as st
import scipy.fftpack as ft
import scipy.integrate as od
import scipy.signal as sig
import scipy.interpolate as intpl
from scikits.odes.odeint import odeint

####################################
####################################
####################################


def update_parameters(tmax_new, scale_new):

    global tmax
    global hbar, kB, gam, B_factor, B0, Bn, inspi, Nsam, cfac, dt, ode_tol
    global N, dwa, dwb, nv, mv, tva, tvb, wva, wvb
    global time_factor, N_extended, tmax_extended, dwa_extended, dwb_extended, nv_extended, mv_extended, tva_extended, tvb_extended, wva_extended, wvb_extended
    global scale, S02, S01, barT2, barT1, y

    # max dynamics time in wL^-1
    # to see short time dynamics: choose small tmax (2*np.pi*9)
    # to calculate steady state magnetisation: choose large tmax (2*np.pi*7200)
    tmax = 2 * np.pi * tmax_new

    hbar = 1.05E-34
    kB = 1.38E-23
    gam = -1.76E11  # electron gyromagnetic ratio

    # external magnetic field
    B_factor = 1
    B0 = 10.0 * B_factor  # magnitude
    Bn = np.array([0.0, 0.0, 1.0])  # direction = z direction

    # initial spin state
    # first 3 entries: initial unit-length spin state,
    #next 6 entries: empty x and p for integration
    inspi = [1, 0, 0, 0, 0, 0, 0, 0, 0]

    # number of samples
    Nsam = 1

    # conversion factor dimensionless temperature -> Kelvin
    cfac = B0 * np.abs(gam) * hbar / kB

    # timestep
    # choose dt =< 0.15
    dt = 0.15

    # tolerance error in the steps of odeint. If code gives
    # error "x_new out of range", reduce this tolerance.
    # NB: doing that, the time to solve the integration increases.
    ode_tol = 10**(-9)

    # number of time-axis steps (defined to be even)
    N = 2 * int(tmax / dt / 2)
    #print(f'N:\n{N}')

    # note: in the original code we had set N and tmax and calculated dt
    # above we now set dt and tmax and calculate N. This results in a 
    # rounding error that may imply a small mismatch of the time-axis of
    # the order of dt. Not to worry for now.
    # print(f'tmax-N*dt:\n{tmax-N*dt}')

    dwa = 2.0 * np.pi / tmax
    dwb = 2.0 * np.pi / (2.0 * tmax)

    # zero-centred number arrays of length N+1 and 2N respectively
    nv = np.arange(-int(N / 2), int(N / 2))
    mv = np.arange(-N, N - 1)

    # timeaxes
    tva = nv * dt
    tvb = mv * dt

    # frequency axes
    wva = nv * dwa + 0.00000000001
    wvb = mv * dwb + 0.00000000001

    # extended time for stochastic field generation
    time_factor = 10000
    N_extended = N + 2 * time_factor
    tmax_extended = N_extended * dt
    dwa_extended = 2.0 * np.pi / tmax_extended
    dwb_extended = 2.0 * np.pi / (2.0 *  tmax_extended)
    nv_extended = np.arange(- int(N / 2) - time_factor, int(N / 2) + time_factor)
    mv_extended = np.arange(- N_extended, N_extended - 1)
    tva_extended = nv_extended * dt
    tvb_extended = mv_extended * dt
    wva_extended = nv_extended * dwa_extended + 0.00000000001
    wvb_extended = mv_extended * dwb_extended + 0.00000000001

    # Spin magnitude values
    scale = scale_new  # scale between two spin magnitudes

    S02 = 1 * hbar / 2.0
    S01 = scale * S02  #S01 is the larger one!!!

    # dimensionless temperature
    barT2 = 0.0745  # use for spin S02
    barT1 = scale * barT2  # use for spin S01 #T1 is the larger one!!!

    # indicate time-point y, relevant for M(T) curve
    # y= starting point of time considered "steady state"
    y = int(N * 3 / 4)

####################################
####################################
####################################

# Lorentzian kernel
def K(w, prm):
    w0, Gamma, A = prm
    return A / (w0**2 - w**2 - 1j * w * Gamma)

# define quantum noise function (coth) to include in psd.
def quantum_noise(w, barT):
    result = np.zeros(len(w))
    if barT == 0.0:
        return np.sign(w)
    else:
        for i in range(len(w)):
            if w[i] == 0:
                result[i] = 0
            else:
                result[i] = 1.0 / np.tanh(w[i] / (2.0 * barT))
        return result

# 2020_09_16: minused np.sign(w) to obtain NO noise at T=0
def barker_noise(w, barT):
    result = np.zeros(len(w))
    if barT == 0.0:
        return 0
    else:
        for i in range(len(w)):
            if w[i] == 0:
                result[i] = 0
            else:
                result[i] = (1.0 / np.tanh(w[i] / (2.0 * barT))) - np.sign(w[i])
        return result

def classical_noise(w, barT):
    result = np.zeros(len(w))
    for i in range(len(w)):
        if w[i] == 0:
            result[i] = 0
        else:
            result[i] = 2.0 * barT / w[i] # 1.0 / (w[i] / (2.0 * barT))
    return result

# stochastic field in time
#  = with coloured noise incl quantum fluctuations
# (blue and red curves)
def b(barTv, S0, prm, Nsam, whatNoise, ns=[]):
    w0, Gamma, A = prm
    barh = hbar / S0
    noise = ns

    if whatNoise == "quantum":
        psdd = barh * np.imag(K(wvb_extended, prm)) * quantum_noise(wvb_extended, barTv)
    elif whatNoise == "barker":
        psdd = barh * np.imag(K(wvb_extended, prm)) * barker_noise(wvb_extended, barTv)
    elif whatNoise == "classical":
        psdd = barh * np.imag(K(wvb_extended, prm)) * classical_noise(wvb_extended, barTv)

    gwv = np.sqrt(psdd)  # c(w)/sqrt(2w)
    gwt = np.real((dwb_extended / (2.0 * np.pi)) * ft.fftshift(ft.fft(gwv)) * np.exp(1j * N_extended * dwb_extended * tvb_extended))
    res = np.zeros((len(noise), Nsam))
    for i in np.arange(Nsam):
        res[:, i] = sig.convolve(noise[:, i], gwt, mode='valid') * dt
    return res

# time integration of dynamics
#  = coloured noise + quantum fluctuations
# (blue and red curves)
def M_sim(barTv, S0, prm, Nsam, nsx, nsy, nsz, whatNoise, matrix):
    w0, Gamma, A = prm
    matrix_c_w_2 = matrix @ np.transpose(matrix)

    bx_int = intpl.interp1d(tva_extended,
                                b(barTv, S0, prm, Nsam, whatNoise, nsx)[:, 0],
                                bounds_error='False',
                                fill_value=0.0)
    by_int = intpl.interp1d(tva_extended,
                                b(barTv, S0, prm, Nsam, whatNoise, nsy)[:, 0],
                                bounds_error='False',
                                fill_value=0.0)
    bz_int = intpl.interp1d(tva_extended,
                                b(barTv, S0, prm, Nsam, whatNoise, nsz)[:, 0],
                                bounds_error='False',
                                fill_value=0.0)

    bn = lambda t: matrix @ np.array([bx_int(t), by_int(t), bz_int(t)])

    def system(t, V, res):
        s = V[0:3]
        p = V[3:6]
        x = V[6:9]
        Beff = np.sign(gam) * (Bn + matrix_c_w_2 @ x + bn(t))
        res[0] = s[1] * Beff[2] - s[2] * Beff[1]
        res[1] = s[2] * Beff[0] - s[0] * Beff[2]
        res[2] = s[0] * Beff[1] - s[1] * Beff[0]
        res[3:6] = -(w0**2) * x - Gamma * p + np.sign(gam) * A * s
        res[6:9] = p

    resa = odeint(system, tva, inspi, rtol=ode_tol, atol=ode_tol)
    M = np.sign(gam) * resa.values.y[:, 0:3]

    return M

# time integration of dynamics
#  = coloured noise + quantum fluctuations
# (blue and red curves)
def M_sim2(barTv, S0, prm, Nsam, nsx, nsy, nsz, whatNoise, matrix):
    w0, Gamma, A = prm
    matrix_c_w_2 = matrix @ np.transpose(matrix)
    diag_2 = np.diag(np.linalg.eig(matrix_c_w_2)[0])
    diag_2[np.abs(diag_2) < np.finfo(float).eps] = 0 # chop elements smaller than eps
    o_transf = np.linalg.eig(matrix_c_w_2)[1]

    bx_int = intpl.interp1d(tva_extended,
                                b(barTv, S0, prm, Nsam, whatNoise, nsx)[:, 0],
                                bounds_error='False',
                                fill_value=0.0)
    by_int = intpl.interp1d(tva_extended,
                                b(barTv, S0, prm, Nsam, whatNoise, nsy)[:, 0],
                                bounds_error='False',
                                fill_value=0.0)
    bz_int = intpl.interp1d(tva_extended,
                                b(barTv, S0, prm, Nsam, whatNoise, nsz)[:, 0],
                                bounds_error='False',
                                fill_value=0.0)

    bn = lambda t: o_transf @ np.sqrt(diag_2) @ np.array([bx_int(t), by_int(t), bz_int(t)])

    def system(t, V, res):
        s = V[0:3]
        p = V[3:6]
        x = V[6:9]
        Beff = np.sign(gam) * (Bn + matrix_c_w_2 @ x + bn(t))
        res[0] = s[1] * Beff[2] - s[2] * Beff[1]
        res[1] = s[2] * Beff[0] - s[0] * Beff[2]
        res[2] = s[0] * Beff[1] - s[1] * Beff[0]
        res[3:6] = -(w0**2) * x - Gamma * p + np.sign(gam) * A * s
        res[6:9] = p

    resa = odeint(system, tva, inspi, rtol=ode_tol, atol=ode_tol)
    M = np.sign(gam) * resa.values.y[:, 0:3]

    return M

# Analytic function for Average magnetization at fixed T
# computed using classical statistical mechanics
# (black curve)
def M_analytical(barTv, S0):
    t1 = np.abs(gam * S0) / np.tanh(np.abs(S0) / (hbar * barTv))
    t2 = hbar * np.abs(gam) * barTv
    result = (t1 - t2) / np.abs(gam * S0)
    return result

# the average over the entire trace from y to Nn=length of trace
# outputs a number
def avg(mag_array):
    n = len(mag_array)
    result = np.sum(mag_array[y:n]) / (n - y)
    return result

# average dynamics over n_runs realizations of the stochastic field.
# the n_runs stochastic fields are defined in the main notebook so to
# run the dynamics over the same samples of stochastic fields.
def avg_dynamics(barTv, S0, prm, Nsam, nsx, nsy, nsz, whatNoise, matrix, n_runs):
    dyn = np.zeros((N, 3, n_runs))
    dyn_sum = np.zeros((N, 3))

    for i in range(n_runs):
        dyn[:, :, i] = M_sim(barTv, S0, prm, Nsam, nsx[:, :, i], nsy[:, :, i], nsz[:, :, i], whatNoise, matrix)
        dyn_sum += dyn[:, :, i]

    return dyn_sum / n_runs

# test for the fluctuation dissipation relation
def fdr_test(index_a, index_b, barTv, S0, prm, Nsam, whatNoise, nsx, nsy, nsz):
    b_tmp_x = np.array(b(barTv, S0, prm, Nsam, whatNoise, nsx)).flatten()
    b_tmp_y = np.array(b(barTv, S0, prm, Nsam, whatNoise, nsy)).flatten()
    b_tmp_z = np.array(b(barTv, S0, prm, Nsam, whatNoise, nsz)).flatten()
    b_tmp = np.array([b_tmp_x, b_tmp_y, b_tmp_z])
    correlation = np.correlate(b_tmp[index_a], b_tmp[index_b], mode='full')
    # result = np.abs(ft.fftshift(ft.ifft(correlation)) * (len(correlation)) * dt)
    result = np.abs(ft.fftshift(ft.ifft(correlation)))
    return result

# average different data from independent stochastic realizations
def avg_data(*argv):
    n = 0
    result = np.zeros(argv[0].shape)
    for arg in argv:
        n += 1
        result = np.add(result, arg)
    return result / n