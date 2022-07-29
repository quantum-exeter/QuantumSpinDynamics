import numpy as np
import oqupy as oq

from tqdm import tqdm
from time import time

import matplotlib.pyplot as plt

import pickle

##### Constants #####

π = np.pi

ωL = 1.0

σx = oq.operators.sigma("x")
σy = oq.operators.sigma("y")
σz = oq.operators.sigma("z")

##### Spectral density #####

def LorentzianSD(α, ω0, Γ):
  ωcutoff = np.sqrt(ω0**2 + (Γ/2)**2) + 10*Γ
  def Jlor(ω):
    return (α*Γ/π)*ω/((ω**2 - ω0**2)**2 + (ω*Γ)**2)
  return Jlor, ωcutoff

def LorentzianCorrelations(α, ω0, Γ, T):
  Jlor, ωcutoff = LorentzianSD(α, ω0, Γ)
  correlations = oq.CustomSD(Jlor, cutoff=ωcutoff, cutoff_type='gaussian',
                             temperature=T, name="Lorentizna correlations")
  # endtime = 200*(ω0/Γ)
  # endtime = 120*np.sqrt(10/T)*(ω0/Γ)
  endtime = 7.5
  return correlations, endtime

##### Tempo setup #####

def compute_bath_pts(correlations, parameters, endtime, print_debug=False):
  bathx = oq.Bath(σx, correlations[0])
  bathy = oq.Bath(σy, correlations[1])
  bathz = oq.Bath(σz, correlations[2])
  if print_debug:
    fig, ax = plt.subplots(1,1)
    oq.helpers.plot_correlations_with_parameters(correlations[0], parameters, ax=ax)
    plt.savefig("debug_correls.pdf")
  ptx = oq.pt_tempo_compute(bathx,
                            start_time=0.0, end_time=endtime,
                            parameters=parameters)
  pty = oq.pt_tempo_compute(bathy,
                            start_time=0.0, end_time=endtime,
                            parameters=parameters)
  ptz = oq.pt_tempo_compute(bathz,
                            start_time=0.0, end_time=endtime,
                            parameters=parameters)
  return ptx, pty, ptz

def compute_dynamics(pts, initial_state=oq.operators.spin_dm("z-")):
  dynamics = oq.compute_dynamics(system=oq.System(σz), process_tensor=pts,
                                 initial_state=initial_state)
  _, sx = dynamics.expectations(σx, real=True)
  _, sy = dynamics.expectations(σy, real=True)
  t, sz = dynamics.expectations(σz, real=True)
  return t, sx, sy, sz

##### Tempo sweeps #####

def tempo_Tsweep_scratch(α, ω0, Γ, T, dt, dkmax, epsrel, ID):
  parameters = oq.TempoParameters(dt=dt, dkmax=dkmax, epsrel=epsrel)

  sss = np.zeros((len(T), 3))
  for n in range(len(T)):
    correlationsx, endtimex = LorentzianCorrelations(α[0], ω0[0], Γ[0], T[n])
    correlationsy, endtimey = LorentzianCorrelations(α[1], ω0[1], Γ[1], T[n])
    correlationsz, endtimez = LorentzianCorrelations(α[2], ω0[2], Γ[2], T[n])
    ptx, pty, ptz = compute_bath_pts([correlationsx, correlationsy, correlationsz],
                                      parameters, np.max([endtimex, endtimey, endtimez]),
                                      print_debug=True)
    t, sx, sy, sz = compute_dynamics([ptx, pty, ptz])
    sss[n,0] = np.mean(sx[-10:])
    sss[n,1] = np.mean(sy[-10:])
    sss[n,2] = np.mean(sz[-10:])

    np.savez("outputs/run_3d_T{}_{}.npz".format(T[n], ID),
             α=α, ω0=ω0, Γ=Γ, dt=dt, dkmax=dkmax, epsrel=epsrel,
             t=t, sx=sx, sy=sy, sz=sz, T=T[n])

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(12,5))
    ax1.plot(t, sx)
    ax2.plot(t, sy)
    ax3.plot(t, sz)
    fig.tight_layout()
    plt.savefig("outputs/run_3d_T{}_{}.pdf".format(T[n], ID))
    plt.close()

    print("done: ", n+1, "/", len(T))

  np.savez("outputs/ss_run_3d_{}.npz".format(ID),
           α=α, ω0=ω0, Γ=Γ, dt=dt, dkmax=dkmax, epsrel=epsrel,
           T=T, sss=sss)

  fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(12,5))
  ax1.semilogx(T, sss[:,0], '--o')
  ax2.semilogx(T, sss[:,1], '--o')
  ax3.semilogx(T, sss[:,2], '--o')
  fig.tight_layout()
  plt.savefig("outputs/ss_run_3d_{}.pdf".format(ID))

def tempo_Tsweep_precomputepts(α, ω0, Γ, T, dt, dkmax, epsrel, ID):
  parameters = oq.TempoParameters(dt=dt, dkmax=dkmax, epsrel=epsrel)

  for n in range(len(T)):
    correlationsx, endtimex = LorentzianCorrelations(α[0], ω0[0], Γ[0], T[n])
    correlationsy, endtimey = LorentzianCorrelations(α[1], ω0[1], Γ[1], T[n])
    correlationsz, endtimez = LorentzianCorrelations(α[2], ω0[2], Γ[2], T[n])
    ptx, pty, ptz = compute_bath_pts([correlationsx, correlationsy, correlationsz],
                                     parameters, np.max([endtimex, endtimey, endtimez]),
                                     print_debug=True)

    pickle.dump({"ptx": ptx, "pty": pty, "ptz": ptz},
                open("outputs/pts_3d_T{}_{}.p".format(T[n], ID), "wb"))

    print("done: ", n+1, "/", len(T))

def tempo_Tsweep_dynfrompts(T, ID):
  sss = np.zeros((len(T), 3))
  for n in range(len(T)):
    pts_struct = pickle.load(open("outputs/pts_3d_T{}_{}.p".format(T[n], ID), "rb"))
    t, sx, sy, sz = compute_dynamics([pts_struct["ptx"], pts_struct["pty"], pts_struct["ptz"]])
    sss[n,0] = np.mean(sx[-10:])
    sss[n,1] = np.mean(sy[-10:])
    sss[n,2] = np.mean(sz[-10:])

    np.savez("outputs/run_3d_T{}_{}.npz".format(T[n], ID),
             α=α, ω0=ω0, Γ=Γ, dt=dt, dkmax=dkmax, epsrel=epsrel,
             t=t, sx=sx, sy=sy, sz=sz, T=T[n])

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(12,5))
    ax1.plot(t, sx)
    ax2.plot(t, sy)
    ax3.plot(t, sz)
    fig.tight_layout()
    plt.savefig("outputs/run_3d_T{}_{}.pdf".format(T[n], ID))
    plt.close()

    print("done: ", n+1, "/", len(T))

  np.savez("outputs/ss_run_3d_{}.npz".format(ID),
           α=α, ω0=ω0, Γ=Γ, dt=dt, dkmax=dkmax, epsrel=epsrel,
           T=T, sss=sss)

  fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(12,5))
  ax1.semilogx(T, sss[:,0], '--o')
  ax2.semilogx(T, sss[:,1], '--o')
  ax3.semilogx(T, sss[:,2], '--o')
  fig.tight_layout()
  plt.savefig("outputs/ss_run_3d_{}.pdf".format(ID))

##### Tempo run #####

ID = int(time())

# α, ω0, Γ = 10.0, 2.0, 0.001 # charlie

# # charlie prma[1.0] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [1.0]*3
# dt, dkmax, epsrel = 0.10, 260, 5.0e-5
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 15.0
# # result ID: 1658822870 [crashed after ~1.5]

# # charlie prma[1.0] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [1.0]*3
# dt, dkmax, epsrel = 0.10, 260, 5.0e-5
# T = np.geomspace(1e-2, 1e+1, 3*4)
# # endtime = 12.0
# # result ID: 1658857049 [crashed after ~1.5]

# # charlie prma[0.9] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.9]*3
# dt, dkmax, epsrel = 0.10, 200, 5.0e-5
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 12.5
# # result ID: 1658879444{pts} [stopped after ~1.5]

# # charlie prma[0.8] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.8]*3
# dt, dkmax, epsrel = 0.12, 200, 8.0e-5
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 14.0
# # result ID: 1658882371{pts} [crashed after ~1.5]

# # charlie prma[0.8] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.8]*3
# dt, dkmax, epsrel = 0.05, 285, 1.0e-4
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 14.0
# # result ID: 1658991560 [crashed after ~2.8]

# # charlie prma[0.6] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.6]*3
# dt, dkmax, epsrel = 0.04, 320, 1.0e-4
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 12.5
# # result ID: 1659002443 [stopped after ~0.0187]

# # charlie prma[0.6] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.6]*3
# dt, dkmax, epsrel = 0.05, 310, 1.0e-4
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 15.0
# # result ID: 1659003670 [crashed after ~]

# # charlie prma[0.8] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.8]*3
# dt, dkmax, epsrel = 0.05, 310, 5.0e-5
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 15.0
# # result ID: 1659011510 [stopped after ~2.8]

# # charlie prma[0.7] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.7]*3
# dt, dkmax, epsrel = 0.04, 260, 5.0e-5
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 10.0
# # result ID: 1659022107 [stopped after ~2.8]

# # charlie prma[0.6] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.6]*3
# dt, dkmax, epsrel = 0.02, 260, 5.0e-5
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 10.0
# # result ID: 1659029233 [stopped after ~2.8]

# # charlie prma[0.6] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.6]*3
# dt, dkmax, epsrel = 0.02, 360, 5.0e-5
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 7.0
# # result ID: 1659035468 [stopped after ~5.3]

# # charlie prma[0.7] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.7]*3
# dt, dkmax, epsrel = 0.04, 260, 1.0e-5
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 7.5
# # result ID: 1659044181 [stopped after ~2.8]

# # charlie prma[0.7] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.7]*3
# dt, dkmax, epsrel = 0.04, 260, 5.0e-6
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 7.0
# # result ID: 1659047906 [crashed after ~2.8]

# # GOOD one
# # charlie prma[0.6] prt1
α, ω0, Γ = [10.0]*3, [2.0]*3, [0.6]*3
dt, dkmax, epsrel = 0.04, 260, 5.0e-6
T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 7.5
# # result ID: 1659053378 [stopped after ~2.8]

# # charlie prma[0.55] prt1
# α, ω0, Γ = [10.0]*3, [2.0]*3, [0.55]*3
# dt, dkmax, epsrel = 0.04, 260, 3.0e-6
# T = np.geomspace(1e-2, 1e1, 3*4)
# # endtime = 7.8
# # result ID: 1659059120 []

# tempo_Tsweep_scratch(α, ω0, Γ, T, dt, dkmax, epsrel, ID)
tempo_Tsweep_precomputepts(α, ω0, Γ, T, dt, dkmax, epsrel, ID)
# tempo_Tsweep_dynfrompts(T, '')