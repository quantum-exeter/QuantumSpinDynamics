# QuantumSpinDynamics
Julia code to solve the open-system of a single quantum spin coupled to multiple heat baths, both statically and dynamically. This is the code used to generate the quantum curves in the [Enhanced entanglement in multi-bath spin-boson models](https://arxiv.org/abs/2306.11036) paper.

### Repo structure

* **lib**: contains the code to calculate the Gibbs, mean force (MF) and dynamic states, as well as the corresponding spin expectation values. Stuctured as two Julia modules: ```Statics.jl``` and ```Dynamics.jl```
* **run**: contains code from ```lib``` to generate the paper data
* **paper_data**: contains the data generated for the paper
* **plot**: contains Jupyter notebooks with code to plot this data
* **tempo**: contains code to run TEMPO as a verification of the RC results - see [Exact Dynamics of Nonadditive Environments in Non-Markovian Open Quantum Systems](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.3.010321) and the [OQuPy package](https://github.com/tempoCollaboration/OQuPy) 
