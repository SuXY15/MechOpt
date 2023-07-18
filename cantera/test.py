import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

# mech = "../mechanism/nordin_34s121r.yaml"
# props ={
#     'T': 1000.,                # Temperature
#     'P': 10.,                  # Pressure
#     'phi': 1.0,                # equivalence ratio
#     'air':'O2:1.0,N2:3.76',    # air components
#     'fuel':'C7H16:1.0',        # fuel components
#     'tot': .1,                 # total time
# }

mech = "h2o2.yaml"
props ={
    'T': 1000.,                # Temperature
    'P': 1.,                   # Pressure
    'phi': 1.0,                # equivalence ratio
    'air':'O2:1.0,AR:3.76',    # air components
    'fuel':'H2:1.0',           # fuel components
    'tot': 1e-3,               # total time
}

gas = ct.Solution(mech)

print(gas.derivative_settings)

gas.set_equivalence_ratio(props['phi'], props['fuel'], props['air'])
gas.TP = props['T'], props['P']*ct.one_atm
r = ct.IdealGasReactor(gas)
sim = ct.ReactorNet([r])

states = ct.SolutionArray(gas, extra=['t', 'ddT'])
while sim.time < props['tot']:
    sim.step()
    states.append( r.thermo.state, t=sim.time, ddT=gas.net_production_rates_ddT)
data = {'t':states.t, 'T':states.T, 'ddT': states.ddT}

ax = plt.subplot()
ax.plot(data['t']*1000, data['T'])

bx = ax.twinx()
bx.plot(data['t']*1000, data['ddT'], '--')

plt.show()