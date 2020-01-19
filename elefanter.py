import numpy as np
import elephant
import neo
import quantities as pq
#ai=np.arange(0,20,)
ai = np.linspace(0,20,7)
print(ai.astype('int'))
num_tsteps = 100
num_elecs = 10
dt = 0.1
z = np.arange(num_elecs) * 100
t = np.arange(num_tsteps) * dt

coords = np.array([z]).T * pq.um

laminar_LFP = np.random.normal(size=(num_elecs, num_tsteps))
pos_neo = z * pq.um
lfp_neo = neo.AnalogSignal(laminar_LFP, units="uV", sampling_rate=1 / dt* 1000 * pq.Hz)
csd_estimate = elephant.current_source_density.estimate_csd(lfp_neo.T, coords, diam=100*pq.um, #h=25*pq.um,
                                                            method="SplineiCSD")
print(csd_estimate.units)