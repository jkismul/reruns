from LFPy import NetworkCell
import numpy as np
import neuron
import scipy.stats as st

# relative path for simulation output:
OUTPUTPATH = 'example_network_output'

tstart = -200.
tstop = 40.  # bryuns-H has 20 and 25 for eeg, 40 for lfp and icsd

soma_d = 0
top_d = soma_d + 800
bot_d = soma_d - 200
stellate_d = 1000
# class NetworkCell parameters:z
cellParameters = [dict() for x in range(2)]  # NOTE: 2 is the number of different cell/pop types

# L5 Excitatory
cellParameters[0] = {  # set params as usual, just include template stuff
    'morphology': 'morphologies/ball_and_2_sticks.hoc',  # load cell morphology
    'templatefile': 'morphologies/ball_and_2_sticks_Template.hoc',  # load template of cell, to create network
    'templatename': 'ball_and_2_sticks_Template',  # template file can hold several templates

    # 'morphology': 'morphologies/ball_and_2_sticks_pas.hoc',  # load cell morphology
    # 'templatefile': 'morphologies/ball_and_2_sticks_pas_Template.hoc',  # load template of cell, to create network
    # 'templatename': 'ball_and_2_sticks_pas_Template',  # template file can hold several templates

    'cm': 1.0,  # membrane capacitance
    'Ra': 150,  # axial resistance
    'tstart': tstart,  # start time
    'passive': False,  # switch on passive mechs
    'nsegs_method': 'lambda_f',  # method for setting number of segments,
    'lambda_f': 100,  # segments are isopotential at this frequency
    #   'passive_parameters': {'g_pas': 0.0002, 'e_pas': -70.},  # passive params
    'tstop': tstop,  # stop time
    'templateargs': None,  # JFK: no idea. Parameters provided to template-definition
    'delete_sections': False,  # JFK: no idea. delete pre-existing section-references.
}

# L4 inhibitory
cellParameters[1] = {  # set params as usual, just include template stuff
    'morphology': 'morphologies/Stellate.hoc',  # load cell morphology
    'templatefile': 'morphologies/Stellate_Template.hoc',  # load template of cell, to create network
    'templatename': 'Stellate_Template',  # template file can hold several templates
    'cm': 1.0,  # membrane capacitance
    'Ra': 150,  # axial resistance
    'tstart': tstart,  # start time
    'passive': False,  # switch on passive mechs
    'nsegs_method': 'lambda_f',  # method for setting number of segments,
    'lambda_f': 100,  # segments are isopotential at this frequency
    #  'passive_parameters': {'g_pas': 0.0002, 'e_pas': -70.},  # passive params ## used 65 a lot!
    'tstop': tstop,  # stop time
    'templateargs': None,  # JFK: no idea. Parameters provided to template-definition
    'delete_sections': False,  # JFK: no idea. delete pre-existing section-references.
}

# class Population parameters
populationParameters = [dict() for x in range(2)]

# L5 excitatory pop
populationParameters[0] = (
    dict(
        Cell=NetworkCell,  # to create network
        cell_args=cellParameters[0],  # load params from above
        pop_args=dict(  # parameters for population
            radius=100.,  # place within this radius, micrometer
            loc=soma_d,  # with this avg
            scale=20.),  # and this std
        rotation_args=dict(x=0., y=0.),  # Cells are randomly rotated around z-axis using the Cell.set_rotation method.
    )
)

# L4 inhibitory pop
populationParameters[1] = (
    dict(
        Cell=NetworkCell,  # to create network
        cell_args=cellParameters[1],  # load params from above
        pop_args=dict(  # parameters for population
            radius=100.,  # place within this radius, micrometer
            loc=stellate_d,
            # 350 == (525\2+190\2) from layer thickness microcircuit                                     # with this avg
            scale=20.),  # and this std
        rotation_args=dict(x=0., y=0.),  # Cells are randomly rotated around z-axis using the Cell.set_rotation method.
    )
)

# class Network parameters:
networkParameters = dict(
    dt=2 ** -5,
    tstart=tstart,
    tstop=tstop,
    v_init=-64.,
    celsius=37.0,  # JFK: Mainen was done at 23C, using a Qfactor to upscale to 37.
    OUTPUTPATH=OUTPUTPATH
)
# class RecExtElectrode parameters:
electrodeParameters = dict(  # creates 13 electrodes at x,y origin, every 100micrometer
    x=np.zeros(13),
    y=np.zeros(13),
    # z=np.linspace(top_d, bot_d, 13),
    z=np.linspace(bot_d, top_d, 13),
    N=np.array([[0., 1., 0.] for _ in range(13)]),
    # Normal vectors [x, y, z] of each circular electrode contact surface, default None
    r=5.,  # radius of each contact surface, default None
    n=50,
    # if N is not None and r > 0, the number of discrete points used to compute the n-point average potential on each circular contact point.
    sigma=0.3,  # extracellular conductivity in units of (S/m).
    method="soma_as_point",
)

# method Network.simulate() parameters:
networkSimulationArguments = dict(
    rec_current_dipole_moment=True,  # If True, compute and record current-dipole moment from transmembrane currents
    rec_pop_contributions=True,
    # If True, compute and return single-population contributions to the extracellular potential during simulation time
    to_memory=True,
    to_file=True,  # Creates OUTPUT.h5
    rec_vmem=True,
    rec_imem=True
)

# population names, sizes and connection probability:
num_cells = 200#60
population_names = ['E', 'I']
E_percentage = 0.80
I_percentage = 1.0 - E_percentage

population_sizes = [int(E_percentage * num_cells), int(I_percentage * num_cells)]
connectionProbability = [[0.005, 0.05], [0.1, 0.1]]  # 0.02 er ganske flott for alle
# connectionProbability = [[0.0, 0.], [0., 0.]]  # 0.02 er ganske flott for alle

# synapse model. All corresponding parameters for weights,
# connection delays, multapses and layerwise positions are
# set up as shape (2, 2) nested lists for each possible
# connection on the form:
# [["E:E", "E:I"],
#  ["I:E", "I:I"]].
weighttrain = 0.007
synapseModel = neuron.h.Exp2Syn
# synapse parameters
synapseParameters = [[dict(tau1=0.1, tau2=1.0, e=0.),  # AMPA (example 7) t1=1,t2=3
                      dict(tau1=0.1, tau2=1.0, e=0.)],  # AMPA
                     [dict(tau1=1.0, tau2=12.0, e=-80.),  # GABA_A (testihng 70, was 80)
                      dict(tau1=1.0, tau2=12.0, e=-80.)]]  # GABA_A
# synapse max. conductance (function, mean, st.dev., min.):
weightFunction = np.random.normal
EEw = .105 # was 0.002
EIw = 0.01  # was 0.002
IEw = 0.01  # was 0.01
IIw = 0.01  # was 0.01
weightArguments = [[dict(loc=EEw, scale=EEw * 0.1),
                    dict(loc=EIw, scale=EIw * 0.1)],  # 0.003, 0.0002 i gamle modell
                   [dict(loc=IEw, scale=IEw * 0.1),  # 0.01, 0.001 i gamle modell
                    dict(loc=IIw, scale=IIw * 0.1)]]  # 0.0035,0.004 har virket fint her og.
minweight = 0.
# conduction delay (function, mean, st.dev., min.):
delayFunction = np.random.normal  # prøve poisson?
EEd = 7.6  # was 4.0, scale 1.0
EId = 5.3
IEd = 5.3
delayArguments = [[dict(loc=EEd, scale=EEd * .2),  # 0.3
                   dict(loc=EId, scale=EId*0.3)],
                  [dict(loc=IEd, scale=IEd*0.3),
                   dict(loc=2.0, scale=0.3)]]
mindelay = .1  # 3.5
maxdelay = 1000.
multapseFunction = np.random.normal
multapseArguments = [[dict(loc=2., scale=.5), dict(loc=2., scale=.5)],
                     [dict(loc=5., scale=1.), dict(loc=5., scale=1.)]]
# method NetworkCell.get_rand_idx_area_and_distribution_norm
# parameters for layerwise synapse positions:
# stimSynapsePos = dict(section='allsec',fun=st.norm, funargs=dict(loc=0, scale=500.))
######################################################################################3
stimSynapsePos = dict(section='allsec', fun=st.loggamma, funargs=dict(c=0.1, loc=400, scale=120.))
# stimSynapsePos2 = dict(section='allsec',fun=st.loggamma, funargs=dict(c = 0.1,loc=400, scale=120.))
stimSynapsePos2 = dict(section='dend', fun=st.loggamma, funargs=dict(c=0.1, loc=400, scale=120.))
stimSynapsePos3 = dict(section='apic', fun=st.loggamma, funargs=dict(c=0.1, loc=400, scale=120.))

# stimSynapsePos2 = dict(section='allsec',fun=st.gamma, funargs=dict(a = 3.,loc=400, scale=200.))
# stimSynapsePos2 = dict(section='allsec',fun=st.gamma, funargs=dict(a = 3.,loc=-100, scale=200.))

# stimSynapsePos2 = dict(section='apic',fun=st.norm, funargs=dict(loc=800, scale=100.))
# stimSynapsePos2 = dict(section='dend',fun=st.norm, funargs=dict(loc=bot_d+50, scale=20.))

# stimSynapsePos3 = dict(section='apic',fun=st.gamma, funargs=dict(a = 2.5,loc=400, scale=100.))
# stimSynapsePos3 = dict(section='allsec',fun=st.loggamma, funargs=dict(c = 0.1,loc=400, scale=120.))
# stimSynapsePos3 = dict(section='apic',fun=st.norm,funargs=dict(loc=700,scale=10))
##########################################################################################
# stimSynapsePos = dict(section='soma',fun=st.norm, funargs=dict(loc=0, scale=500.))


# stimSynapsePos = dict(section='apic', fun=st.uniform, funargs=dict(loc=600, scale=300.))
# stimSynapsePos = dict(section='apic', fun=st.uniform, funargs=dict(loc=-600, scale=800.))


distr3 = st.gamma(a=1.2, loc=6, scale=1)
distr4 = st.gamma(a=1.2, loc=6, scale=2)
distr5 = st.gamma(a=1.2, loc=6, scale=4)
distr6 = st.gamma(a=1.2, loc=6, scale=6)

stimSynTim2 = distr3.rvs(150)
stimSynTim3 = np.concatenate((distr4.rvs(5), distr5.rvs(5), distr6.rvs(5)))
stimSynapseNum = 30#300
stimSynapseNum2 = 30
stimSynapseNum3 = 30

stimSynapseTime = np.random.normal(loc=2, scale=2)
synapsePositionArguments = [  # [dict(section=['soma', 'apic'],
    #  fun=[st.norm, st.norm], ### soma greia lager vel bare krøll her, skal ikke koble E til soma
    #  funargs=[dict(loc=78000.,scale=100.),#loc=77500., scale=0.0001),#750ish,
    #           dict(loc=78000., scale=100.)],
    #  funweights=[0.5, 1.]
    # ) for _ in range(2)],

    [dict(section='apic',
          #  fun=st.norm,  ### soma greia lager vel bare krøll her, skal ikke koble E til soma
          fun=st.gamma,  # loc=400 var fin
          funargs=dict(a=3., loc=400, scale=200.),  # loc=700 var fin  # loc=77500., scale=0.0001),#750ish,
          # funargs=dict(a=3., loc=50, scale=200.),  # loc=700 var fin  # loc=77500., scale=0.0001),#750ish,

          # funargs=dict(a=3.,loc=top_d-200, scale=100.),  #loc=700 var fin  # loc=77500., scale=0.0001),#750ish,
          # loc=top_d,scale=100
          funweights=0.5,
          ),
     dict(section='soma',
          fun=st.norm,  ### soma greia lager vel bare krøll her, skal ikke koble E til soma
          # funargs=dict(loc=top_d, scale=100.),  # loc=77500., scale=0.0001),#750ish,
          funargs=dict(loc=stellate_d, scale=100.),  # loc=77500., scale=0.0001),#750ish,

          funweights=1.0,
          )],
    [dict(section='dend',  # soma',#'allsec',
          fun=st.gamma,
          funargs=dict(a=3., loc=bot_d, scale=200.),  # loc=77500., scale=0.0001),#750ish, loc satt till stellate loc?

          # fun=st.norm,  ### soma greia lager vel bare krøll her, skal ikke koble E til soma
          # funargs=dict(loc=bot_d, scale=500.),  # loc=77500., scale=0.0001),#750ish,

          funweights=1.0,
          ),
     dict(section='soma',
          fun=st.norm,  ### soma greia lager vel bare krøll her, skal ikke koble E til soma
          funargs=dict(loc=top_d, scale=100.),  # loc=77500., scale=0.0001),#750ish,

          funweights=0.5,
          )]
]

spike_t = 2.
spike_std = 1.
np.random.seed(12)
distr_t = np.random.normal(spike_t, spike_std, size=population_sizes[0] + population_sizes[1])
