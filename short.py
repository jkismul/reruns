def conn(network, GLOBALSEED, loop_n):
    import numpy as np
    from LFPy import Network
    from param import maxdelay, connectionProbability, population_names, synapseModel, synapseParameters, \
        weightArguments, \
        minweight, delayFunction, delayArguments, mindelay, multapseFunction, multapseArguments, \
        synapsePositionArguments

    # maxdelay=1000


    if loop_n == 1:
        maxdelay = 1000
        # connectionProbability = [[0.105, 0.], [0., 0.]]  # 0.02 er ganske flott for alle
        connectionProbability[0][1]=0
        connectionProbability[1][0] = 0
        connectionProbability[1][1] = 0
        # synapseParameters = [[dict(tau1=0.1, tau2=3.0, e=0.),  # AMPA (example 7) t1=1,t2=3
        #                       dict(tau1=0.1, tau2=1.0, e=0.)],  # AMPA
        #                      [dict(tau1=1.0, tau2=12.0, e=-80.),  # GABA_A (testihng 70, was 80)
        #                       dict(tau1=1.0, tau2=12.0, e=-80.)]]  # GABA_A
    for i, pre in enumerate(population_names):
        # if loop_n == 1:
        #     connectionProbability = [[0.2, 0.], [0., 0.]]  # 0.02 er ganske flott for alle
        #     print(delayArguments)
        # EEd = 4.  # was 4.0, scale 1.0
        # delayArguments = [[dict(loc=EEd, scale=EEd * 1.25),  # 0.3
        #                    dict(loc=4.0, scale=0.3)],
        #                   [dict(loc=4.0, scale=0.3),
        #                    dict(loc=4.0, scale=0.3)]]
        # mindelay = 1.0  # 3.5
        for j, post in enumerate(population_names):
            # boolean connectivity matrix between pre- and post-synaptic neurons
            # in each population (postsynaptic on this RANK)
            np.random.seed(GLOBALSEED)

            connectivity = network.get_connectivity_rand(
                pre=pre, post=post,
                connprob=connectionProbability[i][j]
            )
            # connect network:
            np.random.seed(GLOBALSEED)

            (conncount, syncount) = network.connect(
                pre=pre, post=post,
                connectivity=connectivity,
                syntype=synapseModel,
                synparams=synapseParameters[i][j],
                weightfun=np.random.normal,
                weightargs=weightArguments[i][j],
                minweight=minweight,
                delayfun=delayFunction,
                delayargs=delayArguments[i][j],
                mindelay=mindelay,
                maxdelay=maxdelay,
                multapsefun=multapseFunction,
                multapseargs=multapseArguments[i][j],
                syn_pos_args=synapsePositionArguments[i][j],
                # save_connections=False,  # Creates synapse_positions.h5
                save_connections=True,  # Creates synapse_positions.h5

            )


def clean_up(network):
    import neuron
    network.pc.gid_clear()  # allows assigning new gids to threads
    electrode = None
    syn = None
    synapseModel = None
    for population in network.populations.values():
        for cell in population.cells:
            cell = None
        population.cells = None
    population = None
    pop = None
    network = None
    neuron.h('forall delete_section()')
