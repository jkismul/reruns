def H1(name, network, GLOBALSEED, k):
    """
    Inhib in apic
    """
    import numpy as np
    from LFPy import Synapse
    from param import synapseParameters, stimSynapsePos, stimSynapseNum, distr_t
    import pickle

    # syns = np.zeros((len(network.populations[name].cells),100))
    # syns=np.zeros(len(network.populations[name].cells))
    # idxs =[]
    np.random.seed(GLOBALSEED)
    if name == 'E':
        syns = np.zeros((len(network.populations[name].cells), stimSynapseNum))

        for j, cell in enumerate(network.populations[name].cells):
            np.random.seed(GLOBALSEED)
            if j % 2 == 0:
                np.random.seed(GLOBALSEED)
                wt_E = 0.5  # martin: 0.05, 0.02
                weighttrain = np.random.normal(wt_E, wt_E * 0.1, size=stimSynapseNum)
                np.random.seed(j)
                idx = cell.get_rand_idx_area_and_distribution_norm(nidx=stimSynapseNum, **stimSynapsePos)
                syns[j] = idx
                np.random.seed(GLOBALSEED)
                for n, i in enumerate(idx):
                    syn = Synapse(cell=cell, idx=i, syntype='Exp2Syn',
                                  weight=weighttrain[n],
                                  **dict(synapseParameters[0][0]))
                    syn.set_spike_times(np.array([np.random.normal(loc=5, scale=1.)]))

            if j % 5 == 0:  # and k ==1:
                np.random.seed(GLOBALSEED)
                wt_E = 0.05
                weighttrain = np.random.normal(wt_E, wt_E * 0.1)
                np.random.seed(GLOBALSEED)

                idx = cell.get_closest_idx(x=0, y=0, z=800)

                syn = Synapse(cell=cell, idx=idx, syntype='Exp2Syn',
                              weight=weighttrain,
                              **dict(synapseParameters[0][0]))
                np.random.seed(GLOBALSEED)
                # syn.set_spike_times(np.array([10.5]))
                syn.set_spike_times(np.array([np.random.normal(loc=10, scale=3.)]))

                # syn.set_spike_times(np.array([distr_t[j]+7.5]))
            # SECOND EXCITATORY, APICAL
            if j % 10 == 0 and k == 0:
                np.random.seed(GLOBALSEED)
                wt_E = 0.05
                weighttrain = np.random.normal(wt_E, wt_E * 0.1)
                np.random.seed(GLOBALSEED)

                idx = cell.get_closest_idx(x=0, y=0, z=800)

                syn = Synapse(cell=cell, idx=idx, syntype='Exp2Syn',
                              weight=weighttrain,
                              **dict(synapseParameters[1][0]))
                np.random.seed(GLOBALSEED + 1)
                # syn.set_spike_times(np.array([10.5]))
                syn.set_spike_times(np.array([np.random.normal(loc=10, scale=1.)]))

                # syn.set_spike_times(np.array([distr_t[j]+7.5]))
            fi = open("example_network_output/stimsyn.pkl", "wb")
            # pickle.dump(idxs,fi)
            pickle.dump(syns, fi)
            fi.close()

    np.random.seed(GLOBALSEED)
    # print(idxs)
    # fi = open("example_network_output/stimsyn.pkl","wb")
    # # pickle.dump(idxs,fi)
    # pickle.dump(syns,fi)
    # fi.close()
    # print(np.shape(syns))


def H2(name, network, GLOBALSEED, k):
    """
    Long-Short pulse
    """
    import numpy as np
    import scipy.stats as st
    import pickle
    import LFPy
    from LFPy import Synapse
    from param import distr_t, synapseParameters, stimSynapseNum, stimSynapsePos,\
        stimSynapsePos2,stimSynapseNum2,stimSynapsePos3,stimSynapseNum3, \
        stimSynTim2, stimSynTim3
    np.random.seed(GLOBALSEED)
    # distrrr = st.gamma(a=1.2, loc=6, scale=2.0)
    # ff = distrrr.rvs(size=200)
    if name == 'E':
        syns = np.zeros((len(network.populations[name].cells), stimSynapseNum))
        # syns2 = np.zeros((len(network.populations[name].cells), stimSynapseNum2))

        for j, cell in enumerate(network.populations[name].cells):
            np.random.seed(GLOBALSEED)
            if j % 2 == 0:# and k==0: #%5 for non-network
                # np.random.seed(GLOBALSEED)
                np.random.seed(j)

                wt_E = 0.001#.001 med 300 synapser  # martin: 0.05, 0.02
                # weighttrain = np.random.normal(wt_E, wt_E * 0.1, size=stimSynapseNum)
                np.random.seed(j)
                idx = cell.get_rand_idx_area_and_distribution_norm(nidx=stimSynapseNum, **stimSynapsePos)
                # idx1 = cell.get_rand_idx_area_and_distribution_norm(nidx=int(stimSynapseNum/2), **stimSynapsePos2)
                # idx2 = cell.get_rand_idx_area_and_distribution_norm(nidx=int(stimSynapseNum/2), **stimSynapsePos3)
                # idx = np.hstack((idx1,idx2))

                lower=2
                upper = 4
                mu , sigma=2, 2.5
                X = st.truncnorm(
                    (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
                # N = stats.norm(loc=mu, scale=sigma)
                syns[j] = idx
                # np.random.seed(GLOBALSEED)
                for n, i in enumerate(idx):
                    weit = np.random.normal(wt_E, wt_E * 0.3, size=1)

                    # LFPy.StimIntElectrode(cell=cell,idx=i,pptype='IClamp',amp=9992,dur=0.3,delay=2)
                    syn = Synapse(cell=cell, idx=i, syntype='Exp2Syn',
                                  weight=weit,
                                  # weight=weighttrain[n],
                                  **dict(synapseParameters[0][1]))
                    # syn.set_spike_times(np.array([np.random.normal(loc=8, scale=5.5)]))
                    syn.set_spike_times(np.array(X.rvs(1)+j*0.01))

            # if j % 4 == 0 and k==1: #%5 for non-network
            #     # np.random.seed(GLOBALSEED)
            #     np.random.seed(j)
            #
            #     wt_E = .01#.001  # martin: 0.05, 0.02
            #     weighttrain = np.random.normal(wt_E, wt_E * 0.1, size=stimSynapseNum)
            #     np.random.seed(j)
            #     idx = cell.get_rand_idx_area_and_distribution_norm(nidx=stimSynapseNum, **stimSynapsePos2)
            #     # idx1 = cell.get_rand_idx_area_and_distribution_norm(nidx=int(stimSynapseNum/2), **stimSynapsePos2)
            #     # idx2 = cell.get_rand_idx_area_and_distribution_norm(nidx=int(stimSynapseNum/2), **stimSynapsePos3)
            #     # idx = np.hstack((idx1,idx2))
            #     syns[j] = idx
            #     # np.random.seed(GLOBALSEED)
            #     for n, i in enumerate(idx):
            #         # LFPy.StimIntElectrode(cell=cell,idx=i,pptype='IClamp',amp=9992,dur=0.3,delay=2)
            #         syn = Synapse(cell=cell, idx=i, syntype='Exp2Syn',
            #                       weight=weighttrain[n],
            #                       **dict(synapseParameters[0][0]))
            #         # print(np.array([np.random.normal(loc=8, scale=5.5)]))
            #         syn.set_spike_times(np.array([np.random.normal(loc=8, scale=2.5)]))#4, 1.5
            #         # syn.set_spike_times(np.array([np.random.normal(loc=4, scale=1.0)]))
            #         # syn.set_spike_times(np.array([np.random.normal(loc=5, scale=1.0)]))
            #         # syn.set_spike_times(np.array([np.random.normal(loc=6, scale=1.0)]))

            # if k==0 and j%4 == 0:
            #     np.random.seed(GLOBALSEED)
            #     wt_E = 0.02  # martin: 0.05, 0.02
            #     weighttrain = np.random.normal(wt_E, wt_E * 0.1, size=stimSynapseNum2)
            #     np.random.seed(j)
            #     idx = cell.get_rand_idx_area_and_distribution_norm(nidx=stimSynapseNum2, **stimSynapsePos2)
            #     syns2[j] = idx
            #     np.random.seed(GLOBALSEED+j)
            #     syn_t_idx = np.random.normal(loc=6,scale=5,size=5)
            #     syn_t_idx = syn_t_idx[syn_t_idx>5]
            #     syn_t_idx = syn_t_idx[syn_t_idx<10]
            #     for n, i in enumerate(idx):
            #         syn = Synapse(cell=cell, idx=i, syntype='Exp2Syn',
            #                       weight=weighttrain[n],
            #                       **dict(synapseParameters[0][0]))
            #         syn.set_spike_times(syn_t_idx[::1])


            # if k == 1 and j%4 ==0:
            #     np.random.seed(GLOBALSEED)
            #     wt_E = 0.02  # martin: 0.05, 0.02
            #     weighttrain = np.random.normal(wt_E, wt_E * 0.1, size=stimSynapseNum2)
            #     np.random.seed(j)
            #     idx = cell.get_rand_idx_area_and_distribution_norm(nidx=stimSynapseNum2, **stimSynapsePos2)
            #     syns2[j] = idx
            #     np.random.seed(GLOBALSEED+j)
            #     syn_t_idx = np.random.normal(loc=6,scale=5,size=5)
            #     syn_t_idx = syn_t_idx[syn_t_idx>5]
            #     # syn_t_idx = syn_t_idx[syn_t_idx<10]
            #     for n, i in enumerate(idx):
            #         syn = Synapse(cell=cell, idx=i, syntype='Exp2Syn',
            #                       weight=weighttrain[n],
            #                       **dict(synapseParameters[0][0]))
            #         syn.set_spike_times(syn_t_idx[::1])

            fi = open("example_network_output/stimsyn.pkl", "wb")
            # pickle.dump(idxs,fi)
            pickle.dump(syns, fi)
            fi.close()
    np.random.seed(GLOBALSEED)
