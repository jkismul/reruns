def laminar_LFP(OUTPUT0, OUTPUT1,hypothesis):
    import numpy as np
    import matplotlib.pyplot as plt
    from param import tstop
    from scipy import signal

    plt.figure()
    for n_put, output in enumerate([OUTPUT0, OUTPUT1]):
        dt = 2 ** -5  # 0.1
        # lavp 0-2-300, hp:3-500-none
        laminar_LFP = output[0]['E']
        fs = 1 / dt * 1000
        fc = 300
        w = fc / (fs / 2.)

        b, a = signal.butter(5, w, 'lowpass')

        laminar_LFP = signal.filtfilt(b, a, laminar_LFP)

        timses = np.linspace(0, tstop, len(OUTPUT0[0]['E'][0]))

        for i in reversed(range(13)):
            if n_put == 0:
                dz = np.max(laminar_LFP[i]) - np.min(laminar_LFP[i])
                # dz=0.01
                if i == 0:
                    plt.plot(timses, laminar_LFP[i] / dz + (i - 13), 'k', label='with inhib')  # (13 - i)
                else:
                    plt.plot(timses, laminar_LFP[i] / dz + (i - 13), 'k')
            else:
                dz2 = np.max(laminar_LFP[i]) - np.min(laminar_LFP[i])
                # dz2=0.01
                if i == 0:
                    plt.plot(timses, laminar_LFP[i] / dz2 + (i - 13), color='gray', linestyle='dashed', label='blocked')
                else:
                    plt.plot(timses, laminar_LFP[i] / dz2 + (i - 13), color='gray', linestyle='dashed')

    plt.legend()
    plt.xlabel('time [ms]')
    plt.ylabel('depth')
    plt.savefig('plots/laminar_{}.pdf'.format(hypothesis))

    # for i in reversed(range(13)):
    #     dz = np.max(OUTPUT0[0]['E'][i]) - np.min(OUTPUT0[0]['E'][i])
    #     if i == 0:
    #         plt.plot(timses, OUTPUT0[0]['E'][i] / dz + (i-13), 'k', label='with inhib')#(13 - i)
    #     else:
    #         plt.plot(timses, OUTPUT0[0]['E'][i] / dz + (i-13), 'k')
    #     dz2 = np.max(OUTPUT1[0]['E'][i]) - np.min(OUTPUT1[0]['E'][i])
    #     if i == 0:
    #         plt.plot(timses, OUTPUT1[0]['E'][i] / dz2 + (i-13), color='gray', linestyle='dashed', label='blocked')
    #     else:
    #         plt.plot(timses, OUTPUT1[0]['E'][i] / dz2 + (i-13), color='gray', linestyle='dashed')
    #
    #
    # plt.legend()
    # plt.xlabel('time [ms]')
    # plt.ylabel('depth')
    # plt.savefig('plots/laminar_{}.pdf'.format(hypothesis))

def MUA_fig(OUTPUT0, OUTPUT1,fig,label='e'):
        import numpy as np
        import elephant
        import neo
        import quantities as pq
        from param import tstop
        import matplotlib.pyplot as plt
        from scipy import signal
        csd = []

        num_xlabels = 5
        num_ylabels = 11

        for output in [OUTPUT0, OUTPUT1]:
            num_tsteps = output[0]['E'].shape[1]
            num_elecs = output[0]['E'].shape[0]
            dt = 2 ** -5  # 0.1
            z = np.arange(num_elecs) * 100
            t = np.arange(num_tsteps) * dt
            # lavp 0-2-300, hp:3-500-none
            coords = np.array([z]).T * pq.um
            laminar_LFP = output[0]['E']
            fs = 1 / dt * 1000
            # fc=300
            fc=500
            w= fc/(fs/2.)
            # # b,a = signal.butter(5,w,'lowpass')
            b,a = signal.butter(5,w,'highpass')
            #
            laminar_LFP = signal.filtfilt(b,a,laminar_LFP)
            #
            laminar_LFP = np.abs(laminar_LFP)

            fc = 300
            fc=100
            w = fc / (fs / 2.)

            b, a = signal.butter(5, w, 'lowpass')

            laminar_LFP = signal.filtfilt(b, a, laminar_LFP)

            # plt.figure()
            # plt.imshow(laminar_LFP,aspect='equal')
            # plt.show()
            laminar_LFP = laminar_LFP/np.max(laminar_LFP)
            # print(np.max(laminar_LFP),np.min(laminar_LFP))
            lfp_neo = neo.AnalogSignal(laminar_LFP, units="uV", sampling_rate=1 / dt * 1000 * pq.Hz);
            csd.append(elephant.current_source_density.estimate_csd(lfp_neo.T, coords, diam=100 * pq.um,
                                                                    method="SplineiCSD"));

        ax = [fig.add_subplot(211, label='e'), fig.add_subplot(212, label='f')]

        vmax = np.max(np.abs(csd[0]))

        from matplotlib.colors import ListedColormap
        import matplotlib.colors

        norm = matplotlib.colors.Normalize(-1, 1)
        colors = [[norm(-1.0), "deepskyblue"],
                  [norm(-0.6), "blue"],
                  [norm(-0.2), "darkblue"],
                  [norm(0.0), "black"],
                  [norm(0.2), "darkred"],
                  [norm(0.6), "red"],
                  [norm(1.0), "orange"]]

        cmapp = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
        for i, row in enumerate(ax):
            row.imshow(csd[i].T, cmap=cmapp, vmax=vmax, vmin=-vmax, origin='bottom')
            row.axis('auto')
            timses = np.linspace(0, tstop, num_xlabels)
            depthses = np.linspace(-200, 800, num_ylabels)
            yticker = row.get_yticks()
            limsy = np.linspace(0, np.max(yticker), num_ylabels)
            row.set_yticks(limsy)
            row.set_yticklabels(depthses.astype('int'), rotation=45)
            xticker = row.get_xticks()
            limsx = np.linspace(0, np.max(xticker) + np.min(xticker), num_xlabels)
            row.set_xticks(limsx)
            row.set_xticklabels(timses.astype('int'))
            row.set_xlabel('Time [ms]')
            row.set_ylabel('Depth [nm]')
        return ax


def MUA_lines(OUTPUT0, OUTPUT1,hypothesis,label='e'):
        import numpy as np
        from param import tstop
        import matplotlib.pyplot as plt
        from scipy import signal

        plt.figure()

        for n_put, output in enumerate([OUTPUT0, OUTPUT1]):
            # num_tsteps = output[0]['E'].shape[1]
            # num_elecs = output[0]['E'].shape[0]
            dt = 2 ** -5  # 0.1
            # lavp 0-2-300, hp:3-500-none
            laminar_LFP = output[0]['E']
            fs = 1 / dt * 1000
            fc=500
            w= fc/(fs/2.)

            b,a = signal.butter(5,w,'highpass')

            laminar_LFP = signal.filtfilt(b,a,laminar_LFP)

            laminar_LFP = np.abs(laminar_LFP)
            #
            # fc = 300
            # w = fc / (fs / 2.)
            #
            # b, a = signal.butter(5, w, 'lowpass')
            #
            # laminar_LFP = signal.filtfilt(b, a, laminar_LFP)

            timses = np.linspace(0, tstop, len(OUTPUT0[0]['E'][0]))

            for i in reversed(range(13)):
                if n_put ==0:
                    dz = np.max(laminar_LFP[i]) - np.min(laminar_LFP[i])
                    # dz=0.01
                    if i == 0:
                        plt.plot(timses, laminar_LFP[i] / dz + (i-13), 'k', label='with inhib')#(13 - i)
                    else:
                        plt.plot(timses, laminar_LFP[i] / dz + (i-13), 'k')
                else:
                    dz2 = np.max(laminar_LFP[i]) - np.min(laminar_LFP[i])
                    # dz2=0.01
                    if i == 0:
                        plt.plot(timses, laminar_LFP[i] / dz2 + (i-13), color='gray', linestyle='dashed', label='blocked')
                    else:
                        plt.plot(timses, laminar_LFP[i] / dz2 + (i-13), color='gray', linestyle='dashed')


        plt.legend()
        plt.xlabel('time [ms]')
        plt.ylabel('depth')
        plt.savefig('plots/laminar_MUA{}.pdf'.format(hypothesis))


def LFP_top(OUTPUT0, OUTPUT1, fig, hypothesis, label='a',save_ecog=False):
    import numpy as np
    from param import tstop
    import matplotlib.pyplot as plt
    from scipy import signal

    layer_number = 12
    # layer_number= 10
    t_article = 25

##########################################
    #LOWPASS
    #
    # dt = 2 ** -5
    # fs = 1 / dt * 1000
    # fc = 300
    # fc=100
    # w = fc / (fs / 2.)
    #
    # b, a = signal.butter(5, w, 'lowpass')
    #
    # sig0 = signal.filtfilt(b, a, OUTPUT0[0]['E'])
    # sig1 = signal.filtfilt(b, a, OUTPUT1[0]['E'])
    #
    # ratio = t_article / tstop
    # dist = int(len(OUTPUT0[0]['E'][layer_number]) * ratio)
    #
    # ax = fig.add_subplot(111, label=label)
    # timses = np.linspace(0, t_article, dist)
    #
    # ax.plot(timses, sig0[layer_number][:dist], 'k', label='Inhib on')
    # ax.plot(timses, sig1[layer_number][:dist], 'k--', label='Inhib blocked')
#LOWPASS END
    #######################################
    ratio = t_article / tstop
    dist = int(len(OUTPUT0[0]['E'][layer_number]) * ratio)

    ax = fig.add_subplot(111, label=label)
    timses = np.linspace(0, t_article, dist)
    ax.plot(timses, OUTPUT0[0]['E'][layer_number][:dist], 'k', label='Inhib on')
    ax.plot(timses, OUTPUT1[0]['E'][layer_number][:dist], 'k--', label='Inhib blocked')

    # ax.legend()
    ax.set_xlabel('time [ms]')
    ax.set_ylabel('response [units] mV?')

    if save_ecog:
        ax.legend()
        if hypothesis == 'H1':
            ax.set_title('Inhib in apic')
        if hypothesis == 'H2':
            ax.set_title('Short/Long pulse')
        plt.savefig('plots/Ecog_{}'.format(hypothesis))
        ax.get_legend().remove()
        ax.set_title('')

    return ax


def spike_raster(SPIKES0, SPIKES1, fig, label='a'):
    import numpy as np
    from param import population_names
    # import matplotlib.pyplot as plt
    spik = [SPIKES0, SPIKES1]

    ax = [fig.add_subplot(211, label=label), fig.add_subplot(212, label=label)]

    for k in range(2):
        for name, spts, gids in zip(population_names, spik[k]['times'], spik[k]['gids']):
            t = []
            g = []
            for spt, gid in zip(spts, gids):
                t = np.r_[t, spt]
                g = np.r_[g, np.zeros(spt.size) + gid]
            ax[k].plot(t[t >= 0], g[t >= 0], '.', label=name)
        ax[k].set_xticks(np.arange(0,41,10))
        # plt.xticks(range(1, 5))
        # plt.scatter(spik[k]['times'], spik[k]['gids'], '.', label=name)

            # remove_axis_junk(ax[k], lines=['right', 'top'])
    ax[1].set_xlabel('t [ms]')
    ax[0].set_ylabel('gid')
    ax[1].set_ylabel('gid')
    return ax


def remove_axis_junk(ax, lines=['right', 'top']):
    """remove chosen lines from plotting axis"""
    for loc, spine in ax.spines.items():
        if loc in lines:
            spine.set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


def LFP_heatmap(OUTPUT0, fig, label='A'):
    import numpy as np
    from param import tstop, top_d, bot_d

    from matplotlib.colors import ListedColormap
    import matplotlib.colors

    norm = matplotlib.colors.Normalize(-1, 1)
    colors = [[norm(-1.0), "deepskyblue"],
              [norm(-0.6), "blue"],
              [norm(-0.2), "darkblue"],
              [norm(0.0), "black"],
              [norm(0.2), "darkred"],
              [norm(0.6), "red"],
              [norm(1.0), "orange"]]

    cmapp = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

    ax = fig.add_subplot(111, label=label)

    n_height_ticks = 11
    n_time_ticks = 5

    timses = np.linspace(0, tstop, n_time_ticks)
    depthses = np.linspace(bot_d, top_d, n_height_ticks)
    vmax = np.max(np.abs(OUTPUT0[0]['E'])) / 2
    # vmax=OUTPUT0[0]['E'].std()

    ax.imshow(OUTPUT0[0]['E'], cmap=cmapp, vmax=vmax, vmin=-vmax, origin='bottom')
    ax.axis('auto')
    yticker = ax.get_yticks()
    limsy = np.linspace(0, np.max(yticker) + np.min(yticker), n_height_ticks)
    ax.set_yticks(limsy)
    ax.set_yticklabels(depthses.astype('int'), rotation=45)
    xticker = ax.get_xticks()
    limsx = np.linspace(0, np.max(xticker) + np.min(xticker), n_time_ticks)
    ax.set_xticks(limsx)
    ax.set_xticklabels(timses.astype('int'))
    ax.set_xlabel('Time [ms]')
    ax.set_ylabel('Depth [nm]')  # is it nm?
    return ax


def Population():
    import LFPy
    import pickle
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np

    # Seems importing the network instance is clumsy?
    # Recreating network
    GLOBALSEED = 1234
    np.random.seed(GLOBALSEED)
    from param import populationParameters, networkParameters, population_names, population_sizes

    # network = LFPy.Network(**networkParameters)
    # for ii, (name, size) in enumerate(zip(population_names, population_sizes)):
    #     np.random.seed(GLOBALSEED)
    #
    #     network.create_population(
    #             name=name, POP_SIZE=size, **populationParameters[ii])
    #
    # for j, cell in enumerate(network.populations['E'].cells):
    #     if j == 0:
    #         print("segments in network cell: ",cell.totnsegs)

    fig = plt.figure(figsize=(6.4, 4.8 * 2))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=5)
    import h5py
    filename = 'example_network_output/cell_positions_and_rotations.h5'
    f = h5py.File(filename, 'r')

    with open('example_network_output/cells.pkl', 'rb') as handle:
        CellParams = pickle.loads(handle.read())
    with open('example_network_output/stimsyn.pkl', 'rb') as handle:
        stimidx = pickle.loads(handle.read())

    flyndre = LFPy.Cell(morphology=CellParams[0]['morphology'])
    # print("segments in cell: ", flyndre.totnsegs)
    cells = [LFPy.Cell(morphology=CellParams[0]['morphology']), LFPy.Cell(morphology=CellParams[1]['morphology'])]
    # cells = [network.populations['E'].cells, network.populations['I'].cells]

    for i, (name, pop) in enumerate(f.items()):
        if i == 0:
            # for name,pop in f.items():
            pop = len(f[name])
            cell = cells[i]
            # cell = network.populations[name].cells[i]
            for j in range(pop):
                if j == 0:
                    print("segments in network cell: ", cell.totnsegs)

                cell.set_pos(f[name][(j)][(1)], f[name][(j)][(2)], f[name][(j)][(3)])
                cell.set_rotation(0, 0, f[name][(j)][(6)])
                c = 'C0' if name == 'E' else 'C1'
                for idx in range(cell.totnsegs):
                    if idx == 0:  # soma
                        ax.scatter(cell.xmid[idx], cell.ymid[idx], cell.zmid[idx], c=c)
                    else:  # dendrites
                        ax.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]],
                                [cell.zstart[idx], cell.zend[idx]], c)
                    if idx in stimidx[j]:  # synapses
                        ax.scatter(cell.xmid[idx], cell.ymid[idx], cell.zmid[idx], c='y', marker='*', s=15)
    f.close()

    ax.set_xlabel('$x$ ($\mu$m)')
    ax.set_ylabel('$y$ ($\mu$m)')
    ax.set_zlabel('$z$ ($\mu$m)')
    ax.set_title('network populations')
    fig.savefig('plots/population.pdf')
    plt.close(fig)


def iCSD(OUTPUT0, OUTPUT1, fig, label='e'):
    import numpy as np
    import elephant
    import neo
    import quantities as pq
    from param import tstop
    import matplotlib.pyplot as plt
    from scipy import signal
    csd = []

    num_xlabels = 5
    num_ylabels = 11

    for output in [OUTPUT0, OUTPUT1]:
        num_tsteps = output[0]['E'].shape[1]
        num_elecs = output[0]['E'].shape[0]
        dt = 2 ** -5  # 0.1
        z = np.arange(num_elecs) * 100
        t = np.arange(num_tsteps) * dt
    #lavp 0-2-300, hp:3-500-none
        coords = np.array([z]).T * pq.um
        laminar_LFP = output[0]['E']
        fs=1 / dt * 1000
        # fc=300
        # fc=500
        # w= fc/(fs/2.)
        # # b,a = signal.butter(5,w,'lowpass')
        # b,a = signal.butter(5,w,'highpass')
        #
        # laminar_LFP = signal.filtfilt(b,a,laminar_LFP)
        #
        # laminar_LFP = np.abs(laminar_LFP)

        fc=300
        w= fc/(fs/2.)

        b,a = signal.butter(5,w,'lowpass')

        laminar_LFP = signal.filtfilt(b,a,laminar_LFP)

        # laminar_LFP = laminar_LFP/np.max(laminar_LFP)
        # print(np.max(laminar_LFP),np.min(laminar_LFP))
        lfp_neo = neo.AnalogSignal(laminar_LFP, units="uV", sampling_rate=1 / dt * 1000 * pq.Hz);
        csd.append(elephant.current_source_density.estimate_csd(lfp_neo.T, coords, diam=100 * pq.um,
                                                                method="SplineiCSD"));

    ax = [fig.add_subplot(211, label='e'), fig.add_subplot(212, label='f')]

    vmax = np.max(np.abs(csd[0]))

    from matplotlib.colors import ListedColormap
    import matplotlib.colors

    norm = matplotlib.colors.Normalize(-1, 1)
    colors = [[norm(-1.0), "deepskyblue"],
              [norm(-0.6), "blue"],
              [norm(-0.2), "darkblue"],
              [norm(0.0), "black"],
              [norm(0.2), "darkred"],
              [norm(0.6), "red"],
              [norm(1.0), "orange"]]

    cmapp = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    for i, row in enumerate(ax):
        row.imshow(csd[i].T, cmap=cmapp, vmax=vmax, vmin=-vmax, origin='bottom')
        row.axis('auto')
        timses = np.linspace(0, tstop, num_xlabels)
        depthses = np.linspace(-200, 800, num_ylabels)
        yticker = row.get_yticks()
        limsy = np.linspace(0, np.max(yticker), num_ylabels)
        row.set_yticks(limsy)
        row.set_yticklabels(depthses.astype('int'), rotation=45)
        xticker = row.get_xticks()
        limsx = np.linspace(0, np.max(xticker) + np.min(xticker), num_xlabels)
        row.set_xticks(limsx)
        row.set_xticklabels(timses.astype('int'))
        row.set_xlabel('Time [ms]')
        row.set_ylabel('Depth [nm]')
    return ax


def make_plot(OUTPUT0, OUTPUT1, SPIKES0, SPIKES1,hypothesis,save_ecog=False):
    import matplotlib.pyplot as plt

    fig = plt.figure()

    ax2 = LFP_top(OUTPUT0, OUTPUT1, fig, hypothesis, 'c',save_ecog)
    ax2.set_position([0.1, 0.6, 0.2, 0.3])
    ax2.set_title('ECoG')

    ax0 = LFP_heatmap(OUTPUT0, fig, 'a')
    ax1 = LFP_heatmap(OUTPUT1, fig, 'b')
    ax0.set_position([0.4, 0.6, 0.2, 0.3])
    ax0.set_title('LFP inhib on')
    ax1.set_position([0.4, 0.1, 0.2, 0.3])  # left, bottom, width, height
    ax1.set_title('LFP inhib off')

    ax3 = spike_raster(SPIKES0, SPIKES1, fig, 'd')
    ax3[1].set_position([0.1, 0.1, 0.2, 0.15])
    ax3[0].set_position([0.1, 0.25, 0.2, 0.15])
    ax3[0].set_title('Spike Raster')

    ax4 = iCSD(OUTPUT0, OUTPUT1, fig, 'e')
    ax4[1].set_position([0.7, 0.1, 0.2, 0.3])
    ax4[1].set_title('iCSD inhib off')
    ax4[0].set_position([0.7, 0.6, 0.2, 0.3])
    ax4[0].set_title('iCSD inhib on')

    plt.savefig('plots/total_{}.pdf'.format(hypothesis))

# def plot_single_layer(arr,OUTPUT0,OUTPUT1,layer_n):
#     import matplotlib.pyplot as plt
#     from param import tstop
#     import numpy as np
#
#     n = len(arr)
#     fig = plt.figure()
#     ax=[]
#     t_article = 25
#
#     print(np.shape(OUTPUT0[layer_n]['E']))
#     ratio = t_article / tstop
#     dist = int(len(OUTPUT0[layer_n]['E'][0]) * ratio)
#
#     for i in arr:
#         ax.append(fig.add_subplot(2,n,i+1,label='{}'.format(i)))
#         ax[i].plot(-OUTPUT0[layer_n]['E'][i][:dist])
#         ax[i].set_aspect('auto')
#         ax[i].set_xticks([])
#         ax[i].set_yticks([])
#         ax[i].set_title('Cell {}'.format(i))
#     for i in arr:
#         ax.append(fig.add_subplot(2,n,n+i+1,label='{}'.format(n+i)))
#         ax[n+i].plot(-OUTPUT1[layer_n]['E'][i][:dist],c='r')
#         ax[n+i].set_aspect('auto')
#         ax[n+i].set_yticks([])
#         ax[n+i].set_xticks([])
#         ax[n+i].set_xticks([0,int(dist*(10/25)),int(dist*(20/25))])
#         ax[n+i].set_xticklabels([0,10,20])
#     ax[0].set_ylabel('Inhib on')
#     ax[5].set_ylabel('Inhib blocked')
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('plots/sing.png')

def plot_soma(SomaVs,SomaTime,arr):
    import matplotlib.pyplot as plt
    from param import tstop
    import numpy as np

    n = len(arr)
    fig = plt.figure()
    ax=[]
    t_article = 25

    ratio = t_article / tstop
    dist = int(len(SomaVs['inhib'][0]) * ratio)

    for i in arr:
        ax.append(fig.add_subplot(2,n,i+1,label='{}'.format(i)))
        ax[i].plot(SomaVs['inhib'][i][:dist])
        ax[i].set_aspect('auto')
        ax[i].set_xticks([])
        # ax[i].set_yticks([])
        ax[i].set_title('Cell {}'.format(i))
    for i in arr:
        ax.append(fig.add_subplot(2,n,n+i+1,label='{}'.format(n+i)))
        ax[n+i].plot(SomaVs['blocked'][i][:dist],c='r')
        ax[n+i].set_aspect('auto')
        # ax[n+i].set_yticks([])
        ax[n+i].set_xticks([])
        ax[n+i].set_xticks([0,int(dist*(10/25)),int(dist*(20/25))])
        ax[n+i].set_xticklabels([0,10,20])
    ax[0].set_ylabel('Inhib on')
    ax[5].set_ylabel('Inhib blocked')
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('plots/SomaV.png')