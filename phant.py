def phant(OUTPUT0,OUTPUT1):
    import numpy as np
    import elephant
    import neo
    import quantities as pq
    import matplotlib.pyplot as plt
    from param import tstop

    csd=[]

    num_xlabels = 5
    num_ylabels = 11

    for output in [OUTPUT0,OUTPUT1]:
        num_tsteps = output[0]['E'].shape[1]
        num_elecs = output[0]['E'].shape[0]
        dt = 2 ** -5  # 0.1
        z = np.arange(num_elecs) * 100
        t = np.arange(num_tsteps) * dt

        coords = np.array([z]).T * pq.um
        laminar_LFP = output[0]['E']
        lfp_neo = neo.AnalogSignal(laminar_LFP, units="uV", sampling_rate=1 / dt * 1000 * pq.Hz);
        csd.append(elephant.current_source_density.estimate_csd(lfp_neo.T, coords, diam=100 * pq.um,
                                                             method="SplineiCSD"));

    fig,ax=plt.subplots(nrows=2,ncols=1)
    vmax = np.max(np.abs(csd[0]))

    from matplotlib.colors import ListedColormap
    import matplotlib.colors

    norm = matplotlib.colors.Normalize(-1, 1)
    colors = [[norm(-1.0), "deepskyblue"],
              [norm(-0.6),"blue"],
              [norm(-0.2), "darkblue"],
              [norm(0.0),"black"],
              [norm(0.2), "darkred"],
              [norm(0.6),"red"],
              [norm(1.0), "orange"]]

    cmapp = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    for i,row in enumerate(ax):
        row.imshow(csd[i].T,cmap=cmapp,vmax=vmax,vmin=-vmax,origin='bottom')

        im = row.get_images()
        extent = im[0].get_extent()
        row.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2])))
        timses = np.linspace(0, tstop, num_xlabels)
        depthses = np.linspace(-200,800,num_ylabels)
        yticker = row.get_yticks()
        limsy = np.linspace(0, np.max(yticker), num_ylabels)
        row.set_yticks(limsy)
        row.set_yticklabels(depthses.astype('int'), rotation=45)
        xticker = row.get_xticks()
        limsx = np.linspace(0,np.max(xticker)+np.min(xticker),num_xlabels)
        row.set_xticks(limsx)
        row.set_xticklabels(timses.astype('int'), rotation=45)
        row.set_xlabel('Time [ms]')
        row.set_ylabel('Depth [nm]')
    return ax
    # fig.savefig('plots/phant.pdf')