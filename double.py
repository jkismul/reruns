#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import modules:
import os
import sys
import numpy as np
import time
import pickle
# import LFPy
# LFPy.__file__
from LFPy import Network, RecExtElectrode
# import self made modules
from phant import phant
from stims import H1,H2
from short import conn,clean_up
from plotter import laminar_LFP, MUA_lines, LFP_top, spike_raster,LFP_heatmap,Population,make_plot,plot_soma, MUA_fig



from param import cellParameters, OUTPUTPATH, populationParameters, \
                  networkParameters, electrodeParameters, networkSimulationArguments,\
                  population_names, population_sizes

# avoid same sequence of random numbers from numpy and neuron on each RANK,
# e.g., in order to draw unique cell and synapse locations and random synapse
# activation times
GLOBALSEED = 1234
np.random.seed(GLOBALSEED)

#handling of choice of hypothesis
if len(sys.argv)>1:
    hypothesis = sys.argv[1]
else:
    print('No hypothesis chosen. Defaulting to H1.')
    hypothesis = 'H1'
if hypothesis != 'H1' and hypothesis != 'H2':
    print('Hypothesis "{}" not defined. Defaulting to H1.'.format(hypothesis))

fi = open("example_network_output/cells.pkl", "wb")
pickle.dump(cellParameters, fi)
fi.close()

if __name__ == '__main__':
    start = time.time()
    # create directory for output:
    if not os.path.isdir(OUTPUTPATH):
        os.mkdir(OUTPUTPATH)

    SomaVs = {'inhib': [],
              'blocked': []}
    MUA =  {'inhib': [],
              'blocked': []}
    SomaTime = 0
    ############################################################################
    # Main simulation
    ############################################################################
    for k in range(2): #Running for inhib on and off
        # if k==0:
            # cellParameters[0]['passive']= True
            # connectionProbability = [[0.08, 0.05], [0.05, 0.05]]  # 0.02 er ganske flott for alle

        # else:
            # cellParameters[0]['passive']= False
            # connectionProbability = [[0.08, 0.], [0., 0.]]  # 0.02 er ganske flott for alle

        ###     create network
        network = Network(**networkParameters)
        # create E and I populations:
        for ii, (name, size) in enumerate(zip(population_names, population_sizes)):
            np.random.seed(GLOBALSEED)
###         Create population
            network.create_population(
                name=name, POP_SIZE=size, **populationParameters[ii])
            np.random.seed(GLOBALSEED)
###         Create initial stimuli
            if hypothesis == 'H2':
                H2(name,network,GLOBALSEED,k)
            else:
                H1(name,network,GLOBALSEED,k)
            np.random.seed(GLOBALSEED)
###     create connectivity matrices and connect populations:
        conn(network,GLOBALSEED,loop_n=k)
###     set up extracellular recording device:
        electrode = RecExtElectrode(**electrodeParameters)
###     run simulation:
        np.random.seed(GLOBALSEED)
        if k == 0:
            SPIKES0, OUTPUT0, DIPOLEMOMENT0 = network.simulate(
                electrode=electrode,
                **networkSimulationArguments,
            )
        if k ==1:
            SPIKES1, OUTPUT1, DIPOLEMOMENT1 = network.simulate(
                electrode=electrode,
                **networkSimulationArguments,
            )
###     store soma potentials

        # for cell in network.populations['E'].cells:
        if k ==0:
            # zz=0
            for cell in network.populations['E'].cells:
                SomaVs['inhib'].append(cell.somav)
                # MUA['inhib'].append(cell.vmem)
                # zz +=1
                # print(zz,np.shape(cell.vmem))

        if k ==1:
            for cell in network.populations['E'].cells:
                SomaVs['blocked'].append(cell.somav)
                # MUA['inhib'].append(cell.vmem)

###     clean up nrn and LFPy constructs
        clean_up(network)
###########################################################################
# Create plots
###########################################################################
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.imshow(MUA['inhib'][:][0][:])
    # plt.show()
    laminar_LFP(OUTPUT0,OUTPUT1,hypothesis)
    # Population()
    make_plot(OUTPUT0,OUTPUT1,SPIKES0,SPIKES1,hypothesis,save_ecog=True)
    plot_soma(SomaVs,SomaTime,np.arange(0,5))
    import matplotlib.pyplot as plt
    figg = plt.figure()
    ax4 = MUA_fig(OUTPUT0, OUTPUT1, figg, 'e')
    ax4[1].set_position([0.7, 0.1, 0.2, 0.3])
    ax4[1].set_title('MUA inhib off')
    ax4[0].set_position([0.7, 0.6, 0.2, 0.3])
    ax4[0].set_title('MUA inhib on')
    # ax1 = MUA_fig(OUTPUT0, OUTPUT1,figg,'a')
    # ax2 = MUA_fig(OUTPUT0, OUTPUT1,figg,'b')
    plt.savefig('plots/MUA_{}.pdf'.format(hypothesis))

    MUA_lines(OUTPUT0,OUTPUT1,hypothesis)

    testfunk = np.sin(10*np.linspace(0,5,100))+np.sin(1000*np.linspace(0,5,100))+np.sin(100000*np.linspace(0,5,100))+np.sin(70*np.linspace(0,5,100))
    plt.figure()
    from scipy import signal

    fs = 8000
    # fc=300
    fc = 500
    w = fc / (fs / 2.)
    # # b,a = signal.butter(5,w,'lowpass')
    b, a = signal.butter(5, w, 'highpass')
    #
    laminar_LFP = signal.filtfilt(b, a, testfunk)
    #
    laminar_LFP = np.abs(laminar_LFP)

    fc = 2000
    w = fc / (fs / 2.)

    b, a = signal.butter(5, w, 'lowpass')

    laminar_LFP = signal.filtfilt(b, a, laminar_LFP)

    fc = 300
    w = fc / (fs / 2.)

    b, a = signal.butter(5, w, 'lowpass')
    lower=signal.filtfilt(b,a,testfunk)
    plt.plot(testfunk,'k')
    plt.plot(laminar_LFP,'r')
    plt.plot(lower,'b')
    plt.savefig('plots/Filter_test.pdf')


total_time = time.time() - start
print("total runtime:", total_time, "s.")
