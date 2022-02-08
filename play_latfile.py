import matplotlib.pyplot as plt
import os

import wavio
import numpy as np
import pdb
import pandas as pd

def make_soundfile(fname, rate=2205):
    # fname = "/Users/akipnis/Desktop/Createc2_210811.100158.LAT"
    f = open(fname,"rb")
    d = f.readlines()
    header_data = [str(a).strip("\\r\\n'").strip("b'") for a in d[0:626]]

    cols = []
    hdata = []
    for h in header_data:
        if "=" in h:
            b = h.split("=")
            if len(b)<2 or len(b[1])==0 or b[0]==" " or b[0]=="":
                pass
            else:
                cols.append(b[0])
                hdata.append(b[1])

    header_data = pd.DataFrame(np.array([hdata, cols]).T, columns=["value","key"])
    #with pd.option_context('display.max_rows', None): print(header_data)
    d = [a.decode('utf-8').split('\t') for a in d[629:-1]]
    curr_trace = [float(a[3]) for a in d]

    #plot the fourier transform to maybe get an Ag111 lattice constant?
    #assuming the step is constant
    stepsize = 0.02 # angstroms/step
    distances = [stepsize*i for i in range(len(curr_trace))]

    ft = np.fft.fft(curr_trace)/len(curr_trace)     #normalize
    ft = ft[range(int(len(curr_trace)/2))]          #exclude sampling frequency
    tpCount = len(curr_trace)
    values = np.arange(int(tpCount/2))
    period = tpCount * stepsize #steps * (angstroms / step) = angstroms
    frequencies = values / period #1/angstroms

    fig, (ax1,ax2) = plt.subplots(2)

    ax2.plot(frequencies, abs(ft))
    ax2.set_ylim([0,np.max(abs(ft[100:-1]))])
    ax2.set_ylabel("FFT")
    ax2.set_xlabel("1/nm")

    ax1.plot(np.array(distances)/10, curr_trace)
    ax1.set_xlabel("(nm)")
    ax1.set_ylabel("Current")
    plt.suptitle(fname.strip(".LAT"))


    plt.show()
    pdb.set_trace()

    # plt.savefig(os.path.join(dir,fname.strip(".LAT")+"FFT.png"))
    plt.close()

    # Parameters
    rate = rate    # samples per second
    T = 3           # sample duration (seconds)
    f = 440.0       # sound frequency (Hz)
    # Compute waveform samples

    # Write the samples to a file
    newf = fname.strip(".LAT")+("_rate_%d_.wav" %(rate))
    try:
        wavio.write(newf, np.array(curr_trace), 2205, sampwidth=1)
    except:
        print("could not save file")
if __name__=='__main__':
    dir = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/Ag 2021-08-09 2p5 nm radius/LATS"
    files = os.listdir(dir)

    for f in files:
        if ".LAT" in f and ".jpeg" not in f and ".png" not in f:
            make_soundfile(os.path.join(dir,f),)
    plt.show()
    pdb.set_trace()

# plt.plot(curr_trace)
# plt.show()
