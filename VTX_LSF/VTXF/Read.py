import matplotlib.pyplot as plt
import numpy as np
import copy


def loadTrackErr(str, event):
    '''Read data from sample file'''
    file = open(str, 'r+')
    line = file.readline()
    HelixError = np.zeros((5, 5))
    m = 0
    while line:
        a = line.split()
        if(int(a[0]) == event):
            for n in range(5):
                HelixError[m, n] = float(a[n+1])
            m = m+1
        if m == 5:
            break
        line = file.readline()
        # if Nevt == 10:
        #     break
    file.close()
    return HelixError


def loadTrack(str, event):
    '''Read data from sample file'''
    file = open(str, 'r+')
    line = file.readline()
    Helix = []

    Nevt = 0
    while line:
        a = line.split()
        if(Nevt == event):
            for i in range(5):
                Helix.append(float(a[i+1]))
        line = file.readline()
        Nevt = Nevt+1
        # if Nevt == 10:
        #     break

    file.close()
    return np.array(Helix).reshape(-1, 1)


if __name__ == '__main__':
    print(loadTrack('/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackP.dat', 1))
    print(loadTrackErr('/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackEP.dat', 50000))
