import matplotlib.pyplot as plt
import numpy as np
import copy
import Read
import Track

''' only for two tracks '''


class trackPool():

    ''' ----------------- initialization method -------------------- '''

    def __init__(self) -> None:
        self.__trackA = []
        self.__track = []
        self.__trkIdx = []

    ''' ----------------- private method -------------------- '''

    ''' ----------------- public get method -------------------- '''

    def getTrackA(self):
        return self.__trackA

    def getTrack(self):
        return self.__track

    def getTrackA(self, num: int) -> Track.track:
        return self.__trackA[num]

    def getTrack(self, num: int) -> Track.track:
        return self.__track[num]

    def getNTrack(self) -> int:
        return len(self.__trackA)

    def getCharge(self, n: int) -> int:
        return self.__trackA[n].getCharge()

    ''' ----------------- public set method -------------------- '''

    def addTrack(self, trk: Track.track):
        self.__trackA.append(copy.deepcopy(trk))
        self.__track.append(copy.deepcopy(trk))

    def setTrackA(self, num: int, trk: Track.track):
        self.__trackA[num] = copy.deepcopy(trk)

    def setTrack(self, num: int, trk: Track.track):
        self.__track[num] = copy.deepcopy(trk)


if __name__ == '__main__':
    mpi = 0.13957

    ptrackFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackP.dat'
    ptrackErrFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackEP.dat'
    mtrackFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackM.dat'
    mtrackErrFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackEM.dat'

    for i in range(1):
        pa = Read.loadTrack(ptrackFile, i)
        pEa = Read.loadTrackErr(ptrackErrFile, i)
        ptrk = Track.track(mpi, pa, pEa)
        ma = Read.loadTrack(mtrackFile, i)
        mEa = Read.loadTrackErr(mtrackErrFile, i)
        mtrk = Track.track(mpi, ma, mEa)

    a = trackPool()
    a.addTrack(ptrk)
    a.addTrack(mtrk)

    print(a.getNTrack())
    print(a.getCharge(1))
    print(a.getTrack(0).getW().T)
    ll = np.zeros(7).reshape(-1,1)
    # print(ll)
    a.getTrack(0).setW(ll)
    print(a.getTrack(0).getW())


