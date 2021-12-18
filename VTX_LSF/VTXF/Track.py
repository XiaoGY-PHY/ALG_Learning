import matplotlib.pyplot as plt
import numpy as np
import copy
import Read


class track():

    ''' ----------------- initialization method -------------------- '''

    def __init__(self, mass=1, helix=np.ones(5), helixErr=np.ones((5, 5))) -> None:
        self.__mass = mass
        self.__helix = copy.deepcopy(helix)
        self.__helixErr = copy.deepcopy(helixErr)
        ''' w7 = (px, py, pz, E, x0, y0, z0).T '''
        self.__w = np.zeros(7).reshape(-1, 1)
        self.__wErr = np.zeros((7, 7))
        ''' w6 = (px, py, pz, E, x0, y0, z0).T '''
        self.__w6 = np.zeros(6).reshape(-1, 1)
        self.__w6Err = np.zeros((6, 6))
        self.__XT = np.zeros(3).reshape(1, -1)
        self.__Vx = np.zeros((3, 3))
        self.__PT = np.zeros(4).reshape(1, -1)
        self.__Vp = np.zeros((4, 4))
        self.__charge = 0
        self.__makeCharge()
        self.__Helix2W()
        self.__makeXP()
        self.__makeWErr()
        self.__makePErr()
        self.__makeXErr()

    ''' ----------------- private method -------------------- '''

    def __Helix2W(self):
        drho = self.__helix[0]
        phi0 = self.__helix[1]
        kappa = self.__helix[2]
        dz = self.__helix[3]
        Lambda = self.__helix[4]
        sinphi0 = np.sin(phi0)
        cosphi0 = np.cos(phi0)
        q = self.__charge
        self.__w[0] = -sinphi0 * q / kappa
        self.__w[1] = cosphi0 * q / kappa
        self.__w[2] = Lambda * q / kappa
        self.__w[3] = np.sqrt((1+Lambda*Lambda) /
                              (kappa*kappa) + self.__mass * self.__mass)
        self.__w[4] = drho * cosphi0
        self.__w[5] = drho * sinphi0
        self.__w[6] = dz

    def __makeCharge(self):
        if(self.__helix[2] > 0):
            self.__charge = 1
        else:
            self.__charge = -1

    def __makeXP(self):
        self.__PT[0, 0] = self.__w[0]
        self.__PT[0, 1] = self.__w[1]
        self.__PT[0, 2] = self.__w[2]
        self.__PT[0, 3] = self.__w[3]
        self.__XT[0, 0] = self.__w[4]
        self.__XT[0, 1] = self.__w[5]
        self.__XT[0, 2] = self.__w[6]

    def __makeWErr(self):
        m = np.zeros((7, 5))
        drho = self.__helix[0]
        phi0 = self.__helix[1]
        kappa = self.__helix[2]
        kappa2 = kappa * kappa
        kappa3 = kappa * kappa2
        dz = self.__helix[3]
        Lambda = self.__helix[4]
        Lambda2 = Lambda * Lambda
        e = np.sqrt((1+Lambda2) / kappa2 + self.__mass * self.__mass)
        sinphi0 = np.sin(phi0)
        cosphi0 = np.cos(phi0)
        q = self.__charge
        m[0, 1] = -cosphi0 * q / kappa
        m[0, 2] = sinphi0 * q / kappa2
        m[1, 1] = -sinphi0 * q / kappa
        m[1, 2] = -cosphi0 * q / kappa2
        m[2, 2] = -Lambda * q / kappa2
        m[2, 4] = q / kappa
        m[3, 2] = -(1+Lambda2) / (kappa3 * e)
        m[3, 4] = Lambda / (kappa2 * e)
        m[4, 0] = cosphi0
        m[4, 1] = -drho * sinphi0
        m[5, 0] = sinphi0
        m[5, 1] = drho * cosphi0
        m[6, 3] = 1
        self.__wErr = np.dot(np.dot(m, self.__helixErr), m.T)

    def __makePErr(self):
        for i in range(4):
            for j in range(4):
                self.__Vp[i, j] = self.__wErr[i, j]

    def __makeXErr(self):
        for i in range(3):
            for j in range(3):
                self.__Vx[i, j] = self.__wErr[i+4, j+4]

    ''' ----------------- public get method -------------------- '''

    def getHelix(self):
        # print(self.__helix)
        return copy.deepcopy(self.__helix)

    def getHelixErr(self):
        # print(self.__helixErr)
        return copy.deepcopy(self.__helixErr)

    def getCharge(self):
        return copy.deepcopy(self.__charge)

    def getW(self):
        return self.__w

    def getEw(self):
        return copy.deepcopy(self.__wErr)

    def getX(self):
        return copy.deepcopy(self.__XT.T)

    def getP(self):
        return copy.deepcopy(self.__PT.T)

    def getEx(self):
        return copy.deepcopy(self.__Vx)

    def getEp(self):
        return copy.deepcopy(self.__Vp)

    def getMass(self):
        return copy.deepcopy(self.__mass)

    ''' ----------------- public set method -------------------- '''

    def setSomething(self) -> None:
        pass

    def setW(self, w: np.array):
        self.__w = copy.deepcopy(w)
        self.__makeXP()

    def setEw(self, Ew: np.array):
        self.__wErr = copy.deepcopy(Ew)


if __name__ == '__main__':
    mpi = 0.13957
    ptrackFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackP.dat'
    ptrackErrFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackEP.dat'

    pa = Read.loadTrack(ptrackFile, 0)
    pEa = Read.loadTrackErr(ptrackErrFile, 0)
    ptrk = track(mpi, pa, pEa)
    print(ptrk.getW())
    # print(ptrk.getX())
    # print(ptrk.getP())
    # print(ptrk.getEw())
    # print(ptrk.getEp())
    # print(ptrk.getEx())
