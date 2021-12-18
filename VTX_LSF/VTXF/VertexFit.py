import matplotlib.pyplot as plt
import numpy as np
import copy
import Read
import Track
import TrackPool
import VtxParameter
import math


''' --------------------- only for two tracks --------------------- '''


class vtxFit(TrackPool.trackPool):

    m_TRA = np.zeros((6, 7))
    m_TRA[0, 0] = 1.0
    m_TRA[1, 1] = 1.0
    m_TRA[2, 2] = 1.0
    m_TRA[3, 4] = 1.0
    m_TRA[4, 5] = 1.0
    m_TRA[5, 6] = 1.0
    m_TRB = np.zeros((7, 6))
    m_TRB[0, 0] = 1.0
    m_TRB[1, 1] = 1.0
    m_TRB[2, 2] = 1.0
    m_TRB[4, 3] = 1.0
    m_TRB[5, 4] = 1.0
    m_TRB[6, 5] = 1.0

    ''' ----------------- initialization method -------------------- '''

    def __init__(self) -> None:
        super().__init__()
        self.__vpar_A = VtxParameter.vtxParameter()
        self.__vpar = VtxParameter.vtxParameter()
        self.__chisqcut = 1000
        self.__chisq = 99999.0

        ''' the large matrix of track parameters defined here '''
        self.__m_ltkParA = np.array([])
        self.__m_ltkPar = np.array([])
        self.__m_ltkCovParA = np.array([])

        ''' the matrix of vertex parameters defined here '''
        self.__m_xParA = np.zeros(6).reshape(-1, 1)
        self.__m_xPar = np.zeros(6).reshape(-1, 1)
        self.__m_xCovParA = np.zeros((6, 6))
        self.__m_xCovPar = np.zeros((6, 6))
        self.__m_xCovParA_I = np.zeros((6, 6))
        self.__m_xCovPar_I = np.zeros((6, 6))

        ''' the matrix of constraint equition defined here '''
        ''' H = Dda + Edx + d '''
        self.__m_E = np.array([])
        self.__m_D = np.array([])
        self.__m_ET = np.array([])
        self.__m_DT = np.array([])
        self.__m_H = np.array([])
        self.__m_VD = np.array([])
        self.__m_Vx_ET_VD = np.array([])
        self.__m_Va0_DT_VD_E_Vx = np.array([])

    ''' ----------------- public method -------------------- '''

    def addVertex(self, vpar: VtxParameter.vtxParameter):
        vx = np.zeros(3).reshape(-1, 1)
        vx = (super().getTrackA(0).getX()+super().getTrackA(1).getX())/2
        n_vpar = copy.deepcopy(vpar)
        n_vpar.setVertex(vx)
        self.setVpar(n_vpar)
        self.setVparA(n_vpar)

    def Fit(self) -> bool:

        ntrk = super().getNTrack()

        ''' the large matrix of track parameters defined here '''
        self.__m_ltkParA = np.zeros(ntrk*6).reshape(-1, 1)
        self.__m_ltkPar = np.zeros(ntrk*6).reshape(-1, 1)
        self.__m_ltkCovParA = np.zeros((ntrk*6, ntrk*6))
        for itk in range(ntrk):
            ETWT = np.dot(
                np.dot(vtxFit.m_TRA, super().getTrackA(itk).getEw()), vtxFit.m_TRA.T)
            for i in range(6):
                self.__m_ltkParA[i +
                                 (itk*6)] = self.Convert76(super().getTrackA(itk).getW())[i]
                for j in range(6):
                    self.__m_ltkCovParA[i+(itk*6), j+(itk*6)] = ETWT[i, j]
        self.__m_ltkPar = copy.deepcopy(self.__m_ltkParA)

        ''' the matrix of vertex parameters initialized here '''
        self.__m_xParA = copy.deepcopy(self.getVparA().getVertex())
        self.__m_xCovParA = copy.deepcopy(self.getVparA().getVtxErr())
        self.__m_xCovParA_I = np.linalg.inv(self.__m_xCovParA)
        self.__m_xPar = copy.deepcopy(self.__m_xParA)

        ''' the matrix of constraint equition defined here '''
        ''' H = Dda + Edx + d '''
        self.__m_E = np.zeros((2*ntrk, 3))
        self.__m_D = np.zeros((2*ntrk, 6*ntrk))
        self.__m_ET = np.zeros((3, 2*ntrk))
        self.__m_DT = np.zeros((6*ntrk, 2*ntrk))
        self.__m_H = np.zeros(2*ntrk).reshape(-1, 1)
        self.__m_VD = np.zeros((2*ntrk, 2*ntrk))
        self.__m_Vx_ET_VD = np.zeros((3, 2*ntrk))
        self.__m_Va0_DT_VD_E_Vx = np.zeros((6*ntrk, 3))

        ''' iteration loop '''
        chisq = []
        for it in range(10):
            self.updateConstraints()
            self.fitVertex()
            chisq.append(self.__chisq)
            if it > 0:
                delchisq = chisq[it]-chisq[it-1]
                if math.fabs(delchisq) < 1e-3:
                    break
        if self.__chisq >= 100:
            return False
        else:
            self.__vertexCovMatrix()
            return True

    def fitVertex(self):
        ntrk = super().getNTrack()
        '''calculate new Vx'''
        self.__m_xCovPar_I = copy.deepcopy(self.__m_xCovParA_I)
        for i in range(ntrk):
            vd = np.zeros((2, 2))
            ed = np.zeros((2, 3))
            for m in range(2):
                for n in range(2):
                    vd[m, n] = self.__m_VD[i*2+m, i*2+n]
            for m in range(2):
                for n in range(3):
                    ed[m, n] = self.__m_E[i*2+m, n]
            self.__m_xCovPar_I += np.dot(np.dot(ed.T, vd), ed)
        self.__m_xCovPar = np.linalg.inv(self.__m_xCovPar_I)

        '''calculate Vx_ET_VD and Va0_DT_VD_E_Vx'''
        for i in range(ntrk):
            ETWT = np.dot(
                np.dot(vtxFit.m_TRA, super().getTrackA(i).getEw()), vtxFit.m_TRA.T)
            vd = np.zeros((2, 2))
            ed = np.zeros((2, 3))
            kq = np.zeros((3, 2))
            eq = np.zeros((6, 3))
            dq = np.zeros((2, 6))
            '''vd(i)'''
            for m in range(2):
                for n in range(2):
                    vd[m, n] = self.__m_VD[i*2+m, i*2+n]
            '''ed(i)'''
            for m in range(2):
                for n in range(3):
                    ed[m, n] = self.__m_E[i*2+m, n]
            '''kq(i)'''
            kq = np.dot(np.dot(self.__m_xCovPar, ed.T), vd)
            for m in range(3):
                for n in range(2):
                    self.__m_Vx_ET_VD[m, i*2 + n] = kq[m, n]
            '''dq(i)'''
            for m in range(2):
                for n in range(6):
                    dq[m, n] = self.__m_D[m+i*2, n+i*6]
            '''eq(i)'''
            eq = np.dot(np.dot(-ETWT, dq.T), kq.T)
            for m in range(6):
                for n in range(3):
                    self.__m_Va0_DT_VD_E_Vx[i*6+m, n] = eq[m, n]

        '''update vertex position'''
        self.__m_xPar = copy.deepcopy(self.__m_xParA)
        for i in range(ntrk):
            vd = np.zeros((2, 2))
            ed = np.zeros((2, 3))
            kq = np.zeros((3, 2))
            hc = np.zeros(2).reshape(-1, 1)
            '''vd(i)'''
            for m in range(2):
                for n in range(2):
                    vd[m, n] = self.__m_VD[i*2+m, i*2+n]
            '''ed(i)'''
            for m in range(2):
                for n in range(3):
                    ed[m, n] = self.__m_E[i*2+m, n]
            '''kq(i)'''
            kq = np.dot(np.dot(self.__m_xCovPar, ed.T), vd)
            for m in range(3):
                for n in range(2):
                    self.__m_Vx_ET_VD[m, i*2 + n] = kq[m, n]
            '''hc(i)'''
            for m in range(2):
                hc[m] = self.__m_H[i*2+m]
            self.__m_xPar -= np.dot(kq, hc)

        self.getVpar().setVertex(self.__m_xPar)
        # print('The vertex', self.getVpar().getVertex().T)
        # print('The new vertex is:', end=' ')
        # print(self.__m_xPar.T[0])
        # print(self.__m_H.T[0])

        '''update track parameter'''
        dq0q = np.zeros(3).reshape(-1, 1)
        dq0q = np.dot(self.__m_xCovPar_I, (self.__m_xParA-self.__m_xPar))
        for i in range(ntrk):
            ETWT = np.dot(
                np.dot(vtxFit.m_TRA, super().getTrackA(i).getEw()), vtxFit.m_TRA.T)
            alpha = np.zeros(6).reshape(-1, 1)
            alpha0 = self.Convert76(super().getTrackA(i).getW())
            dq = np.zeros((2, 6))
            vd = np.zeros((2, 2))
            ed = np.zeros((2, 3))
            kq = np.zeros((3, 2))
            hc = np.zeros(2).reshape(-1, 1)
            '''dq(i)'''
            for m in range(2):
                for n in range(6):
                    dq[m, n] = self.__m_D[m+i*2, n+i*6]
            '''vd(i)'''
            for m in range(2):
                for n in range(2):
                    vd[m, n] = self.__m_VD[i*2+m, i*2+n]
            '''ed(i)'''
            for m in range(2):
                for n in range(3):
                    ed[m, n] = self.__m_E[i*2+m, n]
            '''kq(i)'''
            kq = np.dot(np.dot(self.__m_xCovPar, ed.T), vd)
            for m in range(3):
                for n in range(2):
                    self.__m_Vx_ET_VD[m, i*2 + n] = kq[m, n]
            '''hc(i)'''
            for m in range(2):
                hc[m] = self.__m_H[i*2+m]

            # print('track', i, 'check')
            # print('alpha0', alpha0.T[0])
            # print('ETWT', ETWT)
            # print('dq', dq.T)
            # print('vd', vd)
            # print('hc', hc)
            # print('kq', kq)
            # print('dq0q', dq0q)
            alpha = alpha0-np.dot(np.dot(ETWT, dq.T),
                                  (np.dot(vd, hc)-np.dot(kq.T, dq0q)))
            mass = super().getTrack(i).getMass()
            # print('-----------------new----------------------',i)
            newalpha = copy.deepcopy(self.Convert67(mass, alpha))
            super().getTrack(i).setW(newalpha)
            # print(newalpha.T)
            # print(alpha.T)
            # print(super().getTrack(i).getW().T)

        '''get chisquare value'''
        part1 = np.dot(np.dot(self.__m_H.T, self.__m_VD), self.__m_H)
        evg = np.dot(np.dot(self.__m_E.T, self.__m_VD), self.__m_H)
        part2 = np.dot(np.dot(evg.T, self.__m_xCovPar), evg)
        self.__chisq = part1-part2
        # print(self.__chisq)

        # '''get chisquare value'''
        # part1 = np.dot(np.dot(self.__m_H.T, self.__m_VD), self.__m_H)[0, 0]
        # print("part1", part1)
        # part2 = np.dot(np.dot((self.__m_xPar-self.__m_xParA).T,
        #                self.__m_xCovPar_I), (self.__m_xPar-self.__m_xParA))[0, 0]
        # print("part2", part2)
        # print('---------------')
        # print(super().getTrack(0).getW().T[0])
        # print(super().getTrack(1).getW().T[0])
        # print('')
        # self.__chisq = part1-part2
        # print('')
        # print('----------new iteration--------')
        # print('')
        # print('x = ', self.__m_xPar.T[0])

        # print('track 0 parameters:', super().getTrack(0).getW().T[0])

        # print('track 1 parameters:', super().getTrack(1).getW().T[0])

        # print('chisq = ', self.__chisq[0])
        # print('')

    def updateConstraints(self):
        ''' calculate d、E、D、VD '''
        ntrk = super().getNTrack()
        for i in range(ntrk):
            alpha = np.zeros(6).reshape(-1, 1)
            mass = 0
            e = 0
            p = np.zeros(4).reshape(-1, 1)
            x = np.zeros(3).reshape(-1, 1)  # poca of track
            vx = np.zeros(3).reshape(-1, 1)  # vertex point

            alpha = self.Convert76(super().getTrack(i).getW())
            # print('aaalpha', alpha)
            mass = super().getTrackA(i).getMass()
            e = math.sqrt(
                mass*mass+alpha[0]*alpha[0]+alpha[1]*alpha[1]+alpha[2]*alpha[2])
            p[0] = alpha[0]
            p[1] = alpha[1]
            p[2] = alpha[2]
            p[3] = e
            x[0] = alpha[3]
            x[1] = alpha[4]
            x[2] = alpha[5]
            vx = copy.deepcopy(self.__m_xPar)
            delx = vx-x

            ''' ---------- '''
            a = 0.002925 * super().getTrackA(i).getCharge()
            J = np.float(a*(delx[0]*p[0]+delx[1]*p[1])/(p[0]*p[0]+p[1]*p[1]))
            J = min(J, 1-1e-4)
            J = max(J, -1+1e-4)
            Rx = delx[0] - 2*p[0]*(delx[0]*p[0] + delx[1]
                                   * p[1])/(p[0]*p[0]+p[1]*p[1])
            Ry = delx[1] - 2*p[1]*(delx[0]*p[0] + delx[1]
                                   * p[1])/(p[0]*p[0]+p[1]*p[1])
            S = 1.0 / math.sqrt(1-J*J) / (p[0]*p[0]+p[1]*p[1])

            '''dc'''
            dc = np.zeros(2).reshape(-1, 1)
            dc[0] = delx[1] * p[0] - delx[0] * p[1] - 0.5 * \
                a * (delx[0] * delx[0] + delx[1] * delx[1])
            dc[1] = delx[2] - p[2] / a*np.arcsin(J)

            '''Ec'''
            Ec = np.zeros((2, 3))
            Ec[0, 0] = -p[1] - a * delx[0]
            Ec[0, 1] = p[0] - a * delx[1]
            Ec[0, 2] = 0
            Ec[1, 0] = -p[0] * p[2] * S
            Ec[1, 1] = -p[1] * p[2] * S
            Ec[1, 2] = 1.0
            for m in range(2):
                for n in range(3):
                    self.__m_E[i*2+m, n] = Ec[m, n]
            self.__m_ET = copy.deepcopy(self.__m_E.T)

            '''Dc'''
            Dc = np.zeros((2, 6))
            Dc[0, 0] = delx[1]
            Dc[0, 1] = -delx[0]
            Dc[0, 2] = 0
            Dc[0, 3] = p[1] + a*delx[0]
            Dc[0, 4] = -p[0] + a*delx[1]
            Dc[0, 5] = 0
            Dc[1, 0] = -p[2]*S*Rx
            Dc[1, 1] = -p[2]*S*Ry
            Dc[1, 2] = -np.arcsin(J)/a
            Dc[1, 3] = p[0]*p[2]*S
            Dc[1, 4] = p[1]*p[2]*S
            Dc[1, 5] = -1.0
            for m in range(2):
                for n in range(6):
                    self.__m_D[i*2+m, i*6+n] = Dc[m, n]
            self.__m_DT = copy.deepcopy(self.__m_D.T)

            '''VD'''
            vd = np.zeros((2, 2))
            vd = np.linalg.inv(
                np.dot(np.dot(Dc, np.dot(
                    np.dot(vtxFit.m_TRA, super().getTrackA(i).getEw()), vtxFit.m_TRA.T)), Dc.T))
            for m in range(2):
                for n in range(2):
                    self.__m_VD[i*2+m, i*2+n] = vd[m, n]
            # print('track', i)
            # print(vd)
            '''Hc'''
            hc = np.zeros(2).reshape(-1, 1)
            hc = np.dot(Dc, (self.Convert76(super().getTrackA(
                i).getW())-self.Convert76(super().getTrack(i).getW())))+np.dot(Ec, (self.getVparA().getVertex()-self.getVpar().getVertex()))+dc
            #print('x-x', (self.getVparA().getVertex()-self.getVpar().getVertex()).T)
            for n in range(2):
                self.__m_H[i*2+n] = hc[n]

    def __vertexCovMatrix(self):
        ntrk = super().getNTrack()
        '''calculate new Vx'''
        self.__m_xCovPar_I = copy.deepcopy(self.__m_xCovParA_I)
        for i in range(ntrk):
            vd = np.zeros((2, 2))
            ed = np.zeros((2, 3))
            for m in range(2):
                for n in range(2):
                    vd[m, n] = self.__m_VD[i*2+m, i*2+n]
            for m in range(2):
                for n in range(3):
                    ed[m, n] = self.__m_E[i*2+m, n]
            self.__m_xCovPar_I += np.dot(np.dot(ed.T, vd), ed)
        self.__m_xCovPar = np.linalg.inv(self.__m_xCovPar_I)

    def swim(self):
        ntrk = super().getNTrack()
        for i in range(ntrk):
            a = 0.002925
            q = super().getCharge(i)
            x0 = super().getTrack(i).getX()[0, 0]
            y0 = super().getTrack(i).getX()[1, 0]
            px0 = super().getTrack(i).getP()[0, 0]
            py0 = super().getTrack(i).getP()[1, 0]
            pz0 = super().getTrack(i).getP()[2, 0]
            e = super().getTrack(i).getP()[3, 0]
            x=self.__m_xPar[0,0]
            y=self.__m_xPar[1,0]
            z=self.__m_xPar[2,0]

            new_W = np.zeros(7).reshape(-1, 1)
            new_W[0,0]=px0-a*q*(y-y0)
            new_W[1,0]=py0+a*q*(x-x0)
            new_W[2,0]=pz0
            new_W[3,0]=e
            new_W[4,0]=x
            new_W[5,0]=y
            new_W[6,0]=z

            super().getTrack(i).setW(new_W)     

        # print('')
        # print('----------swim to vertex--------')
        # print('')
        # print('x = ', self.__m_xPar.T[0])

        # print('track 0 parameters:', super().getTrack(0).getW().T[0])

        # print('track 1 parameters:', super().getTrack(1).getW().T[0])

        # print('chisq = ', self.__chisq[0])
        # print('')       

    def Convert76(self, p: np.array):
        m = np.zeros(6).reshape(-1, 1)
        m[0] = p[0]
        m[1] = p[1]
        m[2] = p[2]
        m[3] = p[4]
        m[4] = p[5]
        m[5] = p[6]
        return copy.deepcopy(m)

    def Convert67(self, mass: float, p: np.array):
        m = np.zeros(7).reshape(-1, 1)
        m[0] = p[0]
        m[1] = p[1]
        m[2] = p[2]
        m[4] = p[3]
        m[5] = p[4]
        m[6] = p[5]
        m[3] = np.sqrt(mass*mass + p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
        return copy.deepcopy(m)

    ''' ----------------- public get method -------------------- '''

    def getVparA(self):
        return self.__vpar_A

    def getVpar(self):
        return self.__vpar

    def getP4(self, n: int):
        return copy.deepcopy(super().getTrack(n).getP())

    ''' ----------------- public set method -------------------- '''

    def setVparA(self, vparA: VtxParameter.vtxParameter):
        self.__vpar_A = copy.deepcopy(vparA)

    def setVpar(self, vpar: VtxParameter.vtxParameter):
        self.__vpar = copy.deepcopy(vpar)


if __name__ == '__main__':
    mpi = 0.13957
    '''
    ptrackFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackP.dat'
    ptrackErrFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackEP.dat'
    mtrackFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackM.dat'
    mtrackErrFile = '/mnt/f/ALG_Learning/VTX_LSF/TrackSample/TrackEM.dat'

    for i in range(100):
        pa = Read.loadTrack(ptrackFile, i)
        pEa = Read.loadTrackErr(ptrackErrFile, i)
        ptrk = Track.track(mpi, pa, pEa)
        ma = Read.loadTrack(mtrackFile, i)
        mEa = Read.loadTrackErr(mtrackErrFile, i)
        mtrk = Track.track(mpi, ma, mEa)
        print('------------------', i, '-------------------')
        print(ptrk.getPiont().T)
        print(mtrk.getPiont().T)
    '''
