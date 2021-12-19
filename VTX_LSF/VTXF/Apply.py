import matplotlib.pyplot as plt
import numpy as np
import copy
import Read
import Track
import VtxParameter
import VertexFit
import math

np.set_printoptions(
    infstr='inf',
    nanstr='nan',
    formatter=None,
    precision=6,    # 精度，保留小数点后几位
    threshold=500,
    # 最多可现实的Array元素个数
    # 限制的是基本元素个数，如3*5的矩阵，限制的是15而非3（行）
    # 如果超过就采用缩略显示
    edgeitems=3,
    # 在缩率显示时在起始和默认显示的元素个数
    linewidth=150,  # 每行最多显示的字符数，默认80，超过则换行显示
    suppress=False   # 浮点显示（不用科学计数法）
)

font1 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 23
         }

mpi = 0.13957

if __name__ == '__main__':
    ptrackFile = '../TrackSample/TrackP.dat'
    ptrackErrFile = '../TrackSample/TrackEP.dat'
    mtrackFile = '../TrackSample/TrackM.dat'
    mtrackErrFile = '../TrackSample/TrackEM.dat'

    list_of_mass = []
    Nevent = 5000
    flag = 1
    print('------------------ begin vertex fit -------------------')
    for i in range(int(Nevent)):
        if ((i+1)/Nevent*100) % 5 == 0:
            print('------------------ fit vertex',
                  5*flag, '% -------------------')
            flag += 1
        ''' load tracks, here is pion '''
        pa = Read.loadTrack(ptrackFile, i)
        pEa = Read.loadTrackErr(ptrackErrFile, i)
        ma = Read.loadTrack(mtrackFile, i)
        mEa = Read.loadTrackErr(mtrackErrFile, i)

        ptrk = Track.track(mpi, pa, pEa)
        mtrk = Track.track(mpi, ma, mEa)

        ''' define the initial vertex, usually the orgin '''
        vx = np.array([0.0, 0.0, 0.0]).reshape(-1, 1)
        Evx = np.zeros((3, 3))
        bx = by = bz = 1e12
        Evx[0, 0] = bx*bx
        Evx[1, 1] = by*by
        Evx[2, 2] = bz*bz

        ''' set the vertex parameters '''
        vxpar = VtxParameter.vtxParameter()
        vxpar.setVertex(vx)
        vxpar.setVtxErr(Evx)

        ''' initialized vertex fit '''
        vtxfit = VertexFit.vtxFit()
        vtxfit.addTrack(ptrk)
        vtxfit.addTrack(mtrk)
        vtxfit.addVertex(vxpar)

        p4pip = np.zeros(4).reshape(-1, 1)
        p4pim = np.zeros(4).reshape(-1, 1)
        Ks0 = np.zeros(4).reshape(-1, 1)
        fitWell = vtxfit.Fit()
        if fitWell:
            vtxfit.swim()
            p4pip = copy.deepcopy(vtxfit.getP4(0))
            p4pim = copy.deepcopy(vtxfit.getP4(1))
            Ks0 = p4pip+p4pim
            mKs = math.sqrt(Ks0[3][0]*Ks0[3][0]-(Ks0[0][0]*Ks0[0][0]+Ks0[1][0]*Ks0[1][0] +
                                                 Ks0[2][0]*Ks0[2][0]))
            if mKs > 0.6:
                continue
            if mKs < 0.4:
                continue
            list_of_mass.append(mKs)

    ''' plt '''
    plt.hist(list_of_mass, bins=48, range=(0.487, 0.511), rwidth=0.85)
    aaxis = plt.axis()
    mKs_mean = np.mean(list_of_mass)
    plt.text(0.65*(aaxis[1]-aaxis[0])+aaxis[0], 0.8 *
             aaxis[3], 'mean: {0:.4f} GeV'.format(mKs_mean))
    plt.text(0.65*(aaxis[1]-aaxis[0])+aaxis[0], 0.85 *
             aaxis[3], 'events: {}'.format(len(list_of_mass)))
    plt.xlabel('mKs [GeV]')
    plt.ylabel('events / 0.5 MeV')
    plt.title('mass distribution of Ks')
    plt.grid()
    plt.savefig('mKs_Grid.png')
    plt.show()
