import copy
import matplotlib.pyplot as plt
import numpy as np


class vtxParameter():

    ''' ----------------- initialization method -------------------- '''

    def __init__(self) -> None:
        self.__vtx = np.zeros(3).reshape(-1, 1)
        self.__Evtx = np.zeros((3, 3))

    ''' ----------------- private method -------------------- '''

    ''' ----------------- public get method -------------------- '''

    def getVertex(self):
        return copy.deepcopy(self.__vtx)

    def getVtxErr(self):
        return copy.deepcopy(self.__Evtx)

    ''' ----------------- public set method -------------------- '''

    def setVertex(self, vx: np.array):
        self.__vtx = copy.deepcopy(vx.reshape(-1, 1))

    def setVtxErr(self, Evx: np.array):
        self.__Evtx = copy.deepcopy(Evx)


if __name__ == '__main__':
    vx = np.array([0.0, 0.0, 0.0]).reshape(-1, 1)
    Evx = np.zeros((3, 3))
    bx = by = bz = 1e6
    Evx[0, 0] = bx*bx
    Evx[1, 1] = by*by
    Evx[2, 2] = bz*bz

    vx1 = np.ones(3).reshape(-1,1)
    print(vx1)

    ''' set the vertex parameters '''
    vxpar = vtxParameter()
    vxpar.setVertex(vx1)
    vxpar.setVtxErr(Evx)

    print(vxpar.getVertex())
