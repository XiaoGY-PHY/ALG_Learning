import matplotlib.pyplot as plt
import numpy as np
import copy


def readFile2List(str):
    '''Read data from sample file'''
    file = open(str, 'r+')
    line = file.readline()
    ListX = []
    ListY = []
    ListSigmaY = []
    Nevt = 0
    while line:
        a = line.split()
        if(Nevt % 1 == 0):
            ListX.append(float(a[1]))
            ListY.append(float(a[2]))
            ListSigmaY.append(float(a[3]))
        line = file.readline()
        Nevt = Nevt+1
        # if Nevt == 10:
        #     break

    file.close()
    return Nevt, ListX, ListY, ListSigmaY


def GaussFunction(x, mean, sigma):
    y = (1/(np.sqrt(2*np.pi)*sigma)) * \
        np.exp(-((x-mean)*(x-mean))/(2*(sigma)*(sigma)))
    return y


def Polynomial(x, a=0, b=0, c=0, d=0, e=0, f=0, g=0):
    y = a+b*x+c*x*x+d*x*x*x+e*x*x*x*x+f*x*x*x*x*x+g*x*x*x*x*x*x
    return y


def LinePassPoint(x, a=0, x0=1.96, y0=7.14):
    y = a*(x-x0)+y0
    return y


def PDF(x, Alpha_A):
    PDF_A = 0.0014*((Alpha_A[0] * GaussFunction(x, Alpha_A[2], Alpha_A[3])) + (Alpha_A[1] * LinePassPoint(
        x, Alpha_A[4])))
    return PDF_A


def Amatrix(x, Alpha_A):
    #print('---------------------------Calculate A----------------------------')
    NP = Alpha_A.shape[0]
    NX = x.shape[0]
    AA = np.zeros(shape=(NP, NX))
    h = 1e-6
    for i in range(NP):
        d_alphap = copy.deepcopy(Alpha_A)
        d_alpham = copy.deepcopy(Alpha_A)
        d_alphap[i] = d_alphap[i] + h
        d_alpham[i] = d_alpham[i] - h
        dF = (PDF(x, d_alphap)-PDF(x, d_alpham))/(2*h)
        AA[i] = dF.T
    return copy.deepcopy(AA.T)


if __name__ == '__main__':

    datafile = '/mnt/f/ALG_Learning/HIST_LSF/sample/dataSample.dat'
    Nevt, dataX, dataY, SigmaY = readFile2List(datafile)
    dataX = np.array(dataX).reshape(-1, 1)
    dataY = np.array(dataY).reshape(-1, 1)
    # SigmaY = np.array(SigmaY).reshape(-1, 1)
    # print(dataX)
    # print(dataX.shape[0])

    # '''Print the dataSimple to test'''
    # for idx in range(Nevt):
    #    print(idx, dataX[idx], dataY[idx], SigmaY[idx])

    # '''Plot the GaussFunction to test'''
    # x = np.linspace(1.89, 2.03, 100)
    # y = (GaussFunction(x, 1.97, 0.01))
    # result = sum(y*(2.03-1.89)/100)
    # print(result)
    # plt.plot(x, y)
    # plt.show()

    # '''Plot the Polynomial to test'''
    # x = np.linspace(1.89, 2.03, 100)
    # y = (Polynomial(x, 1))

    '''Fit Start'''
    Nevent = 1656709.0
    '''Nsig,Nbkg,meanGauss,sigmaGauss,a,b,c,d,e,f,g (11parameters)'''
    alpha = (
        np.array([[0.2*Nevent, 0.8*Nevent, 1.97, 0.01, 1]])).reshape(-1, 1)

    '''Vy'''
    vy = np.diag(SigmaY)
    Vy = np.dot(vy, vy)
    Vy_I = np.linalg.inv(Vy)

    '''chi2'''
    chi2_A = 0.

    for D in range(100):
        Alpha_A = copy.deepcopy(alpha)
        '''F(alpha)=F(alpha_A)+ A*eta'''
        F_A = PDF(dataX, Alpha_A)
        # plt.plot(dataX.T[0], F_A.T[0], label='fit function')
        # plt.show()
        # print(D)
        A = Amatrix(dataX, Alpha_A)

        '''dy'''
        dy = dataY-F_A

        '''VA'''
        # print('det(VA_I):', np.linalg.det(np.dot(np.dot(A.T, Vy_I), A)))
        VA_I = np.dot(np.dot(A.T, Vy_I), A)
        # print('VA_I', VA_I)
        VA = np.linalg.inv(VA_I)
        # print('VA:', VA.shape)

        '''eta = alpha - alpha_A'''
        eta = np.dot(np.dot(np.dot(VA, A.T), Vy_I), dy)
        alpha = eta+Alpha_A
        # print(alpha.T)

        '''chi2'''
        chi2 = np.dot(np.dot((dy-np.dot(A, eta)).T, Vy_I),
                      (dy-np.dot(A, eta)))
        if(np.abs(chi2-chi2_A) < 0.00001):
            print('\n----------------------------------------------------')
            print("D: ", D, "chi2: ", float(chi2))
            print('-----')
            print('Alpha: ')
            print("    Nsig      Nbkg       mean       sigma         k")
            print(format(alpha.T[0][0], '.4e'), format(alpha.T[0][1], '.4e'), format(
                alpha.T[0][2], '.4e'), format(alpha.T[0][3], '.4e'), format(alpha.T[0][4], '.4e'))
            print('-----')
            print('VA: ')
            for i in range(VA.shape[0]):
                for j in range(VA.shape[1]):
                    print(format(VA[i][j], '.4e'), end='  ')
                print('')

            break
        print('\n----------------------------------------------------')
        print("D: ", D, "chi2: ", float(chi2))
        print('-----')
        print("    Nsig      Nbkg       mean       sigma         k")
        print(format(alpha.T[0][0], '.4e'), format(alpha.T[0][1], '.4e'), format(
            alpha.T[0][2], '.4e'), format(alpha.T[0][3], '.4e'), format(alpha.T[0][4], '.4e'))

        chi2_A = chi2+0.
        if(chi2 < 100):
            break

    NewF = PDF(dataX, alpha)

    fig, ax = plt.subplots(1, 1, figsize=(9, 6))
    ax.errorbar(dataX.T[0], dataY.T[0], yerr=SigmaY, fmt='.k', label='data')
    ax.plot(dataX.T[0], NewF.T[0], linewidth=2, label='fit function')
    yline = 0.0014*alpha[1]*LinePassPoint(dataX.T[0], alpha[4])
    ax.plot(dataX.T[0], yline, linewidth=2, label='line function')
    ygauss = 0.0014*alpha[0]*GaussFunction(dataX.T[0], Alpha_A[2], Alpha_A[3])
    ax.plot(dataX.T[0], ygauss, linewidth=2, label='Gauss function')
    ax.legend()
    plt.savefig("LSF.png")
    plt.show()
