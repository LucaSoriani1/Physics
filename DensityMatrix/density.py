import numpy as np
import random
import time

# global pos

# def checkpoint(condition, error):
#     if condition:
#         print('ERROR: ', error)
#         print('At position: ', pos)
#         exit()
#     pos += 1

def sep(rows, colo):
    coeff = np.zeros(rows*colo, dtype=np.complex128)
    dmatrix = np.zeros((rows, colo), dtype=np.complex128)
    for ii in range(rows):
        for jj in range(colo):
            dmatrix[ii, jj] = np.complex128(random.random() + 1j*random.random())
    for jj in range(colo):
        norm = 0
        for ii in range(rows):
            norm += dmatrix[ii, jj]*np.conj(dmatrix[ii, jj])
        norm = np.sqrt(norm)
        for ii in range(rows):
            dmatrix[ii, jj] = dmatrix[ii, jj]/(norm*np.sqrt(colo))
    for ii in range(rows):
        for jj in range(colo):
            coeff[jj + (ii-1)*colo] = dmatrix[ii, jj]
    norm = 0
    for ii in range(len(coeff)):
        norm += coeff[ii]*np.conj(coeff[ii])
    norm = np.sqrt(norm)
    # checkpoint(np.abs(norm-1.0)>0.01, 'Something wrong with the norm. It is not 1.')
    return coeff

def nonsep(row, colo):
    lenght = int(colo**row)
    coeff = np.zeros(lenght, dtype=np.complex128)
    for ii in range(lenght):
        coeff[ii] = np.complex128(random.random() + 1j*random.random())
    norm = 0
    for ii in range(lenght):
        norm += coeff[ii]*np.conj(coeff[ii])
    norm = np.sqrt(norm)
    for ii in range(lenght):
        coeff[ii] = coeff[ii]/norm
    norm = 0
    for ii in range(lenght):
        norm += coeff[ii]*np.conj(coeff[ii])
    norm = np.sqrt(norm)
    # checkpoint(np.abs(norm-1.0)>0.01, 'Something wrong with the norm. It is not 1.')
    return coeff

def main():
    dim1 = 2
    dim2 = 2
    time1 = 0
    time2 = 0
    time3 = 0
    tot = 0
    print('===============================================================================')
    print('What system do you want?')
    print('1) Pure SEPARABLE system')
    print('2) Pure NON-SEPARABLE system')
    inte = int(input())
    # checkpoint(inte!=1 and inte!=2, 'You make no choice')
    if inte == 1:
        print('The default setting are:')
        print('Number of components N: ', dim1)
        print('Dimension of sub-system d: ', dim2)
        print('Change them? Type c if you want change or type anything to continue')
        cho = input()
        if cho == 'c':
            print('Number of components N:')
            dim1 = int(input())
            # checkpoint(dim1<=1, 'Please insert a positive integer for N. It must be bigger than 1.')
            print('Dimension of sub-system d:')
            dim2 = int(input())
            # checkpoint(dim2<=1, 'Please insert a positive integer for N. It must be bigger than 1.')
        time1 = time.time()
        a = sep(dim1, dim2)
        densirho = np.zeros((len(a), len(a)), dtype=np.complex128)
        for ii in range(len(a)):
            for jj in range(len(a)):
                densirho[ii, jj] = a[ii]*np.conj(a[jj])
        tracerho = 0
        for ii in range(len(a)):
            tracerho += densirho[ii, ii]
        # checkpoint(np.abs(tracerho-1.0)>0.01, 'The trace is not 1. Please restart')
        # for ii in range(len(a)):
        #     for jj in range(len(a)):
                # checkpoint(np.abs(densirho[ii, jj]-np.conj(densirho[jj, ii]))>0.01, 'The density matrix is not hermitian. Please restart')
        time2 = time.time()
        time1 = time2-time1
        print('------------------DIMENSIONI-----------------')
        print('...............N = ', dim1)
        print('...............d = ', dim2)
        print('...............SEPARABLE!!!!!!')
        print('')
        print('====================================================')
        print('')
        print('The coefficients are:')
        for ii in range(len(a)):
            print(a[ii])
        print('')
        print('====================================================')
        print('')
        print('The density matrix of dimension ', len(a), 'x', len(a), ' is :')
        for ii in range(len(densirho)):
            print(densirho[ii])
        print('')
        print('The trace of the density matrix is : ', tracerho)
        print('')
        print('The density matrix is Hermitian!!!')
        print('')
        print('====================================================')
        print('')
        print('The time spent is', time1)
        print('====================================================')
        print('')
        print('')
        print('====================================================')
        print('====================================================')
        print('====================================================')
        print('====================================================')
        print('')
    else:
        time1 = time.time()
        a = nonsep(dim1, dim2)
        densirho = np.zeros((len(a), len(a)), dtype=np.complex128)
        for ii in range(len(a)):
            for jj in range(len(a)):
                densirho[ii, jj] = a[ii]*np.conj(a[jj])
        tracerho = 0
        for ii in range(len(a)):
            tracerho += densirho[ii, ii]
        # checkpoint(np.abs(tracerho-1.0)>0.01, 'The trace is not 1. Please restart')
        # for ii in range(len(a)):
        #     for jj in range(len(a)):
        #         checkpoint(np.abs(densirho[ii, jj]-np.conj(densirho[jj, ii]))>0.01, 'The density matrix is not hermitian. Please restart')
        dimensione = int(dim2**(dim1-1))
        reduce = np.zeros((dimensione, dimensione), dtype=np.complex128)
        reduce = np.complex128(0)
        print('What system do you want calculate the reduced matrix? Write an integer between 1 and ', dim1)
        kappa = int(input())
        # checkpoint(kappa<=1 and kappa>=dim1, 'Please insert a right number.')
        for ii in range(dimensione):
            for jj in range(dimensione):
                for kk in range(dim2):
                    righe = 1+(np.mod((ii-1), (dim2**(kappa-1)))) + (kk-1)*(dim2**(kappa-1)) + dim2*(ii-1-np.mod((ii-1), dim2**(kappa-1))) + (kk-1)*(2-kappa)
                    colonne = 1+(np.mod((jj-1), dim2**(kappa-1))) + (kk-1)*(dim2**(kappa-1)) + dim2*(jj-1-(np.mod((jj-1), dim2**(kappa-1)))) + (kk-1)*(2-kappa)
                    reduce[ii, jj] = reduce[ii, jj]+densirho[righe, colonne]
        traceredu = 0
        for ii in range(dimensione):
            traceredu += reduce[ii, ii]
        # checkpoint(np.abs(traceredu-1.0)>0.01, 'Something wrong with the computation of the trace of the reduced matrix. Please check it.')
        # for ii in range(dimensione):
            # for jj in range(dimensione):
            #     checkpoint(np.abs(reduce[ii, jj] - np.conj(reduce[jj,ii]))>0.01, 'The density matrix is not Hermitian. WRONG CALCULATION')
        time3 = time.time()
        print('FInish!')
        tot = time3-time1
        print('------------------DIMENSIONI-----------------')
        print('...............N = ', dim1)
        print('...............d = ', dim2)
        print('...............NON-SEPARABLE!!!!!!')
        print('')
        print('====================================================')
        print('')
        print('The coefficients are:')
        for ii in range(len(a)):
            print(a[ii])
        print('')
        print('====================================================')
        print('')
        print('The density matrix of dimension ', len(a), 'x', len(a), ' is :')
        for ii in range(len(densirho)):
            print(densirho[ii])
        print('')
        print('The trace of the density matrix is : ', tracerho)
        print('')
        print('The density matrix is Hermitian!!!')
        print('')
        print('====================================================')
        print('')
        print('The reduced matrix of dimension ', dimensione, 'x', dimensione, ' of the subsystem ', kappa, ' is:')
        for ii in range(len(reduce)):
            print(reduce[ii])
        print('')
        print('The trace of the reduced density matrix is : ', traceredu)
        print('')
        print('The reduced density matrix is Hermitian!!!')
        print('')
        print('====================================================')
        print('The time spent is', tot)
        print('====================================================')
        print('')
        print('')
        print('====================================================')
        print('====================================================')
        print('====================================================')
        print('====================================================')
        print('')

if __name__ == '__main__':
    main()