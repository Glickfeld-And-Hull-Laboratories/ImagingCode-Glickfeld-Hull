import os
import scipy.io as spio
import numpy as np
import sys

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''

    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''

    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''

    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def sbxFrames(filename):
    '''
    Input: filename should be full path excluding .sbx
    '''
    # Check if contains .sbx and if so just truncate
    if '.sbx' in filename:
        filename = filename[:-4]

    # Load info
    info = loadmat(filename + '.mat')['info']
    # print info.keys()

    # Defining number of channels/size factor
    if info['channels'] == 1:
        info['nChan'] = 2;
        factor = 1
    elif info['channels'] == 2:
        info['nChan'] = 1;
        factor = 2
    elif info['channels'] == 3:
        info['nChan'] = 1;
        factor = 2

    # Determine number of frames in whole file
    max_idx = os.path.getsize(filename + '.sbx') / info['recordsPerBuffer'] / info['sz'][1] * factor / 4 - 1

    print(int(max_idx))


sbxFrames(str(sys.argv[1]))



