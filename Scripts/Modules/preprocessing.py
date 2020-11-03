import numpy as np

def remove_artifacts(signal, num_window, coeffRA):
    """
    Inputs:
        - signal: Numeric 2D - array. (Num Samples x Num Channels)  
        - num_window: Number of samples per window.  

    Output:
        - clsignal: Numeric 2D - array (Num Samples x Num Channels)
        It is the clean signal.
    """

    num_samples = np.shape(signal)[0]
    avabs = np.mean(np.abs(signal),axis=1)
    mav = np.convolve(avabs, np.ones((num_window,))/num_window, mode='same')
    threshold = np.median(mav)/coeffRA

    index_artifacts = mav > threshold

    clsignal = np.copy(signal)
    clsignal[index_artifacts,:]=0

    return index_artifacts,clsignal

def referencing(signal,config):
    """
    Inputs:
        - signal: Numeric 2D - array. (Num Samples x Num Channels)  
        - config: struct
            - 'method':     Options: 'car', 'br', 'hard' (default)
                            'car' = Common Average Referencing
                            'br' = Bipolar Referencing
                            'hard' = Hardware Referencing

            - 'electrode':  Numeric 2D-array. Bidimensional structure of electrode 

    Output:
        - refsignal: Numeric 2D - array. (Num Samples x Num Channels Ref)

    Observacion:
    'Bipolar Referencing: La señal se resta hacia abajo. El ultimo contacto está a Hardware.'
    """

    if config.get('method') == None:
        method = 'hard'
    else:
        method = config['method']

    if method != 'hard':
        electrode = config['electrode']
        num_chV, num_chH = np.shape(electrode)
        num_channels =  num_chH * num_chV
         
    if method == 'car':
        car = np.mean(signal,axis=1)
        refsignal = signal - np.kron(np.ones((num_channels,1)) ,car).T

    elif method == 'br':
        refsignal = np.copy(signal)

        for ii in range(num_chV-1):
            for jj in range(num_chH):
                refsignal[:,electrode[ii,jj]] = signal[:,electrode[ii,jj]] - signal[:,electrode[ii+1,jj]]

    else:
        refsignal = np.copy(signal) 
                    
    return refsignal

def fragmentation(signal,indexs,min_len):
    """ 
    Inputs:
        - signal: Numeric array 2D (num_samples x num_channels)
        - indexs: Numeric array 2D (num_periods x 2) [t0,t1]
        - min_len: Numeric value. 
    
    Outputs:
        - List of signals. Each signal corresponds to one period of
        - "list_index". So, length of signal is bigger than "min_length"
        - Index_t.  list_index > min_length
    """
    
    frags_signals = []
    frags_indexs = []
    boolindexs = [0] * np.shape(signal)[0]

    num_periods = len(indexs)

    for ii in range(num_periods):
        I0 = indexs[ii,0]
        If = indexs[ii,1]
        
        if If-I0 > min_len and If < np.shape(signal)[0]:
            frags_signals.append(signal[I0:If,:])
            frags_indexs.append(indexs[ii])
            #boolindexs[I0:If]=1 ver como cambiar

    return frags_signals, frags_indexs, boolindexs