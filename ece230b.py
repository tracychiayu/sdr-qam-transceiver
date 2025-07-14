import numpy as np
import matplotlib.pyplot as plt

# A. Generate Random QAM Symbols
def gen_rand_qam_symbols(N, M=4):
    '''
    Input
        - N: number of random symbols to generate
        - M: order of the QAM constellation
    Output
        - symbols: N randomly selected M-QAM symbols
        - constellation: the full M-QAM constellation
    '''

    A = np.sqrt(1.5/(M-1))
    d = int(np.sqrt(M))     
    constellation = np.zeros((d,d), dtype=complex)
    for im in range(0,d):
        for re in range (0,d):
            a = A * (2 * re + 1 - d)
            b = A * (2 * im + 1 - d)
            constellation[d-1-im][re] = a + b * 1j
    
    constellation = constellation.flatten()  # (sqrt(M), sqrt(M)) -> (M,)
    symbol = np.random.choice(constellation, size=N)
    return symbol, constellation

def qam_symbols_detection(symbols, M):
    '''
    Input
        - symbols: array of received symbols (complex)
        - M: order of the QAM constellation
    Output
        - detected_symbols: array of detected symbols
    '''
    _, constellation=gen_rand_qam_symbols(N=1, M=M)

    detected_symbols = np.zeros_like(symbols, dtype=complex)

    for i,sym in enumerate(symbols):
        dist = np.abs(sym-constellation)
        nearest_idx = np.argmin(dist)
        detected_symbols[i] = constellation[nearest_idx]
    
    return detected_symbols

# B. Create Pulse Train
def create_pulse_train(symbols, sps):
    '''
    Input
        - symbols: QAM symbols
        - sps: sample per symbol
    Output
        - pulse_train: a discrete-time signal where each symbol is separated by (sps - 1) zeros
    '''
    N = len(symbols)
    pulse_train = np.zeros(N * sps, dtype=complex)

    pulse_train[::sps] = symbols

    return pulse_train

def get_rrc_pulse(beta, span, sps):
    '''
    Input
        - beta: rolloff factor (between 0 and 1)
        - span: the integer number of symbol durations spanned by the pulse, not including the symbol at t = 0
        - sps: samples per symbol
    Output
        - pulse: a root raised cosine pulse (normalized such that its peak value is unity), symmetric and centered at t = 0. 
                 The number of zero crossings should be equal to span
    '''

    t = np.linspace(-span/2, span/2, sps * span + 1)
    pulse = np.zeros_like(t)
    T = 1

    if beta == 0.0:
        pulse = np.sinc(t/T) / np.sqrt(T)
        pulse = pulse / np.max(np.abs(pulse)) # normalize to peak = 1
        return pulse

    for i in range(len(t)):
        if t[i] == 0:
            pulse[i] = (1 + beta * (4/np.pi - 1)) / T
        elif t[i] == T / (4*beta) or t[i] == -T / (4*beta):
            pulse[i] = beta * ((1+2/np.pi)*np.sin(np.pi/(4*beta)) + (1-2/np.pi)*np.cos(np.pi/(4*beta))) / (T * np.sqrt(2))
        else:
            pulse[i] = (np.sin(np.pi*t[i]*(1-beta)/T) + 4*beta*t[i]*np.cos(np.pi*t[i]*(1+beta)/T)) / (np.pi*t[i]*(1-(4*beta*t[i]/T)**2))
    
    pulse = pulse / np.max(np.abs(pulse)) # normalize to peak =1
    return pulse

def get_rc_pulse(beta, span, sps):
    '''
    Input
        - beta: rolloff factor (between 0 and 1)
        - span: the integer number of symbol durations spanned by the pulse, not including the symbol at t = 0
        - sps: samples per symbol
    Output
        - pulse: a root raised cosine pulse (normalized such that its peak value is unity), symmetric and centered at t = 0. 
                 The number of zero crossings should be equal to span
    '''

    t = np.linspace(-span/2, span/2, sps * span + 1)
    pulse = np.zeros_like(t)
    T = 1

    if beta == 0.0:
        pulse = np.sinc(t/T) / np.sqrt(T)
        pulse = pulse / np.max(np.abs(pulse)) # normalize to peak = 1
        return pulse

    for i in range(len(t)):
        if t[i] == T / (2*beta) or t[i] == -T / (2*beta):
            pulse[i] = np.pi * np.sinc(1/(2*beta)) / (4*T)
        else:
            pulse[i] = (1/T) * np.sinc(t[i]/T) * np.cos(np.pi*beta*t[i]/T) / (1 - (2*beta*t[i]/T)**2)
    
    pulse = pulse / np.max(np.abs(pulse)) # normalize to peak =1
    return pulse

def zadoff_chu_seq(q, N_zc):
    '''
    Input
        - q: root index
        - N_zc: length of the sequence (should be prime number)
    Output
        - t_q: complex Zadoff-Chu sequence of length N_zc
    '''

    n = np.arange(N_zc)  # t_q[n]: n = 0, ... , N_zc-1
    t_q = np.exp(-1j * np.pi * q * n * (n+1) / N_zc)

    return t_q

def estimate_frame_start(rx_symbols, N_tr):
    '''
    Input:
        - rx_symbols: received symbols extracted after symbol synchronization
        - N_tr: length of the long training sequence (LTF)
    Output:
        - d: estimated starting index of the frame (LTF)
    '''

    N_short = 31
    L_short = 16
    STF_len = N_short * L_short
    corr = np.zeros(len(rx_symbols)-2*N_tr, dtype=complex)
    for i in range(len(rx_symbols)-2*N_tr):
        y1 = rx_symbols[i:i+N_tr]
        y2 = rx_symbols[i+N_tr:i+2*N_tr]
        corr[i] = np.abs(np.sum(y2 * np.conj(y1)))**2 / ((np.sum(np.abs(y1)**2))**2+1e-12)

    # Search for peak correlation after the STF region 
    # This ensures d is large enough to extract a full frame later
    d_relative = np.argmax(corr[STF_len:]) 
    d = d_relative + STF_len   # adjust index back to original scale

    plt.figure(figsize=(12,4))
    plt.plot(corr)
    plt.axvline(x=d, linestyle='--',color='r', label=f"LTF start at index: {d}")
    plt.title('Self-Correlation Output for Frame Synchronization', fontsize=20)
    plt.xlabel('Sample Index', fontsize=14)
    plt.legend(loc='best', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.show()

    return d

def CFO_estimation(training_sequence, zc_len, zc_count, T):
    '''
    Input:
        - training_sequence: the entire STF or LTF
        - zc_len: length of ZC sequence
        - zc_count: number of ZC sequences in a training field
        - T: symbol period
    Output:
        - estimated_CFO (in Hz)
    '''
    y1 = training_sequence[:zc_len*(zc_count-1)]
    y2 = training_sequence[zc_len:zc_len*zc_count]

    corr = np.sum(np.conj(y1) * y2)
    estimated_CFO = np.angle(corr) / (2 * np.pi * zc_len * T)

    return estimated_CFO

def channel_estimation(tx_pilot, rx_pilot):
    '''
    Input:
        - tx_pilot (np.ndarray or complex): transmitted pilot symbol
        - rx_pilot (np.ndarray or complex): received pilot symbol
    Output:
        - h (np.ndarray or complex): estimated channel response
    '''
    h = (np.conj(tx_pilot) * rx_pilot) / (np.conj(tx_pilot) * tx_pilot)

    return h

def custom_corr(x, N, plot=False):
    """
    Computes sliding-window normalized correlations between every
    N-sample segment and the next N-sample segment.

    Parameters:
    - x : np.ndarray
        Input signal (1D array, real or complex)
    - N : int
        Segment length for correlation
    - plot : bool
        Whether to plot correlation values

    Returns:
    - correlations : np.ndarray (complex)
        Array of scalar correlation values (one per sliding window)
    - peak_index : int
        Index (starting sample) where the maximum correlation magnitude occurs
    """
    x = np.asarray(x).flatten()
    M = len(x)
    num_windows = M - 2 * N + 1

    if num_windows < 1:
        raise ValueError("Signal too short for even one sliding N-to-N comparison.")

    correlations = np.zeros(num_windows, dtype=np.complex128)

    for k in range(num_windows):
        seg1 = x[k : k + N]
        seg2 = x[k + N : k + 2 * N]
        correlations[k] = np.vdot(seg1, seg2) / np.sqrt(N) 

    peak_index = np.argmax(np.abs(correlations))  # sample index where max magnitude occurs

    if plot:
        plt.figure(figsize=(8, 4))
        plt.plot(np.abs(correlations), label='|correlation|')
        plt.plot(peak_index, np.abs(correlations[peak_index]), 'ro', label='Peak')
        plt.title('Sliding N-to-N Correlation')
        plt.xlabel('Sliding Window Start Index')
        plt.ylabel('Correlation Magnitude')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    return correlations, peak_index
