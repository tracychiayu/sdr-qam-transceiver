import numpy as np
import matplotlib.pyplot as plt

def plot_tx_signal(tx_signal, config):
    """
    Plot the real part of the transmitted OFDM signal with annotations for different fields.

    Args:
        tx_signal (np.ndarray): Transmitted baseband signal
        config (dict): Dictionary containing the following keys:
            - sps: samples per symbol
            - zc_len_short: length of short ZC sequence
            - zc_count_short: number of short ZC sequences
            - zc_len_long: length of long ZC sequence
            - N_train: number of training symbols (STF + LTF)
            - N_pilots: number of pilot symbols
            - N: number of payload OFDM symbols
            - N_zeros: number of trailing zero symbols
    """
        
    sps = config['sps']
    zc_len_short = config['zc_len_short']
    zc_count_short = config['zc_count_short']
    zc_len_long = config['zc_len_long']
    N_train = config['N_train']
    N_pilots = config['N_pilots']
    N = config['N']
    N_zeros = config['N_zeros']

    plt.figure(figsize=(15, 8))
    plt.subplot(2,1,1)
    plt.plot(np.arange(len(tx_signal)), np.real(tx_signal), zorder=1)
    plt.scatter(np.arange(0, len(tx_signal), sps), np.real(tx_signal[::sps]), s=2, color='red', zorder=2)

    for i in range(1, zc_count_short):
        plt.axvline(x=zc_len_short * i * sps, linestyle='--', linewidth=1, color='gray')

    plt.axvline(x=zc_len_short * zc_count_short * sps, linestyle='--', linewidth=1, color='black')
    plt.axvline(x=(zc_len_short * zc_count_short + zc_len_long) * sps, linestyle='--', linewidth=1, color='gray')
    plt.axvline(x=N_train * sps, linestyle='--', linewidth=1, color='black')
    plt.axvline(x=(N_train + N_pilots) * sps, linestyle='--', linewidth=1, color='black')
    plt.axvline(x=(N_train + N_pilots + N) * sps, linestyle='--', linewidth=1, color='black')

    plt.axvspan(0, zc_len_short * zc_count_short * sps, color='purple', alpha=0.1, label='Short Training Field (STF)')
    plt.axvspan(zc_len_short * zc_count_short * sps, N_train * sps, color='blue', alpha=0.1, label='Long Training Field (LTF)')
    plt.axvspan(N_train * sps, (N_train + N_pilots) * sps, color='red', alpha=0.3, label='Pilot') 
    plt.axvspan((N_train + N_pilots) * sps, (N_train + N_pilots + N) * sps, color='gray', alpha=0.1, label='Payload')  
    plt.axvspan((N_train + N_pilots + N) * sps, (N_train + N_pilots + N + N_zeros) * sps, color='green', alpha=0.1, label='Zero Padding')

    plt.title('Transmit Signal (real)', fontsize=20)
    plt.xlabel('Sample Index', fontsize=16)
    plt.ylabel('Amplitude', fontsize=16)
    plt.xlim(0, (N_train + N_pilots + N + N_zeros) * sps)
    plt.ylim(-3, 3)
    plt.xticks(fontsize=12)
    plt.legend(loc='upper left', bbox_to_anchor=(1,1), prop={'size': 14})
    plt.grid(True)

    plt.subplot(2,1,2)
    plt.plot(np.arange(len(tx_signal)), np.imag(tx_signal), zorder=1)
    plt.scatter(np.arange(0, len(tx_signal), sps), np.imag(tx_signal[::sps]), s=2, color='red', zorder=2)

    for i in range(1, zc_count_short):
        plt.axvline(x=zc_len_short * i * sps, linestyle='--', linewidth=1, color='gray')

    plt.axvline(x=zc_len_short * zc_count_short * sps, linestyle='--', linewidth=1, color='black')
    plt.axvline(x=(zc_len_short * zc_count_short + zc_len_long) * sps, linestyle='--', linewidth=1, color='gray')
    plt.axvline(x=N_train * sps, linestyle='--', linewidth=1, color='black')
    plt.axvline(x=(N_train + N_pilots) * sps, linestyle='--', linewidth=1, color='black')
    plt.axvline(x=(N_train + N_pilots + N) * sps, linestyle='--', linewidth=1, color='black')

    plt.axvspan(0, zc_len_short * zc_count_short * sps, color='purple', alpha=0.1)
    plt.axvspan(zc_len_short * zc_count_short * sps, N_train * sps, color='blue', alpha=0.1)
    plt.axvspan(N_train * sps, (N_train + N_pilots) * sps, color='red', alpha=0.3) 
    plt.axvspan((N_train + N_pilots) * sps, (N_train + N_pilots + N) * sps, color='gray', alpha=0.1)  
    plt.axvspan((N_train + N_pilots + N) * sps, (N_train + N_pilots + N + N_zeros) * sps, color='green', alpha=0.1)

    plt.title('Transmit Signal (real)', fontsize=20)
    plt.xlabel('Sample Index', fontsize=16)
    plt.ylabel('Amplitude', fontsize=16)
    plt.xlim(0, (N_train + N_pilots + N + N_zeros) * sps)
    plt.ylim(-3, 3)
    plt.xticks(fontsize=12)
    plt.grid(True)    

    plt.tight_layout()
    plt.show()

def stem_plot(x, title):
    """
    Plot a stem graph of a 1D array.

    Args:
        x (np.ndarray): 1D array values
        title (str): Title of the plot
    """
    
    plt.figure()
    markerline, stemlines, baseline = plt.stem(np.arange(len(x)), x)
    plt.setp(baseline, visible=False)
    plt.setp(markerline, markersize=3)
    plt.title(title, fontsize=14)
    plt.grid(True)
    plt.show()

def plot_constellation(payload_symbols_equalized, constellation, SER):
    """
    Plot the received symbols on the complex plane along with reference constellation points.

    Args:
        payload_symbols_equalized (np.ndarray): Equalized received symbols
        constellation (np.ndarray): Reference constellation points, e.g., for 16-QAM
        SER (float): Symbol Error Rate, displayed in the plot title
    """

    plt.figure()
    plt.scatter(np.real(payload_symbols_equalized), np.imag(payload_symbols_equalized), s=2, color='red', label='Rx Symbols')
    plt.scatter(np.real(constellation), np.imag(constellation), s=12, color='black', label='Tx Symbols')
    plt.title(f'Received Symbols After Equalization -- OTA SDRs (SER = {SER})', fontsize=16)
    plt.xlabel('Real', fontsize=14)
    plt.ylabel('Imag', fontsize=14)
    plt.legend(loc='upper left', fontsize=13)
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.5,1.5)
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.grid(True)
    plt.show()