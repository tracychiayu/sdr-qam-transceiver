import numpy as np
import matplotlib.pyplot as plt
from ece230b import *
from remoteRF.drivers.adalm_pluto import * 
from scipy.signal import correlate

# transmit gain: -25dB
# carrier frequency: 905~925MHz
# ---------------------------------------------------------------
# Digital communication system parameters.
# ---------------------------------------------------------------
fs = 1e6 # baseband sampling rate (samples per second)
ts = 1 / fs # baseband sampling period (seconds per sample)
sps = 10 # samples per data symbol
T = ts * sps # time between data symbols (seconds per symbol)
NUM_SYMBOLS_TO_SCAN = 40000

# ---------------------------------------------------------------
# Pluto system parameters.
# ---------------------------------------------------------------
sample_rate = fs # sampling rate, between ~600e3 and 61e6
tx_carrier_freq_Hz = 915e6 # transmit carrier frequency, between 325 MHz to 3.8 GHz
rx_carrier_freq_Hz = 915e6 # receive carrier frequency, between 325 MHz to 3.8 GHz
tx_rf_bw_Hz = sample_rate * 1 # transmitter's RF bandwidth, between 200 kHz and 56 MHz
rx_rf_bw_Hz = sample_rate * 1 # receiver's RF bandwidth, between 200 kHz and 56 MHz
tx_gain_dB = -25 # transmit gain (in dB), from -50dB to -25dB
rx_gain_dB = 40 # receive gain (in dB), beteween 0 to 74.5 dB (only set if AGC is 'manual')
rx_agc_mode = 'manual' # receiver's AGC mode: 'manual', 'slow_attack', or 'fast_attack'
rx_buffer_size = 500e3 # receiver's buffer size (in samples), length of data returned by sdr.rx()
tx_cyclic_buffer = True # cyclic nature of transmitter's buffer (True -> continuously repeat transmission)
# ---------------------------------------------------------------
# Initialize Pluto object using issued token.  
# ---------------------------------------------------------------
sdr1 = adi.Pluto(token='f4mJ-J4HLws') # create Pluto object (SDR 3)
sdr1.sample_rate = int(sample_rate) # set baseband sampling rate of Pluto
# ---------------------------------------------------------------
# Setup Pluto's transmitter.
# ---------------------------------------------------------------
sdr1.tx_destroy_buffer() # reset transmit data buffer to be safe
sdr1.tx_rf_bandwidth = int(tx_rf_bw_Hz) # set transmitter RF bandwidth
sdr1.tx_lo = int(tx_carrier_freq_Hz) # set carrier frequency for transmission
sdr1.tx_hardwaregain_chan0 = tx_gain_dB # set the transmit gain
sdr1.tx_cyclic_buffer = tx_cyclic_buffer # set the cyclic nature of the transmit buffer
# ---------------------------------------------------------------
# Setup Pluto's receiver.
# ---------------------------------------------------------------
sdr2 = adi.Pluto(token='bBvnXQSB3XQ') # create Pluto object (SDR 4)
sdr2.sample_rate = int(sample_rate) # set baseband sampling rate of Pluto
sdr2.rx_destroy_buffer() # reset receive data buffer to be safe
sdr2.rx_lo = int(rx_carrier_freq_Hz) # set carrier frequency for reception
sdr2.rx_rf_bandwidth = int(sample_rate) # set receiver RF bandwidth
sdr2.rx_buffer_size = int(rx_buffer_size) # set buffer size of receiver
sdr2.gain_control_mode_chan0 = rx_agc_mode # set gain control mode
sdr2.rx_hardwaregain_chan0 = rx_gain_dB # set gain of receiver

# ---------------------------------------------------------------
# Create transmit signal.
# ---------------------------------------------------------------

# Part I: Generate 1000 random 64-QAM symbols
N = 1000
M = 64
symbols, constellation = gen_rand_qam_symbols(N, M)
sps = 10   # samples per symbol
np.save('tx_symbols.npy', symbols)

zc_len_short = 31      # Length of one short ZC sequence
zc_len_long = 937      # Length of one long ZC sequence
zc_count_short = 16    # Number of short ZC sequences in STF
zc_count_long = 2      # Number of long ZC sequences in LTF
N_train = zc_len_short * zc_count_short + zc_len_long * zc_count_long     # (length: 2370)
N_pilots = 50
N_zeros = int(N/10)
N_frame = N_train + N_pilots + N + N_zeros

# Generate training sequence for frame sync and CFO correction
zc_31 = zadoff_chu_seq(q=1, N_zc=zc_len_short)
zc_937 = zadoff_chu_seq(q=1, N_zc=zc_len_long)
STF = np.tile(zc_31, zc_count_short)  # short training field (length: 496)
LTF = np.tile(zc_937, zc_count_long)  # long training field (length: 1874)

# Generate pilots for channel estimation
pilots = np.ones(N_pilots) * constellation[0]
zero_padding = np.zeros(N_zeros)  # insert 100 zeros after sequence of symbols

tx_symbols = np.concatenate([STF, LTF, pilots, symbols, zero_padding])
pulse_train = create_pulse_train(tx_symbols, sps)


# Pulse shaping with RRC
# T = 1e-5    # Rs = 1e5 [symbols/sec]
beta = 0.5
span = 10
pulse = get_rrc_pulse(beta, span, sps)
tx_signal = np.convolve(pulse_train, pulse, mode='same')

np.save('tx_signal.npy', tx_signal)

# Show the transmit signal
plt.figure(figsize=(16,4))
plt.plot(np.arange(len(tx_signal)), np.real(tx_signal), zorder=1)
plt.scatter(np.arange(0,len(tx_signal),sps), np.real(tx_signal[::sps]), s=2, color='red', zorder=2)
for i in range(1,zc_count_short):
    plt.axvline(x=zc_len_short*i*sps, linestyle='--', linewidth=1, color='gray')
plt.axvline(x=zc_len_short*zc_count_short*sps, linestyle='--', linewidth=1, color='black')
plt.axvline(x=(zc_len_short*zc_count_short+zc_len_long)*sps, linestyle='--', linewidth=1, color='gray')
plt.axvline(x=N_train*sps, linestyle='--', linewidth=1, color='black')
plt.axvline(x=(N_train+N_pilots)*sps, linestyle='--', linewidth=1, color='black')
plt.axvline(x=(N_train+N_pilots+N)*sps, linestyle='--', linewidth=1, color='black')

plt.axvspan(0, zc_len_short*zc_count_short*sps, color='purple', alpha=0.1, label='Short Training Field (STF)')
plt.axvspan(zc_len_short*zc_count_short*sps, N_train*sps, color='blue', alpha=0.1, label='Short Training Field (STF)')
plt.axvspan(N_train*sps, (N_train+N_pilots)*sps, color='red', alpha=0.3, label='Pilot') 
plt.axvspan((N_train+N_pilots)*sps, (N_train+N_pilots+N)*sps, color='gray', alpha=0.1, label='Payload')  
plt.axvspan((N_train+N_pilots+N)*sps, (N_train+N_pilots+N+N_zeros)*sps, color='green', alpha=0.1, label='Zero Padding')
plt.title(f'Transmit Signal (real)', fontsize=20)
plt.xlabel('Sample Index', fontsize=16)
plt.ylabel('Amplitude', fontsize=16)
plt.xlim(0,(N_train+N_pilots+N+N_zeros)*sps)
plt.ylim(-3,3)
plt.xticks(fontsize=12)
plt.legend(loc='lower left', bbox_to_anchor=(-0.02, -1), prop={'size': 14})
plt.grid(True)

plt.show()
# ---------------------------------------------------------------
# Transmit from Pluto!
# ---------------------------------------------------------------
tx_signal_scaled = tx_signal / np.max(np.abs(tx_signal)) * 2**14 # Pluto expects TX samples to be between -2^14 and 2^14 
sdr1.tx(tx_signal_scaled) # will continuously transmit when cyclic buffer set to True

# ---------------------------------------------------------------
# Receive with Pluto!
# ---------------------------------------------------------------
sdr2.rx_destroy_buffer() # reset receive data buffer to be safe
for i in range(1): # clear buffer to be safe
    rx_data_ = sdr2.rx() # toss them out
    
rx_signal = sdr2.rx() # capture raw samples from Pluto
np.save('rx_signal.npy', rx_signal)

# ---------------------------------------------------------------
# Clear buffers to stop transmitting.
# ---------------------------------------------------------------
sdr1.tx_destroy_buffer()                   # reset transmit data buffer to be safe
sdr1.rx_destroy_buffer()                   # reset receive data buffer to be safe
sdr2.tx_destroy_buffer()                   # reset transmit data buffer to be safe
sdr2.rx_destroy_buffer()                   # reset receive data buffer to be safe

# Apply matched filter
filtered_rx_signal = np.convolve(rx_signal, pulse, mode='same')
t_rx = np.arange(len(filtered_rx_signal)) / fs   

# ---------------------------------------------------------------
# 1. Symbol Synchronization
# ---------------------------------------------------------------
# a. Find symbol offset using MOE
offset = np.arange(sps)
E_output = []   # output energy
for tau in offset:
    samples = filtered_rx_signal[tau::sps]
    energy = np.mean(np.abs(samples)**2)
    E_output.append(energy)

sample_offset = np.argmax(E_output)
print(f'Sample offset: {sample_offset}')

# b. Use the estimated sample offset to extract and downsample symbol-rate samples
rx_symbols = filtered_rx_signal[sample_offset::sps]
np.save('rx_symbols.npy', rx_symbols)

# Plot the received symbols
plt.figure()
markerline, stemlines, baseline = plt.stem(np.arange(len(rx_symbols)), rx_symbols)
plt.setp(baseline, visible=False)
plt.setp(markerline, markersize=3)
plt.title('Received Symbols', fontsize=14)
plt.grid(True)
plt.show()

# ---------------------------------------------------------------
# 2. Frame Synchronization: find the starting index of LTF and extract one frame of received symbols
# ---------------------------------------------------------------
d = estimate_frame_start(rx_symbols[:NUM_SYMBOLS_TO_SCAN], zc_len_long)
print(f'start index of LTF, d = {d}')
start_index = d - zc_len_short * zc_count_short
end_index = start_index + N_frame
one_frame_symbols = rx_symbols[start_index:end_index]

# Plot one frame of symbols
plt.figure(figsize=(12,3))
markerline, stemlines, baseline = plt.stem(np.arange(len(one_frame_symbols)), np.real(one_frame_symbols))
plt.setp(baseline, visible=False)
plt.setp(markerline, markersize=3)
for i in range(1,16):
    plt.axvline(x=zc_len_short*i, linestyle='--', linewidth=1, color='gray')
plt.axvline(x=zc_len_short*16, linestyle='--', linewidth=1, color='black')
plt.axvline(x=(zc_len_short*16+zc_len_long), linestyle='--', linewidth=1, color='gray')
plt.axvline(x=N_train, linestyle='--', linewidth=1, color='black')
plt.axvline(x=(N_train+N_pilots), linestyle='--', linewidth=1, color='black')
plt.axvline(x=(N_train+N_pilots+N), linestyle='--', linewidth=1, color='black')

plt.axvspan(0, zc_len_short*16, color='purple', alpha=0.1, label='Short Training Field (STF)')
plt.axvspan(zc_len_short*16, N_train, color='blue', alpha=0.1, label='Short Training Field (STF)')
plt.axvspan(N_train, (N_train+N_pilots), color='red', alpha=0.3, label='Pilot') 
plt.axvspan((N_train+N_pilots), (N_train+N_pilots+N), color='gray', alpha=0.1, label='Payload')  
plt.axvspan((N_train+N_pilots+N), (N_train+N_pilots+N+N_zeros), color='green', alpha=0.1, label='Zero Padding')
plt.title('One Frame of Received Symbols (real)', fontsize=20)
plt.xlabel('Sample Index', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(0,N_frame)
plt.grid(True)
plt.show()

# ---------------------------------------------------------------
# 3. Frequency Synchronization: coarse and fine CFO correction
# ---------------------------------------------------------------
# a. Coarse CFO correction
rx_STF = one_frame_symbols[:zc_len_short*zc_count_short]
coarse_CFO = CFO_estimation(rx_STF, zc_len_short, zc_count_short, T)

# b. Perform coarse CFO correction on the entire received signal & extract LTF symbols
rx_signal_coarse_CFO = filtered_rx_signal * np.exp(-1j * 2 * np.pi * coarse_CFO * t_rx)  
rx_symbols_coarse_CFO = rx_signal_coarse_CFO[sample_offset::sps]
one_frame_symbols_coarse_CFO = rx_symbols_coarse_CFO[start_index:end_index]
LTF_coarse_CFO = one_frame_symbols_coarse_CFO[zc_len_short*zc_count_short:N_train]

# c. Fine CFO correction
fine_CFO = CFO_estimation(LTF_coarse_CFO, zc_len_long, zc_count_long, T)

print(f'coarse CFO: {coarse_CFO} Hz')
print(f'fine CFO: {fine_CFO} Hz')

# d. Perform coarse+fine CFO correction on the received symbols and extract the payload
rx_signal_CFO = filtered_rx_signal * np.exp(-1j * 2 * np.pi * (coarse_CFO+fine_CFO) * t_rx) 
rx_symbols_CFO = rx_signal_CFO[sample_offset::sps]
one_frame_symbols_CFO = rx_symbols_CFO[start_index:end_index]
rx_payload_symbols = one_frame_symbols_CFO[N_train+N_pilots:N_train+N_pilots+N]

# Show pilot symbols after CFO correction
rx_pilots = one_frame_symbols_CFO[N_train:N_train+N_pilots]

plt.figure()
plt.suptitle('Received Pilot Symbols After CFO Correction', fontsize=20)
plt.subplot(2,1,1)
markerline, stemlines, baseline = plt.stem(np.arange(N_pilots), np.real(pilots))
plt.setp(baseline, visible=False)
plt.setp(markerline, markersize=4)
plt.title('Real Component', fontsize=18)
plt.xlabel('Sample Index', fontsize=14)
plt.ylabel('Amplitude', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(0,N_pilots)
plt.grid(True)

plt.subplot(2,1,2)
markerline, stemlines, baseline = plt.stem(np.arange(N_pilots), np.imag(pilots))
plt.setp(baseline, visible=False)
plt.setp(markerline, markersize=4)
plt.title('Imaginary Component', fontsize=18)
plt.xlabel('Sample Index', fontsize=14)
plt.ylabel('Amplitude', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(0,N_pilots)
plt.grid(True)
plt.tight_layout()
plt.show()

# ---------------------------------------------------------------
# Channel Equalization
# ---------------------------------------------------------------
rx_pilot = np.mean(rx_pilots) 
tx_pilot = constellation[0]
estimated_h = channel_estimation(tx_pilot=tx_pilot, rx_pilot=rx_pilot)

print(f"Transmitted pilot symbol: {tx_pilot}")
print(f"Received pilot symbol: {rx_pilot}")
print(f"Estimated channel h: {estimated_h}")

rx_symbols_equalized = one_frame_symbols_CFO / estimated_h
payload_symbols_equalized = rx_symbols_equalized[N_train+N_pilots:N_train+N_pilots+N]

# Calculate SER
detected_symbols = qam_symbols_detection(payload_symbols_equalized, M)
SER = np.sum(detected_symbols != symbols) / len(symbols)
print(f"Symbol Error Rate (SER): {SER:.4f}")

np.save('detected_symbols.npy', detected_symbols)

# Visualize received symbols on complex plane
plt.figure()
plt.scatter(np.real(payload_symbols_equalized), np.imag(payload_symbols_equalized), s=2, color='red', label='Rx Symbols')
plt.scatter(np.real(constellation), np.imag(constellation), s=12, color='black', label='Tx Symbols')
# plt.title(f'Received Symbols After Equalization -- loopback SDR (SER = {SER})', fontsize=16)
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