# SDR QAM Transceiver
This project implements a complete digital communication system using Software Defined Radio (SDR). It transmits a QAM-modulated signal over the air and performs all essential baseband processing steps on the receiver side, including:
- Matched filtering
- Symbol synchronization
- Frame synchronization
- Coarse and fine Carrier Frequency Offset (CFO) correction
- Channel estimation and equalization
- Symbol detection and Symbol Error Rate (SER) evaluation

## File Description
- `main.py`: The main script that runs the full transmit-receive pipeline. It sets up the Pluto SDR, performs all signal processing steps, and evaluates SER.
- `ece230b.py`: Contains all utility functions used in the project.
