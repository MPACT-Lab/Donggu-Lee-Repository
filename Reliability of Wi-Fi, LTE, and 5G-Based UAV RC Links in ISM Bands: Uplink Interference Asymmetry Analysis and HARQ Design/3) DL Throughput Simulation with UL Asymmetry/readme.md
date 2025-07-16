## **DL Throughput Simulation with UL Asymmetry**



### **Description of File Structure**

This subfolder contains two MATLAB scripts for simulating DL throughput under UL SINR asymmetry scenarios:

1. LTE PDSCH/PUCCH  
2. 5G PDSCH/PUCCH  

All simulations assume **single-input single-output (SISO)** operation and support both AWGN and Rayleigh fading channels. 

---

### **Instructions**

Before running the scripts, configure the following simulation parameters:

- `dB_div`: Indicates asymmetric SINR in UL in dB scale. For instance, -5 dB asymmetry at 10 dB DL SINR leads to 5 dB SINR in the UL.

- `is_awgn`: Indicates whether AWGN channel is used. Set to `1` for AWGN. Otherwise, a Rayleigh fading channel with EVA delay profile and 5 Hz maximum Doppler frequency is used.

- `NFrames`: Number of frames simulated per SINR point.
- `SNRIn`: SINR range in dB scale.
- `coding_rate_arr`: Defines the coding rate(s) used for the simulation.
- `simulationParameters.PDSCH.Modulation`: Defines the modulation scheme (e.g., `'QPSK'`, `'16QAM'`).

After setting these variables, run the script. To visualize the throughput in Mbps over the SINR range, use the following command:

```matlab
plot(SNRIn, (throughput_arr/(NFrames*10*10^-3))/10^6)
```



**More Information**

TBD (for publication or arXiv)

https://www.mathworks.com/help/lte/ug/pdsch-throughput-conformance-test-for-single-antenna-tm1-transmit-diversity-tm2-open-loop-tm3-and-closed-loop-tm4-6-spatial-multiplexing.html

https://www.mathworks.com/help/lte/ug/pucch1a-ack-missed-detection-probability-conformance-test.html

https://www.mathworks.com/help/5g/ug/nr-pdsch-throughput.html

https://www.mathworks.com/help/5g/ug/nr-pucch-block-error-rate.html



**Required/Recommended Software**

MATLAB (the scripts are based on MATLAB R2024a version) including LTE, 5G, and WLAN toolboxes.
