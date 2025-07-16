## **DL Throughput Simulation with Perfect ACK/NACK**



### **Description of File Structure**

This subfolder contains four MATLAB scripts for simulating downlink (DL) throughput under different hybrid automatic repeat request (HARQ) options:

1. HARQ Type-III with Incremental Redundancy (IR)  
2. HARQ Type-I without Chase Combining (CC)  
3. HARQ Type-I with CC  
4. Burst Transmission with CC (i.e., four consecutive subframes of 1 ms each)

All simulations assume **single-input single-output (SISO)** operation and support both AWGN and Rayleigh fading channels. The scripts focus solely on DL throughput under the assumption of **perfect ACK/NACK feedback** on the uplink (UL).

---

### **Instructions**

Before running the scripts, configure the following simulation parameters:

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



**Required/Recommended Software**

MATLAB (the scripts are based on MATLAB R2024a version) including LTE, 5G, and WLAN toolboxes.
