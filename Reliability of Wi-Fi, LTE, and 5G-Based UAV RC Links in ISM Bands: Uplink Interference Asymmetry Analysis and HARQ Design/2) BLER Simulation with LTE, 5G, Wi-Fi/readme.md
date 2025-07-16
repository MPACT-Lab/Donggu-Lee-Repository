## **BLER Simulation with LTE, 5G, and Wi-Fi**



### **Description of File Structure**

This subfolder contains six MATLAB scripts for simulating block-error-rate (BLER) under different wireless communication standards:

1. Wi-Fi ACK Frame  
2. Wi-Fi Data Frame  
3. LTE physical uplink control channel (PUCCH) Format1a (ACK/NACK) 
4. 5G PUCCH Format1 (ACK/NACK)
5. 5G physical downlink shared channel (PDSCH, data carrying)

All simulations assume **single-input single-output (SISO)** operation and support AWGN channel. Note that BLER of LTE PDSCH can be simulated with `1)DL Throughput Simulation with Perfect ACK NACK/HARQ_Type_III_with_IR.m` 

---

### **Instructions**

Before running the scripts, configure the following simulation parameters:

- `NFrames`: Number of frames simulated per SINR point.
- `SNRIn`: SINR range in dB scale.

After setting these variables, run the script. To visualize the BLER over the SINR range, use the following command:

```matlab
semilogy(SNRIn, BLER)
```



**More Information**

TBD (for publication or arXiv)

https://www.mathworks.com/help/wlan/ug/802-11ax-packet-error-rate-simulation-for-single-user-format.html

https://www.mathworks.com/help/lte/ug/pucch1a-ack-missed-detection-probability-conformance-test.html

https://www.mathworks.com/help/5g/ug/nr-pdsch-throughput.html

https://www.mathworks.com/help/5g/ug/nr-pucch-block-error-rate.html



**Required/Recommended Software**

MATLAB (the scripts are based on MATLAB R2024a version) including LTE, 5G, and WLAN toolboxes.
