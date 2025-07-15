%% wifi DL data simulation
% reference example link: https://www.mathworks.com/help/wlan/ug/802-11ax-packet-error-rate-simulation-for-single-user-format.html

clear;
clc;

cfgHE = wlanHESUConfig;
cfgHE.ChannelBandwidth = 'CBW20';  
cfgHE.NumSpaceTimeStreams = 1;     
cfgHE.NumTransmitAntennas = 1;     
cfgHE.APEPLength = 1e3;            
cfgHE.ExtendedRange = false;       
cfgHE.Upper106ToneRU = false;      
cfgHE.PreHESpatialMapping = false; 
cfgHE.GuardInterval = 3.2;         
cfgHE.HELTFType = 4;               
cfgHE.ChannelCoding = 'LDPC';      
cfgHE.MCS = 1;                     

chanBW = cfgHE.ChannelBandwidth;
tgaxChannel = wlanTGaxChannel;
tgaxChannel.DelayProfile = 'Model-A';
tgaxChannel.NumTransmitAntennas = cfgHE.NumTransmitAntennas;
tgaxChannel.NumReceiveAntennas = 1;
tgaxChannel.TransmitReceiveDistance = 5; % Distance in meters for NLOS
tgaxChannel.ChannelBandwidth = chanBW;
tgaxChannel.LargeScaleFadingEffect = 'None';
tgaxChannel.NormalizeChannelOutputs = false;
fs = wlanSampleRate(cfgHE);
tgaxChannel.SampleRate = fs;

snr = linspace(0, 15, 30);

numSNR = numel(snr); 
throughput_sim = zeros(1,numSNR);

ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgHE);

ind = wlanFieldIndices(cfgHE);

num_packet = 10^5; 

error_count_sim = [];

for isnr = 1:numSNR
    packetSNR = convertSNR(snr(isnr),"snrsc","snr",...
        FFTLength=ofdmInfo.FFTLength,...
        NumActiveSubcarriers=ofdmInfo.NumTones);

    numPacketErrors = 0;

    for pkt = 1:num_packet
        disp(isnr)
        disp(pkt)

        psduLength = getPSDULength(cfgHE); 
        txPSDU = randi([0 1],psduLength*8,1);
        tx = wlanWaveformGenerator(txPSDU,cfgHE);

        txPad = [tx; zeros(50,cfgHE.NumTransmitAntennas)];
       
        reset(tgaxChannel); 
        rx = txPad; % AWGN

        rx = awgn(rx,packetSNR);

        coarsePktOffset = wlanPacketDetect(rx,chanBW);
        if isempty(coarsePktOffset) 
            numPacketErrors = numPacketErrors + 1;
            continue; 
        end

        lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
        rx = frequencyOffset(rx,fs,-coarseFreqOff);

        nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);

        pktOffset = coarsePktOffset+finePktOffset;

        if pktOffset>50        
            continue; 
        end

        rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
        fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
        rx = frequencyOffset(rx,fs,-fineFreqOff);

        rxHELTF = rx(pktOffset+(ind.HELTF(1):ind.HELTF(2)),:);
        heltfDemod = wlanHEDemodulate(rxHELTF,'HE-LTF',cfgHE);
        [chanEst,pilotEst] = wlanHELTFChannelEstimate(heltfDemod,cfgHE);

        rxData = rx(pktOffset+(ind.HEData(1):ind.HEData(2)),:);
        demodSym = wlanHEDemodulate(rxData,'HE-Data',cfgHE);

        demodSym = wlanHETrackPilotError(demodSym,chanEst,cfgHE,'HE-Data');

        nVarEst = wlanHEDataNoiseEstimate(demodSym(ofdmInfo.PilotIndices,:,:),pilotEst,cfgHE);

        demodDataSym = demodSym(ofdmInfo.DataIndices,:,:);
        chanEstData = chanEst(ofdmInfo.DataIndices,:,:);

        [eqDataSym,csi] = wlanHEEqualize(demodDataSym,chanEstData,nVarEst,cfgHE,'HE-Data');

        rxPSDU = wlanHEDataBitRecover(eqDataSym,nVarEst,csi,cfgHE,'LDPCDecodingMethod','norm-min-sum');
        
        if isequal(txPSDU, rxPSDU)        
            throughput_sim(1, isnr) = throughput_sim(1, isnr) + sum(txPSDU == rxPSDU);
        else
            numPacketErrors = numPacketErrors + 1;
        end

    end

    error_count_sim = [error_count_sim, numPacketErrors];

end

