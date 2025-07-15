clc;
clear;

cfgMAC = wlanMACFrameConfig(FrameType='ACK');
cfgMAC.Address1 = '0012345678AB'; % Receiver MAC address (dummy)

[ackFrameHex , ackBits] = wlanMACFrame(cfgMAC, 'Bits'); 

ackFrameBytes = hex2dec(ackFrameHex); 
ackFrameBits = reshape(de2bi(ackFrameBytes, 8, 'left-msb')', [], 1); 


cfgNonHT = wlanNonHTConfig('Modulation', 'OFDM'); 
cfgNonHT.PSDULength = numel(ackFrameBits) / 8; 

chanBW = cfgNonHT.ChannelBandwidth;
tgaxChannel = wlanTGaxChannel;
tgaxChannel.DelayProfile = 'Model-A';
tgaxChannel.NumTransmitAntennas = cfgNonHT.NumTransmitAntennas;
tgaxChannel.NumReceiveAntennas = 1;
tgaxChannel.TransmitReceiveDistance = 5; 
tgaxChannel.ChannelBandwidth = chanBW;
tgaxChannel.LargeScaleFadingEffect = 'None';
tgaxChannel.NormalizeChannelOutputs = false;
fs = wlanSampleRate(cfgNonHT);
tgaxChannel.SampleRate = fs;

snr = linspace(0, 10, 40);
numSNR = numel(snr); 
ACK_prob_sim = zeros(1, numSNR);

ofdmInfo = wlanNonHTOFDMInfo('NonHT-Data');

ind = wlanFieldIndices(cfgNonHT);

num_packet = 10^5; 

for isnr = 1:numSNR
    packetSNR = convertSNR(snr(isnr), "snrsc", "snr", ...
        FFTLength=ofdmInfo.FFTLength, ...
        NumActiveSubcarriers=ofdmInfo.NumTones);

    num_decode_correct = 0;

    for pkt = 1:num_packet
        disp('--------------------------------')
        disp(isnr)
        disp(pkt)
        disp('--------------------------------')

        txWaveform = wlanWaveformGenerator(ackFrameBits, cfgNonHT); % Generate PHY waveform

        txPad = [txWaveform; zeros(50, cfgNonHT.NumTransmitAntennas)];

        reset(tgaxChannel);

        rxWaveform = txPad;

        rxWaveform = awgn(rxWaveform, packetSNR);

        coarsePktOffset = wlanPacketDetect(rxWaveform, chanBW);
        if isempty(coarsePktOffset)
            continue;
        end

        lstf = rxWaveform(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)), :);
        coarseFreqOff = wlanCoarseCFOEstimate(lstf, chanBW);
        rxWaveform = frequencyOffset(rxWaveform, fs, -coarseFreqOff);

        nonhtfields = rxWaveform(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)), :);
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields, chanBW);

        pktOffset = coarsePktOffset + finePktOffset;
        if pktOffset > 50
            continue;
        end

        rxLLTF = rxWaveform(pktOffset+(ind.LLTF(1):ind.LLTF(2)), :);
        fineFreqOff = wlanFineCFOEstimate(rxLLTF, chanBW);
        rxWaveform = frequencyOffset(rxWaveform, fs, -fineFreqOff);

        rxLLTF = rxWaveform(pktOffset+(ind.LLTF(1):ind.LLTF(2)), :);
        lltfDemod = wlanLLTFDemodulate(rxLLTF, cfgNonHT, 1);
        chanEst = wlanLLTFChannelEstimate(lltfDemod, cfgNonHT);
        nVar = wlanLLTFNoiseEstimate(lltfDemod);

        rxPSDU = wlanNonHTDataRecover(rxWaveform(pktOffset+(ind.NonHTData(1):ind.NonHTData(2)), :), ...
                                      chanEst, nVar, cfgNonHT);

        if isequal(ackFrameBits, rxPSDU)
            num_decode_correct = num_decode_correct + 1;
        end

    end

    ACK_prob_sim(1, isnr) = num_decode_correct / num_packet;

end

% Convert ACK probability to Block Error Rate (BLER)
BLER = 1 - ACK_prob_sim;
BLER(BLER == 0) = NaN; % Avoid log(0) issues


