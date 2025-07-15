clear;
clc;

throughput_arr_sim = [];
harq_error_sim = [];
harq_UL_error_sim = [];

noise_DL_sim = [];
noise_UL_sim = [];

dB_div = [-15; -10; -5; 0];

is_awgn = 1;

for div = 1:length(dB_div)
    error_sim = [];
    
    NFrames = 50;
    SNRIn = linspace(-15, 20, 70);

    throughput_arr = [];
    
    code_rate_arr = [];
    mod_order_arr = [];
    
    harq_error_count_arr = [];

    noise_DL_arr = [];
    noise_UL_arr = [];
    
    cqi_arr = [];
    sinr_arr = [];
        
    txMode = 'TM1';
    
    count_arr = [];
    error_count_arr = [];

    harq_UL_error_count_arr = [];

    block_size_arr = [];
    harq_error_indicator_arr = [];
    
    for snrIdx = 1:numel(SNRIn)

        count = 0;
        error_count = 0;
        harq_UL_error_count = 0;
    
        harq_count = 0;
    
        simulationParameters = [];
        simulationParameters.NDLRB = 50;    
        simulationParameters.PDSCH.PRBSet = (0:49)';
    
        simulationParameters.PDSCH.TargetCodeRate = 0.5;
        simulationParameters.PDSCH.Modulation = {'QPSK'};
        
        simulationParameters.PDSCH.TxScheme = 'Port0';
        simulationParameters.PDSCH.DCIFormat = 'Format1';
        simulationParameters.CellRefP = 1;
        simulationParameters.PDSCH.NSoftbits = 10^10;
    
        simulationParameters.TotSubframes = 1;
        simulationParameters.PDSCH.CSI = 'On';
        
        enb = lteRMCDL(simulationParameters);
        
        rvSequence = enb.PDSCH.RVSeq;
        trBlkSizes = enb.PDSCH.TrBlkSizes;
        
        ncw = length(string(enb.PDSCH.Modulation));
        
        pmiDelay = 8;
        
        hDisplayENBParameterSummary(enb, txMode);
        
        channel.Seed = 6;
        channel.NRxAnts = 1;

        if is_awgn == 1
            channel.DelayProfile = 'Off';        
            channel.DopplerFreq = 0;             
        else
            channel.DelayProfile = 'EVA';        
            channel.DopplerFreq = 5;             
        end

        channel.MIMOCorrelation = 'Low';
        channel.NTerms = 16;
        channel.ModelType = 'GMEDS';
        channel.InitPhase = 'Random';
        channel.NormalizePathGains = 'On';
        channel.NormalizeTxAnts = 'On';
        
        ofdmInfo = lteOFDMInfo(enb);
        channel.SamplingRate = ofdmInfo.SamplingRate;
        
        channel.InitTime = 0;
        [~,chInfo] = lteFadingChannel(channel,0);
        maxChDelay = ceil(max(chInfo.PathSampleDelays)) + chInfo.ChannelFilterDelay;
        
        perfectChanEstimator = false;
        
        cec.PilotAverage = 'UserDefined';
        cec.FreqWindow = 41;
        cec.TimeWindow = 27;
        cec.InterpType = 'Cubic';
        cec.InterpWindow = 'Centered';
        cec.InterpWinSize = 1;

        cec_UL = struct;
        cec_UL.PilotAverage = 'UserDefined';
        cec_UL.FreqWindow = 12;
        cec_UL.TimeWindow = 1;
        cec_UL.InterpType = 'cubic';
        
        displaySimulationInformation = true;
        
        dims = lteDLResourceGridSize(enb);
        P = dims(3);
        
        maxThroughput = zeros(length(SNRIn),1);
        simThroughput = zeros(length(SNRIn),1);
        
        [~,~,enbOut] = lteRMCDLTool(enb, []);
        harqProcessSequence = enbOut.PDSCH.HARQProcessSequence;
        
        legendString = ['Throughput: ' char(enb.PDSCH.TxScheme)];
        allRvSeqPtrHistory = cell(1,numel(SNRIn));
        nFFT = ofdmInfo.Nfft;
        
        code_rate_arr = [code_rate_arr; enb.PDSCH.ActualCodeRate];
    
        SNRdB = SNRIn(snrIdx);
        fprintf('\nSimulating at %g dB SNR for %d Frame(s)\n' ,SNRdB, NFrames);
    
        offsets = 0;
        offset = 0;
        blkCRC = [];
        bitTput = [];
        txedTrBlkSizes = [];
        pmiIdx = 0;
        rvSeqPtrHistory = NaN(ncw, NFrames*10);
        harqProcesses = hNewHARQProcess(enb);
        pmidims = ltePMIInfo(enb,enb.PDSCH);
        txPMIs = randi([0 pmidims.MaxPMI], pmidims.NSubbands, pmiDelay);
        offsetused = 0;

        for subframeNo = 0:(NFrames*10-1)

            enb.NSubframe = subframeNo;

            ue1.NULRB = 6;
            ue1.NSubframe = enb.NSubframe;
            ue1.NCellID = 9;
            ue1.RNTI = 61;
            ue1.CyclicPrefixUL = 'Normal';
            ue1.Hopping = 'Off';
            ue1.Shortened = 0;
            ue1.NTxAnts = 1;

            pucch1.ResourceIdx = 0;
            pucch1.DeltaShift = 1;
            pucch1.CyclicShifts = 1;
            pucch1.ResourceSize = 0;

            harqID = harqProcessSequence(mod(subframeNo, length(harqProcessSequence))+1);

            if harqID == 0
                continue;
            end

            harqProcesses(harqID) = hHARQScheduling(harqProcesses(harqID), subframeNo, rvSequence);

            trBlk = trBlkSizes(:, mod(subframeNo, 10)+1).';

            if displaySimulationInformation
                disp(' ');
                disp(['Subframe: ' num2str(subframeNo) '. HARQ process ID: ' num2str(harqID)]);
            end

            rvSeqPtrHistory(:,subframeNo+1) = harqProcesses(harqID).txConfig.RVIdx.';

            enb.PDSCH = harqProcesses(harqID).txConfig;
            data = harqProcesses(harqID).data;

            if strcmpi(enb.PDSCH.TxScheme,'SpatialMux')
                pmiIdx = mod(subframeNo, pmiDelay);
                enb.PDSCH.PMISet = txPMIs(:, pmiIdx+1);
            end          

            txWaveform = lteRMCDLTool(enb, data);
            txWaveform =  [txWaveform; zeros(maxChDelay, P)];

            harqProcessSequence = enbOut.PDSCH.HARQProcessSequence;

            channel.InitTime = subframeNo/1000;

            if is_awgn == 1
                rxWaveform = txWaveform;
            else
                rxWaveform = lteFadingChannel(channel, txWaveform);
            end   

            SNR = 10^((SNRdB-enb.PDSCH.Rho)/20);
            N0 = 1/(sqrt(2.0*enb.CellRefP*double(nFFT))*SNR);
            noise = N0*complex(randn(size(rxWaveform)), randn(size(rxWaveform)));
            rxWaveform = rxWaveform + noise;

            if (mod(subframeNo,10) == 0)
                offset = lteDLFrameOffset(enb, rxWaveform);
                if (offset > maxChDelay)
                    offset = offsets(end);
                end
                offsets = [offsets offset];
            end

            rxWaveform = rxWaveform(1+offset:end, :);
            rxSubframe = lteOFDMDemodulate(enb, rxWaveform);

            if(perfectChanEstimator)
                estChannelGrid = lteDLPerfectChannelEstimate(enb, channel, offset);
                noiseGrid = lteOFDMDemodulate(enb, noise(1+offset:end ,:));
                noiseEst = var(noiseGrid(:));
            else
                [estChannelGrid, noiseEst] = lteDLChannelEstimate(enb, enb.PDSCH, cec, rxSubframe);
            end

            pdschIndices = ltePDSCHIndices(enb, enb.PDSCH, enb.PDSCH.PRBSet);
            [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxSubframe*(10^(-enb.PDSCH.Rho/20)), estChannelGrid);
            [dlschBits, ~] = ltePDSCHDecode(enb, enb.PDSCH, pdschRx, pdschHest, noiseEst);

            [decbits, harqProcesses(harqID).blkerr,harqProcesses(harqID).decState] = lteDLSCHDecode(enb, enb.PDSCH, trBlk, dlschBits, harqProcesses(harqID).decState);

            count = count + 1;

            info = lteSCFDMAInfo(ue1);
            ue_message = ~harqProcesses(harqID).blkerr;
            SNRdB_UL = SNRdB+dB_div(div, 1);
            SNR_UL = 10^(SNRdB_UL/20);
            N_UL = 1/(SNR_UL*sqrt(double(info.Nfft)))/sqrt(2.0);
            pucch1Sym = ltePUCCH1(ue1, pucch1, ue_message);    
            pucch1DRSSym = ltePUCCH1DRS(ue1, pucch1);
            pucch1Indices = ltePUCCH1Indices(ue1, pucch1);                    
            pucch1DRSIndices = ltePUCCH1DRSIndices(ue1, pucch1);                  
            grid1 = lteULResourceGrid(ue1);
            grid1(pucch1Indices) = pucch1Sym;
            grid1(pucch1DRSIndices) = pucch1DRSSym;
            txwave_UL = lteSCFDMAModulate(ue1, grid1);

            if is_awgn == 1
                rxwave_UL = txwave_UL; 
            else
                rxwave_UL = lteFadingChannel(channel,[txwave_UL; zeros(25,1)]);
            end

            noise_UL = N_UL*complex(randn(size(rxwave_UL)), randn(size(rxwave_UL)));
            rxwave_UL = rxwave_UL + noise_UL;

            if (is_awgn) ~= 1
                offset1 = lteULFrameOffsetPUCCH1(ue1, pucch1, rxwave_UL);
                if (offset1 < 25)
                    offsetused = offset1;
                end
            end
            
            rxgrid1 = lteSCFDMADemodulate(ue1,rxwave_UL(1+offsetused:end,:));
            [H1, n0] = lteULChannelEstimatePUCCH1(ue1, pucch1, cec_UL, rxgrid1);
            [pucchrx1, pucchH1] = lteExtractResources(pucch1Indices, rxgrid1, H1);
            eqgrid1 = lteULResourceGrid(ue1);   
            eqgrid1(pucch1Indices) = lteEqualizeMMSE(pucchrx1, pucchH1, n0);
            ue_message_decoded = ltePUCCH1Decode_custom(ue1, pucch1, length(ue_message), eqgrid1(pucch1Indices));

            disp("--------------------")
            disp(snrIdx)
            disp("--------------------")

            if any(ue_message ~= ue_message_decoded) || isempty(ue_message_decoded)
                harq_UL_error_count = harq_UL_error_count + 1;
                harqProcesses(harqID).blkerr = 1;
                disp("HARQ Error")
            end

            if displaySimulationInformation
                if any(harqProcesses(harqID).blkerr) 
                    disp(['Block error. RV index: ' num2str(harqProcesses(harqID).txConfig.RVIdx) ', CRC: ' num2str(harqProcesses(harqID).blkerr)])
                    error_count = error_count + 1;
                    harq_count = harq_count + 1;
                else
                    disp(['No error. RV index: ' num2str(harqProcesses(harqID).txConfig.RVIdx) ', CRC: ' num2str(harqProcesses(harqID).blkerr)])
                end
            end

            if any(trBlk)
                blkCRC = [blkCRC harqProcesses(harqID).blkerr];
                bitTput = [bitTput trBlk.*(1- harqProcesses(harqID).blkerr)]; 
                txedTrBlkSizes = [txedTrBlkSizes trBlk];
            end

            if strcmpi(enb.PDSCH.TxScheme,'SpatialMux')
                PMI = ltePMISelect(enb, enb.PDSCH, estChannelGrid, noiseEst);
                txPMIs(:, pmiIdx+1) = PMI;
            end
            block_size_arr = [block_size_arr, trBlk];
            harq_error_indicator_arr = [harq_error_indicator_arr, harqProcesses(harqID).blkerr];
        end
        
        count_arr = [count_arr, count];
        error_count_arr = [error_count_arr, error_count];
        mod_order_arr = [mod_order_arr; enb.PDSCH.Modulation];
        harq_UL_error_count_arr = [harq_UL_error_count_arr, harq_UL_error_count];
        harq_error_count_arr = [harq_error_count_arr, harq_count];
        maxThroughput(snrIdx) = sum(txedTrBlkSizes);
        simThroughput(snrIdx) = sum(bitTput,2);
        fprintf([['\nThroughput(Mbps) for ', num2str(NFrames) ' Frame(s) '],'= %.4f\n'], 1e-6*simThroughput(snrIdx)/(NFrames*10e-3));
        allRvSeqPtrHistory{snrIdx} = rvSeqPtrHistory;
        throughput_arr = [throughput_arr, simThroughput(snrIdx)];
    end

    throughput_arr_sim = [throughput_arr_sim; throughput_arr];
    harq_error_sim = [harq_error_sim; harq_error_count_arr];
    harq_UL_error_sim = [harq_UL_error_sim; harq_UL_error_count_arr];
    error_sim = [error_sim; error_count_arr];
end
