% Modified from https://www.mathworks.com/help/lte/ug/pdsch-throughput-conformance-test-for-single-antenna-tm1-transmit-diversity-tm2-open-loop-tm3-and-closed-loop-tm4-6-spatial-multiplexing.html

clear;
clc;

is_awgn = 1; % 1 for AWGN 

for MCS=1:1 %%%%%%%%%%% testing
    
    throughput_arr_sim = [];
    harq_error_sim = [];

    repeat_counter_arr_sim = [];

    NFrames = 50;                % Number of frames, 1 Frame (10 ms) = 10 Subframe
    SNRIn = linspace(-10, 10, 40);   % SNR range in dB
    
    throughput_arr = [];
    
    code_rate_arr = [];
    mod_order_arr = [];
    
    harq_error_count_arr = [];

    repeat_counter_arr = [];
    
    cqi_arr = [];
    sinr_arr = [];
    
    txMode = 'TM1';
    
    count_arr = [];
    error_count_arr = [];

    block_size_arr = [];
    harq_error_indicator_arr = [];

    max_num_retrans = 4;
    coding_rate_arr = [0.25; 0.5; 0.75];
    
    for snrIdx = 1:numel(SNRIn)

        repeat_counter = [];

        count = 0;
        error_count = 0;
        retrans_count = 0;

        data_buffer = [];
    
        harq_count = 0;
        harq_buffer = {0, 0, 0, 0, 0, 0, 0, 0};
    
        simulationParameters = []; 
        simulationParameters.NDLRB = 50;    
        simulationParameters.PDSCH.PRBSet = (0:49)';
        simulationParameters.PDSCH.TargetCodeRate = coding_rate_arr(MCS, 1);
        simulationParameters.PDSCH.Modulation = {'QPSK'};
        
        simulationParameters.PDSCH.TxScheme = 'Port0';
        simulationParameters.PDSCH.DCIFormat = 'Format1';
        simulationParameters.CellRefP = 1;

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

        pmidims = ltePMIInfo(enb,enb.PDSCH);
        txPMIs = randi([0 pmidims.MaxPMI], pmidims.NSubbands, pmiDelay);
    
        for subframeNo = 0:(NFrames*10-1)
            enb.NSubframe = subframeNo;        
    
            trBlk = trBlkSizes(:, mod(subframeNo, 10)+1).';
            data = randi([0 1], trBlk, 1);

            if ~isempty(data_buffer)
                data = data_buffer;
            end

            if displaySimulationInformation
                disp(' ');
                disp(['Subframe: ' num2str(subframeNo)]);
            end        

            txWaveform = lteRMCDLTool(enb, data);
            txWaveform =  [txWaveform; zeros(maxChDelay, P)]; %#ok<AGROW>

            channel.InitTime = subframeNo/1000;
    
            if is_awgn == 1
                rxWaveform = txWaveform; % AWGN
            else
                rxWaveform = lteFadingChannel(channel, txWaveform);
            end   

            SNR = 10^((SNRdB-enb.PDSCH.Rho)/20);
            N0 = 1/(sqrt(2.0*enb.CellRefP*double(nFFT))*SNR);
    
            noise = N0*complex(randn(size(rxWaveform)), ...
                                randn(size(rxWaveform)));

            rxWaveform = rxWaveform + noise;    

            if (mod(subframeNo,10) == 0)
                offset = lteDLFrameOffset(enb, rxWaveform);
                if (offset > maxChDelay)
                    offset = offsets(end);
                end
                offsets = [offsets offset]; %#ok
            end

            rxWaveform = rxWaveform(1+offset:end, :);

            rxSubframe = lteOFDMDemodulate(enb, rxWaveform);
    
            if(perfectChanEstimator)
                estChannelGrid = lteDLPerfectChannelEstimate(enb, channel, offset); %#ok
                noiseGrid = lteOFDMDemodulate(enb, noise(1+offset:end ,:));
                noiseEst = var(noiseGrid(:));
            else
                [estChannelGrid, noiseEst] = lteDLChannelEstimate( ...
                    enb, enb.PDSCH, cec, rxSubframe);
            end
    
            pdschIndices = ltePDSCHIndices(enb, enb.PDSCH, enb.PDSCH.PRBSet);

            [pdschRx, pdschHest] = lteExtractResources(pdschIndices, ...
                rxSubframe*(10^(-enb.PDSCH.Rho/20)), estChannelGrid);
    
            [dlschBits, ~] = ltePDSCHDecode(...
                                 enb, enb.PDSCH, pdschRx, pdschHest, noiseEst);
    
            [decbits, blkerr,~] = ...
                lteDLSCHDecode(enb, enb.PDSCH, trBlk, dlschBits);
    
            count = count + 1;

            if displaySimulationInformation
                if any(blkerr) 
                    disp(['Block error.' ', CRC: ' num2str(blkerr)])
                    error_count = error_count + 1;
                    retrans_count = retrans_count + 1;

                    data_buffer = data;
    
                    harq_count = harq_count + 1;
                else
                    disp('No error.')
                    repeat_counter = [repeat_counter, retrans_count];
                    retrans_count = 0;
                    data_buffer = [];
                end
            end

            if retrans_count == max_num_retrans
                repeat_counter = [repeat_counter, retrans_count];
                retrans_count = 0;
                data_buffer = [];
                disp("MAX RETRANSMISSION")
            end

            if any(trBlk)
                blkCRC = [blkCRC blkerr]; %#ok<AGROW>
                bitTput = [bitTput trBlk.*(1-blkerr)]; %#ok<AGROW>
                txedTrBlkSizes = [txedTrBlkSizes trBlk]; %#ok<AGROW>
            end
    
            block_size_arr = [block_size_arr, trBlk];
            harq_error_indicator_arr = [harq_error_indicator_arr, blkerr];
        end
        
        count_arr = [count_arr, count];
        error_count_arr = [error_count_arr, error_count];
    
        mod_order_arr = [mod_order_arr; enb.PDSCH.Modulation];
    
        harq_error_count_arr = [harq_error_count_arr, harq_count];

        repeat_counter_arr = [repeat_counter_arr, {repeat_counter}];

        maxThroughput(snrIdx) = sum(txedTrBlkSizes); 
        simThroughput(snrIdx) = sum(bitTput,2);      

        fprintf([['\nThroughput(Mbps) for ', num2str(NFrames) ' Frame(s) '],...
            '= %.4f\n'], 1e-6*simThroughput(snrIdx)/(NFrames*10e-3));
    
        allRvSeqPtrHistory{snrIdx} = rvSeqPtrHistory;
    
        throughput_arr = [throughput_arr, simThroughput(snrIdx)];            
    end

    throughput_arr_sim = [throughput_arr_sim; throughput_arr];
    harq_error_sim = [harq_error_sim; harq_error_count_arr];

    repeat_counter_arr_sim = [repeat_counter_arr_sim; repeat_counter_arr];

    
end

