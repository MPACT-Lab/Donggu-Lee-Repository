clear;
clc;

throughput_arr_sim = [];
harq_error_sim = [];
harq_UL_error_sim = [];

dB_div = [-15; -10; -5; 0];

is_awgn = 1;

for div = 1:length(dB_div)
    error_sim = [];

    NFrames = 50;                

    SNRIn = linspace(-10, 20, 60);   % SNR range in dB

    
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

    simParameters = struct();       
    simParameters.NFrames = NFrames;
    simParameters.SNRIn = SNRIn; % SNR range (dB)
    simParameters.PerfectChannelEstimator = false;
    simParameters.DisplaySimulationInformation = true;
    simParameters.DisplayDiagnostics = false;


    simParameters.Carrier = nrCarrierConfig;         
    simParameters.Carrier.NSizeGrid = 52;            
    simParameters.Carrier.SubcarrierSpacing = 15;    
    simParameters.Carrier.CyclicPrefix = 'Normal';   
    simParameters.Carrier.NCellID = 9;               

    simParameters.PDSCH = nrPDSCHConfig;      
    simParameters.PDSCHExtension = struct();  

    
    simParameters.PDSCH.PRBSet = 0:simParameters.Carrier.NSizeGrid-1;                 
    simParameters.PDSCH.SymbolAllocation = [0,simParameters.Carrier.SymbolsPerSlot];  
    simParameters.PDSCH.MappingType = 'A';     

    simParameters.PDSCH.NID = simParameters.Carrier.NCellID;
    simParameters.PDSCH.RNTI = 1;
    
    simParameters.PDSCH.VRBToPRBInterleaving = 0; 
    simParameters.PDSCH.VRBBundleSize = 4;
    
    simParameters.PDSCH.NumLayers = 1;            

    simParameters.PDSCH.Modulation = 'QPSK';     
    simParameters.PDSCHExtension.TargetCodeRate = 0.5;

    simParameters.PDSCH.DMRS.DMRSPortSet = 0:simParameters.PDSCH.NumLayers-1; 
    simParameters.PDSCH.DMRS.DMRSTypeAPosition = 2;      
    simParameters.PDSCH.DMRS.DMRSLength = 1;             
    simParameters.PDSCH.DMRS.DMRSAdditionalPosition = 2; 
    simParameters.PDSCH.DMRS.DMRSConfigurationType = 1;  
    simParameters.PDSCH.DMRS.NumCDMGroupsWithoutData = 1;
    simParameters.PDSCH.DMRS.NIDNSCID = 1;               
    simParameters.PDSCH.DMRS.NSCID = 0;                  

    simParameters.PDSCH.EnablePTRS = 0;                 
    simParameters.PDSCH.PTRS.TimeDensity = 1;            
    simParameters.PDSCH.PTRS.FrequencyDensity = 2;      
    simParameters.PDSCH.PTRS.REOffset = '00';            
    simParameters.PDSCH.PTRS.PTRSPortSet = [];           

    simParameters.PDSCH.ReservedPRB{1}.SymbolSet = [];   
    simParameters.PDSCH.ReservedPRB{1}.PRBSet = [];      
    simParameters.PDSCH.ReservedPRB{1}.Period = [];      

    simParameters.PDSCHExtension.PRGBundleSize = [];     

    simParameters.PDSCHExtension.XOverhead = 6*simParameters.PDSCH.EnablePTRS; 
    simParameters.PDSCHExtension.NHARQProcesses = 8;    
    simParameters.PDSCHExtension.EnableHARQ = true;     

    simParameters.PDSCHExtension.LDPCDecodingAlgorithm = 'Normalized min-sum';
    simParameters.PDSCHExtension.MaximumLDPCIterationCount = 6;

    simParameters.NTxAnts = 1;
    simParameters.NRxAnts = 1;

    simParameters.DataType = 'single';

    waveformInfo = nrOFDMInfo(simParameters.Carrier); 

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

    channel.SamplingRate = waveformInfo.SampleRate;
    
    channel.InitTime = 0;
    [~,chInfo] = lteFadingChannel(channel,0); 
    maxChDelay = ceil(max(chInfo.PathSampleDelays)) + chInfo.ChannelFilterDelay;

    maxThroughput = zeros(length(simParameters.SNRIn),1);
        simThroughput = zeros(length(simParameters.SNRIn),1);
    
    if simParameters.PDSCHExtension.EnableHARQ
        rvSeq = [0 2 3 1];
    else
        rvSeq = 0;
    end

    encodeDLSCH = nrDLSCH;
    encodeDLSCH.MultipleHARQProcesses = true;
    encodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;

    decodeDLSCH = nrDLSCHDecoder;
    decodeDLSCH.MultipleHARQProcesses = true;
    decodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;
    decodeDLSCH.LDPCDecodingAlgorithm = simParameters.PDSCHExtension.LDPCDecodingAlgorithm;
    decodeDLSCH.MaximumLDPCIterationCount = simParameters.PDSCHExtension.MaximumLDPCIterationCount;

    simParameters_UL = struct;

    carrier_UL = nrCarrierConfig;
    carrier_UL.NCellID = simParameters.Carrier.NCellID;
    carrier_UL.SubcarrierSpacing = simParameters.Carrier.SubcarrierSpacing;
    carrier_UL.CyclicPrefix = simParameters.Carrier.CyclicPrefix;

    carrier_UL.NSizeGrid = 6;
    carrier_UL.NStartGrid = 0;

    pucch = nrPUCCH1Config;
    pucch.PRBSet = 0;
    pucch.SymbolAllocation = [0 14];
    pucch.FrequencyHopping = 'intraSlot';

    pucch.SecondHopStartPRB = (carrier_UL.NSizeGrid-1) - (numel(pucch.PRBSet)-1);
    pucch.GroupHopping = 'neither';
    pucch.HoppingID = 0;

    simParameters_UL.NTxAnts = 1;
    simParameters_UL.NRxAnts = 1;
    
    perfectChannelEstimator = false;

    waveformInfo_UL = nrOFDMInfo(carrier_UL);
    
    for snrIdx = 1:numel(SNRIn)
        count = 0;
        error_count = 0;
        harq_UL_error_count = 0;
    
        harq_count = 0;

        rng('default');

        SNRdB = SNRIn(snrIdx);
        fprintf('\nSimulating transmission scheme 1 (%dx%d) and SCS=%dkHz at %gdB SNR for %d 10ms frame(s)\n', ...
            simParameters.NTxAnts,simParameters.NRxAnts,simParameters.Carrier.SubcarrierSpacing, ...
            SNRdB,simParameters.NFrames);

        harqSequence = 0:simParameters.PDSCHExtension.NHARQProcesses-1;

        harqEntity = HARQEntity(harqSequence,rvSeq,simParameters.PDSCH.NumCodewords);

        NSlots = simParameters.NFrames * simParameters.Carrier.SlotsPerFrame;

        offset = 0;
        offset_UL = 0;
 
        SNR = 10^(SNRdB/10);
        N0 = 1/sqrt(simParameters.NRxAnts*double(waveformInfo.Nfft)*SNR);
        nPowerPerRE = N0^2*double(waveformInfo.Nfft);

        for nslot = 0:NSlots-1
            channel.InitTime = nslot/1000;

            carrier.NSlot = nslot;
    
            [pdschIndices,pdschIndicesInfo] = nrPDSCHIndices(simParameters.Carrier, simParameters.PDSCH);
            trBlkSizes = nrTBS(simParameters.PDSCH.Modulation,simParameters.PDSCH.NumLayers,numel(simParameters.PDSCH.PRBSet), ...
                pdschIndicesInfo.NREPerPRB,simParameters.PDSCHExtension.TargetCodeRate,simParameters.PDSCHExtension.XOverhead);

            for cwIdx = 1:simParameters.PDSCH.NumCodewords

                if harqEntity.NewData(cwIdx)
                    trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                    setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);

                    if harqEntity.SequenceTimeout(cwIdx)
                        resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                    end
                end
            end

            codedTrBlocks = encodeDLSCH(simParameters.PDSCH.Modulation,simParameters.PDSCH.NumLayers, ...
                pdschIndicesInfo.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

            pdschGrid = nrResourceGrid(simParameters.Carrier,simParameters.NTxAnts,OutputDataType=simParameters.DataType);

            pdschSymbols = nrPDSCH(simParameters.Carrier,simParameters.PDSCH,codedTrBlocks);
            pdschGrid(pdschIndices) = pdschSymbols;

            dmrsSymbols = nrPDSCHDMRS(simParameters.Carrier,simParameters.PDSCH);
            dmrsIndices = nrPDSCHDMRSIndices(simParameters.Carrier,simParameters.PDSCH);

            pdschGrid(dmrsIndices) = dmrsSymbols;
    

            ptrsSymbols = nrPDSCHPTRS(simParameters.Carrier,simParameters.PDSCH);
            ptrsIndices = nrPDSCHPTRSIndices(simParameters.Carrier,simParameters.PDSCH);

            pdschGrid(ptrsIndices) = ptrsSymbols;
    
            txWaveform = nrOFDMModulate(simParameters.Carrier,pdschGrid);

            if is_awgn == 1
                rxWaveform = txWaveform;
            else
                rxWaveform = lteFadingChannel(channel, txWaveform);
            end   

            noise = N0*randn(size(rxWaveform),"like",rxWaveform);
            rxWaveform = rxWaveform + noise;

            if (simParameters.PerfectChannelEstimator)
                offset = timingOffset;
            else
                [t,mag] = nrTimingEstimate(simParameters.Carrier,rxWaveform,dmrsIndices,dmrsSymbols);
                offset = hSkipWeakTimingOffset(offset,t,mag);
            end

            rxWaveform = rxWaveform(1+offset:end,:);

            rxGrid = nrOFDMDemodulate(simParameters.Carrier,rxWaveform);
            [K,L,R] = size(rxGrid);
            if (L < simParameters.Carrier.SymbolsPerSlot)
                rxGrid = cat(2,rxGrid,zeros(K,simParameters.Carrier.SymbolsPerSlot-L,R));
            end
    
            if (simParameters.PerfectChannelEstimator)
                estChannelGridAnts = ofdmResponse;

                noiseEst = nPowerPerRE;

                [pdschRx,pdschHest,~,pdschHestIndices] = nrExtractResources(pdschIndices,rxGrid,estChannelGridAnts);

                pdschHest = nrPDSCHPrecode(carrier,pdschHest,pdschHestIndices,permute(wtx,[2 1 3]));
            else

                [estChannelGridPorts,noiseEst] = hSubbandChannelEstimate(simParameters.Carrier,rxGrid,dmrsIndices,dmrsSymbols,...
                    simParameters.PDSCHExtension.PRGBundleSize,'CDMLengths',simParameters.PDSCH.DMRS.CDMLengths);

                noiseEst = mean(noiseEst,'all');
    
                [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChannelGridPorts);

            end

            [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

            if ~isempty(ptrsIndices)
                tempGrid = nrResourceGrid(carrier,pdsch.NumLayers);

                [ptrsRx,ptrsHest,~,~,ptrsHestIndices,ptrsLayerIndices] = nrExtractResources(ptrsIndices,rxGrid,estChannelGridAnts,tempGrid);
                ptrsHest = nrPDSCHPrecode(carrier,ptrsHest,ptrsHestIndices,permute(wtx,[2 1 3]));

                ptrsEq = nrEqualizeMMSE(ptrsRx,ptrsHest,noiseEst);
                tempGrid(ptrsLayerIndices) = ptrsEq;

                cpe = nrChannelEstimate(tempGrid,ptrsIndices,ptrsSymbols);

                cpe = angle(sum(cpe,[1 3 4]));

                tempGrid(pdschIndices) = pdschEq;
    
                symLoc = pdschIndicesInfo.PTRSSymbolSet(1)+1:pdschIndicesInfo.PTRSSymbolSet(end)+1;
                tempGrid(:,symLoc,:) = tempGrid(:,symLoc,:).*exp(-1i*cpe(symLoc));

                pdschEq = tempGrid(pdschIndices);
            end

            [dlschLLRs,rxSymbols] = nrPDSCHDecode(simParameters.Carrier,simParameters.PDSCH,pdschEq,noiseEst);

            if (simParameters.DisplayDiagnostics)
                plotLayerEVM(NSlots,nslot,pdsch,size(pdschGrid),pdschIndices,pdschSymbols,pdschEq);
            end

            csi = nrLayerDemap(csi); 
            for cwIdx = 1:simParameters.PDSCH.NumCodewords
                Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); 
                csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 
                dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);  
            end

            decodeDLSCH.TransportBlockLength = trBlkSizes;
            [decbits,blkerr] = decodeDLSCH(dlschLLRs,simParameters.PDSCH.Modulation,simParameters.PDSCH.NumLayers,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    
            procstatus = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschIndicesInfo.G);
            if (simParameters.DisplaySimulationInformation)
                fprintf('\n(%3.2f%%) NSlot=%d, %s',100*(nslot+1)/NSlots,nslot,procstatus);
            end

            SNRdB_UL = SNRdB + dB_div(div);

            carrier_UL.NSlot = mod(nslot, 160);  

            [pucchIndices,pucchIndicesInfo] = nrPUCCHIndices(carrier_UL,pucch);
            dmrsIndices = nrPUCCHDMRSIndices(carrier_UL,pucch);
            dmrsSymbols = nrPUCCHDMRS(carrier_UL,pucch);

            harq_message = ~blkerr;

            sr = [];

            pucchSymbols = nrPUCCH(carrier_UL,pucch,harq_message);

            pucchGrid = nrResourceGrid(carrier_UL,simParameters_UL.NTxAnts);

            pucchGrid(pucchIndices) = pucchSymbols;
            pucchGrid(dmrsIndices) = dmrsSymbols;

            txWaveform_UL = nrOFDMModulate(carrier_UL,pucchGrid);

            if is_awgn == 1
                rxWaveform_UL = txWaveform_UL;
            else
                rxWaveform_UL = lteFadingChannel(channel, txWaveform_UL);
            end           

            SNR_UL = 10^(SNRdB_UL/20);
            N0_UL = 1/(sqrt(2.0*simParameters_UL.NRxAnts*waveformInfo_UL.Nfft)*SNR_UL);
            noise_UL = N0_UL*complex(randn(size(rxWaveform_UL)),randn(size(rxWaveform_UL)));
            rxWaveform_UL = rxWaveform_UL + noise_UL;

            if perfectChannelEstimator == 1
                pathFilters = getPathFilters(channel);
                [offset_UL,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
            else
                [t,mag] = nrTimingEstimate(carrier_UL,rxWaveform_UL,dmrsIndices,dmrsSymbols);
                offset_UL = hSkipWeakTimingOffset(offset_UL,t,mag);
            end

            rxWaveform_UL = rxWaveform_UL(1+offset:end,:);

            rxGrid_UL = nrOFDMDemodulate(carrier_UL,rxWaveform_UL);
            [K_UL,L_UL,R_UL] = size(rxGrid_UL);
            if (L_UL < carrier_UL.SymbolsPerSlot)
                rxGrid_UL = cat(2,rxGrid_UL,zeros(K_UL,carrier_UL.SymbolsPerSlot-L_UL,R_UL));
            end

            if perfectChannelEstimator == 1
                estChannelGrid = nrPerfectChannelEstimate(carrier_UL,pathGains,pathFilters,offset,sampleTimes);

                noiseGrid = nrOFDMDemodulate(carrier_UL,noise(1+offset:end,:));
                noiseEst = var(noiseGrid(:));

                K = size(estChannelGrid,1);
                estChannelGrid = reshape(estChannelGrid,K*symbolsPerSlot*nRxAnts,nTxAnts);
                estChannelGrid = estChannelGrid*F.';
                estChannelGrid = reshape(estChannelGrid,K,symbolsPerSlot,nRxAnts,[]);
            else
                [estChannelGrid_UL,noiseEst_UL] = nrChannelEstimate(carrier_UL,rxGrid_UL,dmrsIndices,dmrsSymbols);                
            end

            [pucchRx,pucchHest] = nrExtractResources(pucchIndices,rxGrid_UL,estChannelGrid_UL);

            [pucchEq,csi] = nrEqualizeMMSE(pucchRx,pucchHest,noiseEst_UL);

            decodedBits = nrPUCCHDecode(carrier_UL, pucch, 1, pucchEq, noiseEst_UL);

            decodedVec = cell2mat(decodedBits);  

            if isempty(decodedVec) || any(decodedVec ~= harq_message)
                blkerr = 1;         
                disp('----------------PUCCH ACK/NACK decoding FAILED--------------------');
                harq_UL_error_count = harq_UL_error_count + 1;
            end

            simThroughput(snrIdx) = simThroughput(snrIdx) + sum(~blkerr .* trBlkSizes);
            maxThroughput(snrIdx) = maxThroughput(snrIdx) + sum(trBlkSizes);

        end
        harq_UL_error_sim = [harq_UL_error_sim, harq_UL_error_count];
    end

    
end