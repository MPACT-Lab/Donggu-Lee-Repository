clear;
clc;

NFrames = 10^4;                
SNRIn = linspace(-10, 5, 30);   
BLER = zeros(length(SNRIn), 1);

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

for snrIdx = 1:numel(SNRIn)
    count = 0;
    error_count = 0;
    harq_UL_error_count = 0;

    harq_count = 0;

    rng('default');

    simLocal = simParameters;
    waveinfoLocal = waveformInfo;

    carrier = simLocal.Carrier;
    pdsch = simLocal.PDSCH;
    pdschextra = simLocal.PDSCHExtension;
    decodeDLSCH = decodeDLSCH;  
    decodeDLSCH.reset();        
    pathFilters = [];

    SNRdB = simLocal.SNRIn(snrIdx);

    harqSequence = 0:pdschextra.NHARQProcesses-1;

    harqEntity = HARQEntity(harqSequence,rvSeq,pdsch.NumCodewords);

    NSlots = simLocal.NFrames * carrier.SlotsPerFrame;

    offset = 0;
    offset_UL = 0;

    SNR = 10^(SNRdB/10);
    N0 = 1/sqrt(simLocal.NRxAnts*double(waveinfoLocal.Nfft)*SNR);
    nPowerPerRE = N0^2*double(waveinfoLocal.Nfft);

    for nslot = 0:NSlots-1
        
        disp("--------------------------------")
        disp(snrIdx)
        disp(nslot)
        disp("--------------------------------")

        carrier.NSlot = nslot;

        [pdschIndices,pdschIndicesInfo] = nrPDSCHIndices(carrier,pdsch);
        trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschIndicesInfo.NREPerPRB,pdschextra.TargetCodeRate,pdschextra.XOverhead);

        for cwIdx = 1:pdsch.NumCodewords            
            if harqEntity.NewData(cwIdx)
                trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
                
                if harqEntity.SequenceTimeout(cwIdx)
                    resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                end
            end
        end

        codedTrBlocks = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers, ...
            pdschIndicesInfo.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

        pdschGrid = nrResourceGrid(carrier,simLocal.NTxAnts,OutputDataType=simLocal.DataType);

        pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlocks);
        pdschGrid(pdschIndices) = pdschSymbols;

        dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
        dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);
        pdschGrid(dmrsIndices) = dmrsSymbols;

        ptrsSymbols = nrPDSCHPTRS(carrier,pdsch);
        ptrsIndices = nrPDSCHPTRSIndices(carrier,pdsch);
        pdschGrid(ptrsIndices) = ptrsSymbols;

        txWaveform = nrOFDMModulate(carrier,pdschGrid);

        rxWaveform = txWaveform; % AWGN

        noise = N0*randn(size(rxWaveform),"like",rxWaveform);
        rxWaveform = rxWaveform + noise;

        if (simLocal.PerfectChannelEstimator)
            offset = timingOffset;
        else
            [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
            offset = hSkipWeakTimingOffset(offset,t,mag);
        end

        rxWaveform = rxWaveform(1+offset:end,:);

        rxGrid = nrOFDMDemodulate(carrier,rxWaveform);
        [K,L,R] = size(rxGrid);
        if (L < carrier.SymbolsPerSlot)
            rxGrid = cat(2,rxGrid,zeros(K,carrier.SymbolsPerSlot-L,R));
        end

        if (simLocal.PerfectChannelEstimator)
            estChannelGridAnts = ofdmResponse;

            noiseEst = nPowerPerRE;

            [pdschRx,pdschHest,~,pdschHestIndices] = nrExtractResources(pdschIndices,rxGrid,estChannelGridAnts);

            pdschHest = nrPDSCHPrecode(carrier,pdschHest,pdschHestIndices,permute(wtx,[2 1 3]));
        else

            [estChannelGridPorts,noiseEst] = hSubbandChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,pdschextra.PRGBundleSize,'CDMLengths',pdsch.DMRS.CDMLengths);

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

        [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);

        if (simLocal.DisplayDiagnostics)
            plotLayerEVM(NSlots,nslot,pdsch,size(pdschGrid),pdschIndices,pdschSymbols,pdschEq);
        end

        csi = nrLayerDemap(csi); 
        for cwIdx = 1:pdsch.NumCodewords
            Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); 
            csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 
            dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   
        end

        decodeDLSCH.TransportBlockLength = trBlkSizes;
        [decbits,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

        procstatus = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschIndicesInfo.G);

        if blkerr == 1
            error_count = error_count + 1;
        end
    end
    BLER(snrIdx) = error_count / NSlots;
end
