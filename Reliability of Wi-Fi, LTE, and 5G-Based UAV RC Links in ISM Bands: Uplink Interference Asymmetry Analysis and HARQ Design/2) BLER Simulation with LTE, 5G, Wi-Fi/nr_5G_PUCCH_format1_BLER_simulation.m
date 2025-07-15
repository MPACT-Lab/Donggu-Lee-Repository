clear;
clc;

SNRIn = linspace(10, 20, 60); 

NSubFrames = 10^5;
BLER = zeros(length(SNRIn), 1);

simParameters_UL = struct;

carrier_UL = nrCarrierConfig;
carrier_UL.NCellID = 9;
carrier_UL.SubcarrierSpacing = 15;
carrier_UL.CyclicPrefix = 'Normal';

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
    error_count = 0;
    rng('default');

    SNRdB = SNRIn(snrIdx);

    for nslot = 1:NSubFrames
        disp("--------------------------------")
        disp(snrIdx)
        disp(nslot)
        disp("--------------------------------")

        carrier_UL.NSlot = mod(nslot, 160);  % For PUCCH sequence generation

        [pucchIndices,pucchIndicesInfo] = nrPUCCHIndices(carrier_UL,pucch);
        dmrsIndices = nrPUCCHDMRSIndices(carrier_UL,pucch);
        dmrsSymbols = nrPUCCHDMRS(carrier_UL,pucch);
        
        harq_message = randi([0 1]);

        sr = [];

        pucchSymbols = nrPUCCH(carrier_UL,pucch,harq_message);

        pucchGrid = nrResourceGrid(carrier_UL,simParameters_UL.NTxAnts);

        pucchGrid(pucchIndices) = pucchSymbols;
        pucchGrid(dmrsIndices) = dmrsSymbols;

        txWaveform_UL = nrOFDMModulate(carrier_UL,pucchGrid);

        rxWaveform_UL = txWaveform_UL;

        SNR = 10^(SNRdB/20);
        N0 = 1/(sqrt(2.0*simParameters_UL.NRxAnts*waveformInfo_UL.Nfft)*SNR);
        noise = N0*complex(randn(size(rxWaveform_UL)),randn(size(rxWaveform_UL)));
        rxWaveform_UL = rxWaveform_UL + noise;

        rxGrid_UL = nrOFDMDemodulate(carrier_UL,rxWaveform_UL);
        [K,L,R] = size(rxGrid_UL);
        if (L < carrier_UL.SymbolsPerSlot)
            rxGrid_UL = cat(2,rxGrid_UL,zeros(K,carrier_UL.SymbolsPerSlot-L,R));
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
            error_count = error_count+1; 
        end
    end
    BLER(snrIdx) = error_count / NSubFrames;
end
