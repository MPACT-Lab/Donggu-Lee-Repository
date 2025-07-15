clear; clc;

SNRdB = -20:0.5:10;              % SNR range
NSubframes = 10^3;              % Number of slots/subframes per SNR
nmse_lte = zeros(size(SNRdB));
nmse_5g  = zeros(size(SNRdB));

nTxAnts = 1; nRxAnts = 1;
nSizeGrid = 6;                 
nCellID = 9;
cpType = 'Normal';
scSpacing = 15;

for idx = 1:length(SNRdB)
    snr_db = SNRdB(idx);
    snr = 10^(snr_db/10);
    
    nmse_lte_temp = zeros(NSubframes,1);
    nmse_5g_temp  = zeros(NSubframes,1);

    rmse_5g_numerator = 0;
    rmse_lte_numerator = 0;

    for n = 1:NSubframes
        disp('-------------------------')
        disp(idx)
        disp(n)
        disp('-------------------------')

        carrier = nrCarrierConfig;
        carrier.NCellID = nCellID;
        carrier.NSizeGrid = nSizeGrid;
        carrier.SubcarrierSpacing = scSpacing;
        carrier.CyclicPrefix = cpType;
        carrier.NSlot = mod(n, 160);
        
        pucch = nrPUCCH1Config;
        pucch.PRBSet = 0;
        pucch.SymbolAllocation = [0 14];
        pucch.FrequencyHopping = 'intraSlot';
        pucch.GroupHopping = 'neither';
        pucch.HoppingID = 0;
        pucch.SecondHopStartPRB = 5;

        [pucchInd,pucchInfo] = nrPUCCHIndices(carrier, pucch);
        dmrsInd = nrPUCCHDMRSIndices(carrier, pucch);
        dmrsSym = nrPUCCHDMRS(carrier, pucch);
        bits = randi([0 1],1);
        symbols = nrPUCCH(carrier, pucch, bits);

        nrgrid = nrResourceGrid(carrier, nTxAnts);
        nrgrid(pucchInd) = symbols;
        nrgrid(dmrsInd) = dmrsSym;

        txWave = nrOFDMModulate(carrier, nrgrid);
        waveformInfo = nrOFDMInfo(carrier);
        N0 = 1 / sqrt(nRxAnts * waveformInfo.Nfft * snr);
        noise = N0 * complex(randn(size(txWave)), randn(size(txWave)));
        rxWave = txWave + noise;

        rxGrid = nrOFDMDemodulate(carrier, rxWave);
        [Hest, ~] = nrChannelEstimate(carrier, rxGrid, dmrsInd, dmrsSym);
        Htrue = ones(size(Hest));  

        validInd = unique([pucchInd(:); dmrsInd(:)]);

        err5g = abs(Hest(validInd) - Htrue(validInd)).^2;
        rmse_5g_numerator = rmse_5g_numerator + sum(err5g(:));

        ue = struct;
        ue.NULRB = nSizeGrid;
        ue.NTxAnts = nTxAnts;
        ue.CyclicPrefixUL = cpType;
        ue.NCellID = nCellID;
        ue.NSubframe = mod(n-1, 10);

        pucch_lte = struct;
        pucch_lte.ResourceSize = 0;
        pucch_lte.DeltaShift = 1;
        pucch_lte.CyclicShifts = 0;
        pucch_lte.ResourceIdx = 0;

        grid_lte = lteULResourceGrid(ue);
        ACK = 1;
        pucchSym = ltePUCCH1(ue, pucch_lte, ACK);
        drsSym = ltePUCCH1DRS(ue, pucch_lte);
        pucchInd_lte = ltePUCCH1Indices(ue, pucch_lte);
        drsInd_lte = ltePUCCH1DRSIndices(ue, pucch_lte);
        grid_lte(pucchInd_lte) = pucchSym;
        grid_lte(drsInd_lte) = drsSym;
        txWave_lte = lteSCFDMAModulate(ue, grid_lte);

        scfdmaInfo = lteSCFDMAInfo(ue);
        N0 = 1 / (snr * sqrt(double(scfdmaInfo.Nfft))) / sqrt(2 * nTxAnts);
        noise = N0 * complex(randn(size(txWave_lte)), randn(size(txWave_lte)));
        rxWave_lte = txWave_lte + noise;
        rxGrid = lteSCFDMADemodulate(ue, rxWave_lte);

        cec = struct;
        cec.PilotAverage = 'UserDefined';
        cec.FreqWindow = 12;
        cec.TimeWindow = 1;
        cec.InterpType = 'cubic';
        [Hest_lte, ~] = lteULChannelEstimatePUCCH1(ue, pucch_lte, cec, rxGrid);
        Htrue_lte = ones(size(Hest_lte));

        validInd_lte = unique([pucchInd_lte(:); drsInd_lte(:)]);

        errLte = abs(Hest_lte(validInd_lte) - Htrue_lte(validInd_lte)).^2;
        rmse_lte_numerator = rmse_lte_numerator + sum(errLte(:));

        count_5g = numel(err5g) * NSubframes;
        count_lte = numel(errLte) * NSubframes;
    end

    nmse_5g(idx) = sqrt(rmse_5g_numerator / count_5g);
    nmse_lte(idx) = sqrt(rmse_lte_numerator / count_lte);
end
