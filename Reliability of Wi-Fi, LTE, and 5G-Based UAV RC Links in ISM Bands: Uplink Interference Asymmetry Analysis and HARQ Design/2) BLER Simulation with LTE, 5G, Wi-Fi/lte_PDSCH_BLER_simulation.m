%% PUCCH Format 1a baseline
% reference example link: https://www.mathworks.com/help/lte/ug/pucch1a-ack-missed-detection-probability-conformance-test.html
% Note: You can modify the threshold in the decoding function by stepping
% into ltePUCCH1Decode function. The initial setting is 0.97585. 

clear;
clc;

NSubframes = 10^5;                          
SNRIn = linspace(-30, -5, 50);         
actual_SNR = zeros(length(SNRIn), 1);
NTxAnts = 1;                                    

ue = struct;                                
ue.NULRB = 6;                              
ue.CyclicPrefixUL = 'Normal';             
ue.Hopping = 'Off';                                            
ue.NCellID =  9;                            
ue.Shortened = 0;                           
ue.NTxAnts = NTxAnts;
                            
pucch = struct;                         
pucch.ResourceSize = 0;
pucch.DeltaShift = 1;                       
pucch.CyclicShifts = 0;
pucch.ResourceIdx = 0:ue.NTxAnts-1;

channel = struct;                       
channel.NRxAnts = 1;                    
channel.DelayProfile = 'Off';           
channel.DopplerFreq = 0;             

channel.MIMOCorrelation = 'Low';        
channel.NTerms = 16;                    
channel.ModelType = 'GMEDS';                 
channel.Seed = 13;                      
channel.InitPhase = 'Random';            
channel.NormalizePathGains = 'On';     
channel.NormalizeTxAnts = 'On';         
    
scfdmaInfo = lteSCFDMAInfo(ue);

channel.SamplingRate = scfdmaInfo.SamplingRate;   

cec = struct;                     
cec.PilotAverage = 'UserDefined'; 
cec.FreqWindow = 12;              
cec.TimeWindow = 1;                                                     
cec.InterpType = 'cubic';        
            
PMISS = zeros(size(SNRIn));
error_arr = [];

for nSNR = 1:length(SNRIn)

    missCount = 0;
    errors = [];

    offsetused = 0;

    SNR = 10^(SNRIn(nSNR)/20);
    N = 1/(SNR*sqrt(double(scfdmaInfo.Nfft)))/sqrt(2.0*ue.NTxAnts);

    for nsf = 1:NSubframes
        ACK = 1;

        disp("---------------")
        disp(nsf)
        disp(nSNR)
        
        ue.NSubframe = mod(nsf-1,10);            
        reGrid = lteULResourceGrid(ue);

        pucch1Sym = ltePUCCH1(ue,pucch,ACK);    
        pucch1DRSSym = ltePUCCH1DRS(ue,pucch); 

        pucch1Indices = ltePUCCH1Indices(ue,pucch);
        pucch1DRSIndices = ltePUCCH1DRSIndices(ue,pucch);

        reGrid(pucch1Indices) = pucch1Sym;
        reGrid(pucch1DRSIndices) = pucch1DRSSym;

        txwave = lteSCFDMAModulate(ue,reGrid);  

        channel.InitTime = (nsf-1)/1000;

        rxwave = txwave; % AWGN

        noise = N * complex(randn(size(rxwave)),randn(size(rxwave)));
        rxwave = rxwave + noise;

        rxgrid = lteSCFDMADemodulate(ue,rxwave(1+offsetused:end,:));

        [H,n0] = lteULChannelEstimatePUCCH1(ue,pucch,cec,rxgrid);

        [pucch1Rx,pucch1H] = lteExtractResources(pucch1Indices,rxgrid,H);

        eqgrid = lteULResourceGrid(ue);    
        eqgrid(pucch1Indices) = lteEqualizeMMSE(pucch1Rx,pucch1H,n0);

        rxACK = ltePUCCH1Decode(ue,pucch,length(ACK),eqgrid(pucch1Indices));  % You can modify decoding threshold step-into this line 

        if (isempty(rxACK) || any(rxACK ~= ACK))
            missCount = missCount + 1;  
        end 
        
    end
            
    PMISS(nSNR) = (missCount/NSubframes);
        
end
