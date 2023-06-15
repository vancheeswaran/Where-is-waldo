% Description:
% Configure FMCW Mode with sequential activation of Tx antennas and measure upchirp IF signal
%
% (1) Connect to Radarbook2 with ADF24 Frontend
% (2) Enable Supply
% (3) Configure RX
% (4) Configure TX
% (5) Start Measurements with TX sequence
% (6) Configure signal processing
% (7) Calculate DBF algorithm

clc;
clear;
close all;

pause(5);

filenumber = 10;
folder = "trackingRadarDataJACOBSHALL_day2/";
mkdir(folder);
filename = folder + "radar_test" + filenumber + ".mat";
filename1 = folder + "radar_test" + filenumber + ".h5";
% Configure script
Disp_FrmNr = 1;
Disp_TimSig = 0;      % display time signals
Disp_JOpt = 1;      % display cost function for DBF

%--------------------------------------------------------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
CurPath = pwd();
addpath([CurPath,'/../../Class']);
addpath([CurPath,'/../../PNet'])

%--------------------------------------------------------------------------
% Define Constants
%--------------------------------------------------------------------------
c0 = 299792458; 

Brd = Rbk2Adf24Tx2Rx8('PNet', '192.168.1.1');
% Verify if sampling framework is installed
% Brd.BrdChkSocSysId();

%--------------------------------------------------------------------------
% Reset Board and Enable Power Supply
%--------------------------------------------------------------------------
Brd.BrdRst();
Brd.BrdPwrEna();

%--------------------------------------------------------------------------
% Software Version
%--------------------------------------------------------------------------
% Brd.BrdDispSwVers();

%--------------------------------------------------------------------------
% Status Information
%--------------------------------------------------------------------------
Brd.BrdDispInf();


% Configure Board as Master and set ADC clock to 40 MHz
Brd.BrdSetRole('Ms', 40e6);

%--------------------------------------------------------------------------
% Load Calibration Data for the entire array
%--------------------------------------------------------------------------
CalCfg.Mask = 1;
CalCfg.Len = 32;
CalData = Brd.BrdGetCalData(CalCfg);

%--------------------------------------------------------------------------
% Enable Receive Chips
%--------------------------------------------------------------------------
Brd.RfRxEna();
Brd.RfTxEna(1, 80);

%--------------------------------------------------------------------------
% Configure AFE5801
%--------------------------------------------------------------------------
Brd.Set('AfeLowNoise',0);
% Enable/Disable internal DC coupling
Brd.Set('AfeIntDcCoupling',0);
% Ramp pattern can be used to test communication: sample data is replaced
% by linarly increasing ramp
Brd.Set('AfePatRamp','Off');
% Set AfeGain in dB (-5 - 30 dB); The closest available value is configured
Brd.Set('AfeGaindB', 20);

%--------------------------------------------------------------------------
% Configure Up-Chirp with a sequence
% SPI Command is required to configure TX to different antenna the spi
% command takes < 10 us;
% 
%--------------------------------------------------------------------------
Cfg.fStrt = 24.0e9; 
Cfg.fStop = 24.2e9;
Cfg.TRampUp = 120e-6;   % Chirp Time
Cfg.TInt = 100e-3;      % Inter-frame interval??
Cfg.Tp = 1e-3;        % Inter-chirp time??
Cfg.N = 1000;           % Samples per chirp
Cfg.IniTim = 100e-3;    % Check if delay is consistent               
Cfg.IniEve = 0;         % Start automatically after IniTim
Cfg.TxSeq = [1, 2];     % Activate 1 and then antenna2
Cfg.Np = 50;            % Number of chirps in one frame
NFrames = 200;

%--------------------------------------------------------------------------
% Configure DMA Transfer to copy numel(Cfg.TxSeq)*Np frames simultaneously.
% Required to achiev maximal data transfer between FPGA and Soc 
Brd.Set('DmaMult', numel(Cfg.TxSeq).*Cfg.Np*NFrames);
% Brd.Set('DmaMult', Cfg.Np);

Brd.RfMeas('ExtTrigUp_TxSeq',Cfg);

%--------------------------------------------------------------------------
% Read Settings for N and fs
%--------------------------------------------------------------------------
N = Brd.Get('N');
fs = Brd.Get('fs');
fc = Brd.RfGet('fc');

%--------------------------------------------------------------------------
% Configure Signal Processing
%--------------------------------------------------------------------------
% Processing of range profile
NFFT = 2^12;
NrChn = Brd.Get('NrChn');
Win2D = Brd.hanning(N-1,2*NrChn-1);
ScaWin = sum(Win2D(:,1));
kf = Brd.RfGet('kf');
vRange = [0:NFFT-1].'./NFFT.*fs.*c0/(2.*kf);
NFFTVel = 2^10;

B = Cfg.fStop-Cfg.fStrt;
k = B/Cfg.TRampUp;
range_res = c0/(2*B);
range_max = (fs*c0)/(2*k);

range_axis = (0:NFFT-1)*range_max/NFFT;
WinVel = Brd.hanning(Cfg.Np); %creating hann window for range doppler length Cfg.Np (1 for each chirp)
ScaWinVel = sum(WinVel); %used for scaling range doppler fft
WinVel2D = repmat(WinVel.',numel(range_axis),1);

% Configure range interval to be displayed
RMin = 1;
RMax = 10;

[Val RMinIdx] = min(abs(vRange - RMin));
[Val RMaxIdx] = min(abs(vRange - RMax));
vRangeExt = vRange(RMinIdx:RMaxIdx);

% Window function for receive channels
NFFTAnt = 256;
WinAnt = Brd.hanning(2*NrChn-1);
ScaWinAnt = sum(WinAnt);
WinAnt2D = repmat(WinAnt.',numel(vRangeExt),1);
vAngDeg = asin(2*[-NFFTAnt./2:NFFTAnt./2-1].'./NFFTAnt)./pi.*180;
freqVel_axis = (-NFFTVel/2:NFFTVel/2-1).'./NFFTVel.*(1/Cfg.Tp);
vel_axis = freqVel_axis*c0/(2.*fc); 

% Select virtual channels for processing
% Remove the overlapping channel
AntIdx = [1:8,10:16].';

% Calibration data for the selected channesl
CalData = CalData(AntIdx);
mCalData = repmat(CalData.',N-1,1);

% Positions for polar plot of cost function
vU = linspace(-1,1,NFFTAnt);
[mRange , mU] = ndgrid(vRangeExt,vU);
mX = mRange.*mU;
mY = mRange.*cos(asin(mU));

fprintf("Starting data collection........\n");

t1 = datetime('now','Format',"yyyy-MM-dd HH:mm:ss:SSSSSS")
% Record data for Tx1 and Tx2
Data = Brd.BrdGetData(Cfg.Np*numel(Cfg.TxSeq)*NFrames); 
t2 = datetime('now','Format',"yyyy-MM-dd HH:mm:ss:SSSSSS")
t = [t1 t2];

Brd.BrdRst();
Brd.BrdPwrDi();

fprintf("Data collection done\n");

%% Data collection
% (Add more params for further processing if required, only adding params required for vendor provided processing)
Cfg.NFFT = NFFT;
Cfg.NrChn = NrChn;
Cfg.Win2D = Win2D;
Cfg.ScaWin = ScaWin;
Cfg.kf = kf;
Cfg.vRange = vRange;
Cfg.RMin = RMin;
Cfg.RMax = RMax;
Cfg.RMinIdx = RMinIdx;
Cfg.RMaxIdx = RMaxIdx;
Cfg.vRangeExt = vRangeExt;
Cfg.NFFTAnt = NFFTAnt;
Cfg.WinAnt = WinAnt;
Cfg.ScaWinAnt = ScaWinAnt;
Cfg.WinAnt2D = WinAnt2D;
Cfg.vAngDeg = vAngDeg;
Cfg.AntIdx = AntIdx;
Cfg.CalData = CalData;
Cfg.mCalData = mCalData;
Cfg.radarsca = Brd.FuSca;
Cfg.mX = mX;
Cfg.mY = mY;
Cfg.N = N;
Cfg.fs = fs;
Cfg.NFFTVel = NFFTVel;
Cfg.ScaWinVel = ScaWinVel;
Cfg.WinVel2D = WinVel2D;
Cfg.vel_axis = vel_axis;
Cfg.NFrames = NFrames;
Cfg.B = B;

%% Experiment parameters
fprintf("Saving Meta Data....\n");
tic
save(filename, 'Cfg', 't', '-v6');
toc
fprintf("Saved Meta data\n");

%%
fprintf("Saving Data....\n");
tic
h5create(filename1,"/Data",[Cfg.Np*numel(Cfg.TxSeq)*NFrames*Cfg.N 8]);
h5write(filename1, "/Data", Data);
toc
fprintf("Saved data\n");

%%
clear Brd;
