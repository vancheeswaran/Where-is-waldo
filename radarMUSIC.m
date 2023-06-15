clc;
clear;
close all;

%% Parameters 
% Downsample the number of samples in the chirp to reduce Rxx matrix size
% and therefore the complexity - compromise the range which is already
% quite high; Also helps to perform averaging across different samples of
% the same frequency spacing
freqDownSampling = 10;
% Downsample the number of antenna to allow for averaging across reduced
% number of virtual antenna with same spacing
spaceDownSampling = 5;
% Controls the amount of averaging that can be performed. Should not exceed
% downsampling values. 
freqAvgFactor = 10;
spaceAvgFactor = 5;

c = 299792458;                                         % Speed of light
lambda = 12.434e-3;                                    % Wavelength
thresholdMUSIC = 0.0025;                                 % Threshold for MUSIC
theta_vals = -pi/2:0.01:pi/2;                          % AoA values 
d_vals = 2:0.5:70;                                     % Range values                                  

subfolder = "trackingRadarDataImages_day2/";

%%
filenumber = 4;
subfolderPath = subfolder + filenumber + "/";
[MeasData, Cfg] = loadData(filenumber);
fc = Cfg.fStrt;                                        % Center Frequency
d = lambda / 2;                                        % Antenna Spacing
df = (Cfg.fStop - Cfg.fStrt)/Cfg.N*freqDownSampling;   % Frequency Spacing
n_sub = Cfg.N/freqDownSampling;                        % Number of Subcarriers
n_ant = length(Cfg.AntIdx) - spaceDownSampling;        % Number of Antennas corrected for downsampling
saveFlag = 0;
plotFlag = 1;
frameByFrame = 1;

% Get the 2D steering matrix for the given AoA and Range values
S = getSteeringMatrix2D(theta_vals, d_vals, n_sub, fc, df, d, n_ant);

processMUSIC(MeasData, S, thresholdMUSIC, d_vals, theta_vals, Cfg, n_sub, n_ant, freqDownSampling, spaceDownSampling, freqAvgFactor, spaceAvgFactor, subfolderPath, saveFlag, plotFlag, frameByFrame);

%% Save Images
for filenumber = 1:10
    filenumber
    subfolderPath = subfolder + filenumber + "/";
    mkdir(subfolderPath);
    [MeasData, Cfg] = loadData(filenumber);
    fc = Cfg.fStrt;                                        % Center Frequency
    d = lambda / 2;                                        % Antenna Spacing
    df = (Cfg.fStop - Cfg.fStrt)/Cfg.N*freqDownSampling;   % Frequency Spacing
    n_sub = Cfg.N/freqDownSampling;                        % Number of Subcarriers
    n_ant = length(Cfg.AntIdx) - spaceDownSampling;
    saveFlag = 1;
    plotFlag = 0;
    frameByFrame = 0;

    if filenumber == 1
        % Get the 2D steering matrix for the given AoA and Range values
        S = getSteeringMatrix2D(theta_vals, d_vals, n_sub, fc, df, d, n_ant);
    end
    
    processMUSIC(MeasData, S, thresholdMUSIC, d_vals, theta_vals, Cfg, n_sub, n_ant, freqDownSampling, spaceDownSampling, freqAvgFactor, spaceAvgFactor, subfolderPath, saveFlag, plotFlag, frameByFrame);
end

%%
function [MeasData, Cfg] = loadData(filenumber)
    folder = "trackingRadarDataJACOBSHALL_day2/";
    radar_filename_meta = folder + "radar_test" + filenumber + ".mat";
    radar_filename = folder + "radar_test" + filenumber + ".h5";
    load(radar_filename_meta);
    Data = h5read(radar_filename, "/Data"); 

    frameLevelData = reshape(Data, Cfg.N*numel(Cfg.TxSeq)*Cfg.Np, [], 8);
    frameLevelData = permute(frameLevelData, [1 3 2]);
    clear Data;

    tx1_idx = [];
    tx2_idx = [];
    for i=1:2:Cfg.Np*numel(Cfg.TxSeq)
        tx1_idx = [tx1_idx 1000*(i-1)+1:1000*(i-1)+1000];
        tx2_idx = [tx2_idx 1000*i+1:1000*i+1000];
    end

    chirpLevelData = [frameLevelData(tx1_idx, :, :) frameLevelData(tx2_idx, :, :)];
    clear frameLevelData;
    MeasData = chirpLevelData(:, Cfg.AntIdx, :);
    clear chirpLevelData;
    c0 = 299792458;
end

%%
function A = processMUSIC(MeasData, S, thresholdMUSIC, d_vals, theta_vals, Cfg, n_sub, n_ant, freqDownSampling, spaceDownSampling, freqAvgFactor, spaceAvgFactor, subfolderPath, saveFlag, plotFlag, frameByFrame)
    for frame_num = 1:Cfg.NFrames
        frame_num
        Data = MeasData(:, :, frame_num);
        % Reshape and permute Data to dimension: N_samples x N_antenna x N_chirps
        DataFrame = reshape(Data, 1000, [], 15);
        DataFrame = permute(DataFrame, [1 3 2]);
        rangeData = DataFrame(2:1000, :, :);
        % Obtain the IQ samples for the first chirp 
        rawData = squeeze(rangeData(:, :, 1));
        rawData = [rawData; zeros(1,15)];
        A = getAveragedData(rawData, n_sub, n_ant, freqDownSampling, spaceDownSampling, freqAvgFactor, spaceAvgFactor);
        P = estimateMUSIC(A, S, thresholdMUSIC, d_vals, theta_vals);
        
        if saveFlag || plotFlag
            if plotFlag
                figure(1);
            else 
                f=figure;
                set(f, 'visible', 'off');
            end

            surf(theta_vals*(180/pi), d_vals/2, db(P));
            title('MUSIC Range-AoA', 'FontSize', 15, 'FontWeight','bold');
            xlabel('Angle of Arrival', 'FontSize', 14);
            ylabel('Range (m)', 'FontSize', 14);
            shading flat;
            view(0, 90);
            colorbar;
            
            if frameByFrame
                w = waitforbuttonpress;
            end
       
%             colormap("gray");
%             clim([-63.4 -63.35]);
    
            if saveFlag
                filename = subfolderPath + "rangeAoA_frame" + frame_num + ".png";
                saveas(f, filename);
            end
        end
    end
end


%%
function P = estimateMUSIC(A, S, thresholdMUSIC, d_vals, theta_vals)
    % Compute the covariance matrix
    R = A*A';
    % Find eigenvalues and eigenvectors of the covariance matrix
    [V, D] = eig(R);
    % Extract the eigenvalues
    eig_vals = diag(D);
    % Apply threshold on the eigenvalues to identify the noise space
    idx = find(eig_vals<thresholdMUSIC*max(abs(eig_vals)));
    % Apply MUSIC estimator function
    B = 1./vecnorm((V(:, idx))'*S, 2).^2;
    % Reshape to get Range-AoA profile
    P = reshape(B, length(d_vals), length(theta_vals));
end

%% Function to obtain the 2D steering matrix
function S = getSteeringMatrix2D(theta_vals, d_vals, N, f, df, d, N_Ant)
    % N     -  number of subcarriers
    % f     -  center frequency
    % df    -  subcarrier spacing
    % d     -  antenna spacing
    % N_Ant -  number of antennas

    Phi = zeros(N_Ant, length(theta_vals));
    for i = 1:length(theta_vals)
        phi_theta(:, i) = exp(-1i*2*pi*f/3e8*sin(theta_vals(i))*d*(0:N_Ant-1));
    end

    Omega = zeros(N, length(d_vals));
    for i = 1:length(d_vals)
        omega_t = exp(-1i*2*pi*df*d_vals(i)/3e8);
        Omega(:, i) = omega_t.^((1:N)');
    end

    S = kron(phi_theta, Omega);
end

%% Function to perform averaging on Raw Data
function A = getAveragedData(rawData, n_sub, n_ant, freqDownSampling, spaceDownSampling, freqAvgFactor, spaceAvgFactor)
    assert(freqDownSampling >= freqAvgFactor, "Frequency Average Factor should not exceed the Down-Sampling value");
    assert(spaceDownSampling >= spaceAvgFactor, "Space Average Factor should not exceed the Down-Sampling value");
    A = zeros(n_sub * n_ant, freqAvgFactor * spaceAvgFactor);
    for f = 1:freqAvgFactor
        for s = 1:spaceAvgFactor
            tempA = rawData(1:freqDownSampling:end, s:(n_ant + s - 1)); 
            A(:, f*s) = tempA(:);
        end
    end

end
