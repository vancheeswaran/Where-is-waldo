clc;
clear;
close all;

% Specify the radar data to be processed
filenumber = 10;
folder = "trackingRadarDataJACOBSHALL_day2/";
radar_filename_meta = folder + "radar_test" + filenumber + ".mat";
radar_filename = folder + "radar_test" + filenumber + ".h5";

load(radar_filename_meta);

%% Data Preprocessing

Data = h5read(radar_filename, "/Data"); %load data
% Formatting the data
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

%% Range-AoA
cap = rangeAoACancellation(MeasData, Cfg, c0);

%% Range-AoA dB
cap = rangeAoACancellationdB(MeasData, Cfg, c0);

%% Cancellation Range-AoA 
function [cancelledAoAProfile] = rangeAoACancellation(MeasData, Cfg, c0)
    prevProfile = zeros(204, Cfg.NFFTAnt);
    cancelledAoAProfile = zeros(204, Cfg.NFFTAnt, Cfg.NFrames);
    for i=1:Cfg.NFrames
        Data = MeasData(:, :, i);
        DataFrame = reshape(Data, 1000, [], 15);
        DataFrame = permute(DataFrame, [1 3 2]);
        rangeData = DataFrame(2:1000, :, :);

        rangeProfile = fft(rangeData.*Cfg.Win2D.*Cfg.mCalData, Cfg.NFFT, 1).*Cfg.radarsca/Cfg.ScaWin;
        rangeProfile = permute(rangeProfile, [1 3 2]);
        aoaData = squeeze(rangeProfile(:, 1, :));
        aoaProfile = fftshift(fft(aoaData(1:204, :), Cfg.NFFTAnt, 2), 2);
        aoaProfileC = aoaProfile - prevProfile;
        aoa_db = db(aoaProfile);
        aoa_db(aoa_db<-80) = -80;

        figure(1);
        range_axis = (1:Cfg.NFFT)*(Cfg.N/Cfg.NFFT)*c0/(2*(Cfg.fStop-Cfg.fStrt));
        surf(Cfg.vAngDeg', range_axis(1:204), abs(aoaProfile));
        title('Range-AoA without cancellation');
        shading flat;
        view(0, 90);
        colorbar;

        saveas(gcf, "range_AoA_without_cancellation")

        figure(2);
        range_axis = (1:Cfg.NFFT)*(Cfg.N/Cfg.NFFT)*c0/(2*(Cfg.fStop-Cfg.fStrt));
        surf(Cfg.vAngDeg', range_axis(1:204), abs(aoaProfileC));
        title('Range-AoA with cancellation');
        shading flat;
        view(0, 90);
        colorbar;

        if(mod(i,4))
            prevProfile = aoaProfile;
        else
            prevProfile = prevProfile;
        end
        cancelledAoAProfile(:, :, i) = aoaProfileC;
    end
end

%% Cancellation Range-AoA dB
function [cancelledAoAProfile] = rangeAoACancellationdB(MeasData, Cfg, c0)
    prevProfile = zeros(204, Cfg.NFFTAnt);
    cancelledAoAProfile = zeros(204, Cfg.NFFTAnt, Cfg.NFrames);
    for i=1:Cfg.NFrames
        Data = MeasData(:, :, i);
        DataFrame = reshape(Data, 1000, [], 15);
        DataFrame = permute(DataFrame, [1 3 2]);
        rangeData = DataFrame(2:1000, :, :);

        rangeProfile = fft(rangeData.*Cfg.Win2D.*Cfg.mCalData, Cfg.NFFT, 1).*Cfg.radarsca/Cfg.ScaWin;
        rangeProfile = permute(rangeProfile, [1 3 2]);
        aoaData = squeeze(rangeProfile(:, 1, :));
        aoaProfile = fftshift(fft(aoaData(1:204, :), Cfg.NFFTAnt, 2), 2);
        aoaProfileC = aoaProfile - prevProfile;
        aoa_db = db(aoaProfile);
        aoa_db(aoa_db<-80) = -80;

        figure(1);
        range_axis = (1:Cfg.NFFT)*(Cfg.N/Cfg.NFFT)*c0/(2*(Cfg.fStop-Cfg.fStrt));
        surf(Cfg.vAngDeg', range_axis(1:204), aoa_db);
        caxis([-80, -45]);
        title('Range-AoA without cancellation');
        shading flat;
        view(0, 90);
        colorbar;

        saveas(gcf, "range_AoA_without_cancellation")

        figure(2);
        range_axis = (1:Cfg.NFFT)*(Cfg.N/Cfg.NFFT)*c0/(2*(Cfg.fStop-Cfg.fStrt));
        surf(Cfg.vAngDeg', range_axis(1:204), db(aoaProfileC));
        title('Range-AoA with cancellation');
        shading flat;
        view(0, 90);
        colorbar;

        if(mod(i,4))
            prevProfile = aoaProfile;
        else
            prevProfile = prevProfile;
        end
        cancelledAoAProfile(:, :, i) = aoaProfileC;
    end
end
