close all
clear all

targetdir = '/Volumes/PS2Akermanlab/GG/Data';

% folder for each condition
P60_folders = {'P60 control'; 'P60 test'};
P21_folders = {'P21 control/'; 'P21 test/'; 'Chronos/'};
   

% which experiments do i want to analyse?
% 1 = control
% 2 = test
% 3 = chronos

for i = 1
    
    
    dataFolder = [targetdir P60_folders{i}];
    
    myFolder          = dir(dataFolder);
    qremove           = ismember({myFolder.name},{'.','..','.DS_Store'}); % locate where the unwanted names are
    myFolder(qremove) = []; % set these names to 0
    animal_number     = length(myFolder);
    

    
    %% Set up variables
    
    layers = 1:32;
    layers_idx = 1:length(layers);

    
    % how many spontaneous windows either side do i want to normalise to?
    spon_n = 5;
    
    timecount       = linspace(0,3,3000);
    
    response_start  = 1.005;
    response_end    = 1.030;
    response_window = response_end-response_start;
    
    stim_idx        = timecount >= response_start & timecount < response_end; % find spikes which are x msec after stimulus
    burst_idx       = timecount >= 0.975 & timecount < 1;
    
    spon1_idx       = timecount >= 0.025 & timecount < 1;
    spon1_burst_idx = timecount >= 0 & timecount < 0.025;
    
    idx = [burst_idx; spon1_burst_idx];
    
    % pad vectors with NaNs
    stim_response_all = NaN(31,14);
    trial_counter = zeros(31,14);
    spon1_response_all = NaN(31,14);
    
    
    %% Load data
    for k = 1:length(myFolder)
        
    folderName = myFolder(k).name
    fileName        = dir([dataFolder folderName '/*.mat']);
    
    for a = 1:length(fileName)
        
        
        data        = load([dataFolder folderName filesep fileName(a).name]);
        
        %% Whisk response - all layers
        
        spikecount      = histc(data.ephys_data.conditions.spikes(:,:,:),timecount,3);
        single_trial_length = size(spikecount,2);

        %% Calculate spiking response
        
        % spike count for each trial and each channel (e.g. 32x30) in Hz
        shortspike = median(sum(spikecount(layers_idx,:,stim_idx),3));
        spon1_shortspike = median(sum(spikecount(layers_idx,:,spon1_idx),3));
        
     
        % sum over channels per trial (e.g. 1x30)
        stim_response_all(1:single_trial_length,a) = shortspike/response_window;
        spon1_response_all(1:single_trial_length,a) = spon1_shortspike/0.975;
        
        
    end
    
    %% Get spikes into one vector for each animal
    
    % reshape to get one column for each animal
    stim_data = reshape(stim_response_all,[],1);
    
    spon_data = reshape(spon1_response_all,[],1);

    
    %% STATS
    
    [p h] = ranksum(spon_data,stim_data);
    
    signal(k) = nanmean(stim_data) / nanmean(spon_data)
    
    
    if signal > 1.5
        response = 'yes :D'

    else
        response = 'no :('

    end 
    
    %% PLOTTING
    
    figure(1)
    % mean spikes before and after induction
    sponvsstim = bar([nanmean(spon_data) nanmean(stim_data)],'FaceColor',[0.0784313753247261 0.168627455830574 0.549019634723663])
    fixplot
    hold on
    errorbar([nanmean(spon_data) nanmean(stim_data)],[nanstd(spon_data)/sqrt(60) nanstd(stim_data)/sqrt(60)],'LineStyle','none','LineWidth',1.5,'Color',[0.600000023841858 0.600000023841858 0.600000023841858])
    xlabel('Spon spikes             Stim spikes')
    ylabel('Mean spike rate')
    ylim([0 140])
    
%     set(gcf,'PaperUnits','inches','PaperPosition',fig_size)
%  
%     % Optional saving
%     if save_fig
%         save_file   = fullfile(save_dir, ['Spon vs stim spikes ' fig_title]);
%         print(gcf,save_file,fig_format)
%     end
    
    end 
    
end
    