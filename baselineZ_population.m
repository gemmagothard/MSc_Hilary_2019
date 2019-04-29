close all
clear all

targetdir = '/Volumes/PS2Akermanlab/GG/Data';

% folder for each condition
P60_folders = {'P60 control'; 'P60 test'};
P21_folders = {'P21 control/'; 'P21 test/'; 'Chronos/'};

% what trial number does the induction occur at? (30 /120 for P60 and
% 60/180 for P21)
ind = 30;
end_trial = 120;


% Saving options for output figures (Note - will overwrite files of the same name if script is re-run) 
save_fig        = false; % If false, no figures are automatically saved; if true, figures are saved in directory and format specified below
save_dir        = '/Users/gemmagothard/Documents/OXFORD/ROTATION1/RESULT FIGURES/figure 2 bursty/';
fig_format      = '-depsc'; % '-depsc' for vector graphics or '-dpng' for png image file

% which experiments do i want to analyse?
% 1 = P60 control
%      '2019_01_31'
%      '2019_02_01'
%      '2019_02_13'
% 2 = P60 test
%     '2019_01_08'
%     '2019_01_11'
%     '2019_02_08'
%     '2019_02_15'
% 1 = P21 control
%     '2019_03_30'
% 2 = P21 test
%      '2019_03_27'
%      '2019_03_28'
%      '2019_03_29'
%      '2019_04_01'

for i = 1
    
    
    
    dataFolder = [targetdir P60_folders{i}];
    
    myFolder          = dir(dataFolder);
    qremove           = ismember({myFolder.name},{'.','..','.DS_Store'}); % locate where the unwanted names are
    myFolder(qremove) = []; % set these names to 0
    animal_number     = length(myFolder);

    population_bursts = NaN(180,32,length(myFolder));
    plot_bursts = NaN(10,180);
    
    track_burst_spikes = [];
    
    for k = 1:length(myFolder)
        
        folderName = myFolder(k).name
        fig_title = folderName;
        
    %% Set up variables
    
    layers = 1:32;
    layers_idx = 1:length(layers);
    
    
    
    n = 1; % one plotted point is n number of data points
    
    
    % how many spontaneous windows either side do i want to normalise to?
    spon_n = 5;
    
    timecount       = linspace(0,3,3000);
    
    response_start  = 1.005;
    response_end    = 1.030;
    response_window = response_end-response_start;
    
    stim_idx        = timecount >= response_start & timecount < response_end; % find spikes which are x msec after stimulus
    burst_idx       = timecount >= 0.975 & timecount < 1;
    burst_threshold = 2;
    
    spon1_idx       = timecount >= 0.04 & timecount < 0.08;
    spon1_burst_idx = timecount >= 0.0375 & timecount < 0.04;
    
    spon2_idx       = timecount >= 0.54 & timecount < 0.58;
    spon2_burst_idx = timecount >= 0.5375 & timecount < 0.54;
    
    
    idx = [burst_idx; spon1_burst_idx; spon2_burst_idx];
    
    % pad vectors with NaNs
    stim_response_all = NaN(31,14);
    trial_counter = zeros(31,14);
    spon1_response_all = NaN(31,14);
    spon2_response_all = NaN(31,14);
    P2P = NaN(31,14);
    
    all_exp_bursts = NaN(6,32,30);
    burst_tracker = NaN(30,6);
    
    %% Load data
    
    
    fileName        = dir([dataFolder folderName '/*.mat']);
    
    for a = 1:length(fileName)
        
        
        data        = load([dataFolder folderName filesep fileName(a).name]);
        
        %% Whisk response - all layers
        
        
        spikecount      = histc(data.ephys_data.conditions.spikes(layers,:,:),timecount,3);
        LFPtimeframe = data.ephys_data.conditions.LFP_trace(:,:,1000:1150);
        
        single_trial_length = size(spikecount,2);
        trial_length(a)    = single_trial_length;
        trial_counter(1:trial_length,a) = 1;
        
        
        %% QC the data
        
        for m = 1:3
            burst_channels = false(32,single_trial_length);
            spikecount      = histc(data.ephys_data.conditions.spikes(layers,:,:),timecount,3);
            
            for h = layers_idx
                % How many spikes happen 25ms before stimulus?
                burst_spikes     = sum(spikecount(h,:,idx(m,:)),3);
                    
                    if m == 1
                    track_burst_spikes = [track_burst_spikes; burst_spikes(:)];
                    end 
                
                for c = 1:single_trial_length

                    if burst_spikes(c) >= burst_threshold
                        burst_channels(h,c) = true;

                    else
                        burst_channels(h,c) = false;
                
                    end
                    
                    
                    
                end
                
            end
            
        % get the total amount of channels which contain bursts in each trial (1x30)
        all_bursts = sum(burst_channels);
        burst_tracker(1:length(all_bursts),a) = all_bursts;
        
                if m == 1
                    % set any trials with prestim bursts in any channels to NaN
                    spikecount(:,all_bursts>0,:) = NaN;
                    shortspike = sum(spikecount(layers_idx,:,stim_idx),3);
                    LFPtimeframe(:,all_bursts>0,:) = NaN;
                end
                
                if m == 2
                    % set any trials with bursts in any channels to NaN
                    % window to NaN
                    spikecount(:,all_bursts>0,:) = NaN;
                    spon1_shortspike = sum(spikecount(layers_idx,:,spon1_idx),3);
                end
                
                if m == 3
                    % set any trials with bursts in any channels to NaN
                    % window to NaN
                    spikecount(:,all_bursts>0,:) = NaN;
                    spon2_shortspike = sum(spikecount(layers_idx,:,spon2_idx),3);
                end
                
                


      
        end
        
        
        
        %% Calculate spiking response

        % Total spike count for each trial (median over channels)
        stim_response_all(1:single_trial_length,a) = nanmean(shortspike);
        spon1_response_all(1:single_trial_length,a) = nanmean(spon1_shortspike);
        spon2_response_all(1:single_trial_length,a) = nanmean(spon2_shortspike);

    end
    
    
    
    
    %% Get spikes into one vector for each animal
   
    
    all_trials = sum(trial_length);
    
    % all burst trials is how many channels had bursts for each trial (e.g. 1x180)
    all_burst_trials = reshape(burst_tracker,[],1);
    
    % Get rid of the padded NaNs
    all_burst_trials(isnan(all_burst_trials)) = [];
    
    %this tells you what percentage of trials were removed in each animal
    exp_bursts = (sum(all_burst_trials~=0) / all_trials)*100
    
    exp_bursts_all(k) = exp_bursts;

    % this is how many bursts were found in each channel in each trial of the experiment
    all_exp_bursts2 = reshape(all_exp_bursts,[],32);

    population_bursts(1:length(all_exp_bursts2),:,k) = all_exp_bursts2;
    
    plot_bursts(1:length(all_burst_trials),k) = all_burst_trials;
    
    

    
    % reshape to get one column for each animal
    plot_data = reshape(stim_response_all,[],1);
    spon1_data = reshape(spon1_response_all,[],1);
    spon2_data = reshape(spon2_response_all,[],1);
    trial_counter = reshape(trial_counter,[],1);
    
    
    
    %remove padded NaNs
    plot_data(trial_counter==0) = [];
    
    % spontaneous spikes without spontaneous activity removed
    not_sub_spon = plot_data;
    
    %% Sort out spontaneous data into one vector
    
    spon_data = [spon1_data spon2_data];
    
    % mean across the two windows for each trial
    spon_data = nanmean(spon_data,2);
    
    %% Get whisker response relative to spontaneous spikes
    
    for u = 1:length(plot_data)-spon_n
        
        if u <= spon_n
            spon_spikes = nanmedian(spon_data(1:(spon_n*2)+1));
        else
            spon_spikes = nanmedian(spon_data(u-spon_n:u+spon_n));
        end
        
        plot_data(u) = plot_data(u) - spon_spikes;
        
        
    end
    
    baseline = plot_data(1:ind);

    
    % get median and variance for first and last 25% of pre trials
    median_1 = nanmedian(baseline(1:ind/2));
    iqr_1 = iqr(baseline(1:ind/2));
    median_2 = nanmedian(baseline(ind/2+1:end));
    iqr_2 = iqr(baseline(ind/2+1:end));
    
    z(k) = abs((median_1 - median_2) / iqr_1)
    
    if abs(z) > 1
        baseline_flat = 'no :('
    else
        baseline_flat = 'yes :D'
    end 
    
%     if abs(z) > 1
%         baseline_flat = 'no :('
%     else
%         baseline_flat = 'yes :D'
%     end 
%     
    all_spikes = [track_burst_spikes ]
    
    end 

    histogram(track_burst_spikes)
    set(gca,'YScale','log')
    fixplot
    print(gcf,fullfile('/Users/gemmagothard/Documents/OXFORD/ROTATION1/RESULT FIGURES',['bursty hist']),'-depsc','-painters')

    
    % Optional saving
    if save_fig
        save_file   = fullfile(save_dir, ['Spontaneous spikes ' fig_title]);
        print(gcf,save_file,fig_format,'-painters')
    end
    
end

  