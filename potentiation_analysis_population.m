close all
clear all

targetdir = '[target folder]';

% Saving options for output figures (Note - will overwrite files of the same name if script is re-run) 
save_fig        = false; % If false, no figures are automatically saved; if true, figures are saved in directory and format specified below
save_dir        = '[save figure folder]';
fig_format      = '-depsc'; % '-depsc' for vector graphics or '-dpng' for png image file

%% Change these things depending on experiment: 

% folder for each condition
P60_folders = {'Control_100Hz_P60/'; 'Test_100Hz_P60/'};
P21_folders = {'P21 control/'; 'P21 test/'; 'Chronos/'};

% what trial number does the induction occur at? (30 /120 for P60 and 60/180 for P21)
ind = 60;
end_trial = 180;

n = 10; % one plotted point is n number of data points

% which layers do i want to analyse?
layers = [1:6 9:12 15:26];
% which experiments do i want to analyse?
% 1 = control
% 2 = test
% 3 = chronos

fig_title = 'P60 test vs control';

for i = 1:2
    
    
    dataFolder = [targetdir P21_folders{i}];
    
    myFolder          = dir(dataFolder);
    
    qremove           = ismember({myFolder.name},{'.','..','.DS_Store'}); % locate where the unwanted names are
    myFolder(qremove) = []; % set these names to 0
    animal_number     = length(myFolder);
    
    %% Set up variables
    
    layers_idx = 1:length(layers);

    % how many spontaneous windows either side do i want to normalise to?
    spon_n = 5;
    
    timecount       = linspace(0,3,3000);
    
    response_start  = 1.005;
    response_end    = 1.03;
    response_window = response_end-response_start;
    
    % whisker-evoked response window
    stim_idx        = timecount >= response_start & timecount < response_end; 
    burst_idx       = timecount >= 0.975 & timecount < 1;
    burst_threshold = 2;
    
    % spontaneous window 1
    spon1_idx       = timecount >= 0.04 & timecount < 0.08;
    spon1_burst_idx = timecount >= 0.0375 & timecount < 0.04;
    
    % spontaneous window 2
    spon2_idx       = timecount >= 0.54 & timecount < 0.58;
    spon2_burst_idx = timecount >= 0.5375 & timecount < 0.54;

    
    idx = [burst_idx; spon1_burst_idx; spon2_burst_idx];
    
    % pad vectors with NaNs
    stim_response_all = NaN(30,14);    
    trial_counter = zeros(30,14);    
    spon1_response_all = NaN(30,14);
    spon2_response_all = NaN(30,14);    
    P2P = NaN(30,14);
    PP = NaN(30,14);
    NP = NaN(30,14);
    
    
    %% Load data
    for k = 1:length(myFolder)
        
        folderName      = myFolder(k).name;
        
        fileName        = dir([dataFolder folderName '/*.mat']);
        
        plot_data2 = [];
        plot_P2P_2 = [];
        plot_data = [];
        
        for a = 1:length(fileName)
            
            
            data        = load([dataFolder folderName filesep fileName(a).name]);
            
            %% Whisk response - all layers


            spikecount      = histc(data.ephys_data.conditions.spikes(layers,:,:),timecount,3);
            LFPtimeframe = data.ephys_data.conditions.LFP_trace(:,:,1000:1300);

            single_trial_length = size(spikecount,2);
            trial_length(a)    = single_trial_length;
            trial_counter(1:trial_length,a) = 1;
           

         
            %% Remove trials which show burst activity
         
        for m = 1:3
            burst_channels = false(32,single_trial_length);
            spikecount      = histc(data.ephys_data.conditions.spikes(layers,:,:),timecount,3);
            
            for h = layers_idx
                
                % How many spikes happen 25ms before stimulus?
                burst_spikes     = sum(spikecount(h,:,idx(m,:)),3);
                
                
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
        
        

            % Total spike count for each trial (median over channels)
            stim_response_all(1:single_trial_length,a) = nanmean(shortspike);
            spon1_response_all(1:single_trial_length,a) = nanmean(spon1_shortspike);
            spon2_response_all(1:single_trial_length,a) = nanmean(spon2_shortspike);
            
            
       
            
            %% Get the LFP data
    
            % LFP in 30 x 41 format - sum over channels
            LFP = squeeze(nanmean(LFPtimeframe));

            LFP = notch_filt(LFP',1000,50);
            LFP = LFP';
            LFP = notch_filt(LFP',1000,60);
            LFP = LFP';

            % peak to peak for each trial
            P2P(1:length(PP_this_exp),a) = PP_this_exp - NP_this_exp;

            
        end
        
    
    %% Get spikes into one vector for each animal
    
    
    % reshape to get one column for each animal
    plot_data = reshape(stim_response_all,[],1);
    spon1_data = reshape(spon1_response_all,[],1);
    spon2_data = reshape(spon2_response_all,[],1);
    trial_counter = reshape(trial_counter,[],1);
    
    plot_P2P = reshape(P2P,[],1);
    plot_PP = reshape(PP,[],1);
    plot_NP = reshape(NP,[],1);
    
    % get LFP values pre and post
    P2P_pre1(:,k) = nanmedian(plot_P2P(1:ind));
    P2P_post1(:,k) = nanmedian(plot_P2P(ind+1:end_trial));
    
    
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

    pre_data(:,k) = nanmedian(plot_data(1:ind));
    post_data(:,k) = nanmedian(plot_data(ind+1:end_trial));
    

    
    
    %% Normalise to the "pre" trials

        % calculate median of first 60 trials and divide by this
        norm_mean = nanmedian(plot_data(1:ind));
        plot_data = (plot_data./norm_mean)*100;
        
        post_data_mean(:,k) = nanmedian(plot_data(ind+1:end_trial));
        pre_data_mean(:,k) = nanmedian(plot_data(1:ind));
        
        delta(:,k) = nanmedian(plot_data(ind+1:end_trial)) - nanmedian(plot_data(1:ind));
        
        norm_P2P = nanmedian(plot_P2P(1:ind));
        plot_P2P = (plot_P2P./norm_P2P)*100;
        
        % median across every n values
        plot_data2(:,k) = arrayfun(@(ii) nanmedian(plot_data(ii:ii+n-1)),1:n:length(plot_data)-n+1)'; % the averaged vector
        plot_P2P_2(:,k) = arrayfun(@(ii) nanmedian(plot_P2P(ii:ii+n-1)),1:n:length(plot_P2P)-n+1)';

    end

    
    pre_timevec = 2.5:0.25*n:ind*0.25;
    post_timevec = ind*0.25+5:0.25*n:180*0.25;
    timevec = [pre_timevec post_timevec];
    
    
    % mean across animals and get standard error for error bars
    plot_data3 = nanmean(plot_data2,2);
    plot_stderr2 = nanstd(plot_data2,[],2)/sqrt(animal_number);
    
    plot_P2P_3 = nanmean(plot_P2P_2,2);
    plot_P2P_stderr = nanstd(plot_P2P_2,[],2)/sqrt(animal_number);    
    
    %% PLOT SPIKES AND P2P OVER TIME
    
    % control data
    if i == 1
        figure(1)
        errorbar(timevec,plot_data3(1:length(timevec)),plot_stderr2(1:length(timevec)),'s','MarkerEdgeColor',[0.600000023841858 0.600000023841858 0.600000023841858],'MarkerFaceColor',	[0.600000023841858 0.600000023841858 0.600000023841858],'Color',[0.600000023841858 0.600000023841858 0.600000023841858],'LineWidth',1.5,'MarkerSize',10)
        fixplot
        xlabel('Time (min)')
        ylabel('Normalised spike count')
        hold on
        
        figure(2)
        errorbar(timevec,plot_P2P_3(1:length(timevec)),plot_P2P_stderr(1:length(timevec)),'s','MarkerEdgeColor',[0.600000023841858 0.600000023841858 0.600000023841858],'MarkerFaceColor',	[0.600000023841858 0.600000023841858 0.600000023841858],'Color',[0.600000023841858 0.600000023841858 0.600000023841858],'LineWidth',1.5,'MarkerSize',10)
        fixplot
        xlabel('Time (min)')
        ylabel('Normalised LFP peak to peak')
        hold on
        
    end
    
    % test data
    if i == 2
        figure(1)
        errorbar(timevec,plot_data3(1:length(timevec)),plot_stderr2(1:length(timevec)),'o','MarkerEdgeColor',[0.0784313753247261 0.168627455830574 0.549019634723663],'MarkerFaceColor',[0.0784313753247261 0.168627455830574 0.549019634723663],'Color',[0.0784313753247261 0.168627455830574 0.549019634723663],'LineWidth',1.5,'MarkerSize',10)
        fixplot
        xlabel('Time (min)')
        ylabel('Normalised spike count')
        legend('- RWS', '+ RWS')
        hold on
        %ylim([0 600])

        
        fig_size        = [0 0 5 3];
        set(gcf,'PaperUnits','inches','PaperPosition',fig_size)
 
        % Optional saving
        if save_fig
            save_file   = fullfile(save_dir, ['Spikes ' fig_title]);
            print(gcf,save_file,fig_format,'-painters')
        end
        
        figure(2)
        errorbar(timevec,plot_P2P_3(1:length(timevec)),plot_P2P_stderr(1:length(timevec)),'o','MarkerEdgeColor',[0.0784313753247261 0.168627455830574 0.549019634723663],'MarkerFaceColor',[0.0784313753247261 0.168627455830574 0.549019634723663],'Color',[0.0784313753247261 0.168627455830574 0.549019634723663],'LineWidth',1.5,'MarkerSize',10)
        fixplot
        xlabel('Time (min)')
        ylabel('Normalised LFP peak to peak')
        legend('- RWS', '+ RWS')
        hold on
        %ylim([50 160])

        
        fig_size        = [0 0 5 3];
        set(gcf,'PaperUnits','inches','PaperPosition',fig_size)
 
        % Optional saving
        if save_fig
            save_file   = fullfile(save_dir, ['P2P ' fig_title]);
            print(gcf,save_file,fig_format,'-painters')
        end
        
       
    end

    
end   
    
