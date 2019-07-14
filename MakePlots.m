%% Collect Trials from MATs

subjectMatFolder = './SubjectDataMats';

subjectMats = {
                'Subject20Data.mat', ...
                'Subject21Data.mat', ...
                'Subject22Data.mat', ...
                'Subject23Data.mat', ...
                'Subject24Data.mat', ...
                'Subject25Data.mat', ...
                'Subject26Data.mat', ...
                'Subject28Data.mat', ...
                'Subject29Data.mat', ...
                'Subject30Data.mat'
              };

clear trials
trials = [];
for i = 1:length(subjectMats)
    file = fullfile(subjectMatFolder, subjectMats{i});
    temp = load(file, 'trials');
    trials = [trials; temp.trials];
end
          

% Rise Time

for i = 1:length(trials)
    
    x = trials(i).Data.EndEffPos_FromJA;
    
    zero_level = GetZeroLevel(x);
    ss_level = GetSteadyStateLevel(x);
    
    ten_percent_level = zero_level + 0.1*(ss_level - zero_level);
    ninty_percent_level = zero_level + 0.9*(ss_level - zero_level);
    
    if (ss_level < 0)
        inds = x < ten_percent_level;
        ten_percent_time = find(inds, 1, 'first');
        
        inds = x < ninty_percent_level;
        ninty_percent_time = find(inds, 1, 'first');
    else
        inds = x > ten_percent_level;
        ten_percent_time = find(inds, 1, 'first');
        
        inds = x > ninty_percent_level;
        ninty_percent_time = find(inds, 1, 'first');
        
    end
    trials(i).RiseTime = ninty_percent_time - ten_percent_time;
    
end


%% Make Plots Individually

figCount = 0;
for targetDirNum = 1:4
    % Make plots
    figs = [];
    for dampNum = 1:3
        figs = [figs, MakeX_VA_Damp_Plot(trials, targetDirNum, dampNum)];
    end

    % Fix Axes Limits
    
    % Get ylims for x and va subplots on all figures
    x_ylims = [];
    va_ylims = [];
    damp_ylims = [];
    for i = 1:length(figs)
        fig = figs(i);
        axes_all = fig.Children; % There should be 3 axes per fig b/c each fig has 3 subplots

        x_ylims = [x_ylims; ylim(axes_all(3))];
        va_ylims = [va_ylims; ylim(axes_all(2))];
        damp_ylims = [damp_ylims; ylim(axes_all(1))];
    end

    % Get the extreme y lims for x and va
    x_ylim = [min(x_ylims(:,1)), max(x_ylims(:,2))];
    va_ylim = [min(va_ylims(:,1)), max(va_ylims(:,2))];
    damp_ylim = [min(damp_ylims(:,1)), max(damp_ylims(:,2))];
    
    % Get range of limits 
    x_range = x_ylim(2) - x_ylim(1);
    va_range = va_ylim(2) - va_ylim(1);
    damp_range = damp_ylim(2) - damp_ylim(1);

    % Make limits 10% wider
    x_ylim(1) = x_ylim(1) - 0.1*x_range;
    x_ylim(2) = x_ylim(2) + 0.1*x_range;
    
    va_ylim(1) = va_ylim(1) - 0.1*va_range;
    va_ylim(2) = va_ylim(2) + 0.1*va_range;
    
    damp_ylim(1) = damp_ylim(1) - 0.1*damp_range;
    damp_ylim(2) = damp_ylim(2) + 0.1*damp_range;
    
    % Set ylims
    for i = 1:length(figs)
        fig = figs(i);
        axes_all = fig.Children;

        ylim(axes_all(1), damp_ylim);
        ylim(axes_all(2), va_ylim);
        ylim(axes_all(3), x_ylim);
    end
        
end

%% Make X, VA, Damp Plots Combined by Target Location (Left, Right, Down, Up)

locations = {'Left', 'Right', 'Down', 'Up'};
for targetDirNum = 1:4
    fig = Make_X_VA_Damp_Plot_Combined(trials, targetDirNum, 'ShowMovementMatch', true, 'ShowViolations', false);
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    figName = sprintf('CombinedPlot_%s.pdf', locations{targetDirNum});
    
%     print(fig, '-dpdf', figName, '-fillpage')
    
end

%% Make max velocity bar graphs

subject = 20;

% Collect all filtered data
data_processed = [];
for targetNum = 1:4
    for dampNum = 1:3
        temp = FilterAndAnalyzeData(trials, targetNum, dampNum);
        data_processed = [data_processed; temp];
        
        % print data table
        maxVelAll = [];
        maxVelArr = temp.XdotMax;
        fprintf('\n\n')
        fprintf('Direction: %s        Damping: %s\n', temp.TargetDirText, temp.DampingText)    
        fprintf('--------------------------------------------------------------\n');
        fprintf('Subject\tMean[m/s]\tStd[m/s]\tMax[m/s]\tn\n')
        fprintf('--------------------------------------------------------------\n');
        for i =1:length(maxVelArr)
            maxVel = maxVelArr(i);
            fprintf('%d\t%f\t%f\t%f\t%d\n', ...
                    maxVel.SubjectNumber, ...
                    maxVel.max_xdot_mean, ...
                    maxVel.max_xdot_std, ...
                    maxVel.max_xdot_max, ...
                    maxVel.nTrials)
            maxVelAll = [maxVelAll, maxVel.max_xdot];
        end
        fprintf('--------------------------------------------------------------\n');
        fprintf('Total\t%f\t%f\t%f\t%d\n', ...
                mean(maxVelAll), ...
                std(maxVelAll), ...
                max(maxVelAll), ...
                length(maxVelAll));
        
        
        
    end
end

fig = figure;
set(fig,'Color',[1,1,1]);
for targetNum = 1:4
    % Get subset of data of based on target number
    set_inds = [data_processed.TargetDirNum] == targetNum;
    data_processed_set = data_processed(set_inds);
    
    % Get means, std, damp label text
    xdotmax_mean = zeros(1,3);
    xdotmax_std = zeros(1,3);
    nTrials = zeros(1,3);
    for i = 1:3
        xdotmax_data = data_processed_set(i).XdotMax;
        subjectInds = [xdotmax_data.SubjectNumber];
        subjectInd = subjectInds == subject;
        
        xdotmax_subject = xdotmax_data(subjectInd);
        xdotmax_mean(i) = xdotmax_subject.max_xdot_mean;
        xdotmax_std(i) = xdotmax_subject.max_xdot_std;
        nTrials(i) = xdotmax_subject.nTrials;
    end
    
    dampTexts = {};
    dampTextsCat = {};
    for i = 1:length(data_processed_set)
        dText = sprintf('%s (n=%d)', data_processed_set(i).DampingText, ...
                                     nTrials(i));
        dampTexts{i} = dText;
        dampTextsCat{i} = data_processed_set(i).DampingText;
    end
    c = categorical(dampTextsCat);
    
    % Plot 
    ax = subplot(1, 4, targetNum);
    for i = 1:length(xdotmax_mean)
        bar(c(i),xdotmax_mean(i));
        hold on
    end
    
    for i = 1:length(xdotmax_mean)
        er = errorbar(c(i),xdotmax_mean(i), xdotmax_std(i));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
    end
    
    % Make plot pretty
    if (targetNum ~= 2)
        title(data_processed_set(1).TargetDirText);
    else
        title(sprintf('Subject %d\n\n%s', subject, data_processed_set(1).TargetDirText))
    end
    if (targetNum == 1)
        ylabel('$\dot{x}$ [m/s]', 'Interpreter', 'latex');
    end
    box('off')
    set(ax, 'ticklength', [0,0]);
    legend(dampTexts);
    legend('boxoff')
    hold off
end

% Set Y Limits equal
child = fig.Children;
axes_arr = [];
for i = length(child):-1:1
    if (~isa(child(i), 'matlab.graphics.axis.Axes'))
        child(i) = [];
    end
end
SetYLimsEqual(child, 'BottomPadPercent', 0)


%% Make percent overshoot bar graphs
subject = 20;

po_rm_meas = zeros(10,4,4);


% Collect all filtered data
data_processed = [];
for targetNum = 1:4
    for dampNum = 1:3
        temp = FilterAndAnalyzeData(trials, targetNum, dampNum);
        data_processed = [data_processed; temp];
        
        po_rm_means_temp = zeros(10,2);
        
        % print data table
        poAll = [];
        poArr = temp.PercentOvershoot;
        fprintf('\n\n')
        fprintf('Direction: %s        Damping: %s\n', temp.TargetDirText, temp.DampingText)    
        fprintf('--------------------------------------------------------------\n');
        fprintf('Percent Overshoot\n')
        fprintf('Subject\tMean\tStd\tMax\tn\n')
        fprintf('--------------------------------------------------------------\n');
        for i =1:length(poArr)
            po = poArr(i);
            fprintf('%d\t%f\t%f\t%f\t%d\n', ...
                    po.SubjectNumber, ...
                    po.po_mean, ...
                    po.po_std, ...
                    po.po_max, ...
                    po.nTrials)
            poAll = [poAll, po.po];
            
            % put data into temp repeated measures anova mean array
            po_rm_means_temp(i,1) = po.SubjectNumber;
            po_rm_means_temp(i,2) = po.po_mean;
        end
        fprintf('--------------------------------------------------------------\n');
        fprintf('Total\t%f\t%f\t%f\t%d\n', ...
                mean(poAll), ...
                std(poAll), ...
                max(poAll), ...
                length(poAll));
            
        % put temp data into repeated measures anova mean array 
        po_rm_means_temp = sortrows(po_rm_means_temp);
        po_rm_meas(:,1,targetNum) = po_rm_means_temp(:,1);
        po_rm_meas(:,dampNum+1,targetNum) = po_rm_means_temp(:,2);
            
    end
end

fig = figure;
set(fig,'Color',[1,1,1]);
for targetNum = 1:4
    % Get subset of data of based on target number
    set_inds = [data_processed.TargetDirNum] == targetNum;
    data_processed_set = data_processed(set_inds);
    
    % Get means, std, damp label text
    po_mean = zeros(1,3);
    po_std = zeros(1,3);
    nTrials = zeros(1,3);
    for i = 1:3
        po_data = data_processed_set(i).PercentOvershoot;
        subjectInds = [po_data.SubjectNumber];
        subjectInd = subjectInds == subject;
        
        po_subject = po_data(subjectInd);
        po_mean(i) = po_subject.po_mean;
        po_std(i) = po_subject.po_std;
        nTrials(i) = po_subject.nTrials;
    end
    
    dampTexts = {};
    dampTextsCat = {};
    for i = 1:length(data_processed_set)
        dText = sprintf('%s (n=%d)', data_processed_set(i).DampingText, ...
                                     nTrials(i));
        dampTexts{i} = dText;
        dampTextsCat{i} = data_processed_set(i).DampingText;
    end
    c = categorical(dampTextsCat);
    
    % Plot 
    ax = subplot(1, 4, targetNum);
    for i = 1:length(po_mean)
        bar(c(i),po_mean(i));
        hold on
    end
    
    for i = 1:length(po_mean)
        er = errorbar(c(i),po_mean(i), po_std(i));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
    end
    
    % Make plot pretty
    if (targetNum ~= 2)
        title(data_processed_set(1).TargetDirText);
    else
        title(sprintf('Subject %d\n\n%s', subject, data_processed_set(1).TargetDirText))
    end
    if (targetNum == 1)
        ylabel('Percent Overshoot', 'Interpreter', 'latex');
    end
    box('off')
    set(ax, 'ticklength', [0,0]);
    legend(dampTexts);
    legend('boxoff')
    hold off
end

% Set Y Limits equal
child = fig.Children;
axes_arr = [];
for i = length(child):-1:1
    if (~isa(child(i), 'matlab.graphics.axis.Axes'))
        child(i) = [];
    end
end
SetYLimsEqual(child, 'BottomPadPercent', 0)


% Repeated Measures Anova - Percent Overshoot
fprintf('RM ANOVA - Percent Overshoot - All Damping Groups')
% Left
fprintf('\n\n\nLeft\n')
CalcRepeatedAnova(po_rm_meas(:,2:4,1));


% Left
fprintf('\n\n\nRight\n')
CalcRepeatedAnova(po_rm_meas(:,2:4,2));


% Left
fprintf('\n\n\nDown\n')
CalcRepeatedAnova(po_rm_meas(:,2:4,3));

% Left
fprintf('\n\n\nUp\n')
CalcRepeatedAnova(po_rm_meas(:,2:4,4));


%% Repeated Measures Anova - Percent Overshoot - Positive vs Negative
% Left
po_rm_results = [];
targetDirTexts = {'Left', 'Right', 'Down', 'Up'};
for targetNum = 1:4
    fprintf('\n\n\n%s\n', targetDirTexts{targetNum});
    fprintf('Positive vs Negative\n');
    temp = CalcRepeatedAnova(po_rm_meas(:,[2,3],targetNum));
    temp.Direction = targetDirTexts{targetNum};
    temp.DampGroup1 = 'Positive';
    temp.DampGroup2 = 'Negative';
    po_rm_results = [po_rm_results, temp];
    
    fprintf('\nPositive vs Variable\n');
    CalcRepeatedAnova(po_rm_meas(:,[2,4],targetNum));
    temp = CalcRepeatedAnova(po_rm_meas(:,[2,4],targetNum));
    temp.Direction = targetDirTexts{targetNum};
    temp.DampGroup1 = 'Positive';
    temp.DampGroup2 = 'Variable';
    po_rm_results = [po_rm_results, temp];

    fprintf('\nNegative vs Variable\n');
    temp = CalcRepeatedAnova(po_rm_meas(:,[3,4],targetNum));
    temp.Direction = targetDirTexts{targetNum};
    temp.DampGroup1 = 'Negative';
    temp.DampGroup2 = 'Variable';
    po_rm_results = [po_rm_results, temp];
end

%% Make Force Plots During Movement


% Collect all filtered data
data_processed = [];
for targetNum = 1:4
    for dampNum = 1:3
        data_processed = [data_processed; FilterAndAnalyzeData(trials, targetNum, dampNum)];
    end
end

% Line Color RGB vals
lightGreyColor = [179, 179, 179]/256;

fprintf('\n\n\n')
fprintf('--------------------------------------------------------------\n');
fprintf('Direction\tDamping\tSubject\tMeanMean[N]\tMeanStd[N]\tMaxMean[N]\tMaxStd[N]\tn\n')
fprintf('--------------------------------------------------------------\n');

x_all = {[],[],[];
         [],[],[];
         [],[],[];
         [],[],[]};
     
force_all = {[],[],[];
             [],[],[];
             [],[],[];
             [],[],[]};

force_rm_max = zeros(10,4,4);
force_rm_mean = zeros(10,4,4);
         
% Print force data
for targetNum = 1:4
    % Get subset of data of based on target number
    set_inds = [data_processed.TargetDirNum] == targetNum;
    data_processed_set = data_processed(set_inds);

    % Get all trials without violations with this target number
    trials_wov = [];
%     va_mean = {};
    va_mean_movement_length = [];
    for i = 1:length(data_processed_set)
        trials_wov = [trials_wov; data_processed_set(i).wov.trials];
        va_mean = data_processed_set(i).wov.va.mean_;
        va_mean_movement_length(end+1) = GetLastMoveInd(va_mean) - GetFirstMoveInd(va_mean);
    end
    
    max_movement_length = max(va_mean_movement_length);
    

    
    subjectNumbersAll = [trials_wov.SubjectNumber];
    subjectNumbers = unique(subjectNumbersAll);
    for i = 1:length(subjectNumbers)
        subjectNumber = subjectNumbers(i);
        subjectInds = subjectNumbersAll == subjectNumber;
        trials_subject = trials_wov(subjectInds);
        
        
 
        
        for dampNum = 1:3
            dampNumAll = [trials_subject.DampingNumber];
            dampNumInds = dampNumAll == dampNum;
            trials_subject_damp = trials_subject(dampNumInds);
            dampText = trials_subject_damp(i).DampingText;
            tarDirText = trials_subject_damp(i).TargetDirText;
            
            force = [];
            x = [];
            for j = 1:length(trials_subject_damp)
                firstMoveInd = trials_subject_damp(j).FirstMoveInd;
                lastMoveInd = firstMoveInd + max_movement_length;
                force = [force, trials_subject_damp(j).Data.Force(firstMoveInd:lastMoveInd, 1)];
                x = [x, trials_subject_damp(j).Data.EndEffPos_FromJA(firstMoveInd:lastMoveInd)];
            end
            
            % add raw data to arrays
            force_all{targetNum, dampNum} = [force_all{targetNum, dampNum}, force];
            x_all{targetNum, dampNum} = [x_all{targetNum, dampNum}, x];
            
            % calc mean, std, max of rms of the mean
            force_rms = sqrt(force.^2);
            force_rms_mean = mean(force_rms, 1);
            force_rms_max = max(force_rms, [], 1);
            
            force_rms_mean_mean = mean(force_rms_mean);
            force_rms_mean_std = std(force_rms_mean);
            force_rms_max_mean = mean(force_rms_max);
            force_rms_max_std = std(force_rms_max);
            
            % put data into repeated measures anova data array
            force_rm_mean(i,1,targetNum) = subjectNumber;
            force_rm_mean(i,dampNum+1,targetNum) = force_rms_mean_mean;
            force_rm_max(i,1,targetNum) = subjectNumber;
            force_rm_max(i,dampNum+1,targetNum) = force_rms_max_mean;
            
            % print data table
            fprintf('%s\t%s\t%d\t%f\t%f\t%f\t%f\t%d\n', ...
                    tarDirText, ...
                    dampText, ...
                    subjectNumber, ...
                    force_rms_mean_mean, ...
                    force_rms_mean_std, ...
                    force_rms_max_mean, ...
                    force_rms_max_std, ...
                    length(trials_subject_damp))


        end
    end
end

dampText = {'Positive', 'Negative', 'Variable'};
targetDirText = {'Left', 'Right', 'Down', 'Up'};

for targetNum = 1:4
    fig = figure;
    set(fig,'Color',[1,1,1]);
    for dampNum = 1:3
        x_temp = GetDataSetStats(x_all{targetNum, dampNum});
        force_temp = GetDataSetStats(force_all{targetNum, dampNum});
        
        % Plot x
        subplot(2,3,dampNum)
        plot(x_temp.mean_, 'Color', 'k');
        hold on
        plot(x_temp.ub_, 'Color', 'k', 'LineStyle', '--');
        plot(x_temp.lb_, 'Color', 'k', 'LineStyle', '--');
        hold off
        if (dampNum == 2)
            title(sprintf('%s\n\n%s', targetDirText{targetNum}, dampText{dampNum}));
        else
            title(sprintf('%s', dampText{dampNum}));
        end
        if (dampNum == 1)
            ylabel('x [m]')
        end
        
        
        % Plot force
        subplot(2,3,dampNum+3)
        force_rms = sqrt(force_temp.mean_.^2);
        plot(force_temp.mean_, 'k');
        hold on
        plot(force_temp.ub_, 'k--');
        plot(force_temp.lb_, 'k--', 'HandleVisibility', 'off');
        plot(mean(force_rms)*ones(size(force_rms)), 'r')
        plot(max(force_rms)*ones(size(force_rms)), 'b')
        hold on
        legend('Nominal', '+/-3std', 'RMS Average', '|Max|');
        legend('boxoff');
        if (dampNum == 1)
            ylabel('Force [N]')
        end
        xlabel('Time [ms]')
    end
    
    axes_arr = [subplot(2,3,1), subplot(2,3,2), subplot(2,3,3)];
    SetYLimsEqual(axes_arr);
    
    axes_arr = [subplot(2,3,4), subplot(2,3,5), subplot(2,3,6)];
    SetYLimsEqual(axes_arr);
    
    
end


fprintf('\n\n\n')
fprintf('--------------------------------------------------------------\n');
fprintf('Direction\tDamping\tMean[N]\tMax[N]\tn\n')
fprintf('--------------------------------------------------------------\n');

for targetNum = 1:4
    for dampNum = 1:3
        force_temp = GetDataSetStats(force_all{targetNum, dampNum});
        force_rms = sqrt(force_temp.mean_.^2);
        fprintf('%s\t%s\t%f\t%f\t%d\n', ...
                targetDirText{targetNum}, ...
                dampText{dampNum}, ...
                mean(force_rms), ...
                max(force_rms), ...
                size(force_temp.data_, 2));
                
        
    end
end

%%
% Repeated Measures Anova - Force Mean
% Left
fprintf('\n\n\nLeft - Force Mean\n')
CalcRepeatedAnova(force_rm_mean(:,2:4,1));

fprintf('\n\n\nLeft - Force Mean\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_mean(i,1,1), ...
                                force_rm_mean(i,2,1), ...
                                force_rm_mean(i,3,1), ...
                                force_rm_mean(i,4,1))
end

% Right
fprintf('\n\n\nRight - Force Mean\n')
CalcRepeatedAnova(force_rm_mean(:,2:4,2));

fprintf('\n\n\nRight - Force Mean\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_mean(i,1,2), ...
                                force_rm_mean(i,2,2), ...
                                force_rm_mean(i,3,2), ...
                                force_rm_mean(i,4,2))
end

% Down
fprintf('\n\n\nDown - Force Mean\n')
CalcRepeatedAnova(force_rm_mean(:,2:4,3));

fprintf('\n\n\nDown - Force Mean\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_mean(i,1,3), ...
                                force_rm_mean(i,2,3), ...
                                force_rm_mean(i,3,3), ...
                                force_rm_mean(i,4,3))
end

% Up
fprintf('\n\n\nUp - Force Mean\n')
CalcRepeatedAnova(force_rm_mean(:,2:4,4));

fprintf('\n\n\nUp - Force Mean\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_mean(i,1,4), ...
                                force_rm_mean(i,2,4), ...
                                force_rm_mean(i,3,4), ...
                                force_rm_mean(i,4,4))
end

% Repeated Measures Anova - Force Max
% Left
fprintf('\n\n\nLeft - Force Max\n')
CalcRepeatedAnova(force_rm_max(:,2:4,1));

fprintf('\n\n\nLeft - Force Max\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_max(i,1,1), ...
                                force_rm_max(i,2,1), ...
                                force_rm_max(i,3,1), ...
                                force_rm_max(i,4,1))
end

% Right
fprintf('\n\n\nRight - Force Max\n')
CalcRepeatedAnova(force_rm_max(:,2:4,2));

fprintf('\n\n\nRight - Force Max\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_max(i,1,2), ...
                                force_rm_max(i,2,2), ...
                                force_rm_max(i,3,2), ...
                                force_rm_max(i,4,2))
end

% Down
fprintf('\n\n\nDown - Force Max\n')
CalcRepeatedAnova(force_rm_max(:,2:4,3));

fprintf('\n\n\nDown - Force Max\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_max(i,1,3), ...
                                force_rm_max(i,2,3), ...
                                force_rm_max(i,3,3), ...
                                force_rm_max(i,4,3))
end

% Up
fprintf('\n\n\nUp - Force Max\n')
CalcRepeatedAnova(force_rm_max(:,2:4,4));

fprintf('\n\n\nUp - Force Max\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_max(i,1,4), ...
                                force_rm_max(i,2,4), ...
                                force_rm_max(i,3,4), ...
                                force_rm_max(i,4,4))
end

%% Rise Time


% Collect all trials without violations
trials_wov = [];
for targetNum = 1:4
    for dampNum = 1:3
        temp = FilterAndAnalyzeData(trials, targetNum, dampNum);
        trials_wov = [trials_wov; temp.wov.trials];
    end
end


rise_time_all = {[],[],[];
                [],[],[];
                [],[],[];
                [],[],[]};
dampText = {'Positive', 'Negative', 'Variable'};
targetDirText = {'Left', 'Right', 'Down', 'Up'};

rt_rm_meas = zeros(10, 4, 4);

fprintf('Rise Time\n')
fprintf('Direction\tDamping\tSubject\tMean[ms]\tStd[ms]\tn\n')

subjectNumbersAll = [trials_wov.SubjectNumber];
subjectNumbers = unique(subjectNumbersAll);
    
for iSubject = 1:length(subjectNumbers)
    subjectNumber = subjectNumbers(iSubject);

    subjectInds = subjectNumbersAll == subjectNumber;
    trials_subject = trials_wov(subjectInds);

    targetNumbers = [trials_subject.TargetDirNum];
    dampNumbers = [trials_subject.DampingNumber];

    rt_rm_meas(iSubject, 1, :) = subjectNumber;
    
    
    for targetNum = 1:4
        for dampNum = 1:3
            targetInds = targetNumbers == targetNum;
            dampInds = dampNumbers == dampNum;

            setInds = targetInds & dampInds;
            trials_set = trials_subject(setInds);

            rise_time_subject = [trials_set.RiseTime];
            
            if (dampNum == 3)
                mer = 1;
            end
            
            rt_rm_meas(iSubject, dampNum+1, targetNum) = mean(rise_time_subject);


            fprintf('%s\t%s\t%d\t%f\t%f\t%d\n', ...
                    targetDirText{targetNum}, ...
                    dampText{dampNum}, ...
                    subjectNumber, ...
                    mean(rise_time_subject), ...
                    std(rise_time_subject), ...
                    length(rise_time_subject));

            rise_time_all{targetNum, dampNum} = [rise_time_all{targetNum, dampNum}, rise_time_subject];

        end
    end

end

fprintf('\n\n')
fprintf('Rise Time\n')
fprintf('Direction\tDamping\tMean[ms]\tStd[ms]\tn\n')
for targetNum = 1:4
    for dampNum = 1:3

        rise_time = rise_time_all{targetNum, dampNum};

        fprintf('%s\t%s\t%f\t%f\t%d\n', ...
                targetDirText{targetNum}, ...
                dampText{dampNum}, ...
                mean(rise_time), ...
                std(rise_time), ...
                length(rise_time));

        rise_time_all{targetNum, dampNum} = [rise_time_all{targetNum, dampNum}, ];

    end
end

% Make plot
fig = figure;
set(fig,'Color',[1,1,1]);

c = categorical(dampText);
for targetNum = 1:4
    rise_time_means = zeros(1,3);
    rise_time_stds = zeros(1,3);
    for dampNum = 1:3
        rise_time_means(dampNum) = mean(rise_time_all{targetNum, dampNum});
        rise_time_stds(dampNum) = std(rise_time_all{targetNum, dampNum});
    end
    
    % Plot 
    ax = subplot(1, 4, targetNum);
    for i = 1:3
        bar(c(i),rise_time_means(i));
        hold on
    end
    
    for i = 1:3
        er = errorbar(c(i),rise_time_means(i), rise_time_stds(i));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
    end
    
    % Make plot pretty
    title(targetDirText{targetNum});
    if (targetNum == 1)
        ylabel('Rise Time[ms]', 'Interpreter', 'latex');
    end
    box('off')
    set(ax, 'ticklength', [0,0]);
    legend(dampText);
    legend('boxoff')
    hold off
end


% Rise Time Repeated Measures Anova

% Left
fprintf('\n\n\nLeft\n')
CalcRepeatedAnova(rt_rm_meas(:,2:4,1));


% Left
fprintf('\n\n\nRight\n')
CalcRepeatedAnova(rt_rm_meas(:,2:4,2));


% Left
fprintf('\n\n\nDown\n')
CalcRepeatedAnova(rt_rm_meas(:,2:4,3));

% Left
fprintf('\n\n\nUp\n')
CalcRepeatedAnova(rt_rm_meas(:,2:4,4));

