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
    fig = Make_X_VA_Damp_Plot_Combined(trials, targetDirNum, 'ShowMovementMatch', true);
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    figName = sprintf('CombinedPlot_%s.pdf', locations{targetDirNum});
    
%     print(fig, '-dpdf', figName, '-fillpage')
    
end

%% Make max velocity bar graphs

% Collect all filtered data
data_processed = [];
for targetNum = 1:4
    for dampNum = 1:3
        data_processed = [data_processed; FilterAndAnalyzeData(trials, targetNum, dampNum)];
    end
end

fig = figure;
set(fig,'Color',[1,1,1]);
for targetNum = 1:4
    % Get subset of data of based on target number
    set_inds = [data_processed.TargetDirNum] == targetNum;
    data_processed_set = data_processed(set_inds);
    
    % Get means, std, damp label text
    xdotmax_means = [data_processed_set.XdotMax_data_wov_mean];
    xdotmax_std = [data_processed_set.XdotMax_data_wov_std];
    dampTexts = {};
    dampTextsCat = {};
    for i = 1:length(data_processed_set)
        dText = sprintf('%s (n=%d)', data_processed_set(i).DampingText, ...
                                     length(data_processed_set(i).PercentOvershoot_data_wov));
        dampTexts{i} = dText;
        dampTextsCat{i} = data_processed_set(i).DampingText;
    end
    c = categorical(dampTextsCat);
    
    % Plot 
    ax = subplot(1, 4, targetNum);
    for i = 1:length(xdotmax_means)
        bar(c(i),xdotmax_means(i));
        hold on
    end
    
    for i = 1:length(xdotmax_means)
        er = errorbar(c(i),xdotmax_means(i), xdotmax_std(i));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
    end
    
    % Make plot pretty
    title(data_processed_set(1).TargetDirText);
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

% Collect all filtered data
data_processed = [];
for targetNum = 1:4
    for dampNum = 1:3
        data_processed = [data_processed; FilterAndAnalyzeData(trials, targetNum, dampNum)];
    end
end

fig = figure;
set(fig,'Color',[1,1,1]);
for targetNum = 1:4
    % Get subset of data of based on target number
    set_inds = [data_processed.TargetDirNum] == targetNum;
    data_processed_set = data_processed(set_inds);
    
    % Get means, std, damp label text
    po_means = [data_processed_set.PercentOvershoot_data_wov_mean];
    po_std = [data_processed_set.PercentOvershoot_data_wov_std];
    dampTexts = {};
    dampTextsCat = {};
    for i = 1:length(data_processed_set)
        dText = sprintf('%s (n=%d)', data_processed_set(i).DampingText, ...
                                     length(data_processed_set(i).PercentOvershoot_data_wov));
        dampTexts{i} = dText;
        dampTextsCat{i} = data_processed_set(i).DampingText;
    end
    c = categorical(dampTextsCat);
    
    % Plot 
    ax = subplot(1, 4, targetNum);
    for i = 1:length(po_means)
        bar(c(i),po_means(i));
        hold on
    end
    
    for i = 1:length(po_means)
        er = errorbar(c(i),po_means(i), po_std(i));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
    end
    
    % Make plot pretty
    title(data_processed_set(1).TargetDirText);
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

for targetNum = 1:4
    % Get subset of data of based on target number
    set_inds = [data_processed.TargetDirNum] == targetNum;
    data_processed_set = data_processed(set_inds);

    lastMoveInds = [];
    for i = 1:length(data_processed_set)
        va_mm_mean = mean(data_processed_set(i).va_mm, 2);
        inds = abs(va_mm_mean) > 0.001;
        lastMoveInds = [lastMoveInds, find(inds, 1, 'last')];        
    end
    lastMoveInd = max(lastMoveInds);

    fig = figure;
    set(fig,'Color',[1,1,1]);
    for i = 1:length(data_processed_set)
        dampNum = data_processed_set(i).DampingNum;
        dampText = data_processed_set(i).DampingText;
        tarDirText = data_processed_set(i).TargetDirText;
        
        pltcols = 3;
        pltrows = 2;
        
        % x movement matched data
        x_mm = data_processed_set(i).x_mm_wov(1:lastMoveInd, :);
        x_mm_mean = mean(x_mm, 2);
        x_mm_std = std(x_mm, [], 2);
        x_mm_ub = x_mm_mean + 3*x_mm_std;
        x_mm_lb = x_mm_mean - 3*x_mm_std;
        
        force_mm = abs(data_processed_set(i).force_mm_wov(1:lastMoveInd,:));
        force_mm_mean = mean(force_mm, 2);
        
        subplot(pltrows, pltcols, dampNum)
        plot(x_mm, 'Color', lightGreyColor);
        hold on
        plot(x_mm_mean, 'Color', 'k')
        plot(x_mm_ub, 'Color', 'k', 'LineStyle', '--');
        plot(x_mm_lb, 'Color', 'k', 'LineStyle', '--');
        hold off
        if (dampNum == 1)
            ylabel('x [mm]')
            title(dampText);
        elseif (dampNum == 2)
            title(sprintf('Direction: %s\n\n%s', tarDirText, dampText));
        elseif (dampNum == 3)
            title(dampText)
        end
            
        
        subplot(pltrows, pltcols, pltcols*1+dampNum)
        plot(force_mm, 'Color', lightGreyColor)
        hold on
        plot(force_mm_mean, 'Color', 'k')
        plot(mean(force_mm_mean)*ones(size(force_mm_mean)), 'Color', 'r', 'LineWidth', 1.25);
        plot(max(force_mm_mean)*ones(size(force_mm_mean)), 'Color', 'b', 'LineWidth', 1.25);
        hold on
        if (dampNum == 1)
            ylabel('|Force| [N]')
            xlabel('Time [ms]')
        elseif (dampNum == 2)
            xlabel('Time [ms]')
        elseif (dampNum == 3)
            xlabel('Time [ms]')
        end

    end
    
    ax1 = subplot(pltrows, pltcols, 1);
    ax2 = subplot(pltrows, pltcols, 2);
    ax3 = subplot(pltrows, pltcols, 3);
    
    axes_arr = [ax1, ax2, ax3];
    SetYLimsEqual(axes_arr, 'TopPadPercent', 0.05, 'BottomPadPercent', 0.05)
    
    ax1 = subplot(pltrows, pltcols, 4);
    ax2 = subplot(pltrows, pltcols, 5);
    ax3 = subplot(pltrows, pltcols, 6);
    
    axes_arr = [ax1, ax2, ax3];
    SetYLimsEqual(axes_arr, 'TopPadPercent', 0, 'BottomPadPercent', 0)
    
end
