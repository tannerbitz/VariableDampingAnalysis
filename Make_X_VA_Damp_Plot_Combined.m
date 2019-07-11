function res = Make_X_VA_Damp_Plot_Combined(trials, targetDirNum, varargin)
    
    % Input Parser to add 'ShowViolations' option
    p = inputParser;
    argName = 'ShowViolations';
    defaultVal = false;
    addParameter(p, argName, defaultVal);
    addParameter(p, 'ShowMovementMatch', false);
    
    
    parse(p,varargin{:});
    showViolations = p.Results.ShowViolations;
    showMovementMatch = p.Results.ShowMovementMatch;

    % Plot
    fig1 = figure;
    set(fig1,'Color',[1,1,1]);
    
    if (showMovementMatch)
        fig2 = figure;
        set(fig2, 'Color', [1,1,1]);
    end
    
    dampColors = {'b', 'r', 'g'};
    dampTexts = {};
    
    % Make Plots
    for dampNum = 1:3

        data = FilterAndAnalyzeData(trials, targetDirNum, dampNum);
    
        % Line Color RGB vals
        lightGreyColor = [179, 179, 179]/256;
        lightGreenColor = [66, 245, 69]/256;
        
        % set fig1 as current figure
        set(0, 'CurrentFigure', fig1);

        % Position
        subplot(3,4,1+(dampNum-1));
        plot(data.x_data_wov, 'Color', lightGreyColor)
        hold on
        if (showViolations)
            plot(data.x_data_wv, 'Color', lightGreenColor)
        end
        plot(data.x_data_wov_mean, 'Color', 'k')
        plot(data.x_data_wov_ub, 'Color', 'k', 'LineStyle', '--');
        plot(data.x_data_wov_lb, 'Color', 'k', 'LineStyle', '--');
        hold off
        if (dampNum == 1)
            ylabel('x [m]', 'Interpreter', 'latex');
        end
        if (dampNum == 2)
            title(sprintf('Direction: %s\n\nDamping: %s', ...
                        data.TargetDirText, ...
                        data.DampingText));
        else
            title(sprintf('Damping: %s', ...
                          data.DampingText));
        end

        xlim([0, data.nSamples]);
        box('off');
                
        % VA
        subplot(3,4,5+(dampNum-1));
        plot(data.va_data_wov, 'Color', lightGreyColor)
        hold on 
        if (showViolations)
            plot(data.va_data_wv, 'Color', lightGreenColor)
        end
        plot(data.va_data_wov_mean, 'Color', 'k')
        plot(data.va_data_wov_ub, 'Color', 'k', 'LineStyle', '--');
        plot(data.va_data_wov_lb, 'Color', 'k', 'LineStyle', '--');
        hold off
        if (dampNum == 1)
            ylabel('$\dot{x} \ddot{x}$', 'Interpreter', 'latex');
        end
        xlim([0, data.nSamples]);
        box('off');

        % Damping
        subplot(3,4,9+(dampNum-1));
        plot(data.damp_data_wov, 'Color', lightGreyColor)
        hold on
        if (showViolations)
            plot(data.damp_data_wv, 'Color', lightGreenColor)
        end
        plot(data.damp_data_wov_mean, 'Color', 'k')
        plot(data.damp_data_wov_ub, 'Color', 'k', 'LineStyle', '--');
        plot(data.damp_data_wov_lb, 'Color', 'k', 'LineStyle', '--');
        hold off
        if (dampNum == 1)
            ylabel('Damping [Ns/m]', 'Interpreter', 'latex');
        end
        xlabel('Time [ms]', 'Interpreter', 'latex');
        xlim([0, data.nSamples]);
        box('off');
        
        % Means - x
        subplot(3,4,4)
        if (dampNum > 1)
            hold on
        end
        plot(data.x_data_wov_mean, dampColors{dampNum});
        hold off
        box('off');

        % Means - va
        subplot(3,4,8)
        if (dampNum > 1)
            hold on
        end
        plot(data.va_data_wov_mean, dampColors{dampNum});
        hold off
        box('off');

        % Means - damp
        subplot(3,4,12)
        if (dampNum > 1)
            hold on
        end
        plot(data.damp_data_wov_mean, dampColors{dampNum});
        hold off
        box('off');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%            MOVEMENT MATCHED PLOT        %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (showMovementMatch)
            % set fig2 as current figure
            set(0, 'CurrentFigure', fig2);

            % Position
            subplot(3,4,1+(dampNum-1));
            plot(data.x_mm_wov, 'Color', lightGreyColor)
            hold on
            if (showViolations)
                plot(data.x_mm_wv, 'Color', lightGreenColor)
            end
            plot(data.x_mm_wov_mean, 'Color', 'k')
            plot(data.x_mm_wov_ub, 'Color', 'k', 'LineStyle', '--');
            plot(data.x_mm_wov_lb, 'Color', 'k', 'LineStyle', '--');
            hold off
            if (dampNum == 1)
                ylabel('x [m]', 'Interpreter', 'latex');
            end
            if (dampNum == 2)
                title(sprintf('Direction: %s\n\nDamping: %s', ...
                            data.TargetDirText, ...
                            data.DampingText));
            elseif (dampNum == 3)
                title(sprintf('%s\n\n%s', ...
                              'DELAY PORTIONS REMOVED', ...
                              data.DampingText));
            else
                title(sprintf('Damping: %s', ...
                              data.DampingText));
            end

            xlim([0, size(data.x_mm, 1)]);
            box('off');

            % VA
            subplot(3,4,5+(dampNum-1));
            plot(data.va_mm_wov, 'Color', lightGreyColor)
            hold on 
            if (showViolations)
                plot(data.va_mm_wv, 'Color', lightGreenColor)
            end
            plot(data.va_mm_wov_mean, 'Color', 'k')
            plot(data.va_mm_wov_ub, 'Color', 'k', 'LineStyle', '--');
            plot(data.va_mm_wov_lb, 'Color', 'k', 'LineStyle', '--');
            hold off
            if (dampNum == 1)
                ylabel('$\dot{x} \ddot{x}$', 'Interpreter', 'latex');
            end
            xlim([0, size(data.x_mm, 1)]);
            box('off');

            % Damping
            subplot(3,4,9+(dampNum-1));
            plot(data.damp_mm_wov, 'Color', lightGreyColor)
            hold on
            if (showViolations)
                plot(data.damp_mm_wv, 'Color', lightGreenColor)
            end
            plot(data.damp_mm_wov_mean, 'Color', 'k')
            plot(data.damp_mm_wov_ub, 'Color', 'k', 'LineStyle', '--');
            plot(data.damp_mm_wov_lb, 'Color', 'k', 'LineStyle', '--');
            hold off
            if (dampNum == 1)
                ylabel('Damping [Ns/m]', 'Interpreter', 'latex');
            end
            xlabel('Time [ms]', 'Interpreter', 'latex');
            xlim([0, size(data.x_mm, 1)]);
            box('off');

            % Means - x
            subplot(3,4,4)
            if (dampNum > 1)
                hold on
            end
            plot(data.x_mm_wov_mean, dampColors{dampNum});
            hold off
            box('off');

            % Means - va
            subplot(3,4,8)
            if (dampNum > 1)
                hold on
            end
            plot(data.va_mm_wov_mean, dampColors{dampNum});
            hold off
            box('off');

            % Means - damp
            subplot(3,4,12)
            if (dampNum > 1)
                hold on
            end
            plot(data.damp_mm_wov_mean, dampColors{dampNum});
            hold off
            box('off');
        end
        
        % Add to dampText
        dampText{dampNum} = data.DampingText;
        
    end
    
    % Set fig1 as current figure
    set(0, 'CurrentFigure', fig1)
    
    % Apply legends
    subplot(3,4,4)
    legend(dampText)
    legend('boxoff');
    if (targetDirNum == 1 || targetDirNum == 3)
        legend('Location', 'southeast');
    else
        legend('Location', 'northeast');
    end
    
    % Apply legends
    subplot(3,4,8)
    legend(dampText)
    legend('boxoff');
    if (targetDirNum == 1 || targetDirNum == 3)
        legend('Location', 'southeast');
    else
        legend('Location', 'northeast');
    end

    % Apply legends
    subplot(3,4,12)
    legend(dampText)
    legend('boxoff');
    if (targetDirNum == 1 || targetDirNum == 3)
        legend('Location', 'southeast');
    else
        legend('Location', 'northeast');
    end
    xlabel('Time [ms]', 'Interpreter', 'latex');
    
    
    % Set X - YLims Equal
    axes_arr = [subplot(3,4,1), ...
                subplot(3,4,2), ...
                subplot(3,4,3), ...
                subplot(3,4,4)];
    SetYLimsEqual(axes_arr, ...
                 'TopPadPercent', 0.1, ...
                 'BottomPadPercent', 0.1);

    % Set VA - YLims Equal
    axes_arr = [subplot(3,4,5), ...
                subplot(3,4,6), ...
                subplot(3,4,7), ...
                subplot(3,4,8)];
    SetYLimsEqual(axes_arr, ...
                 'TopPadPercent', 0.1, ...
                 'BottomPadPercent', 0.1);
    
    % Set Damping - YLims Equal
    axes_arr = [subplot(3,4,9), ...
                subplot(3,4,10), ...
                subplot(3,4,11), ...
                subplot(3,4,12)];
    SetYLimsEqual(axes_arr, ...
                 'TopPadPercent', 0.05, ...
                 'BottomPadPercent', 0.05);    
             
             
    if (showMovementMatch)
        % Set fig1 as current figure
        set(0, 'CurrentFigure', fig2)

        % Apply legends
        subplot(3,4,4)
        legend(dampText)
        legend('boxoff');
        if (targetDirNum == 1 || targetDirNum == 3)
            legend('Location', 'southeast');
        else
            legend('Location', 'northeast');
        end

        % Apply legends
        subplot(3,4,8)
        legend(dampText)
        legend('boxoff');
        if (targetDirNum == 1 || targetDirNum == 3)
            legend('Location', 'southeast');
        else
            legend('Location', 'northeast');
        end

        % Apply legends
        subplot(3,4,12)
        legend(dampText)
        legend('boxoff');
        if (targetDirNum == 1 || targetDirNum == 3)
            legend('Location', 'southeast');
        else
            legend('Location', 'northeast');
        end
        xlabel('Time [ms]', 'Interpreter', 'latex');


        % Set X - YLims Equal
        axes_arr = [subplot(3,4,1), ...
                    subplot(3,4,2), ...
                    subplot(3,4,3), ...
                    subplot(3,4,4)];
        SetYLimsEqual(axes_arr, ...
                     'TopPadPercent', 0.1, ...
                     'BottomPadPercent', 0.1);

        % Set VA - YLims Equal
        axes_arr = [subplot(3,4,5), ...
                    subplot(3,4,6), ...
                    subplot(3,4,7), ...
                    subplot(3,4,8)];
        SetYLimsEqual(axes_arr, ...
                     'TopPadPercent', 0.1, ...
                     'BottomPadPercent', 0.1);

        % Set Damping - YLims Equal
        axes_arr = [subplot(3,4,9), ...
                    subplot(3,4,10), ...
                    subplot(3,4,11), ...
                    subplot(3,4,12)];
        SetYLimsEqual(axes_arr, ...
                     'TopPadPercent', 0.05, ...
                     'BottomPadPercent', 0.05);    
    end
                 
                 
    res = [fig1];
end