function functionPlotSetup(squareLength, RP_positions, AP_positions, TP_positions)
    persistent figCount
    if isempty(figCount)
        figCount = 1;
    end

    [X, Y] = meshgrid(0:100:squareLength);
    figure; hold on;
    
    % Plot grid lines
    plot(X, Y, ":k");
    plot(Y, X, ":k");
    
    % Extract real and imaginary parts for RP, AP, and TP positions
    RP_positions_X = real(RP_positions);
    RP_positions_Y = imag(RP_positions);
    
    AP_positions_X = real(AP_positions);
    AP_positions_Y = imag(AP_positions);
    
    TP_positions_X = real(TP_positions);
    TP_positions_Y = imag(TP_positions);
    
    % Plot RP positions (blue 'x')
    hRP = plot(RP_positions_X, RP_positions_Y, "xb", 'MarkerSize', 10, 'LineWidth', 1.5);
    
    % Plot AP positions (red triangle)
    hAP = plot(AP_positions_X, AP_positions_Y, "^r", 'MarkerFaceColor', "r", 'MarkerSize', 9, 'LineWidth', 2);
    
    % Plot TP positions (green cross)
    hTP = plot(TP_positions_X, TP_positions_Y, "+g", 'MarkerSize', 11.5, 'LineWidth', 2.5);
    
    % Add title and labels
    title('Positioning Simulation Setup - randomized RP');
    xlabel('Horizontal Distance (m)');
    ylabel('Vertical Distance (m)');
    
    % Add legend with markers and labels
    lgd = legend([hRP(1), hAP(1), hTP(1)], {'RP', 'AP', 'TP'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
    
    % Increase space between legend items
    lgd.ItemTokenSize = [40, 10]; % Adjust the size as needed
    
    % Adjust figure size
    set(gcf, 'Position', [100, 100, 800, 600]); % Width x Height in pixels
end
