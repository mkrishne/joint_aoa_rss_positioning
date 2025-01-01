rng(1);
diary myDiaryFile_cdf_err_with_diff_M_pow;
noiseVariancedBm = -96;
%center frequency of 2GHz
format long
fc = 2e9;
% M = Number of antennas per AP
M = [8,16,32,64,128];
len_M = length(M);
array_spacing = 0.5;
single_side_angle_spread = 4; %4 degree angular spread
sigma_sf = db2pow(8); %8 shadowing std of 8db 
%Height of BS (in meters)
h_BS = 10;
%Height of UT (in meters)
h_UT = 1.5;
squareLength = 100; %each side of a square simulation area 
nbrOfSetups = 100; %The number of random realizations of the AP/user locations
N = 10; %total number of APs in simulation area
eta = 5;%inter RP distance
RP_positions_per_row = floor(squareLength / eta) + 1; %creates a grid of (RP_positions_per_row x RP_positions_per_row) RP positions
K = RP_positions_per_row^2;
num_users = 20;
num_channel_realisations = 100;  % Number of channel realizations
%Total uplink transmit power per UE (mW)
p = 100;
%JOINT AOA RSS Positioning
num_clusters = 15
Nc_sel = 3
N_max_RP = 4

avg_positioning_err_joint_aoa_rss = zeros(1,len_M);
avg_positioning_err_joint_aoa_rss_aoa_only = zeros(1,len_M);

positioning_err_joint_aoa_rss = zeros(len_M,nbrOfSetups,num_users);
positioning_err_joint_aoa_rss_aoa_only = zeros(len_M,nbrOfSetups,num_users);

for m_idx = 1:len_M
    ant_number = M(m_idx)
    for setup_idx = 1:nbrOfSetups
        setup_idx
        beta_fngprnt = zeros(K,N); %25x9
        %JOINT AOA RSS Positioning
        Joint_AOA_fngprnt = zeros(M(m_idx),N,K);
        
        [RP_positions,TP_positions,AP_positions] = cell_free_layout_setup(RP_positions_per_row,squareLength,N,num_users);
        %functionPlotSetup(squareLength,RP_positions,AP_positions,TP_positions);
        for RP_idx = 1:K
            for AP_idx = 1:N
                %disp(['Running RP' num2str(RP_idx) ' and AP ' num2str(AP_idx)]);
                d_2D = abs(RP_positions(RP_idx) - AP_positions(AP_idx));
                if d_2D < 10
                    PL = -81.2; %PL = Path Loss
                elseif d_2D < 50
                    PL = -61.2 - 20 * log10(d_2D);
                else
                    PL = -35.7 - 35 * log10(d_2D) + sigma_sf * (randn + 1i * randn);
                end                
                beta_fngprnt(RP_idx,AP_idx) = PL;
                nom_azi_angle = rad2deg(angle(RP_positions(RP_idx)- AP_positions(AP_idx)));
                for ch_idx = 1:num_channel_realisations
                    % Generate the channel realization
                    h_nk = functionChannelEstimates(M(m_idx), beta_fngprnt(RP_idx, AP_idx), array_spacing, single_side_angle_spread, nom_azi_angle);
                    
                    % Take the FFT of the channel realization
                    G_k = fft(h_nk);
                    
                    % Compute the squared magnitude and accumulate it
                    if ch_idx == 1
                        squared_magnitudes = abs(G_k).^2;  % Initialize with the first realization
                    else
                        squared_magnitudes = squared_magnitudes + abs(G_k).^2;  % Accumulate for the average
                    end
                end
                Joint_AOA_fngprnt(:, AP_idx, RP_idx) = squared_magnitudes / num_channel_realisations;
            end %for AP_idx = 1:L
        end %for RP_idx = 1:numel(RP_positions)
    
    RSS_fngprnt_mW = M(m_idx)*p*db2pow(beta_fngprnt);
    RSS_fngprnt_dB = 10*log10(RSS_fngprnt_mW/100); 
    
    [joint_aoa_rp_idx, Joint_AOA_fngprnt_centroids] = joint_AOA_RSS_positioning_offline(Joint_AOA_fngprnt,num_clusters);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ONLINE STAGE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta_fngprnt_tp = zeros(num_users,N); %2000x10 matrix
    Joint_AOA_fngprnt_tp = zeros(M(m_idx),N,num_users);

    for TP_idx = 1:num_users
        TPDistances = abs(RP_positions(:) - TP_positions(TP_idx)); %distance of the TP from each of the 25 RPs; 25x1 vector
        for AP_idx = 1:N
            d_2D_TP = abs(TP_positions(TP_idx)- AP_positions(AP_idx));
            if d_2D_TP < 10
                PL = -81.2; %PL = Path Loss
            elseif d_2D_TP < 50
                PL = -61.2 - 20 * log10(d_2D_TP);
            else
                PL = -35.7 - 35 * log10(d_2D_TP) + sigma_sf * (randn + 1i * randn);
            end   
            beta_fngprnt_tp(TP_idx,AP_idx) = PL;
            nom_azi_angle = rad2deg(angle(TP_positions(TP_idx)- AP_positions(AP_idx)));
            for ch_idx = 1:num_channel_realisations
                % Generate the channel realization
                h_tp = functionChannelEstimates(M(m_idx),beta_fngprnt_tp(TP_idx,AP_idx),array_spacing,single_side_angle_spread,nom_azi_angle);
                
                % Take the FFT of the channel realization
                G_tp = fft(h_tp);
                
                % Compute the squared magnitude and accumulate it
                if ch_idx == 1
                    squared_magnitudes = abs(G_tp).^2;  % Initialize with the first realization
                else
                    squared_magnitudes = squared_magnitudes + abs(G_tp).^2;  % Accumulate for the average
                end
            end
            Joint_AOA_fngprnt_tp(:, AP_idx, TP_idx) = squared_magnitudes / num_channel_realisations;
        end %for AP_idx = 1:N
    end %TP_idx = 1:num_users
    
    RSS_tp_mW = M(m_idx)*p*db2pow(beta_fngprnt_tp); %RSS_tp_mW --> 16x9 vector
    RSS_tp_dB = 10*log10(RSS_tp_mW/100); %16x9, with TPs in the centre of the 16 subregions
    
    for TP_idx = 1:num_users
        per_tp_RSS = RSS_tp_dB(TP_idx,:); %per_tp_RSS --> 9x1 vector of the TP
        TRUE_TP_COORDS = [real(TP_positions(TP_idx)) imag(TP_positions(TP_idx))];

        [JOINT_AOA_RSS_PRED,JOINT_AOA_RSS_PRED_AOA_ONLY] = joint_AOA_RSS_positioning_online(Joint_AOA_fngprnt,Joint_AOA_fngprnt_tp(:,:,TP_idx),joint_aoa_rp_idx, Joint_AOA_fngprnt_centroids,num_clusters,Nc_sel,N_max_RP,db2pow(RSS_fngprnt_dB),db2pow(per_tp_RSS),RP_positions);
        dist_coords = [JOINT_AOA_RSS_PRED;TRUE_TP_COORDS];
        positioning_err_joint_aoa_rss(m_idx, setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        dist_coords = [JOINT_AOA_RSS_PRED_AOA_ONLY;TRUE_TP_COORDS];
        positioning_err_joint_aoa_rss_aoa_only(m_idx,setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
    end %TP_idx = 1:num_users
    
    end %setup_idx = 1:nbrOfSetups
    avg_positioning_err_joint_aoa_rss(m_idx) = mean(positioning_err_joint_aoa_rss(m_idx,:));
    avg_positioning_err_joint_aoa_rss_aoa_only(m_idx) = mean(positioning_err_joint_aoa_rss_aoa_only(m_idx,:));

   fprintf("joint aoa gpr : %f %f\n",avg_positioning_err_joint_aoa_rss(m_idx),avg_positioning_err_joint_aoa_rss_aoa_only(m_idx));
 
  end %m_idx = 1:len_M
diary off

% Define the font to be used
font = 'Arial';

% Create a figure and set default font properties for axes and text
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

hold on; 
box on;

% Define improved color scheme for clarity
colors = [
    0, 0.4470, 0.7410;   % Blue
    0.8500, 0.3250, 0.0980; % Orange
    0.9290, 0.6940, 0.1250; % Yellow
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.6350, 0.0780, 0.1840; % Maroon
];

% Line styles for variety
lineStyles = {'--', '-.', '-', '--', '-.', '-'};

% Loop through each M value and plot the CDF
for m_idx = 1:len_M
    % Get the positioning errors for the current M value
    positioning_errors = positioning_err_joint_aoa_rss(m_idx, :, :);
    
    % Flatten the matrix into a vector
    positioning_errors = reshape(positioning_errors, [], 1);

    % Sort the positioning errors in ascending order and calculate the CDF
    sorted_errors = sort(positioning_errors);
    cdf_values = linspace(0, 1, length(sorted_errors));

    % Plot the CDF for the current value of M
    plot(sorted_errors, cdf_values, 'Color', colors(m_idx, :), ...
        'LineStyle', lineStyles{m_idx}, 'LineWidth', 2, ...
        'DisplayName', sprintf('M = %d', M(m_idx)));
end

% Add legend with tex interpreter
legend('show', 'Interpreter', 'tex', 'Location', 'SouthEast', 'FontSize', 14);

% Set axis properties
set(gca, 'FontSize', 15);

% Add labels with tex interpreter
xlabel('Positioning Estimation Error (m)', 'Interpreter', 'tex', 'FontSize', 15, 'FontName', font);
ylabel('CDF', 'Interpreter', 'tex', 'FontSize', 15, 'FontName', font);

% Set axis limits
xlim([0 5]);

% Add grid
grid on;
set(gca, 'GridLineStyle', ':', 'GridColor', 'k', 'GridAlpha', 0.5);

% Force MATLAB to use 'painters' renderer
set(gcf, 'Renderer', 'painters');

% Save the figure in high resolution
print(gcf, 'positioning_error_cdf_diff_M', '-dpng', '-r300'); % Save as PNG

