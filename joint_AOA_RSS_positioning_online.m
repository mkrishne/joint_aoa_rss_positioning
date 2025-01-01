function [estimated_location_joint, estimated_location_aoa_only] = joint_AOA_RSS_positioning_online(AOA_fngprnt, AOA_test_vec, joint_aoa_rp_idx, AOA_fngprnt_centroids, num_clusters, Nc_sel, N_max_RP, RSS_fngprnt, RSS_test, RP_positions)
    dims = size(AOA_fngprnt);
    M = dims(1);
    N = dims(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Step 1: Compute the Angle Similarity Coefficient (AS) %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Reshape AOA_test_vec to match the dimensionality of centroids
    test_vec = reshape(AOA_test_vec, M * N, 1);
    % Initialize the similarity coefficient array
    angle_similarity = zeros(num_clusters, 1);

    % Compute the angle similarity for each centroid
    for c = 1:num_clusters
        centroid_vec = reshape(AOA_fngprnt_centroids(:, :, c), M * N, 1);
        % Using dot product and norm to calculate cosine similarity
        angle_similarity(c) = (test_vec' * centroid_vec) / (norm(test_vec) * norm(centroid_vec));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Step 2: Find Top Nc_sel clusters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select the clusters with the largest angle similarity coefficients:
    [~, sorted_indices] = sort(angle_similarity, 'descend');
    selected_clusters = sorted_indices(1:Nc_sel);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Step 3: Calculate Similarity for RPs in Selected Clusters %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Storage for angle similarities and corresponding RP indices
    global_similarity = [];
    global_RP_indices = [];

    for n = 1:Nc_sel
        cluster_idx = selected_clusters(n);
        RPs_in_cluster = find(joint_aoa_rp_idx == cluster_idx);
        num_RPs = length(RPs_in_cluster);
        local_similarity = zeros(num_RPs, 1);
        
        for i = 1:num_RPs
            RP_vec = reshape(AOA_fngprnt(:, :, RPs_in_cluster(i)), M * N, 1);
            local_similarity(i) = (test_vec' * RP_vec) / (norm(test_vec) * norm(RP_vec));
        end
        
        % Store the local similarities and RP indices globally
        global_similarity = [global_similarity; local_similarity];
        global_RP_indices = [global_RP_indices; RPs_in_cluster];
    end
    
    if length(global_RP_indices) < N_max_RP
        N_max_RP = length(global_RP_indices);  % Update N_max_RP to the length of global_RP_indices
    end

    % Find the global top N_max_RP RPs with the highest angle similarities
    [sorted_similarities, sort_indices] = sort(global_similarity, 'descend');
    top_global_RP_indices = global_RP_indices(sort_indices(1:N_max_RP));
    top_global_similarities = sorted_similarities(1:N_max_RP);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Step 4: Compute Reciprocal of Euclidean Distance for RSS %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    reciprocal_euclidean = zeros(N_max_RP, 1);
    for i = 1:N_max_RP
        RP_index = top_global_RP_indices(i);
        if RP_index == 0
            continue;
        end
        distance = norm(RSS_test - RSS_fngprnt(RP_index, :));
        reciprocal_euclidean(i) = 1 / distance;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% Step 5: Apply WKNN for Location Estimation  %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_RP_positions = real(RP_positions(top_global_RP_indices));
    y_RP_positions = imag(RP_positions(top_global_RP_indices));
    weights = reciprocal_euclidean .* top_global_similarities;  % Combine both metrics
    normalized_weights = weights / sum(weights);  % Normalize weights
    
    x_estimated_location = sum(normalized_weights .* x_RP_positions, 1);  % Weighted sum of coordinates
    y_estimated_location = sum(normalized_weights .* y_RP_positions, 1);  % Weighted sum of coordinates
    estimated_location_joint = [x_estimated_location y_estimated_location];

    weights = top_global_similarities;  % Combine both metrics
    normalized_weights = weights / sum(weights);  % Normalize weights
    
    x_estimated_location = sum(normalized_weights .* x_RP_positions, 1);  % Weighted sum of coordinates
    y_estimated_location = sum(normalized_weights .* y_RP_positions, 1);  % Weighted sum of coordinates
    estimated_location_aoa_only = [x_estimated_location y_estimated_location];
end
