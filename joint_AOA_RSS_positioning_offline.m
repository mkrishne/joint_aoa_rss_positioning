function [rp_idx,AOA_fngprnt_centroids] = joint_AOA_RSS_positioning_offline(AOA_fngprnt,num_clusters)
    dims = size(AOA_fngprnt);
    M = dims(1);
    N = dims(2);
    K = dims(3);
    AOA_fngprnt_reshaped = reshape(AOA_fngprnt,N*M,K).';
    % Perform K-means clustering
    [rp_idx, centroids] = kmeans(AOA_fngprnt_reshaped, num_clusters);
    AOA_fngprnt_centroids = reshape(centroids.', M, N, num_clusters);
end