function [RP_positions,TP_positions,AP_positions] = cell_free_layout_setup(RP_positions_per_row,squareLength,N,num_tp_points)
inter_RP_dist = squareLength/(RP_positions_per_row-1); %1000/20 => each RP spaced apart by 50
RP_positions_x = 1:RP_positions_per_row;
RP_positions_y = 1:RP_positions_per_row;
[RP_positions_X,RP_positions_Y] = meshgrid(RP_positions_x,RP_positions_y);
RP_positions_X = (RP_positions_X-1)*inter_RP_dist;
RP_positions_Y = (RP_positions_Y-1)*inter_RP_dist ;
RP_positions = RP_positions_X + 1i*RP_positions_Y;
RP_positions = RP_positions.'; % A.' will transpose a complex matrix. A' will provide complex conjugate. 

AP_positions = zeros(1,N);
%Random AP locations with uniform distribution
for AP_idx = 1:N
    AP_positions(AP_idx) =  (rand + 1i*rand) * squareLength;
end

%Random TP locations with uniform distribution
TP_positions = zeros(1,num_tp_points);
for TP_idx = 1:num_tp_points
    TP_positions(TP_idx) = (rand + 1i*rand) * squareLength;
end




