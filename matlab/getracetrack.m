function [target_arc,inner,outer,center] = getracetrack(PathX,PathY,width,N)
    n = min(N, length(PathX));
    dx = diff(PathX);
    dy = diff(PathY);
    distances = sqrt(dx.^2 + dy.^2);
    arc_length = [0; cumsum(distances)];
    total_length = arc_length(end);
    normalized_arc_length = arc_length / total_length;
    target_arc = linspace(0, 1, n);
    path_x_normalized = interp1(normalized_arc_length, PathX, target_arc);
    path_y_normalized = interp1(normalized_arc_length, PathY, target_arc);
    normalized_path = [path_x_normalized', path_y_normalized'];

    % Combine pathX and pathY into track_data
    track_data = normalized_path;
    track_xy = track_data;

    % Extend track data to form a loop (p4-p1-p2-p3-p4-p1)
    extended = zeros(size(track_xy, 1) + 2, 2);
    extended(2:end-1, :) = track_xy;
    extended(1, :) = track_xy(end, :);
    extended(end, :) = track_xy(1, :);
    
    % Width of the track
    outer_d = -width / 2;
    inner_d = width / 2;
    
    outer_track = [];
    inner_track = [];

    for it = 2:length(extended)-1
        % Take triplets p1-p2-p3
        p1 = extended(it-1, :);
        p2 = extended(it, :);
        p3 = extended(it+1, :);
        p1d = p1;
        p3d = p3;

        for j = 1:10
            p1d = (p1d + p2) / 2;
            p3d = (p2 + p3d) / 2;
            p2 = (p1d + p3d) / 2;
        end
        
        % Angle made by line perpendicular to p1-p3
        perp_ang = atan2(-(p3d(1) - p1d(1)), (p3d(2) - p1d(2)));
        
        % End points of the perpendicular
        outer_p2d_x = p2(1) + outer_d * cos(perp_ang);
        outer_p2d_y = p2(2) + outer_d * sin(perp_ang);
        outer_track = [outer_track; outer_p2d_x, outer_p2d_y];
        
        inner_p2d_x = p2(1) + inner_d * cos(perp_ang);
        inner_p2d_y = p2(2) + inner_d * sin(perp_ang);
        inner_track = [inner_track; inner_p2d_x, inner_p2d_y];
    end
    inner = inner_track;
    outer = outer_track;
    center = track_data;
end
