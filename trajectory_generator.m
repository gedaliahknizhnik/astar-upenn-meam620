function [desired_state] = trajectory_generator(t, qn, varargin)
% TRAJECTORY_GENERATOR: Turn a Dijkstra or A* path into a trajectory
% t: A scalar, specifying inquiry time
%
% varargin: variable number of input arguments. In the framework,
% this function will first (and only once!) be called like this:
%
% trajectory_generator([],[], 0, path)
%
% i.e. map = varargin{1} and path = varargin{2}.
%
% path: A N x 3 matrix where each row is (x, y, z) coordinate of a
% point in the path. N is the total number of points in the path
%
% This is when you compute and store the trajectory.
%
% Later it will be called with only t and qn as an argument, at
% which point you generate the desired state for point t.
%

% Declare persistent variables
persistent waypoints timesteps coeffsMatrix t_max

% Maximum allowable acceleration (experimentally determined)
maxima.a = 10;

% Output variable
desired_state = [];

% During the first call, calculate the polynomial coefficients
if nargin > 2

    % Store the map
    map = varargin{1};
    % Store the path as waypoints
    waypoints  = varargin{2};
    % Optimize the Dijkstra/A* path by ray-tracing
    waypoints = optimizePath(waypoints, map);
    % Calculate the time-steps for te trajectory.
    [timesteps, t_max] = getTimeSteps(waypoints, maxima);
    fprintf('Planned trajectory will take %f s\n', t_max);

    % Calculate polynomial coefficients to use continuous cubic splines
    [coeffsMatrix] = getCoeffs(waypoints, timesteps);

% On anything but the first call, return the interpolated state
%   variables.
else
    % Limit time to [0,t_max]
    t = max(min(t,t_max),0);

    if (t == 0)
        pos = waypoints(1,:);
        vel = [0,0,0];
        acc = [0,0,0];           
        yaw = 0;
        yawdot = 0;
    elseif t== t_max
        pos = waypoints(end,:);
        vel = [0,0,0];
        acc = [0,0,0];           
        yaw = 0;
        yawdot = 0;
    else
        % Get index of bottom node relevant to t - i.e. the segment
        % relevant to the current time. The max function will return 1 only
        % for t=0, which is handled elsewhere. Otherwise, if we are in the
        % first segment it will return 2, and then we subtract 1 to find
        % the segment.
        [~,currind] = max(timesteps >= t); currind = currind - 1;

        % Plug into the polynomial equations to get the state values.
        [pos,vel,acc,yaw,yawdot] = getVals(t, coeffsMatrix(:,:,currind));

    end

    % Formate the desired state.
    desired_state.pos = pos(:);
    desired_state.vel = vel(:);
    desired_state.acc = acc(:);
    desired_state.yaw = yaw;
    desired_state.yawdot = yawdot;
end

end

%% HELPER FUNCTIONS

function [newpath] = optimizePath(path, map)
% OPTIMIZEPATH Optimize the path produced by 26-connected Dijkstra/A*. 
%   Starting at the start position, we attempt to draw a straight line to
%   the farthest possible point without causing a collision. The
%   intermediate points on the original path are then removed in favor of
%   this new line. The new line is represented as three points (start,
%   middle, end). The process is repeated until the goal is reached. The
%   search is performed starting with both the start and the goal, and the
%   shorter path is selected. While this is not guaranteed to be the most 
%   optimal path, it is certainly better than the original path.
%   I have found experimentally that 3 points provide enough connection to
%   the straight line to avoid collisions while allowing the quad enough
%   speed to be efficient.
%
% INPUTS
%   path - an Nx3 array of waypoints produced by DIJKSTRA
%   map  - a structure containing information about the map, including an
%           obstacle grid.
%   
% OUTPUTS
%   newpath - an Mx3 array of waypoints optimized from path. In the worst
%               case, this will just be path

    % EXPLORE FORWARDS
    
    % Let newpath be the start position
    newpathf = path(1,:);
        
    % Start with the first point
    jj = 0;
    canbebetter = 1;    % Improvement flag
    
    % While we can still improve the path
    while canbebetter
        
        canbebetter = 0;
        % Take the current node
        currnode = path(jj+1,:);
        % Search the rest of the path
        for ii = (jj+1):size(path,1)
            % Select the next node
            newnode = path(ii,:);

            % Linearly interpolate between the two nodes with a density
            % twice that of the grid (to make sure we don't miss an
            % obstacle). While this is expensive, it is not done in
            % real-time.
            list = interp1([1 2],[currnode;newnode],linspace(1,2,2*(ii-jj)));
            % Check for collision on any of the interpolated nodes
            inds = pos2ind(map, list);
            collided = map.occgrid(inds);         
            % If a collision is detected
            if (nnz(collided) == 0)
                % The relevant node is the last non-collision-inducing one.
                max = ii-1;
            else
                % Take the last non-collision-inducing node and exit the
                % search.
                newnode = path(max,:);
                break
            end
        end
        
        % If we got to the goal without hitting an obstacle - that is the
        % end of the ray.
        if ii ~= size(path,1)
            canbebetter = 1;
        end
        
        % Experimentally I have determined one intermediate node to be a
        % good compromise between speed and safety.
        howmany = 3;
        % Interpolate nodes
        list2add = interp1([1 2],[currnode;newnode],linspace(1,2,howmany));

        % Add all the nodes to the new path except the first (which is
        % already on it from the last ray).
        newpathf = [newpathf; list2add(2:end,:)];
        
        % Start with the current last safe node on the next round.
        jj = max;           
    end
    
    % Just in case, prune any repeated nodes (this creates a 0 dt, which
    % makes the coefficient matrix singular). We count backwards to prevent
    % indexing issues when the nodes are deleted.
    for kk = size(newpathf,1):-1:2
        if norm(newpathf(kk,:) - newpathf(kk-1,:)) < 10e-5
            newpathf(kk,:) = [];
        end
    end    
    
    % EXPLORE BACKWARDS
    
    % Let newpathb be the end position
    newpathb = path(end,:);
    % Start with the first point
    jj = size(path,1);
    canbebetter = 1;    % Improvement flag
    
    % While we can still improve the path
    while canbebetter
        
        canbebetter = 0;
        % Take the current node
        currnode = path(jj-1,:);
        % Search the rest of the path
        for ii = (jj-1):-1:1
            % Select the next node
            newnode = path(ii,:);

            % Linearly interpolate between the two nodes with a density
            % twice that of the grid (to make sure we don't miss an
            % obstacle). While this is expensive, it is not done in
            % real-time.
            list = interp1([1 2],[currnode;newnode],linspace(1,2,2*(jj-ii)));
            % Check for collision on any of the interpolated nodes
            inds = pos2ind(map, list);
            collided = map.occgrid(inds);         
            % If a collision is detected
            if (nnz(collided) == 0)
                % The relevant node is the last non-collision-inducing one.
                max = ii+1;
            else
                % Take the last non-collision-inducing node and exit the
                % search.
                newnode = path(max,:);
                break
            end
        end
        
        % If we got to the goal without hitting an obstacle - that is the
        % end of the ray.
        if ii ~= 1
            canbebetter = 1;
        end
        
        % Experimentally I have determined one intermediate node to be a
        % good compromise between speed and safety.
        howmany = 3;
        % Interpolate nodes
        list2add = interp1([1 2],[newnode;currnode],linspace(1,2,howmany));

        % Add all the nodes to the new path except the first (which is
        % already on it from the last ray).
        newpathb = [list2add(1:(end-1),:); newpathb];
        
        % Start with the current last safe node on the next round.
        jj = max;           
    end
    
    % Just in case, prune any repeated nodes (this creates a 0 dt, which
    % makes the coefficient matrix singular). We count backwards to prevent
    % indexing issues when the nodes are deleted.
    for kk = size(newpathb,1):-1:2
        if norm(newpathb(kk,:) - newpathb(kk-1,:)) < 10e-5
            newpathb(kk,:) = [];
        end
    end    
    
    % Compare the lengths of the two paths
    distf = sum(vecnorm(newpathf,2,2));
    distb = sum(vecnorm(newpathb,2,2));
    
    % Choose the shorter one
    if distf <= distb
        newpath = newpathf;
    else
        newpath = newpathb;
    end
end

function [timesteps, t_max] = getTimeSteps(waypoints, maxima)
% GETTIMESTEPS Calculate the timesteps for the given trajectory under the
%   assumption of a constant acceleration/deceleration specified in 
%   maxima. Although we do not require the quadrator to start and stop each
%   segment at v=0, this approach allows each segment time to accelerate
%   at the beginning of the trajectory and decelerate at the end.

    % Initialize timesteps. First one will stay zero.
    timesteps = zeros(size(waypoints,1),1);
    
    % For each waypoint after the initial position
    for ii=2:size(waypoints,1)
        % Find the distance to the next waypoint
        dist = norm(waypoints(ii,:) - waypoints(ii-1,:));
        % Calculate the dt according to d = 1/2*a*dt^2, where we want to
        % get to d/2 at t/2. Thus dt = 2*sqrt(d/a). Credit to Luca Scheuer
        % for this approach.
        dt = 2*sqrt(dist/maxima.a);
        % Find the timestep.
        timesteps(ii) = timesteps(ii-1) + dt;
                            
    end
    
    % Output the max time
    t_max = timesteps(end);
                 
end

function [coeffsMatrix] = getCoeffs(waypoints, timesteps)
% GETCOEFFS Calculate the coefficients for the cubic splines that make up
%   the trajectory. We specify the position and velocity of the first and
%   last points, the position of all intermediate waypoints, and that the
%   velocity and acceleration at intermediate waypoints must be continuous.

    % There is 1 fewer segment than there are waypoints
    numSegments = size(waypoints,1) - 1;

    % Matrix of time variables X*c = V, where c is the coefficient vector
    % we are solving for.
    xMatrix = zeros(4*numSegments,4*numSegments);
    % xVec encodes the contraints.
    xVec = zeros(4*numSegments,3);

    % Row coordinate in matrix
    yloc = 1;

    % Position Constraints for all waypoints
    for ii=1:(size(waypoints,1)-1)
        % Beginning of spline
        subMatrix  = [timesteps(ii)^3, timesteps(ii)^2, timesteps(ii), 1];
        % End of spline
        subMatrix1 = [timesteps(ii+1)^3, timesteps(ii+1)^2, timesteps(ii+1), 1];

        % Column coordinate
        xind = 4*(ii-1)+1;

        % Place the polynomials in the matrix
        xMatrix(yloc,xind:(xind+3)) = subMatrix;
        xMatrix(yloc+1,xind:(xind+3)) = subMatrix1;

        % Place the position contraints in the vectors
        xVec(yloc,:) = waypoints(ii,:);
        xVec(yloc+1,:) = waypoints(ii+1,:);

        % Update the row location
        yloc = yloc+2;
    end

    % Velocity Constraints on the first and last waypoints
    xMatrix(yloc,1:4) = [3*timesteps(1)^2, 2*timesteps(1), 1, 0];
    xMatrix(yloc+1,(end-3):end) = [3*timesteps(end)^2, 2*timesteps(end), 1, 0];

    yloc = yloc + 2;

    % If there are intermediate waypoints
    if (size(waypoints,1) > 2)
        % For each intermediate waypoints
        for ii=1:(size(waypoints,1)-2)
            % Velocity
            subMatrix_vel = [3*timesteps(ii+1)^2, 2*timesteps(ii+1), 1, 0];
            % Acceleration
            subMatrix_acc = [6*timesteps(ii+1), 2, 0, 0];
            
            % Column position for end of current spline
            xind = 4*(ii-1)+1;
            % Column position for beginning of next spline
            xind1 = 4*ii + 1;

            % Place the polynomials in the matrix
            xMatrix(yloc,xind:(xind+3)) = subMatrix_vel;
            xMatrix(yloc,xind1:(xind1+3)) = -subMatrix_vel;
            xMatrix(yloc+1,xind:(xind+3)) = subMatrix_acc;
            xMatrix(yloc+1,xind1:(xind1+3)) = -subMatrix_acc;
            % The vector constraints remain zero since the constraint is
            % that they be equal, not what the equal.
            
            yloc = yloc+2;
        end
    end
    
    % Solve for the coefficients
    xCoeff = xMatrix\xVec;
    
    % Place them in a matrix [xCoeffs;yCoeffs;zCoeffs] with the third index
    % being the segment of the trajectory
    coeffsMatrix = zeros(3,4,numSegments);
    
    % Extract the values from columns (MATLAB's reshape does the wrong
    % ordering)
    for ii = 1:numSegments
        first = 4*(ii-1)+1; last = first + 3;
        coeffsMatrix(1,:,ii) = xCoeff(first:last,1);
        coeffsMatrix(2,:,ii) = xCoeff(first:last,2);
        coeffsMatrix(3,:,ii) = xCoeff(first:last,3);
    end
end

function [p, v, a, y, yd] = getVals(t, coeffs)
% GETVALS Find the state value at the given time by using continuous cubic 
%   splines. Adapted from a similar function by Jessica McWilliams, though
%   cubic splines and matrix multiplication is my own.
% 
% INPUTS
%   t       - the current time in [0,t_max].
%   coeffs  - the a 3x4 matrix containing the coefficients to the cubic
%             in order of descending power. First row is x, second y, third
%             z. 
%
% OUTPUTS
%   p  - 1x3 vector of positions
%   v  - 1x3 vector of velocities
%   a  - 1x3 vector of accelerations
%   y  - 1x1 value of yaw
%   yd - 1x1 value of yawdot

    % Time vectors for position, velocity, acceleration
    ts = [t^3; t^2; t^1; 1];
    dts = [3*t^2; 2*t; 1; 0];
    ddts = [6*t; 2; 0; 0];
    % Extract state by plugging into polynomial
    p = (coeffs*ts)';
    v = (coeffs*dts)'; 
    a = (coeffs*ddts)';
    y = 0;
    yd = 0;
end