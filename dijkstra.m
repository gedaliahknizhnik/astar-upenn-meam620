function [path, num_expanded] = dijkstra(map, start, goal, astar)
% DIJKSTRA Find the shortest path from start to goal.
%   PATH = DIJKSTRA(map, start, goal) returns an mx3 matrix, where each row
%   consists of the (x, y, z) coordinates of a point on the path. The first
%   row is start and the last row is goal. If no path is found, PATH is a
%   0x3 matrix. Consecutive points in PATH should not be farther apart than
%   neighboring voxels in the map (e.g. if 5 consecutive points in PATH are
%   co-linear, don't simplify PATH by removing the 3 intermediate points).
% 
%   PATH = DIJKSTRA(map, start, goal, astar) finds the path using euclidean
%   distance to goal as a heuristic if astar is true.
%
%   [PATH, NUM_EXPANDED] = DIJKSTRA(...) returns the path as well as
%   the number of nodes that were expanded while performing the search.
%   
% paramaters:
%   map     - the map object to plan in
%   start   - 1x3 vector of the starting coordinates [x,y,z]
%   goal:   - 1x3 vector of the goal coordinates [x,y,z]
%   astar   - boolean use astar or dijkstra
%
% AUTHOR - This function is written by Gedaliah Knizhnik, based on the code
%   skeleton created by MEAM620 staff.

%% Prep Code

if nargin < 4
    astar = false;
end

path = [];
num_expanded = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  Algortihm Starts Here             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings

% Save info about the world to the map structure.
map.maxy = size(map.occgrid,1);
map.maxx = size(map.occgrid,2);
map.maxz = size(map.occgrid,3);
map.maxima = [map.maxy;map.maxx;map.maxz];
map.num_nodes = map.maxy*map.maxx*map.maxz;
map.num_free = nnz(~map.occgrid);

% Variables for calculating the cost of travel to neighbors.
map.info.xy = map.res_xyz(1); map.info.z = map.res_xyz(3); 
map.info.ratio = map.info.z/map.info.xy;

map.info.neighborlists = containers.Map('KeyType','int32','ValueType','any');
map.info.neighborcostlists = containers.Map('KeyType','int32','ValueType','any');

%% Check for valid start and goal

valid = validStartGoal(map,start,goal);

if ~valid
    return
end

%% Data Structures

% Initialize obstacles matrix - infinity on obstacles, zero elsewhere.
%   Arranged by index.
world.obstacles = reshape(Inf*(map.occgrid==1),[],1);
world.obstacles(isnan(world.obstacles)) = 0;

% Initial 1D cost matrix - all points initially have infinite cost -
%   and matrix to store previous cost. Arranged by index.
world.cost = Inf*ones(size(world.obstacles));
world.cost_prev = Inf*ones(size(world.cost));

% Store parent information in a nx1 vector where the value corresponds to
%   the index of the parent node (initially all zero since no zero index).
world.parent = zeros(size(world.cost));

% Store status of explored nodes. 1 for not explored and 0 for explored.
world.explored = ones(size(world.cost));

% Calculate heuristics.
world.heuristics = heuristic(map, goal, astar);

%CHANGE CHANGE CHANGE
world.neighbors2search = world.cost;

%% Implement Dijkstra's Algorithm or A*

% Find the index corresponding to the start position.
startind = pos2ind(map, start);

% Find the index of the goal node.
goalind = pos2ind(map, goal);

% Set the cost of the start position to 0 and the parent index of the start
%   node to itself. Add the start position to the path.
world.cost(startind) = 0;
world.cost_prev(startind) = 0;
world.parent(startind) = startind;

% Set the current node as the start node and initialize minimum cost as
%   zero. Initialize array of current indices to look at.
currind = startind;

% Continue expanding nodes as long as we haven't found the goal node.
while (world.parent(goalind) == 0)
     
    % Set current node as explored in both explored and neighbors list.
    world.explored(currind) = 0;
    world.neighbors2search(currind) = Inf;
        
    % Get list of neighbor indices and their associated costs.
    [neighbors, neighborcosts] = getNeighbors(currind, map.maxima, map.info);

    
    % Add distance cost if this will result in a cheaper result but leave
    %   infinity cost if it's an obstacle. 
    world.cost(neighbors) = max(min(world.cost(currind) + neighborcosts, ...
                                    world.cost_prev(neighbors)), ...
                                world.obstacles(neighbors));                        
    
    % Add cost of only those nodes that have changed to the list of nodes 
    %   to be searched.
    world.neighbors2search(neighbors(world.explored(neighbors) ~= 0)) = ...
        world.cost(neighbors(world.explored(neighbors) ~= 0)) ...
        + world.heuristics(neighbors(world.explored(neighbors) ~= 0));
    
    % Change parent of nodes whose cost changed. We retain the original 
    %   parents by elementwise multiplying with the matrix of costs that
    %   stayed the same and then sum with the parent for the costs that
    %   changed, which will have zero for all those that remained. The
    %   first multiplicaton accounts for costs that decreased from already
    %   having a parent.
    world.parent(neighbors) = ...
        world.parent(neighbors).*(world.cost(neighbors) == world.cost_prev(neighbors)) + ...
        currind*(world.cost(neighbors) ~= world.cost_prev(neighbors));

    % Adjust the previous costs vector to maintain continuity.
    world.cost_prev(neighbors) = world.cost(neighbors);
    
    % Get a vector of all nodes with the next minimum cost
    currind = getNextNode(world.neighbors2search);
    
    % Update number of expanded nodes 
    num_expanded = num_expanded + 1;
    
    if (num_expanded > map.num_nodes)
        fprintf('\n\n Explored all nodes. No path found. \n\n');
        return
    end
        
end


%% Generate and plot the path

% If we got this far without returning, a path exists from the start to the
%   goal. We now extract it from the parent vector.

path = getPath(map, start, goal, world.parent);
plot_path(map, path);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS 


function [neighbors, neighborcosts] = getNeighbors(currind, maxima, info)
% GETNEIGHBORS - find the indices and costs of nodes that are neighbors to
%    the current node being analyzed. Accounts for faces, corners, etc. 
% 
% INPUTS:
%   currind: index of current node under observation.
%   maxima:  3x1 vector of max values, defining the dimension of the array.
%       Ordered as (maxy;maxx;maxz)
%   info:    a structure containing distance cost info.
%
% OUTPUTS:
%   neighbors:     a Nx1 matrix of neighbors. At minimum it will be 8 items
%       (corners), at max it will be 27. The node itself is included.
%   neighborcosts: a Nx1 matrix of costs corresponding to the neighbor
%       nodes. Indices match, and cost for the node itself will be zero.

    % Get location of current node
    [J,I,K] = ind2sub(maxima, currind);
    currjik = [J;I;K];
    
    % Determine whether we're at a corner, an edge, a face, or in the
    % space. The edge case where a dimension is reduced is accounted for by
    % the last statement.
    idval = [1,10,100]*(currjik <= [1;1;1]) ...
                + [2,20,200]*(currjik >= maxima) ...
                + [4,40,400]*(currjik > [1;1;1] & currjik < maxima);
        
    if ~isKey(info.neighborlists, idval)
        vals = [0,1].*(currjik <= [1;1;1]) + [-1,0].*(currjik >= maxima) ...
            + [-1,1].*(currjik > [1;1;1] & currjik < maxima) ...
            - [-1,1].*(maxima == [1;1;1;]); 

        % Transform from a 1x2 vector in vals to a 1x2 or 1x3 vector, adding in
        % the cost to move. The cost to move in y is 1 since index numbering
        % goes down the columns. The cost to move in x is maxy, since we have
        % to move a full column over. The cost to move in z is maxx*maxy.
        jvals = vals(1,1):vals(1,2);
        ivals = maxima(1)*(vals(2,1):vals(2,2));
        kvals = maxima(1)*maxima(2)*(vals(3,1):vals(3,2));

        % Use the cost matrices to compute the cost of motion along an axis.
        costjvals = info.xy*(vals(1,1):vals(1,2));
        costivals = info.xy*(vals(2,1):vals(2,2));
        costkvals =  info.z*(vals(3,1):vals(3,2));

        % Compute the neighbors by adding column matrices to row matrices.
        % Credit to Luca Scheuer for showing me this approach.
        flatvals = jvals' + ivals;
        cubevals = flatvals(:) + kvals;
        info.neighborlists(idval) = cubevals(:);

        % Apply the same approach to the cost, but square and square-root to
        % calculate distance.
        costflatvals = ((costjvals').^2 + costivals.^2).^0.5;
        costcubevals = ((costflatvals(:)).^2 + costkvals.^2).^0.5;
        info.neighborcostlists(idval) = costcubevals(:);
    end
    
    neighbors = info.neighborlists(idval) + currind;
    neighborcosts = info.neighborcostlists(idval);
    
end

function [distances] = heuristic(map, goal, astar)
% HEURISTIC - calculates the heuristic to apply to all the nodes to choose
%   the closes one. For Dijkstra this is zero, for A* its the Manhattan
%   distance. Although the Manhattan distance is not technically admissible
%   as a heuristic for a 26-connected grid, it produces paths with fewer
%   turns and significantly reduces run-time. Given the post-processing
%   that happens during the trajectory generation, this is considered ok.
%
% INPUTS:
%   map:   structure containing map info
%   goal:  (x,y,z) of the goal node.
%   astar: flag for whether A* algorithm is used (yes if astar=true).
% OUTPUTS:
%   distances: a Nx1 vector of distances from the indexed node to the goal

    % For all the logical indices.
    indices = [1:map.num_nodes]';

    if (astar)
        % Compute the euclidian distance heuristic.
        %distances = vecnorm(ind2pos(map,indices) - goal,2,2);
        % Calculate the Manhattan distance heuristic. 
        distances = vecnorm(ind2pos(map,indices) - goal,1,2);
    else
        distances = 0*indices;
    end
    
end

function [nextinds] = getNextNode(neighbors2search)
% GETNEXTNODE - find the set of nodes that are next in line to be checked
%   (i.e. the set of unexplored nodes with the lowest cost).
%
% INPUTS:
%   cost:       the Nx1 vector of costs by index.
%   heuristics: the Nx1 vector of heuristics by index.
%   explored:   the Nx1 vector of explored state by index (0 explored, 1 not)
%   currinds:   the previous set of nodes we were checking.
%
% OUTPUTS:
%   nextinds: a vector of indices of the nodes that should be considered
%       next. Only empty if map has been searched completely w/out success.
%   err:      a flag for nextinds being an empty set.
        
        err = 0;
        
        % Find the minimum cost of an unexplored node.
        [~,nextinds] = min(neighbors2search);
                
end

function [path] = getPath(map, start, goal, parent)
% GETPATH - trace the parent of each node back from the goal to the start
%   to find the optimal path located by Dijkstra or A*.
%
% INPUTS:
%   map:    a struct with map info
%   start:  (x,y,z) of the start node
%   goal:   (x,y,z) of the goal node
%   parent: an Nx1 vector of parents for each node.
%
% OUTPUTS:
%   path: an Nx3 array of nodes corresponding to the path from start to
%       goal

    % Add the goal node to the path
    path = goal;
    pathind = pos2ind(map, goal);
    
    % Find the snapped node the goal corresponds to
    pathnode = ind2pos(map, pathind);

    % Until we reach the node that is its own parent (i.e. the snapped 
    %   start node)
    while (parent(pathind) ~= pathind)
        % Add this node to the path.
        path = [pathnode;path];

        % Switch to the parent node
        pathind = parent(pathind);
        pathnode = ind2pos(map, pathind);
    end

    % Add the actual start node to the path
    path = [start; path];
    % Move path nodes from their snapped corners to the middle of the
    % voxels they were in.
    path(2:(end-1),:) = path(2:(end-1),:)+ 1/2*map.res_xyz;

end

function [valid] = validStartGoal(map, start, goal)
% VALIDSTARTGOAL checks for valid start and goal positions, both for
%   positions outside the map and ones in obstacles. Invalid message will
%   be printed and the main loop will exit when the error flag is passed.
% 
% INPUTS:
%   map: a structure containing map details
%   start: (x,y,z) of the start position
%   goal: (x,y,z) of the goal position
%
% OUTPUTS
%   valid: a flag for validity - 1 for valid or 0 for not.
    
    % Assume valid
    valid = 1;
    
    % Find the voxel corresponding to the start position.
    startjik = pos2sub(map, start);
    
    % If the voxel is outside the map, return error
    if nnz((startjik < 1)|(startjik > map.maxima')) > 0
        fprintf("\nThe start position you requested is outside the map.\n");
        fprintf("\nPlease try again.\n");
        
        valid = 0;
        
        return
    end
    
    % Find the voxel of the goal node.
    goaljik = pos2sub(map, goal);

    % If the voxel is outside the map, return error
    if nnz((goaljik < 1)|(goaljik > map.maxima')) > 0

        fprintf("\nThe goal position you requested is outside the map.\n");
        fprintf("\nPlease try again.\n");
        
        valid = 0;
        
        return
    end
    
    % Check if start or end are in an obstacle
    startobstacle = map.occgrid(startjik(1), startjik(2),startjik(3));
    goalobstacle  = map.occgrid(goaljik(1),  goaljik(2), goaljik(3));

    % Print error message
    if (startobstacle || goalobstacle)
        if (startobstacle && ~goalobstacle)
            errorstring = "start position";
            isare = "is";
        elseif (~startobstacle && goalobstacle)
            errorstring = "goal position";
            isare = "is";
        else
            errorstring = "start and goal position";
            isare = "are";
        end

        fprintf("\nThe %s you requested %s located in an obstacle.\n\nPlease try again.\n", errorstring, isare);
        
        valid = 0;
        
        return
    end
end
