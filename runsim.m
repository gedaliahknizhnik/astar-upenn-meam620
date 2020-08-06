close all;
clear all;
clc;
addpath(genpath('./'));

%% Plan path
disp('Planning ...');
map = load_map('map1.txt', 0.1, 1.0, 0.25);
start = {[0.0, -4.9, 0.2]};
stop = {[8.0, 18.0, 3.0]};
% map = load_map('map2.txt', 0.1, 1.0, 0.25);
% start = {[5.0, -4.9, 0.25]};
% stop = {[2.0, 30, 3.0]};
% map = load_map('map3.txt', 0.1, 1.0, 0.25);
% start = {[0, 0, 5.5]};
% stop = {[19.5, 4, 1.0]};
% map = load_map('maze.txt', 0.1, 0.5, 0.1);
% start = {[-1.5, 1, 1]};
% stop = {[2.5, -1.5, 1.5]};
% map = load_map('mymap.txt',0.1,0.5,0.25);
% start = {[5.0, -4.0, 0]};
% stop = {[0, 12, 2.5]};
% stop = {[9.8, 16, 5.2]};

nquad = length(start);
for qn = 1:nquad
    [path{qn}, num] = dijkstra(map, start{qn}, stop{qn}, true);
    num
end

pause (1)

if nquad == 1
    plot_path(map, path{1});
else
    % you could modify your plot_path to handle cell input for multiple robots
end

%% Generate trajectory
disp('Generating Trajectory ...');
trajectory_generator([], [], map, path{1})

%% Run trajectory
trajectory = test_trajectory(start, stop, map, path, true); % with visualization