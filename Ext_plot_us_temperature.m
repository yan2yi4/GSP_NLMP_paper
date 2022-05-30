% The plotting code based on the paper
% Spelta, M.J.M., Martins, W.A.: Normalized 
% LMS algorithm and data-selective strategies
% for adaptive graph signal estimation. Signal
% Processing 167(107326) (2020)

clc;clear;format shortEng;

addpath('/GraSP-master')
grasp_start

clc
K = 7;  
earthRadius = 6360; 
[numbDataTable, textDataTable, rawDataTable] = ...
    xlsread("extracted",2);
numberStations = size(numbDataTable,1);
load('Real_graph_and_sampling.mat')
latitudeVector = numbDataTable(:,1); 
longitudeVector = numbDataTable(:,2); 

coordinates_matrix = [longitudeVector latitudeVector];
load('Ext_Real_graph_and_sampling.mat')

graphTopology = grasp_plane_rnd( length(completeDatasetMatrix(:,1)));
graphTopology.A = A;
graphTopology.layout = [coordinates_matrix(:,1), coordinates_matrix(:,2)];

graphTopology.A_layout = 0;

x_offset = 2;
y_offset = 1;

limits =  min(min( completeDatasetMatrix(:,1) ) ) ;

figure;
graph_signal = completeDatasetMatrix(:,1);
graphTopology.A = A;
grasp_show_graph(gca, graphTopology, 'node_values', graph_signal,'node_display_size',300); 
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])

colorbar
box off
set(gca,'Visible','off')
colorbar
colormap('turbo')
