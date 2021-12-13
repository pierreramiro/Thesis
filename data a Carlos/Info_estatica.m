clear;close all
%% Regeneración mina
load('malla_triangular_sin_trasl')
load('nube_puntos_sin_trasl')
trimesh(T,temp(:,1),temp(:,2),temp(:,3))
load('esfera_mina_sin_trasl')%% data de esfera
scatter3(temp(:,1),temp(:,2),temp(:,3))
%% Regeneracion ideal
load('malla_triang_ideal')
load('esfera_ideal')
trimesh(T,temp(:,1),temp(:,2),temp(:,3))
