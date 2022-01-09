clear;close all;clc
T=readmatrix("Two_donut_fill.csv")+1;
load('csv_data.mat');
for i=1:length(TwoDonutFill)
    if ~isequal(TwoDonutFill(i,:),T(i,:))
        fprintf("Este es el indice donde dejan de ser iguales: %d\n",i)
        break
    end
end