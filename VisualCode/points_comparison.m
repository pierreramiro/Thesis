clear;close all;clc
load('csv_data.mat');
%%
T=readmatrix("Two_donut_fill.csv")+1;
for i=1:length(TwoDonutFill)
    if ~isequal(TwoDonutFill(i,:),T(i,:))
        fprintf("Two:Este es el indice donde dejan de ser iguales: %d\n",i)
        break
    end
end
%%
T=readmatrix("Tri_donut_fill.csv")+1;
for i=1:length(TriDonutFill)
    if ~isequal(TriDonutFill(i,:),T(i,:))
        fprintf("Tri:Este es el indice donde dejan de ser iguales: %d\n",i)
        break
    end
end
%%
T=readmatrix("Mid_donut_fill.csv")+1;
for i=1:length(TriMiddleFill)
    if ~isequal(TriMiddleFill(i,:),T(i,:))
        fprintf("Tri:Este es el indice donde dejan de ser iguales: %d\n",i)
        break
    end
end
%%
fprintf("Si no sale nada, es porque todo salió bien! =D\n")