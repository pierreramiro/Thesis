function y=line_equation(x,parameters)
m=parameters(1);
x_point=parameters(2);
y_point=parameters(3);
y=m*(x-x_point)+y_point;
end