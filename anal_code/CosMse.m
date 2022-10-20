function mse = CosMse(data,params)

a = params(1);
b = params(2);

y_model = a*cos(b*data.x);
y_actual = data.y;
mse = mean((y_actual-y_model).^2);
end