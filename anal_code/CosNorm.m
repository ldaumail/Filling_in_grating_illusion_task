function nor = CosNorm(data,params)

a = params(1);
b = params(2);

y_model = a*cos(b*data.x);
y_actual = data.y;
nor = norm(y_model-y_actual);
end