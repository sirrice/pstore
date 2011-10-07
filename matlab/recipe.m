function [] = recipe()
f = @(x)x*[100, 115]'

%x = [b, vs]

x_L = [0,0]
x_U = [100,100]
A   = [1,1; 6,1; 1,4]
b_L = [0, 34, 44]
b_U = [25, 500, 500]
c   = []
c_L = []
c_U = []

Prob = glcAssign('f', x_L, x_U, 'Name', A, b_L, b_U, c, c_L, c_U)
result = tomRun('minos', Prob, 3)
result.x_k

%b = 6s + 1p
%vs = 1s + 4p
%b+v <= 25
%6b + 1vs > 34
%1b + 4vs > 44

%min 100b + 115vs



// function f = objfun(x)
// f = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1);



x0 = [10,10];     % Make a starting guess at the solution
options = optimset('Algorithm','active-set');
[x,fval] = fmincon(@f,x0,[],[],[],[],x_L,x_U,@confun,options);


25 < x
-25 > -x
-x < -25
-x + 25 < 0
