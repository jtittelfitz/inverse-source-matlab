x = -10:0.1:10;
t = 0:0.1:10;

g = @(x) (cos(pi*x/2) + sin(pi*x/2)).*(x > -9 & x < 9);
%g = @(x) (-1).*(x > -2*pi & x < 0) + (1).*(x > 0 & x < 2*pi );
phi = @(x) (1 -0.7*x + 0.1*sin(5*x) + cos(x)).*(x > -pi & x < pi);
%phi = @(x) (x + cos(x)).*(x > -pi & x < pi);


u = zeros(1,length(x));
v = zeros(1,length(x));
mask = zeros(1,length(x));

for s = 0:1/20:3
    for j = 1:length(x);
        u(j) = integral(g,x(j) - s,x(j) + s);    
        v(j) = integral(phi,x(j) - s,x(j) + s);    
        mask(j) = (abs(0.5*g(x(j) + s) + 0.5*g(x(j)-s)) > 0) & (abs(u(j)) < 0.01);        
    end
    if s == pi; disp('t = pi'); end
    %sum(v.*g(x))
    subplot(1,1,1); plot(x,u,x,mask); axis([-10,10,-10,10]); title(sprintf('u(%2.2f)',s));
    %subplot(1,2,2); plot(x,0.5*g(x + s) + 0.5*g(x-s)); axis([-10,10,-5,5])
    %subplot(1,3,3); plot(x,mask); axis([-10,10,-2,2])
    pause(0.1)
end
hold off
subplot(1,1,1); plot(x,g(x)); axis([-10,10,-10,10])