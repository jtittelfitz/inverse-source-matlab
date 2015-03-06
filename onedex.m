x = linspace(0,10,1000);
h = x(2) - x(1);
k = h/5;
t = 0:k:2;
N = length(x);
c = 1 + 0.1*sin(5*x);%ones(N,1);%

I = speye(N);
e = (1/2)*ones(N - 1,1);
D = sparse(diag(e,1));
D = D - D';

D(1,1) = -3/2;
D(1,2) = 4/2;
D(1,3) = -1/2;
D(N, N - 2) = 1/2;
D(N, N - 1) = -4/2;
D(N, N) = 3/2;

e = ones(N - 1,1);
DR = (sparse(diag(e,1)) - I)/(h);
DL = (I - sparse(diag(e,-1)))/(h);

e = ones(N - 1,1);
DPH = 0.5*(sparse(diag(e,1)) + I);
DMH = 0.5*(sparse(diag(e,-1)) + I);

f = zeros(N,1);
g = zeros(N,1); %((x > 0.4) & (x < 0.6))';

extforce = zeros(N, 100000);

for j = 1:1000
    a = zeros(N,1);
    a(abs(x - 2 - 0.01*(j-1) ) < 0.5) = 1;
    %plot(x,a)
    %drawnow
    extforce(:, j) = a;
end

diff_op = DR*spdiags(c(:),0,N,N)*DL;
diff_op2 = DR*DL;

u = cell(3);

u{1} = f;
u{2} = f + k*g + 0.5*k^2*diff_op*u{1};

last = 0;
prev = 0;
vel = zeros(length(t),1);
disp = zeros(length(t),1);

for n = 1:length(t)
    u{3} = 2*u{2} - u{1} + k^2*diff_op*u{2} + k^2*extforce(:,n);    
    u{1} = u{2};
    u{2} = u{3};
    subplot(1,2,1); plot(x,u{3}); axis([0 10 -1 1]);
    v = (u{3} - u{1})/k;
    [max_v,max_x] = max(v);
    a = zeros(length(v),1);
    x_bar = 0;
    for i = 1:length(v)
        if (abs(v(i) - max_v) < 0.2*max_v)
            x_bar = x_bar + v(i)*x(i);
            a(i) = v(i);
        end        
    end
    x_bar = x_bar / sum(a(a>0));
    vel(n) = (3*x_bar - 4*prev + last)/(2*k);
    disp(n) = x_bar - prev;
    last = prev;
    prev = x_bar;    
     
    subplot(1,2,2); plot(x,v); hold on; plot(x_bar,max_v,'ro'); hold on; plot(x,max_v*extforce(:,n),'g'); axis([0 10 -2 2]); hold off; %plot(2 + 0.01*(n-1),max_v,'*');max_v*extforce(:,n)
    drawnow
end

clf; plot(t,vel); axis([0 2 0 10]);
    