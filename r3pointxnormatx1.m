n=1000;
a=-5;
b=5;
x = linspace(a,b,n);
[X,Y] = meshgrid(x,x);
f = @(x) (1-x.^2).^(1/4); %boundary fun
n = @(x) 1./sqrt(c-2*log(x)); %known normal function
u = @(x,y,j) pi*sin(x)+j*pi*sin(y); %value (needs to be experimented with)
c = -2;%see 16

%for doing derivatives symbolically
syms x y j
u2 = pi*sin(x)+j*pi*sin(y); 
n2 = 1./sqrt(c-2*log(x));

%testing to make sure we have orthagonality and to find good values of c
% c seems to range from about -3 to about 3
% for c = -3:0.2:5
%     pause(3);
%     n = @(x) 1./sqrt(c-2*log(x)); %known function of normal vectors
%     plot(x,n(x),x,f(x));
% end
%plot(x,n(x),x,f(x));

%calculating intersection
syms k
xint = vpasolve((1-k^2)^(1/4) == 1/sqrt(c-2*log(k)),k,[-5 5]);
yint = f(xint);
int = [xint,yint];

%calculating actual normal vector
nprime = diff(n2);
slope = subs(nprime,x,xint);
V = int + [1,slope];
normal = V/norm(V);

%differentiate
Ux = diff(u2,x);
Uy = diff(u2,y);

%substitute intersection values
Ux = subs(Ux,x,xint);
Uy = subs(Uy,y,yint);
Uy = subs(Uy,j,closestj);
jrange = linspace(-2*pi,2*pi,2000);
%Find constant j such that [Ux,Uy] dot norm = 0
%find the closest j in jrange
% closestj = 0;
% err = 1000;
% for i = jrange
%     if(abs(double(Ux*normal(1)+subs(Uy,j,double(i))*normal(2)))<double(err))
%         closestj = i;
%         err = abs(double(Ux*normal(1)+subs(Uy,j,double(i))*normal(2)));
%     end
% end
% double(err)
closestj= -0.4055;
ughost = [];
err = [];
delta = 1;
while(delta>.0001)
    %create meshgrid based around delta/2 distance from target edge point
    center = int - [0,delta/2];
    yt = linspace(center(2)-delta,center(2)+delta,3);
    xt = linspace(center(1)-delta,center(1)+delta,3);
    [XT,YT] = meshgrid(xt,yt);
%     x2 = linspace(a,b,1000);
%     plot(XT,YT,x2,n(x2),x2,f(x2));


    %gather function values at meshgrid points
    U = u(XT,YT,closestj);
    
    %save u12 for error
    u12 = double(U(1,2));

    %calculate ustar
    Ustar = (3*Uy/(8*Ux))*(double(U(2,1))-double(U(2,3)))+3*double(U(2,2))/4 +double(U(3,2))/4;

    %solve for quadratic fit
    Amatrix = [xt(2)^2, xt(2), 1;
               xt(3)^2, xt(3), 1;
              (xt(2)-delta/2)^2, xt(2)-delta/2, 1;
               ];
    rhs = [double(U(2,2)), double(U(3,2)), double(Ustar)]';
    sol = Amatrix\rhs;
    quad = @(x)sol(1)*x.^2 + sol(2)*x +sol(3);

    %use quadratic fit to find Ughost
    U(1,2) = quad(xt(1)); 
    %measure how close we are and calculate error in deltaloop



    %measure how close we are to actual fn value outside of region
   
    err = [err; u12-double(U(1,2))];
     %test1: 1st difference normal, 3ptx2pt
    test1 = [test1;double(Ux*(U(2,2)-U(1,2))/(delta)+Uy*(U(2,1)-U(2,3))/(2*delta))];

    %test2: approximate laplacian
    alap = [alap;-double((double(U(1,2))-2*double(U(2,2))+double(U(3,2)))/delta^2 + (double(U(2,1))-2*double(U(2,2))+double(U(2,3)))/delta^2)];

    %increment delta
    delta = delta/2;
end


order = log(err(1:end-1)./err(2:end))/log(2)
order1 = log(test1(1:end-1)./test1(2:end))/log(2)

%test2: compare lap to alap: approximate laplacian of u @ x1 using 3 point central difference
Uxx = diff(diff(u2,x), x);
Uyy = diff(diff(u2,y), y);
lap = double(subs(Uxx, x, xt(2)) + subs(subs(Uyy, y, yt(2)),j,closestj));
disp("alap - lap =" + string(double(alap(length(alap)-1)-lap)));
disp("with order:");
order2 = log(alap(1:end-1)./alap(2:end))/log(2)
