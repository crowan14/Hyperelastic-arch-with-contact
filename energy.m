function [ obj , dir ] = energy( p )

%hyperelasticity

global L1;
global L2;
global L3;

global C1; global D1;

global pts;
global tx; global dxt; 
global cx; global dxc;
global r2x; global r2y; global dx2; global dy2;
global r3x; global r3y; global dx3; global dy3;
global r4x; global r4y; global dx4; global dy4;

global traction;

global pen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%substitute parameters into symbolic functions
u = @(x1,x2) ( N( [x1;x2] , p ) );
jj = @(x1,x2) ( J( [x1;x2] , p ) );
i1 = @(x1,x2) ( I1( [x1;x2] , p ) );

bulk = @(x1,x2) ( C1*( i1(x1,x2) - 2 - 2 * log( jj(x1,x2) ) ) + D1*( jj(x1,x2) - 1 )^2 );
work = @(x1) ( traction(x1,L2)' * u(x1,L2) );

gap = @(x1) ( L2 - L3 + [0,1]*u(x1,L2-L3) );

obj = 0;

%integration of traction boundary
for i=1:(length(tx)-1)
    obj = obj - dxt * work( tx(i)+dxt/2 );
end

%integration of contact boundary
for i=1:(length(cx)-1)
    obj = obj + dxc * pen * I( gap( cx(i)+dxc/2 ) )^2;
end

%integration of region 2
for i=1:(length(r2x)-1)
    for j=1:(length(r2y)-1)
        obj = obj + dx2*dy2*bulk( r2x(i)+dx2/2 , r2y(j)+dy2/2 );
    end
end


%integration of region 3
for i=1:(length(r3x)-1)
    for j=1:(length(r3y)-1)
        obj = obj + dx3*dy3*bulk( r3x(i)+dx3/2 , r3y(j)+dy3/2 );
    end
end

%integration of region 4
for i=1:(length(r4x)-1)
    for j=1:(length(r4y)-1)
        obj = obj + dx4*dy4*bulk( r4x(i)+dx4/2 , r4y(j)+dy4/2 );
    end
end

% %print objective value
% disp(obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%substitute parameters into symbolic functions
j_p = @(x1,x2) ( J_p( [x1;x2] , p ) );
i1_p = @(x1,x2) ( I1_p( [x1;x2] , p ) );
displacement_grad = @(x1,x2) ( u_p( [x1;x2] , p ) );

bulk_grad = @(x1,x2) ( C1*( i1_p(x1,x2) - ( 2/jj(x1,x2) )*j_p(x1,x2) ) + D1*2*( jj(x1,x2)-1 )*j_p(x1,x2) );
work_grad = @(x1) ( traction(x1,L2)' * displacement_grad(x1,L2) );

dir = zeros( length(p) , 1 );

%integration of traction boundary
for i=1:(length(tx)-1)
    dir = dir - dxt * work_grad( tx(i)+dxt/2 )';
end

%integration of contact boundary
for i=1:(length(cx)-1)
    du_dp = displacement_grad( cx(i)+dxc/2 , L2-L3 );
    dir = dir + dxc * 2 * pen * I( gap( cx(i)+dxc/2 ) ) * dI( gap( cx(i)+dxc/2 ) ) * du_dp(2,:)';
end

%integration of region 2
for i=1:(length(r2x)-1)
    for j=1:(length(r2y)-1)
        dir = dir + dx2*dy2*bulk_grad( r2x(i)+dx2/2 , r2y(j)+dy2/2 );
    end
end

%integration of region 3
for i=1:(length(r3x)-1)
    for j=1:(length(r3y)-1)
        dir = dir + dx3*dy3*bulk_grad( r3x(i)+dx3/2 , r3y(j)+dy3/2 );
    end
end

%integration of region 4
for i=1:(length(r4x)-1)
    for j=1:(length(r4y)-1)
        dir = dir + dx4*dy4*bulk_grad( r4x(i)+dx4/2 , r4y(j)+dy4/2 );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plotting

% pgridx = linspace(0,L1,pts); pgridy = linspace(0,L2,pts);
% f1 = zeros(pts,pts); f2 = zeros(pts,pts); 
% 
% [ X , Y ] = meshgrid( pgridx,pgridy );
% 
% for i=1:pts
%     for j=1:pts
%         if or( or( X(i,j) < L3 , X(i,j) > (L1-L3) ) ,  Y(i,j) > (L2-L3) )
%             entry = u( X(i,j) , Y(i,j) );
% 
%             f1(i,j) = entry(1);
%             f2(i,j) = entry(2);
%         end
%     end
% end
% 
% f1(f1==0) = nan; 
% f2(f2==0) = nan; 
% 
% subplot(1,2,1)
% surf( X , Y , f1 )
% title('U1')
% grid on
% xlabel('x1')
% ylabel('x2')
% drawnow;
% 
% subplot(1,2,2)
% surf( X , Y , f2 )
% title('U2')
% grid on
% xlabel('x1')
% ylabel('x2')
% drawnow;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%local functions

function val = I( g )
    eps = 1E-4;
    val = 0.5 * ( ( g^2 + eps^2 )^0.5 - g );
end

function val = dI( g )
    eps = 1E-4;
    val = 0.5 * ( g * ( g^2 + eps^2 )^(-0.5) - 1 );
end




