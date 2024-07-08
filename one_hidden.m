function [ num_params ] = one_hidden(iflag)

global N1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%construct network

x = sym( 'x' , [ 2 , 1 ] ); %input variables

num_layer1 = 2*N1 + N1;
num_output = 2*N1;

num_params = num_layer1 + num_output; %total parameter count

fprintf('number of parameters: %d\n',num_params)

if iflag
    
    p = sym( 'p' , [ num_params , 1 ] ); assume(p,'real'); %symbolic parameter vector
    
    %first layer
    w1 = reshape( p( 1:(2*N1) ) , [N1,2] );
    b1 = p( (2*N1+1):3*N1 );

    %output layer
    w2 = reshape( p( (3*N1+1):(3*N1+2*N1) ) , [2,N1] );
    
    %define forward pass
    first = act( w1*x + b1 );
    basis = tanh( x(2) ) * first;

    output = w2 * basis;
    
    %build solution with linear parameter 
    N = output;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('begin symbolic computations...')

    %compute strain vector
    
    u_xx = diff( N(1) , x(1) , 1 );
    u_xy = diff( N(1) , x(2) , 1 );
    u_yx = diff( N(2) , x(1) , 1 );
    u_yy = diff( N(2) , x(2) , 1 );
    
    grad = sym( zeros( 2 , 2 ) );

    %displacement gradient
    grad(1,1) = u_xx; grad(1,2) = u_xy; 
    grad(2,1) = u_yx; grad(2,2) = u_yy;

    %deformation gradient
    F = eye(2,2) + grad;

    %determinant of deformation gradient
    J = F(1,1)*F(2,2) - F(1,2)*F(2,1);

    %right cauchy-green tensor
    C = transpose( F ) * ( F );

    %first invariant
    I1 = C(1,1) + C(2,2);

    %parameter derivatives of J and I1
    J_p = sym( zeros(num_params,1) );
    I1_p = sym( zeros(num_params,1) );

    for i=1:num_params
        J_p(i) = diff( J , p(i) , 1 );
        I1_p(i) = diff( I1 , p(i) , 1 );
    end
    
    %parameter gradient of displacement vector
    u_p = sym( zeros( 2 , num_params ) );

    for i=1:2
        for j=1:num_params
            u_p(i,j) = diff( N(i) , p(j) , 1 );
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %convert to numerical functions
    
    disp('writing functions to files...')

    %function approximation
    matlabFunction( N ,'vars' , { x , p } , 'File' , 'N' , 'Optimize' , true );

    %determinant of deformation gradient
    matlabFunction( J ,'vars' , { x , p } , 'File' , 'J' , 'Optimize' , true );

    %first invariant of cauchy green
    matlabFunction( I1  ,'vars' , { x , p } , 'File' , 'I1' , 'Optimize' , true );

    %parameter gradient of deformation gradient determinant
    matlabFunction( J_p ,'vars' , { x , p } , 'File' , 'J_p' , 'Optimize' , true );

    %parameter gradient of first invariant
    matlabFunction( I1_p ,'vars' , { x , p } , 'File' , 'I1_p' , 'Optimize' , true );

    %parameter gradient of displacement
    matlabFunction( u_p ,'vars' , { x , p } , 'File' , 'u_p' , 'Optimize' , true );

    disp('done')

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%local functions


%nonlinear activation function
function vec = act( layer )
%vec = tanh(  layer );
%vec = log( 1 + exp( layer ) );
vec = 1 ./ ( 1 + exp( -layer ) );
% eps = 0.05;
% vec = ( ( layer.^2 + eps^2 ).^0.5 - eps + 1.1*layer)/2;
end






