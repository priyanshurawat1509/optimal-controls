function [G,F,S,Cc,W] = mpcsolver(A,B,C,Q,P,R,n,m,p,N,umax,umin,ymax,ymin)

%% mpcsolver.m - Coverts the matrices of the MPC control problem to a quadratic program which `quadprog` solver

% INPUTS:
    % A     :     System matrix
    % B     :     Input matrix
    % C     :     Output matrix 
    % Q     :     States weights matrix 
    % P     :     % States weights matrix
    % R     :     Input penalty matrix 
    % n     :     Number of states
    % m     :     Number of inputs
    % p     :     Number of outputs
    % N     :     Length of prediction horizon
    % umax  :     Vector defining the upper input limits
    % umin  :     Vector defining the lower input limits
    % ymax  :     Vector defining the upper output limits
    % ymin  :     Vector defining the lower input limits

% OUTPUTS:
    % G, F, S, Cc, W     
   

%% Matrix formulation (cost function)

% Φ matrix formulation (let's call it `phi`)
for i=1:N
    phi((i-1)*n+1:i*n,:) = A^i;
end

% Γ matrix formulation (let's call it `gamma`)
for i=1:N
    for j = 1:N
        if i>j
            gamma((i-1)*n+1:i*n,(j-1)*m+1:j*m)=A^(i-j)*B;
        elseif i==j
            gamma((i-1)*n+1:i*n,(j-1)*m+1:j*m)=B;
        else
            gamma((i-1)*n+1:i*n,(j-1)*m+1:j*m)=zeros(n,m);
        end
    end
end 

% Ω matrix formulation (let's call it `omega`)   
omega=Q;
for i=1:N-1
  omega=blkdiag(omega,Q);
end
if sum(P(:)) ~= 0
    omega=blkdiag(omega,P);
end

% Ψ matrix formualtion (let's call it `psi`)
psi=R;
for i=1:N-1
  psi=blkdiag(psi,R);
end

% G & F matrices
G=2*(psi+(gamma'*omega*gamma));
F=2*(gamma'*omega*phi);


%% Matrix formulation (constraints)

%% If constraints are applicable only on the input
if isempty(umin)==0 && isempty(ymin)==1
    M =[];
    for i = 1:N-1
        Mi = [zeros(m,n); zeros(m,n); zeros(p,n); zeros(p,n)];
        M=blkdiag(M,Mi);
    end
    Mn = [zeros(p,n); zeros(p,n)]; %for the last stage
    M = blkdiag(M,Mn);
    M = [zeros(2*(p+m),n*N) ; M]; %initial zeros
    
    % D matrix
    D = zeros(size(M,1),n);
    D(1:2*(m+p),1:n) = [zeros(m,n); zeros(m,n); zeros(p,n); zeros(p,n)];
    
    % Σ matrix (let's call it `E`)
    E = [];
    for i = 1:N
        Ei = [-eye(m,m); eye(m,m); zeros(p,m); zeros(p,m)]; %since no constrains on 'u'
        E=blkdiag(E,Ei);
    end
    E = [E; zeros(size(M,1)-size(E,1),m*N)]; %appending last row zeros
    
    % Ç matrix (let's call it `Cc`)
    Cc=[];
    for i = 1:N
        bi = [-umin; umax; zeros(p,1); zeros(p,1)];
        Cc = [Cc;bi];
    end
    
    bn = [zeros(p,1); zeros(p,1)]; %for the last stage
    Cc = [Cc;bn];
    
    % S & W matrices
    S = M*gamma+E;
    W = -(M*phi+D);
end

%% If constraints are applicable only on the output
if isempty(umin)==1 && isempty(ymin)==0
    % M matrix
    M =[];
    for i = 1:N-1
        Mi = [zeros(m,n); zeros(m,n); -C; +C];
        M=blkdiag(M,Mi);
    end
    Mn = [-C; +C]; %for the last stage
    M = blkdiag(M,Mn);
    M = [zeros(2*(p+m),n*N) ; M]; %initial zeros
    
    % D matrix
    D = zeros(size(M,1),n);
    D(1:2*(m+p),1:n) = [zeros(m,n); zeros(m,n); -C; +C];
    
    % Σ matrix (let's call it `E`)
    E = [];
    for i = 1:N
        Ei = [zeros(m,m); zeros(m,m); zeros(p,m); zeros(p,m)]; % since no constrains on 'u'
        E=blkdiag(E,Ei);
    end
    E = [E ; zeros(size(M,1)-size(E,1),m*N)]; %appending last row zeros
    
    
    % Ç matrix (let's call it `Cc`)
    Cc=[];
    for i = 1:N
        bi = [zeros(m,1); zeros(m,1); -ymin; ymax];
        Cc = [Cc;bi];
    end
    
    bn = [-ymin; +ymax]; %for the last stage
    Cc = [Cc;bn];
    
    % S & W matrices
    S = M*gamma+E;
    W = -(M*phi+D);
end

