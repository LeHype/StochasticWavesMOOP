

% Round 2
A_guess = A;
B_guess = B;
Z_guess = Z;
ocp_old = ocp.copy();
ocp.set_initial([ocp.x(:); ocp.lam_g(:)], ocp.value([ocp.x(:); ocp.lam_g(:)]]));
load ("tp_implicit_1_0.mat")
load ("solX_implicit.mat")
sol = sol2;
l = sol.x(6,:);
% load ('tp_explicit_1_0.mat')
% load ("solX_explicit.mat")
% l = sol.x(12,:);
l_star = repmat(tp.ell, 1, (size(sol.x,2)-1)/length(tp.ell));

x_rep = repmat(tp.x, 1, (size(sol.x,2)-1)/length(tp.ell)+1);
u_rep = repmat(tp.u, 1, (length(sol.u)-1)/length(tp.ell)+1);
options = struct;
options.ipopt.max_iter = 50000;
ocp.solver('ipopt', options);

% A*x + B*u 
A = ocp.variable(50,5);
B = ocp.variable(5,5*50);
% C = ocp.variable(1,1);
% D = ocp.variable(1,1);
% E = ocp.variable(1,5);
% F = ocp.variable(5,5);
% % G = ocp.variable(1,1);
% H = ocp.variable(3,1);
Z = ocp.variable(1,1);



% syms X [5,1];
% syms G [125,1];
% 
% 
% 
% cubix_term = 0;
% m = 1;
% Zor i = 1:5
%    for j = 1:5
%        for k = 1:5
% cubix_term = cubix_term +  X(i)*X(j)*X(k)*G(m) ;
% m = m+1;
%        end
%    end 
% end
% cubix_term = matlabFunction(cubix_term);
% G = ocp.variable(125,1);
% G1 = num2cell(G);


% A = rand(1,5);
% B = rand(1);
% C = rand(5,5);
% D = rand(1)

% Storagefunction = @(x,u) A*x  + x'*B*x + C(1)*sin(x(1)*C(2)+C(3))+ D(1)*sin(x(2)*D(2)+D(3)) + E(1)*exp(x(1)*E(2)+E(3))+ F(1)*exp(x(2)*F(2)+F(3));

% Storagefunction = @(x,u) A*x  + x'*B*x + cubix_term(G1{:},x(1),x(2),x(3),x(4),x(5)) +  C(1)*sin(x(1)*C(2)+C(3)) + D(1)*sin(x(2)*D(2)+D(3))+ E(1)*sin(x(3)*E(2)+E(3))+ F(1)*sin(x(4)*F(2)+F(3))+ H(1)*sin(x(5)*H(2)+H(3)) ;
% Storagefunction = @(x,u) A*x  + x'*B*x + C(1)*exp(x(1)*C(2)+C(3))+ D(1)*exp(x(2)*D(2)+D(3))+ E(1)*exp((x(3)+x(4)+x(5))*E(2)+E(3));

Storagefunction = @(x,u,A,B) A*x+   x'*B*x ;

for j = 1:50


    for i = j:50:(length(sol.x)-1)
        %     SF1 = Storagefunction(sol.x(1:5,i+1) , sol.u(i+1) );
        %     SF2 = Storagefunction(sol.x(1:5,i), sol.u(i) );
        %       SR  = l(i+1) -l(i);
        if (j == 50)
            SF1 = Storagefunction(sol.x(1:5,i+1) - x_rep(:,1), sol.u(1) - u_rep(:,i+1),A(1,:),B(:,1:5));

        else
            SF1 = Storagefunction(sol.x(1:5,i+1) - x_rep(:,i+1), sol.u(i+1) - u_rep(:,i+1),A(j+1,:),B(:,j*5+1:(j+1)*5));
        end

        SF2 = Storagefunction(sol.x(1:5,i) - x_rep(:,i), sol.u(i) - u_rep(:,i),A(j,:),B(:,(j-1)*5+1:(j)*5));

        SR  = l(i+1) - l(i)-tp.ell(j);
        ocp.subject_to (SF1-SF2-SR-Z <= 0);



    end
end

% ocp.minimize(sum(abs(A))+sum(abs(B))+sum(sum(abs(C)))+sum(abs(D)));
ocp.minimize(sum(abs(A(:)))+sum(sum(abs(B))));
ocp.subject_to(Z<= 1E-9)
ocp.set_initial(A, A_guess);
ocp.set_initial(B, B_guess);
% ocp.set_initial(C, zeros(1,1));
% ocp.set_initial(D, zeros(1,1));
% ocp.set_initial(E, zeros(1,5));
% ocp.set_initial(F, rand(5,5));
% ocp.set_initial(G, rand(125,1)); 
% ocp.set_initial(H, rand(3,1));
ocp.set_initial(Z, Z_guess);


ocp.solve()


A = ocp.debug.value(A);
B = ocp.debug.value(B);
% C = ocp.debug.value(C);
% D = ocp.debug.value(D);
% E = ocp.debug.value(E);
% F = ocp.debug.value(F);
% G1 = num2cell(ocp.debug.value(G));
% H = ocp.debug.value(H);
Z = ocp.debug.value(Z);

Storagefunction = @(x,u,A,B) A*x+   x'*B*x ;
plooting = [];
for j = 1:50


    for i = j:50:(length(sol.x)-1)
        %     SF1 = Storagefunction(sol.x(1:5,i+1) , sol.u(i+1) );
        %     SF2 = Storagefunction(sol.x(1:5,i), sol.u(i) );
        %       SR  = l(i+1) -l(i);
        if (j == 50)
            SF1 = Storagefunction(sol.x(1:5,i+1) - x_rep(:,1), sol.u(1) - u_rep(:,i+1),A(1,:),B(:,1:5));

        else
            SF1 = Storagefunction(sol.x(1:5,i+1) - x_rep(:,i+1), sol.u(i+1) - u_rep(:,i+1),A(j+1,:),B(:,j*5+1:(j+1)*5));
        end

        SF2 = Storagefunction(sol.x(1:5,i) - x_rep(:,i), sol.u(i) - u_rep(:,i),A(j,:),B(:,(j-1)*5+1:(j)*5));

        SR  = l(i+1) - l(i)-tp.ell(j);
        plooting = [plooting (SF1-SF2 -SR)];



    end
end


    
plot(plooting)
xlabel('time')
ylabel('Verletzung der diss ungleichung')
EGFixFigure
%      ocp.subject_to (SF1-SF2-SR(i)<=0)











