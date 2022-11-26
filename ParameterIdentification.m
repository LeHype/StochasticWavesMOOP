ocp = casadi.Opti();
load ("tp_implicit_1_0.mat")
load ("solX_implicit.mat")
 
l = sol.x(6,:);
% load ('tp_explicit_1_0.mat')
% load ("solX_explicit.mat")
% l = sol.x(12,:);
l_star = repmat(tp.ell, 1, (size(sol.x,2)-1)/length(tp.ell));

x_rep = repmat(tp.x, 1, (size(sol.x,2)-1)/length(tp.ell)+1);
u_rep = repmat(tp.u, 1, (length(sol.u)-1)/length(tp.ell)+1);
x_rep = x_rep;
options = struct;
options.ipopt.max_iter = 50000;
ocp.solver('ipopt', options);

% A*x + B*u 
A = ocp.variable(50,5);
B = ocp.variable(5,5*50);
C = ocp.variable(50,1);
D = ocp.variable(50,5);
E = ocp.variable(5,5*50);
F = ocp.variable(5,5*50);

Z = ocp.variable(1,1);

Storagefunction = @(x,u,xe,A,B,C,D,E,F) A*x+   x'*B*x +C + D*xe+   xe'*E*xe +  x'*F*x*u;

for j = 1:50


    for i = j:50:(length(sol.x)-1)
        %     SF1 = Storagefunction(sol.x(1:5,i+1) , sol.u(i+1) );
        %     SF2 = Storagefunction(sol.x(1:5,i), sol.u(i) );
        %       SR  = l(i+1) -l(i);
        if (j == 50)
            SF1 = Storagefunction(sol.x(1:5,i+1), sol.u(1) - u_rep(:,i+1),sol.x(1:5,i+1) - x_rep(:,1),A(1,:),B(:,1:5),C(1,1),D(1,:),E(:,1:5),F(:,1:5));

        else
            SF1 = Storagefunction(sol.x(1:5,i+1), sol.u(i+1) - u_rep(:,i+1),sol.x(1:5,i+1) - x_rep(:,i+1),A(j+1,:),B(:,j*5+1:(j+1)*5),C(j+1,1),D(j+1,:),E(:,j*5+1:(j+1)*5),F(:,j*5+1:(j+1)*5));
        end

        SF2 = Storagefunction(sol.x(1:5,i), sol.u(i) - u_rep(:,i),sol.x(1:5,i) - x_rep(:,i),A(j,:),B(:,(j-1)*5+1:(j)*5),C(j,1),D(j,:),E(:,(j-1)*5+1:(j)*5),F(:,(j-1)*5+1:(j)*5));

        SR  = l(i+1) - l(i)-tp.ell(j);
        ocp.subject_to (SF1-SF2-SR-Z <= 0);



    end
end
%
ocp.minimize(sum(abs(A(:)))+sum(abs(B(:)))+sum(sum(abs(C(:))))+sum(abs(D(:)))+sum(abs(Z(:))));
% ocp.minimize(Z);
ocp.subject_to(Z<=0);
ocp.set_initial(A, randn(50,5));
ocp.set_initial(B, randn(5,5*50));
ocp.set_initial(C, randn(50,1));
ocp.set_initial(D, randn(50,5));
ocp.set_initial(E, randn(5,5*50));
ocp.set_initial(F, randn(5,5*50));
% ocp.set_initial(G, rand(125,1)); 
% ocp.set_initial(H, rand(3,1));
ocp.set_initial(Z, zeros(1,1));


ocp.solve()


A = ocp.debug.value(A);
B = ocp.debug.value(B);
C = ocp.debug.value(C);
D = ocp.debug.value(D);
E = ocp.debug.value(E);
F = ocp.debug.value(F);
% G1 = num2cell(ocp.debug.value(G));
% H = ocp.debug.value(H);
Z = ocp.debug.value(Z);
plooting = [];
Storagefunction = @(x,u,xe,A,B,C,D,E,F) A*x+   x'*B*x +C + D*xe+   xe'*E*xe +  x'*F*x*u;

for j = 1:50


    for i = j:50:(length(sol.x)-1)
        %     SF1 = Storagefunction(sol.x(1:5,i+1) , sol.u(i+1) );
        %     SF2 = Storagefunction(sol.x(1:5,i), sol.u(i) );
        %       SR  = l(i+1) -l(i);
        if (j == 50)
            SF1 = Storagefunction(sol.x(1:5,i+1), sol.u(1) - u_rep(:,i+1),sol.x(1:5,i+1) - x_rep(:,1),A(1,:),B(:,1:5),C(1,1),D(1,:),E(:,1:5),F(:,1:5));

        else
            SF1 = Storagefunction(sol.x(1:5,i+1), sol.u(i+1) - u_rep(:,i+1),sol.x(1:5,i+1) - x_rep(:,i+1),A(j+1,:),B(:,j*5+1:(j+1)*5),C(j+1,1),D(j+1,:),E(:,j*5+1:(j+1)*5),F(:,j*5+1:(j+1)*5));
        end

        SF2 = Storagefunction(sol.x(1:5,i), sol.u(i) - u_rep(:,i),sol.x(1:5,i) - x_rep(:,i),A(j,:),B(:,(j-1)*5+1:(j)*5),C(j,1),D(j,:),E(:,(j-1)*5+1:(j)*5),F(:,(j-1)*5+1:(j)*5));

        SR  = l(i+1) - l(i)-tp.ell(j);
        plooting = [plooting (SF1-SF2 -SR)];



    end
end


    
plot(plooting)
xlabel('time')
ylabel('Verletzung der diss ungleichung')
EGFixFigure
%      ocp.subject_to (SF1-SF2-SR(i)<=0)

