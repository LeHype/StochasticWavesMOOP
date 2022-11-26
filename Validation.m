
load ("tp_implicit_1_0.mat")
% load ("solX_implicit.mat")
sol = sol2;
l = sol.x(6,:);
% load ('tp_explicit_1_0.mat')
% load ("solX_explicit.mat")
% l = sol.x(12,:);
l_star = repmat(tp.ell, 1, (size(sol.x,2)-1)/length(tp.ell));

x_rep = repmat(tp.x, 1, (size(sol.x,2)-1)/length(tp.ell)+1);
u_rep = repmat(tp.u, 1, (length(sol.u)-1)/length(tp.ell)+1);


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