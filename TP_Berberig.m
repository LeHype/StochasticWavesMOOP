global PathToParameters
PathToParameters= 'src/PolySurge_inputs.mat';
load(PathToParameters);             
%% HARD CODED !!
tp_start = 100/timestep+1;
tp_end   = tp_start +10/timestep-1;

% tp = struct;
% tp.x = sol_old.x(1:5,tp_start:1:tp_end);
% tp.u = sol_old.u(tp_start:1:tp_end);
% tp.time = sol_old.time(tp_start:1:tp_end);
% tp.ell = sol_old.x(6,tp_start:1:tp_end)-sol_old.x(6,tp_start);
% tp.ell_pi = mean(tp.ell);

tp = struct;
tp.x = x_ana(1:5,tp_start:1:tp_end);
tp.u = u_ana(tp_start:1:tp_end)';
tp.time = t_ana(tp_start:1:tp_end);
tp.ell = x_ana(6,tp_start:1:tp_end)-sol_old.x(6,tp_start);
tp.ell_pi = mean(tp.ell);

AbstandzuTP = AbstandzuTP_vec(x_ana(1:5,1:end-1), u_ana(1:end-1)', tp);


 
function Abstand = AbstandzuTP_vec(x,u,tp)

tp_x = reshape(tp.x, size(tp.x,1), 1, size(tp.x,2));
tp_u = reshape(tp.u, size(tp.u,1), 1, size(tp.u,2));
Abstand_mat = vecnorm(x-tp_x,2,1)+vecnorm(u-tp_u,2,1);
Abstand = min(Abstand_mat, [], 3);

end




