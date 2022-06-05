function [pareto_params,up,np, varargout] = scalarize_moocp(ocp, costfun, args)
arguments
   ocp             (1,1) casadi.Opti 
   costfun         (1,:) casadi.MX
   args.method     {mustBeMember(args.method, ["ws", "wmm", "ps", "nbi"])}  = "wmm"
   args.normalize  {mustBeMember(args.normalize, ["none", "fix", "adapt"])} = "fix"
   args.mwa_weight (1,1) {mustBeNumeric} = 1e-4
   args.simplify   (1,1) logical = true
end

k = size( costfun, 2 );

if ismember(args.normalize, ["fix", "adapt"])
    if isequal(args.normalize, "fix")
        for i = 1:k
            w_vec    = args.mwa_weight*ones(k,1);
            w_vec(i) = 1;
            ocp.minimize( costfun*w_vec )

            sol(i) = ocp.solve();
            ep(i,:) = sol(i).value(costfun);
        end
    
        up = min( ep, [], 1 );
        np = max( ep, [], 1 );
        costfun = (costfun - up)./(np - up);
        ep_norm = (ep - up)./(np - up);
        varargout{1} = ep_norm;
        varargout{2} = costfun;
        varargout{3} = ep;
        varargout{4} = sol;

    elseif isequal(args.normalize, "adapt")
        ep_problem = ocp.copy();
        w_ep       = ep_problem.parameter(k,1);
        ep_problem.minimize(costfun*w_ep);
        varargout{4} = ep_problem;
        varargout{5} = w_ep;

        up = ocp.parameter(1,k);
        np = ocp.parameter(1,k);
        costfun = (costfun - up)./(np - up);
        varargout{1} = up;
        varargout{2} = np;
        varargout{3} = costfun;
    end
   
end

switch args.method
    case "ws"
        w = ocp.parameter(k,1);
        ocp.minimize( costfun*w )
        
        pareto_params = w;

    case "wmm"
        w = ocp.parameter(1,k);
        max_cost = ocp.variable(1,1);
        ocp.subject_to( max_cost >=  w.*costfun )
        ocp.minimize( max_cost )
        
        pareto_params = w;

    case {"ps", "nbi"}
        l   = ocp.variable( 1 );
        pt  = ocp.parameter( 1, k );
        dir = ocp.parameter( 1, k );
        ocp.subject_to( pt + l*dir >= costfun )
        
        if args.simplify
            ocp.minimize( sum(pt + l*dir) )
        else
            ocp.minimize( -l )
        end
        
        pareto_params = {pt, dir};

end

end