function rocket()

    vele = 3000; % Exhaust velocity in m/s
    masspay = 1000; % Mass of payload in kg
    leov = 7800; % LEO velocity in m/s
    maxstages = 10;
    
    % Initial D.V. values
    x0 = zeros(1, maxstages);
    
    % Setup for linear equality constraints
    Aeq = zeros(1, maxstages);
    beq = [0];
    
    % Setup for upper bounds
    ub = ones(1, maxstages).*inf;
        
    % Setup for lower bounds
    lb = zeros(maxstages);
    
    % Setup for linear inequality constraints
    A = -(ones(1, maxstages)); 
    b = [0];

    % Nonlinear constraint function
    function [c ceq] = nlincon(x)
        totalv = 0;
        mstr = x.*.25;
        for i = 1:maxstages
            mstage = mstr(i)+x(i);
            mrest = sum(x(i+1:end)) + sum(mstr(i+1:end)) + masspay;
            deltav = vele*log((mstage + mrest)/(mstr(i)+mrest));
            totalv = totalv + deltav;
        end

        veq = totalv - leov;

        ceq = [veq];
        c = [0];

    end

    % Objective function
    function mtot = objfun(x)
        mtot = sum(x);
    end


    function veq = checkRes(masses)
        totalv = 0;
        mstr = x.*.25;
        for i = 1:maxstages
            mstage = mstr(i)+x(i);
            mrest = sum(x(i+1:end)) + sum(mstr(i+1:end)) + masspay;
            deltav = vele*log((mstage + mrest)/(mstr(i)+mrest));
            totalv = totalv + deltav;
        end
        veq = totalv - leov;
    end

    x = fmincon(@objfun, x0, A, b, Aeq, beq, lb, ub, @nlincon);
    
    checkRes(x)

end