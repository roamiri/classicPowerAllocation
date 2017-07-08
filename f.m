
function y=f(x,FBS,i, pmax, pmin, eta, I_th, num)
    [lambda, gamma, mu] = FBS{i}.getLagrangeVars();
    K = size(FBS,2);
    switch num
        case 1
            X = lambda - x * (pmax - FBS{i}.P);
        case 2
            X = gamma - x * (FBS{i}.P - pmin);
        case 3
            X = mu - x * (FBS{i}.C_FUE - FBS{i}.R_FUE);
        case 4
            dum = 0.0;
            for kk=1:K
                dum = dum + FBS{kk}.P * FBS{kk}.gmf;
            end
            X = eta - x * (I_th - dum);
    end
    
    feval = 0.0;
    for j=1:K
        fbs = FBS{j};
        [ll, gg, mm] = fbs.getLagrangeVars();
        if i==j
            switch num
                case 1
                    ll = x;
                case 2
                    gg = x;
                case 3
                    mm = x;
                case 4
                    eta = x;
            end
        end        
        feval = feval + (1+mm)*fbs.C_FUE + (gg-ll)*fbs.P - eta*fbs.P*fbs.gmf + ll*pmax - gg*pmin - mm*fbs.R_FUE + eta*I_th;
    end
    y = feval;
end