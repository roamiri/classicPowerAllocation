
function beta=gss(FBS,i, pmax, pmin, eta, I_th, num) % maximizing, so v is +gradient of f
    plotting = 0;
    if plotting
        hold on;
    end
    [lambda ,gamma ,mu] = FBS{i}.getLagrangeVars();
    switch num
        case 1
            a = min(lambda, lambda + 100 * (-pmax + FBS{i}.P));
            b = max(lambda, lambda + 100 * (-pmax + FBS{i}.P));
        case 2
            a = min(gamma, gamma + 100 * (-FBS{i}.P + pmin));
            b = max(gamma, gamma + 100 * (-FBS{i}.P + pmin));
        case 3
            a = min(mu, mu + 100 * (-FBS{i}.C_FUE + FBS{i}.R_FUE));
            b = max(mu, mu + 100 * (-FBS{i}.C_FUE + FBS{i}.R_FUE));
        case 4
            dum = 0.0;
            for kk=1:size(FBS,2)
                dum = dum + FBS{kk}.P * FBS{kk}.gmf;
            end
            a = min(eta, eta + 100 * (-I_th + dum));
            b = max(eta, eta + 100 * (-I_th + dum));            
    end
% a=0;                            % start of interval
% b=2000;                            % end of interval
epsilon=0.000001;               % accuracy value
iter= 50;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations


x1=a+(1-tau)*(b-a);             % computing x values
x2=a+tau*(b-a);

f_x1=f(x1, FBS, i, pmax, pmin, eta, I_th, num);                     % computing values in x points
f_x2=f(x2, FBS, i, pmax, pmin, eta, I_th, num);

if plotting
plot(x1,f_x1,'rx');           % plotting x
plot(x2,f_x2,'rx');
end
while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1>f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
        
        f_x1=f(x1,FBS,i, pmax, pmin, eta, I_th, num);
        f_x2=f(x2,FBS,i, pmax, pmin, eta, I_th, num);
        if plotting
        plot(x1,f_x1,'rx');
        end
    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
        
        f_x1=f(x1,FBS,i, pmax, pmin, eta, I_th, num);
        f_x2=f(x2,FBS,i, pmax, pmin, eta, I_th, num);
        if plotting
        plot(x2,f_x2,'rx');
        end
    end
    
    k=k+1;
end


% chooses minimum point
if(f_x1>f_x2)
    beta = x1;
%     sprintf('x_min=%f', x1)
%     sprintf('f(x_min)=%f ', f_x1)
if plotting
    plot(x1,f_x1,'ro');
end
else
    beta = x2;
%     sprintf('x_min=%f', x2)
%     sprintf('f(x_min)=%f ', f_x2)
if plotting
    plot(x2,f_x2,'ro')
end
end


end
