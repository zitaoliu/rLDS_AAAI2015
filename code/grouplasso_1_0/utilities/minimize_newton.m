function [x,exit_parameters,f0,fx] = minimize_newton(x,fun,optparam,varargin)
x000 = x;

tol = optparam.tol;
kmax = optparam.kmax;
display = optparam.display;
n = size(x,1);
dx=2*tol;
k=1;
go_on = 1;
if isfield(optparam,'kmax_ls'), kmax_ls = optparam.kmax_ls;
else, kmax_ls = 40 ; end
if isfield(optparam,'rho_GOLDSTEIN'), rho_GOLDSTEIN = optparam.rho_GOLDSTEIN;
else, rho_GOLDSTEIN = .25 ; end
if isfield(optparam,'t_GOLDSTEIN'), kmax_ls = optparam.t_GOLDSTEIN;
else, t_GOLDSTEIN = 2 ; end
if isfield(optparam,'tol_df'), tol_df = optparam.tol_df;
else, tol_df = sqrt(tol)/100; end
df = 2*tol_df;
if isfield(optparam,'m1_WOLFE'), m1_WOLFE = optparam.m1_WOLFE;
else, m1_WOLFE = .2 ; end
if isfield(optparam,'m1_WOLFE'), m2_WOLFE = optparam.m2_WOLFE;
else, m2_WOLFE = .9 ; end

fxold = Inf;
xold = x;
nevals = 0;
lambda2old = 0;
alpha0 = 1;
nevals = 0;
kls = 0;
descent=zeros(n,1);
while go_on & k<kmax & df > tol_df;
    % compute gradient
    [fx,gradient,hessian,lambda2, descent] = feval(fun,x,varargin{:});
    if isinf(fx), exit_parameters.type='infinite_f'; fprintf('infinite_f -> exit\n'); return; end
    if k==1, f0=fx; end
    if k==1, f000=fx; end

    nevals = nevals + 1;
    if 0,
        'derivatives'
        %check gradients
        ddx = randn(size(x));
        dr =1e-8;
        dfx  = feval(fun,x+dr*ddx,varargin{:});
        abs((dfx-fx)-dr * gradient'*ddx)/abs(dfx-fx)

        %          ddx = randn(size(x));
        dr =1e-10;
        dfx = feval(fun,x+dr*ddx,varargin{:});
        abs((dfx-fx)-dr * gradient'*ddx)/abs(dfx-fx)

        'hessian'

        % check hessian
        ddx = randn(size(x));
        dr =1e-8;
        [dfx,dgradient] = feval(fun,x+dr*ddx,varargin{:});
        D1= dgradient - gradient ;
        D2= dr * hessian * ddx;
        norm( D1-D2 ) / norm(D1 )
   
        ddx = randn(size(x));
        dr =1e-10;
        [dfx,dgradient] = feval(fun,x+dr*ddx,varargin{:});
        D1= dgradient - gradient ;
        D2= dr * hessian * ddx;
        norm( D1-D2 ) / norm(D1 )
 
        return
    end


    df = fxold-fx;
    if norm(gradient)<1e-14, exit_parameters.type='zero grad'; return;  end
    %if norm(descent)<1e-14,  exit_parameters.type='zero descent'; return;  end



    beta = descent'*gradient / norm(descent) / norm(gradient);
    if beta > -.1e-10,
        descent = - gradient;
        beta = descent'*gradient / norm(descent) / norm(gradient); end

    lambda2 = - gradient' * descent ;
    if isnan(lambda2), H = eye(n); descent = - gradient; lambda2 = - gradient' * descent ; end


    lambda2crit = lambda2;
    lambda2old = lambda2;
    if ~isinf(display)
        if display==1,
            fprintf('k=%d - lambda2=%e - lambda2=%e - beta = %f - f=%f - df=%e - nev = %d, dx=%e\n',k,gradient' * gradient,lambda2,...
                beta ,fx,df,kls,norm(xold-x));
        end
        if display>1,
            if mod(k,display)==1
                fprintf('k=%d - lambda2=%e - lambda2=%e - beta = %f - f=%f - df=%e - nev = %d, dx=%e\n',k,gradient' * gradient,lambda2,...
                    beta ,fx,df,kls,norm(xold-x));
            end
        end
    end
    xold=x;
    if lambda2crit < 2 * tol
        % usual end of iteration
        go_on = 0;
        exit_parameters.type = 'small_lambda2';
        exit_parameters.lambda2 = lambda2;
        exit_parameters.k = k;
        exit_parameters.nevals = nevals;

    else

        % Wolfe line search
        a = 0;
        b = Inf;

        kls = 0;
        alpha = alpha0;
        alphas = [];
        falphas = [];
        while kls < kmax_ls

            [falpha,gradalpha] = feval(fun,x + alpha * descent,varargin{:});
            alphas = [ alphas alpha ];
            falphas= [ falphas falpha ];

            nevals = nevals + 1;



            if (falpha > fx - m1_WOLFE * lambda2 * alpha)
                % alpha is too large, reduce it
                b = alpha;
                alpha = (a+b)/2;
                %   if display, fprintf('-'); end
            else
                if gradalpha'*descent < m2_WOLFE * gradient' * descent
                    % alpha is too small
                    a = alpha;
                    alpha = 2 * alpha;
                    %  if display,  fprintf('+'); end
                else
                    % alpha is correct
                    break;
                end
            end
            kls = kls + 1;
        end


        %         % Goldstein line search
        %         a = 0;
        %         b = Inf;
        %         kls = 0;
        %         alpha = alpha0;
        %         alphas = [];
        %         falphas = [];
        %         already_one = 0;
        %            while kls < kmax_ls
        %
        %             falpha = feval(fun,x + alpha * descent,varargin{:});
        %             alphas = [ alphas alpha ];
        %             falphas= [ falphas falpha ];
        %
        %             nevals = nevals + 1;
        %             if (falpha <= fx - rho_GOLDSTEIN * lambda2 * alpha)
        %                 fprintf('+');
        %                 already_one = 1;
        %                 alpha = alpha * t_GOLDSTEIN * 2;
        %             else
        %                 fprintf('-');
        %                 if kls>0,
        %                     if falphas(end)>falphas(end-1), break; end
        %                 end
        %                 if already_one, break;
        %                 else alpha = alpha / t_GOLDSTEIN; end
        %             end
        %             kls = kls + 1;
        %         end
        %         if kls==kmax_ls, keyboard; end
        [a,b]=min(falphas);
        alpha = alphas(b);
        alpha0 = alpha;
        gradold = gradient;
        fxold = fx;
        xold = x;
        x = x + alpha * descent;

        % save x x
        k=k+1;
    end
end
if k==kmax,  exit_parameters.type = 'max_iterations';  end
if df<=tol_df, exit_parameters.type = 'df_small'; exit_parameters.df =df; end
exit_parameters.nevals = nevals;
exit_parameters.k=k;
if isinf(display)
    fprintf('k=%d - lambda2=%e - f=%f - df=%e - dx=%e\n',k,lambda2crit,...
        fx,f000-fx,norm(x000-x));
end