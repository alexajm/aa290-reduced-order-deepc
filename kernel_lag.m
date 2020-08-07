function kernel_lag(w,m,p,lmax)
    %T = length(w) / (m + p)
    tol = 1e-5;
    for lag = 2:lmax
        % Run CVX to get feasible kernel representation
        cvx_begin
        variables R(1,(m+p)*lag)
        %minimize 0
        minimize sum(R(end-(m+p):end))
        subject to
            % Kernel difference equation
%             for j = 0:lmax-lag
%                 0 == R * w((m+p)*j+1:(m+p)*(j+lag));
%             end
            for j = 0:1
                0 == R * w((m+p)*j+1:(m+p)*(j+lag));
            end

            % Full rank condition
            % NOTE: already satisfied for p = 1

            % Establish nonzero element of R to prevent trivial solution
            % NOTE: perhaps this too restrictive, as it failed at lag = 20
            % R(1) == 1;
            
            % Constrain lag to current value
%             for j = 0:lag-1
%                 sum(R((m+p)*j+1:(m+p)*(j+1))) >= tol;
%             end
        cvx_end
        
        if ~strcmp(cvx_status, 'Solved')
            disp(['failed for lag = ', num2str(lag)]);
        end
    end
end