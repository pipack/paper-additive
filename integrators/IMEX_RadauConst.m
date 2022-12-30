classdef IMEX_RadauConst < IntegratorConst & IMEXConst
    
    properties
        graph_line_style = {};
        extrapolate_initial_guess = false;   % if true coefficients A_extrapolate and B_extrapolate will be used to form initial guess for any nonlinear systems
        linearly_implicit                    % if set to true, then only linear solver will be used for implicit component
        linearize                            % if set to true implicit component will be \frac{\partial J}{\partial x}
        kappa                                % number of iterations
    end
    
    properties(SetAccess = protected)
        Bp
        Bi
        A_extrapolate
        nz_rhs_inds
        description = '';
        order = [];
        star                                 % if true, IMEX-Radau*. If false IMEX-Radau
        starting_times = [0];
    end
    
    properties(SetAccess = protected, Dependent)
        name;
    end
    
    methods
        
        function this = IMEX_RadauConst(options)
            if(nargin == 0)
            	options = struct();
            end
            default_field_values = {{'q', 3}, {'kappa', 0}, {'star', true}};
            options = setDefaultOptions(options, default_field_values);
            this = this@IMEXConst(options);
            this = this@IntegratorConst(options);
            this.kappa = options.kappa;
            this.setQ(options.q, options.star);
        end
        
        function name = get.name(this)
            q = size(this.Bp,1);
            if(this.star)
                tag = '$^*$';
            else
                tag = '';
            end
            
            if(this.linearize)
                name = sprintf('LI-Radau%s(%i,%i)', tag, q, this.kappa);
            else
                name = sprintf('SI-Radau%s(%i,%i)', tag, q, this.kappa);
            end
        end
        
        function setQ(this, q, star)
            
            if(nargin == 2)
                star = true;
            end

            % -- coefficients for method ( derive in r,alpha, then convert to h using alpha = 2 ------------------------
            this.Bp = zeros(q);
            this.Bi = zeros(q);            
            z  = this.nodes(q);
            
            % propagator
            if(star)
                ab = [z(q) * ones(q,1), z + 2];
                this.Bp(:,1:q) = double(iPolyCoefficients(z(1:end), ab, 'vpa')) / 2; % convert from r to h by dividing with alpha
                this.nz_rhs_inds = 1:q;
            else
                ab = [z(q) * ones(q,1), z + 2];
                this.Bp(:,2:q) = double(iPolyCoefficients(z(2:end), ab, 'vpa')) / 2; % convert from r to h by dividing with alpha
                this.nz_rhs_inds = 2:q;
            end
            
            % iterator
            ab = [z(1) * ones(q,1), z];
            this.Bi(:,2:q) =  double(iPolyCoefficients(z(2:end), ab, 'vpa')) / 2; % convert from r to h by dividing with alpha
            
            % extrapolator
            this.A_extrapolate = double(polyCoefficients(z,1:q,[],z + 2));
            
            % set-order
            max_radau_order = 2*(q-1) - 1;
            if(star)
                this.order = min(max_radau_order, q + this.kappa);
            else
                this.order = min(max_radau_order, q - 1 + this.kappa);
            end  
            
            this.star = star;
            
        end

    end
    
     methods(Access = protected)
        
        function [step_struct, y_in] = initStepStruct(this, t_n, y_in, problem)
            
            % -- define remainder function -----------------------------------------------------------------------------
            q = size(this.Bi, 1);
            
            step_struct = struct( ... 
                'Fe', repmat(this.evalF(problem, y_in, y_in, 2), 1, q),  ...
                'BiT', transpose(this.Bi), ...
                'BpT', transpose(this.Bp), ...
                'AextrapT', transpose(this.A_extrapolate), ...
                'q', q, ...
                'ode_dim', length(y_in) ...
            );
        
            y_in = repmat(y_in, 1, q);            
            
            % iterate for initial condition
            for i = 1 : q + this.kappa
                [~, y_in, step_struct] = this.iterate(t_n, y_in, step_struct, problem);
            end
            
            % -- if linearize, then re-evaluate linearized RHS at y_q --------------------------------------------------
            if(this.linearize)
                for j = this.nz_rhs_inds
                    step_struct.Fe(:, j) = this.evalF(problem, y_in(:,j), y_in(:,q), 2);
                end                
            end
            
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            step_start_time = tic;
            
            [t_out, y_out, step_struct] = this.propagate(t_in, y_in, step_struct, problem);
            
            % iterate for initial condition
            for i = 1 : this.kappa
                [t_out, y_out, step_struct] = this.iterate(t_out, y_out, step_struct, problem);
            end
            
            if(final_step)
                y_out = y_out(:, 1);
            end
            
            % Update stats
            this.step_stats.recordStep(toc(step_start_time));
            
        end
        
        function  [t_out, y_out, step_struct] = iterate(this, t_n, y_in, step_struct, problem)
            
            h = this.h;
            q = step_struct.q;
            ode_dim = step_struct.ode_dim;
            nz_inds = this.nz_rhs_inds;
            
            % -- evalulate explicit RHS --------------------------------------------------------------------------------
            for j = nz_inds
                step_struct.Fe(:, j) = this.evalF(problem, y_in(:,j), y_in(:,1), 2);
            end
            
            % -- solve fully implicit system ---------------------------------------------------------------------------
            b = bsxfun(@plus, y_in(:,1), step_struct.Fe(:, nz_inds) * h * step_struct.BiT(nz_inds, :));            
            y_out = reshape( ...
                this.solveBCKronF(problem, b(:), h * this.Bi, y_in(:), y_in(:, 1)), ... % linearize about first input
                [ode_dim, q] ...
            );
            t_out = t_n;
            
        end
        
        function  [t_out, y_out, step_struct] = propagate(this, t_n, y_in, step_struct, problem)
            
            % -- read params -------------------------------------------------------------------------------------------
            h = this.h;
            q = step_struct.q;
            ode_dim = step_struct.ode_dim;
            nz_inds = this.nz_rhs_inds;
            
            % -- if necessary, re-evaluate linearized RHS --------------------------------------------------------------
            for j = nz_inds
                step_struct.Fe(:, j) = this.evalF(problem, y_in(:,j), y_in(:,q), 2);
            end                
            
            % -- solve fully implicit system ---------------------------------------------------------------------------
            b = bsxfun(@plus, y_in(:,q), step_struct.Fe(:, nz_inds) * h * step_struct.BpT(nz_inds, :)); 
            if(this.extrapolate_initial_guess)
                y_guess = y_in * step_struct.AextrapT;
            else
                y_guess = repmat(y_in(:,q), q, 1);
            end
            
            y_out = reshape( ...
                this.solveBCKronF(problem, b(:), h * this.Bi, y_guess(:), y_in(:, q)), ... % linearize about first input
                [ode_dim, q] ...
            );
        	t_out = t_n + h;
        
        end
        
        % -- additional functions --------------------------------------------------------------------------------------
        
        function z = nodes(this, q)

            if(q > 8) % use chebfun if q > 8
                if(isempty(which('radaupts')))
                    error('q > 8 requires Chebfun (https://www.chebfun.org/)');
                end
                z = [-1; -flip(radaupts(q-1))];
                return;
            end

            nodes = {
                [1]'
                [-0.333333333333333,  1]'
                [-0.689897948556636,  0.289897948556636, 1]'
                [-0.822824080974592, -0.181066271118531,  0.575318923521694,  1]'
                [-0.885791607770965, -0.446313972723752,  0.167180864737834,  0.720480271312439,  1]'
                [-0.920380285897062, -0.603973164252784, -0.124050379505228,  0.390928546707272,  0.802929828402347,  1]'
                [-0.941367145680430, -0.703842800663031, -0.326030619437691,  0.117343037543100,  0.538467724060109,  0.853891342639482,  1]'
            }; 
            z = [-1; nodes{q-1}];

        end
        
     end
     
end