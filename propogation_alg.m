function M = propogation_alg(M, varargin)
% Propagation based on closest points and length of connections

%% Input parsing

f = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
    function out = read_pair(name, default)
        ind = find(f(name));
        if ind, out = varargin{ind + 1};
        else out = default;
        end
    end

if isfield(M, 'P') && ~any(f('rewrite')), return, end

init_insult = read_pair('init_insult', [0 0]);
meduim_speed = read_pair('meduim_speed', .02);
bundle_speed = read_pair('bundle_speed', .02);
neroun_speed = read_pair('neroun_speed', 1);
all_edge_delay = read_pair('edge_delay', 2);

neuron_speed_formula = @(R) 2./R; % 2/R, For  2/R^2  write this code   ->>  2./R.^2

%% init insult GUI

    function btn_Down(~,~)
        cp = get(gca,'currentpoint');
        init_insult = cp(1,1:2);
        set(ip, 'XData', init_insult(1) ,'YData', init_insult(2));
    end
ok = 0;
    function OK(~,~)
        delete(fig_ui);
        ok = 1;
        drawnow;
    end

if any(f('GUI'))
    fig_ui = figure('units', 'normalized', 'toolbar','figure', 'Name', 'Bundle Settings', 'outerposition', [0.2 0.15 .6 .8]);
    ax_ui = axes('units', 'normalized', 'position', [0.1 0.1 .7 .8]);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .4 .12 .05], 'String', 'Done/Continue', 'Callback', @OK);
    
    M.plot.model();
    hc = allchild(ax_ui);
    for k = 1:length(hc)
        set(hc(k), 'HitTest', 'off');
    end
    ip = plot(init_insult(1), init_insult(2), 'r.', 'markerSize', 20, 'HitTest', 'off');
    set(gca, 'ButtonDownFcn', @btn_Down)
    waitfor(fig_ui);
    if ~ok, return; end
end

%% Find each domain speed

bt = M.csg.bt;
nb_domains = length(bt);
% Find mapping between domains and geometry elements

dom_map = [find(sum(bt, 2) == 1);... % Nerve
    find(sum(bt, 2) == 2);... % Bundles
    find(sum(bt, 2) == 3)]; % Neurons

subdomain.meduim = dom_map(1);
subdomain.bundles = dom_map(1 + (1:length(M.neuron)));
all_neurons_g = [M.neuron{:}];
subdomain.neuron = dom_map(1 + length(M.neuron) + (1:length(all_neurons_g)));

% assigning domain speed
domain_speed = zeros(nb_domains,1);
domain_speed(subdomain.meduim) = meduim_speed;
domain_speed(subdomain.bundles) = bundle_speed;
domain_speed(subdomain.neuron) = neroun_speed * neuron_speed_formula(all_neurons_g(3,:));

p = M.mesh.p; e = M.mesh.e; t = M.mesh.t;
% Setting propagation speed to each point based on their domain speed
propagation_speed = ones(1, length(p)); % speed of each point
edge_delay = zeros(1, length(p)); % additional imagenary distance infection has to go in order to pass the boundaries

for k = 1:nb_domains
    [domain_inter, domain_edge] = pdesdp(p,e,t,k); % find each domain edge and interior points
    propagation_speed(domain_edge) = meduim_speed;
    propagation_speed(domain_inter) = domain_speed(k);
end
Edge = e(1,:); % set of points on any boundary (slower propagation)
edge_delay(Edge) = all_edge_delay;

%% Init

Conn = t(1:3,:); % Connections / Connectivity matrix
active_cols_mask = true(length(Conn),1); % columns of Conn matrix that have a dead point
newly_dead_col_mask = false(length(Conn),1);

newly_dead = [];
[~,id] = min(((p(1,:) - init_insult(1)).^2 + (p(2,:) - init_insult(2)).^2)); % id: Closest mesh point to the starting point
newly_infected = [id; 0];
infected = []; % first row is the infected point index, second row is the death timer
death_time = zeros(length(p),1);
Now = 0; % current time

%% Propogation algorithm
% ************************************************************
% Point states are: live, infected, dead. Those infected have a death timer
% which will kill them after a certain time. While in infected state A
% single point can get infected from several sorrounding points, each time,
% death timer will be updated to the minimum death time assigned.
% ************************************************************

fprintf('Simulation... ');

del = @(x) fprintf(repmat('\b',1,x));
nDispChar = 0; last_disp_num = 0;

while any(active_cols_mask) || ~isempty(infected)
    
    new_disp_num = round(100*(1 - sum(active_cols_mask)/length(Conn)));
    if new_disp_num ~= last_disp_num
        del(nDispChar); nDispChar = fprintf('%.0f%%\t', new_disp_num);
        last_disp_num = new_disp_num;
    end
    if ~isempty(newly_dead)
        % ************************************************************
        % finds all connected points to the newly dead point and based on
        % the length of connection assigns a death timer to those points
        % results go to newly_infected (first row is the infected point
        % index, second row is the death timer)
        % ************************************************************
        newly_dead_col_idx = ceil(find(Conn == newly_dead)/3);
        
        newly_dead_col_mask(:) = 0;
        newly_dead_col_mask(newly_dead_col_idx) = 1;
        x = Conn(:,active_cols_mask & newly_dead_col_mask); % all surrounding connected points (repetitive and including current point idx)
        y = unique(x(:)); % unique points (including current point idx)
        cpi = y(y ~= newly_dead)'; % connected point index
        
        dp = sum((p(:,newly_dead)*ones(1,length(cpi)) - p(:,cpi)).^2); % distance of connected points to current dead point
        death_timer = dp./propagation_speed(cpi) + edge_delay(cpi); % death timer is set based on distance, current meduim speed and if that point is on a boundary
        newly_infected = [cpi; death_timer]; % first row is the infected point index, second row is the death timer
        active_cols_mask(newly_dead_col_mask) = 0; % update active connection set
        
    end
    if ~isempty(infected)
        % ************************************************************
        % Finds points already infected and checks if the new infection
        % death timer is less than the death time already assigned, if so
        % updates the death timer.
        % ************************************************************
        rm = []; % cols to be removed from newly_infected (those already infected)
        for k = 1:size(newly_infected,2)
            idx = find(infected(1,:) == newly_infected(1,k)); % find intersection of two sets
            if ~isempty(idx)
                infected(2,idx) = min(infected(2,idx), newly_infected(2,k)); % choose minimum death time
                rm = [rm k];
            end
        end
        newly_infected(:,rm) = []; % remove already-infected points from the set
    end
    
    infected = [infected newly_infected]; % add non-repetitive newly-infected points to the infected set
    
    [Dt, idx] = min(infected(2,:)); % find the next infected point waiting for death and its timer value
    Now = Now + Dt; % update time to that tragic time
    newly_dead = infected(1,idx); % mark that point as newly_dead for next loop iteration
    death_time(newly_dead) = Now; % record its absolute death time
    infected(:,idx) = []; % remove point from infected set
    infected(2,:) = infected(2,:) - Dt; % update all death timers
end
fprintf('\n');

end_time = Now;

P.init_insult = init_insult;
P.death_time = death_time;
P.end_time = end_time;
P.anim = @anim;
M.P = P;

%% Play (and record) Animation
    function anim(varargin)
        f = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
        function out = read_pair(name,def)
            ind = find(f(name));
            if ind, out = varargin{ind + 1};
            else out = def;
            end
        end
        movie_file = read_pair('movie_file', ''); % Record animation name (records if not empty)
        dt = read_pair('dt', 0.2);
        show_cont = any(f('show_cont')); % show infection contour
        show_dots = any(f('show_dots')); % show infected mesh points
        
        record_flag = ~isempty(movie_file); % record flag
        
        if record_flag
            wO = VideoWriter(movie_file);
            wO.FrameRate = 20;
            wO.Quality = 100;
            wO.open;
        end
        
        fig = figure('units', 'normalized', 'toolbar','figure', 'Name', 'Bundle Settings', 'outerposition', [0.2 0.15 .6 .8]);
        ax = axes('units', 'normalized', 'position', [0.1 0.1 .8 .8]);
        sld = uicontrol('style','slider', 'Units','Normalized', 'position', [.1 .95 .8 .05]);
        M.plot.model();
        plot(init_insult(1), init_insult(2), 'b.', 'markersize', 20);
        
        if show_cont
            upt = Contour(p, t, death_time);
        end
        
        if show_dots
            h_new = plot(inf, inf, 'g.', 'markersize', 10);
            h_inf = plot(inf, inf, 'r.', 'markersize', 10);
            setXY = @(h, p) set(h, 'Xdata', p(1,:), 'YData', p(2,:));
        end
        sld.Value = .51;
        tim = 0;
        while true
            try
                dt = ((sld.Value-.5)*end_time)/20;
                tim = tim + dt;
                tim = max(tim, dt);
                tim = min(tim, end_time);
                if show_dots
                    setXY(h_new, p(:,tim - dt < death_time & death_time < tim));
                    setXY(h_inf, p(:,tim  > death_time));
                end
                if show_cont, upt(tim); end
                drawnow, if record_flag, wO.writeVideo(getframe); end
            catch
                break;
            end
        end
        if record_flag, wO.close; end
    end
end

%% Plots the flooding contour
function upt = Contour(p, t, death_time)
xmin = min(p(1, t)); xmax = max(p(1, t));
ymin = min(p(2, t)); ymax = max(p(2, t));
nt = size(t, 2);
nxy = ceil(sqrt(nt/2)) + 1;
x = linspace(xmin, xmax, nxy);
y = linspace(ymin, ymax, nxy);
[xx, yy] = meshgrid(x, y);
zz=tri2grid(p, t, death_time, x, y);

[~, h] = contourf(xx, yy, max(death_time) - zz);
colormap(hot)
% caxis([-20,20])

upt = @Upt;
upt(0);
    function Upt(tim)
        h.LevelList = max(death_time) - tim;
    end
end
