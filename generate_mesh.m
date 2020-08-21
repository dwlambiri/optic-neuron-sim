function M = generate_mesh(M, varargin)
% GENERATE_MESH
% parameters:
%   Type, Value
%       refine
%   Binary
%       rewrite
%

f = @(g) (cellfun(@(x) ischar(x) && strcmp(x,g), varargin));
    function out = read_pair(name, default)
        ind = find(f(name));
        if ind
            out = varargin{ind + 1};
        else
            out = default;
        end
    end
refine = read_pair('refine', 0);

if isfield(M, 'mesh') && ~any(f('rewrite')), return, end

dl = M.csg.dl;

fprintf('Generating mesh... ');
[p,e,t] = initmesh(dl);
fprintf('DONE\n');

if refine > 0
    fprintf('Refining mesh... ');
    for k = 1:refine
        [p,e,t] = refinemesh(dl, p,e,t);
    end
    fprintf('DONE\n');
end

p = jigglemesh(p,e,t);

M.mesh.p = p; M.mesh.e = e; M.mesh.t = t;

M.plot.mesh = @plot_mesh;

    function plot_mesh
        fprintf('Plotting... ');
        h = pdemesh(p,e,t); axis equal
        h(2).Visible = 'off';
        % set(h(2), 'LineWidth', 2, 'Color', [.9 .2 .5]); % Boundaries linewidth and color
        hold on, M.plot.model('LineWidth', 2);
        fprintf('DONE\n');
        drawnow
    end

end