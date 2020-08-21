function lhon_sim(varargin)

f = @(g) (cellfun(@(x) ischar(x) && strcmp(x,g), varargin));

    function out = read_pair(name,def)
        ind = find(f(name));
        if ind
            out = varargin{ind + 1};
        else
            out = def;
        end
    end

neuron_scale = read_pair('neuron_scale', 30);
%% Create Model
M = generate_model('neuron_scale', neuron_scale, 'file', 'mdl2', 'rewrite', 'GUI', 'refine', 0, 'init_insult', [-0 -20]); % Create an eye nerve model

%% Playback/record

% if any(f('show_dots'))
%     if any(f('plot_mesh'))
%         M.P.anim('plot_mesh', 'show_dots', 'dt', 200); % playback
%     else
%         M.P.anim('show_dots', 'dt', 200); % playback
%     end
% else
%     if any(f('plot_mesh'))
%         M.P.anim('plot_mesh', 'show_cont', 'dt', 200); % playback
%     else
%         M.P.anim('show_cont', 'dt', 200); % playback
%     end
% end
% % M.P.anim('show_dots', 'dt', 200 , 'movie_file', 'nerve1.avi'); % Record into avi file

end