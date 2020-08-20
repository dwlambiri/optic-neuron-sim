
%% Create Model
M = generate_model('file', 'mdl2', 'GUI', 'rewrite'); % Create an eye nerve model

M.plot.model();
% M.plot.histogram();

%% Process Model
M = M.create_csg(M, 'rewrite'); % Create Constructive Solid Geometry Model
M = generate_mesh(M, 'refine', 0, 'plot_mesh', 'rewrite'); % Generate Mesh

%% Run Simulation
M = propogation_alg(M, 'init_insult', [-0 -20], 'GUI', 'rewrite'); % Run Propogation Algorithm

% M.save(M); % Save the model?

%% Playback/record

M.P.anim('show_cont', 'dt', 200); % playback
% P.anim('show_cont', 'dt', 10 , 'movie_file', 'nerve1.avi'); % Record into avi file

