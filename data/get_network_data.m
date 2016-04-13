function [net_name,net_matfile,net_vertices,net_edges,net_unique_edges,...
    net_type,net_weight] =  get_network_data(num_edg_low,num_edg_up,sort_by,filter_by_type,filter_by_weight)
% Loads the meta data for koblenz konect networks dataset
%    * undirected network are loaded as only upper triangular
%
% Input:
%     num_edg_low - lower limit of networks with number of edges to include
%     num_edg_up  - upper limit of networks with number of edges to include
%     sort_by - controls how to sort the list of networks 
%                 possilbe values: 'edges' (default) ,'vertices','unique_edges'
%     filter_by_type - cell array to filter by the network type
%                 possible values: 'undirected','directed','bipart'
%                 default value: {'undirected','directed'}
%     filter_by_weight -  cell array to filter by the network type. 
%                 possible value: 'simple_edge','multiple_unweighted','positive_weights',...
%                    'signed','unknown','ratings','mutiple_ratings','dynamic';
%                 default value: {'simple_edge','multiple_unweighted','positive_weights',...
%                    'unknown','ratings','mutiple_ratings','dynamic'};
% 
%   Output:
%     net_name         - the names of the networks
%     net_matfile      - the file name of each network
%     net_vertices     - the number of vertrices for each network
%     net_edges        - the number of edges for each networks
%     net_unique_edges - the number of unique edges for each network 
%                 (for dynamic networks many edges repeat multiple times)
%     net_type         - the type of edges - directed, undirected or bipart
%     net_weight       - the type of weight for the edges
        
if ~exist('sort_by','var')
    sort_by = 'edges';
end
if ~exist('filter_by_type','var')
    filter_by_type = {'undirected','directed'}; %,'bipart'};
end
if ~exist('filter_by_weight','var')
    filter_by_weight = {'simple_edge','multiple_unweighted','positive_weights',...
            'unknown','ratings','mutiple_ratings','dynamic'}; % ,'signed'};
end



net_meta = textscan(fopen('konect_meta.csv'),'%q %q %q %q %f %f %f %f %f %f','delimiter',',','Headerlines',1);
net_url = net_meta{1};
% konect_code = konect_meta{2};
net_name = net_meta{3};
net_type = net_meta{5};
net_weight = net_meta{6};
net_vertices = net_meta{8};
net_edges = net_meta{9};
net_unique_edges = net_meta{10};

net_types_dict = {'undirected','directed','bipart'};
net_type = net_types_dict(net_type);

net_weight_dict = {'simple_edge','multiple_unweighted','positive_weights',...
    'signed','unknown','ratings','mutiple_ratings','dynamic'};
net_weight = net_weight_dict(net_weight);


net_matfile = cell(size(net_url));
for i = 1:length(net_url);
    if strcmp(net_url{i} ,'')
        net_matfile{i} = '';
    else
        filename = strsplit(net_url{i},'/');
        filename = strsplit(filename{end},'.');
        filename = filename{1};
        net_matfile{i} = sprintf('net_mat/%s.mat',filename);
    end
end

% if do_filter
    filter_net = true(size(net_edges));
    filter_net = filter_net & num_edg_low <= net_edges;
    filter_net = filter_net & net_edges <= num_edg_up;
    filter_net = filter_net & ismember(net_type, filter_by_type)';
    filter_net = filter_net & ismember(net_weight,filter_by_weight)';
    filter_net = filter_net & ~strcmp(net_matfile, '');

    net_edges = net_edges(filter_net);
    net_matfile = net_matfile(filter_net);
    net_name = net_name(filter_net);
    net_vertices = net_vertices(filter_net);
    net_type = net_type(filter_net);
    net_weight = net_weight(filter_net);
    net_unique_edges = net_unique_edges(filter_net);
% end

switch sort_by
    case 'edges'
        [~, sort_ind] = sort(net_edges);
    case 'vertices'
        [~, sort_ind] = sort(net_vertices);
    case 'unique_edges'
        [~, sort_ind] = sort(net_unique_edges);
    otherwise
        sort_ind = 1:length(net_edges);
end
    
net_matfile = net_matfile(sort_ind);
net_name = net_name(sort_ind);
net_vertices = net_vertices(sort_ind);
net_type = net_type(sort_ind);
net_weight = net_weight(sort_ind);
net_edges = net_edges(sort_ind);
net_unique_edges = net_unique_edges(sort_ind);

end