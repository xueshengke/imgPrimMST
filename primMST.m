%% thie script generate a minimum spanning tree from a downsampled image
% adjustable parameters:
%     scale = 4, 2, 1
%     neighbor = 8, 4 
%     begin = highest, center

clear all;
% close all;
clc;

scale = 1;      % 4, 2 or 1 dowbsample ratio
neighbor = 8;   % 8 or 4 neighbors
begin = 'highest';
fprintf('Set parameters: scale = %d, neighbor = %d, begin = %s\n', scale, neighbor, begin);

imgFile = './8_1.jpg';  % load a gray image
img = imread(imgFile);
if ndims(img) > 2
    img = img(:, :, 1);
end
img = single(img); 
dims = size(img);

figure;
subplot(2, 2, 1);

imagesc(img); title(['original image, ', num2str(dims(1)), ' x ', num2str(dims(2))]); axis off;
fprintf('Using image file: %s\n', imgFile);
%% calculate neighbor costs in MST model
if neighbor == 8    % define the 8-neighbor set
    dn_area = [1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1;  1 -1];
elseif neighbor == 4    % define the 4-neighbor set
    dn_area = [1 0; 0 1; -1 0; 0 -1];
end
nodes = floor(dims ./ scale);  % nodes after scaling the original image in graph model
edges = prod(nodes) * neighbor; % edges for all nodes (each one has 8 or 4)
edgeCost = zeros(nodes(1), nodes(2), neighbor); % cost for each edge

tic;
intensity = conv2(img, ones(scale), 'valid');   % intensity matrix
intensity = intensity(1:scale:end, 1:scale:end);% each element is the sum of a group of pixels

for i = 1 : neighbor    % calculate the costs matrix for each neighbor
    edgeCost(:, :, i) = abs(intensity - imgshift(intensity, dn_area(i, :)));
end
%% generate MST model
nodeStatus = zeros(nodes);  % status = 1, in tree; status = 0, outside tree
if strcmp(begin, 'highest')
    [root_x, root_y] = find(intensity == max(intensity(:))); % take the maximum intensity as the root node
elseif strcmp(begin, 'center')
    root_x = floor(nodes(1) / 2); root_y = floor(nodes(2) / 2);
end
nodeStatus(root_x(1), root_y(1)) = true; % root node status = 1
treeNodes = [root_x(1), root_y(1)]; % add the root node in tree
treeEdges = [];
edgeCandidate = [];
minEdge = [0, 0, 0, 0, Inf]; 
newNode = [];

for i = 1 : prod(nodes)  % each node in tree
    pos= treeNodes(i, :);   % get current node's position
    pos_n = repmat(pos, [neighbor 1]) + fliplr(dn_area); % neighbor position ...
%     (dimension and coordinate have different origins, that is why fliplr is utilized)
    pos_n(:, end+1) = 1:neighbor;   % record edge's order
    [out1, ~] = find(pos_n(:, 1:end-1) < 1);    % find positions out of boundary
    [out2, ~] = find(pos_n(:, 1:end-1) > repmat(nodes, [neighbor 1]));
    outIndex = union(out1, out2);
    pos_n(outIndex, :) = [];    % remove positions outside intensity matrix
    inTreeIndex = nodeStatus((pos_n(:, 2)-1)*nodes(1) + pos_n(:, 1)) == 1;
    pos_n(inTreeIndex, :) = []; % remove positions already in tree

    list = [repmat(pos, [size(pos_n,1) 1]), pos_n(:, 1:end-1), ... % construct an edge list
        reshape(edgeCost(pos(1), pos(2), pos_n(:, end)), [size(pos_n,1) 1])];
    edgeCandidate = [edgeCandidate; list]; % merge list with previous one

    [edgeValue, index] = min(edgeCandidate(:,end));     % find the minimum edge
    minEdge = edgeCandidate(index, :);      % read min edge information
    newNode = edgeCandidate(index, 3:4);        % read new node
    edgeCandidate(index, :) = [];       % remove the min edge
    unecessIndex = ...       % find some unecessary edges may occur
        prod(edgeCandidate(:, 3:4) == repmat(newNode, [size(edgeCandidate, 1) 1]), 2) == 1;
    edgeCandidate(unecessIndex, :) = [];    % remove edges

if ~isempty(newNode)    % find a new node with minimum edge
    treeNodes = [treeNodes; newNode];   % add the node into tree
    nodeStatus(newNode(1), newNode(2)) = true;  % update tree status
end
if ~isempty(minEdge)   % find a valid edge
    treeEdges = [treeEdges; minEdge];   % add edge into edge set
end

end

fprintf('MST calculation completes. ');
toc;
%% illustration for MST model
subplot(2, 2, 2);
imagesc(intensity); title(['downsampled image, ', num2str(nodes(1)), ' x ', num2str(nodes(2))]); axis off;

subplot(2, 2, 3);
drawMST(treeNodes, treeEdges, nodes); title('generated MST'); 
