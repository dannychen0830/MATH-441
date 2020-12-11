d = 0.8; % density
n = 200; % grid dimension 
N = d*n^2; % number of trees

% initialize trees
tree_map = zeros(n,n);
perm = randperm(n^2);
for i = 1:N
    [x,y] = ind2sub([n,n],perm(i));
    tree_map(x,y) = 1;
end

trees = find(tree_map == 1);
perm = randperm(N);
count = 5;
for i = 1:count
    [x,y] = ind2sub([n,n],trees(perm(i)));
    tree_map(x,y) = 2;
end

temp_map = zeros(n,n);
oxy_map = zeros(n,n);

burning = find(tree_map == 2);
for i = 1:n
    for j = 1:n
        temp = 0;
        oxy = 1;
        for k = 1:length(burning)
            [x,y] = ind2sub([n,n],burning(k));
            temp = temp + 600*exp(-1*dist([x,y],[i,j])/8);
            oxy = oxy*(1 - 0.5*exp(-1*dist([x,y],[i,j])/32));
        end
        temp_map(i,j) = temp;
        oxy_map(i,j) = oxy;
    end
end

for t = 1:200
    
    for i = 1:n
        for j = 1:n
            if tree_map(i,j) == 1 && temp_map(i,j) > 300 && oxy_map(i,j) > 0.2
                tree_map(i,j) = 2;
            end 
            if tree_map(i,j) == 2 && rand < 0.7
                tree_map(i,j) = 0;
            end 
        end
    end
    
    figure(1)
    imagesc(tree_map);
    title('Tree Map')
    colormap parula 
    axis('square')
    
    burning = find(tree_map == 2);
    for i = 1:n
        for j = 1:n
            temp = 0;
            oxy = 1;
            for k = 1:length(burning)
                [x,y] = ind2sub([n,n],burning(k));
                temp = temp + 600*exp(-1*dist([x,y],[i,j])/8);
                oxy = oxy*(1 - 0.5*exp(-1*dist([x,y],[i,j])/32));
            end
            temp_map(i,j) = temp;
            oxy_map(i,j) = oxy;
        end
    end
    
    figure(2)
    imagesc(temp_map)
    title('Temperature Map')
    colormap hot
    axis('square')

    figure(3)
    imagesc(oxy_map)
    title('Oxygen Map')
    colormap gray
    axis('square')
    
    if t == 1
        w = waitforbuttonpress;
    end
end

function d = dist(u,v)
    d = (u(1)-v(1))^2 + (u(2)-v(2))^2;
end
