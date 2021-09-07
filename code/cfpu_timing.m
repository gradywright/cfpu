% Evaluation grid size
mm = 384-7;

% Regularization parameters:
reginfo.lambda = 1e-3;
reginfo.regularization = 0;
reginfo.schurcmplmnt = 0;
reginfo.exactinterp = 1;
reginfo.regularizationi = 0;
reginfo.lambdai = 1e-4;

% Kernels to use
% RBFs to use
kernelinfo.phi = @(r) -r;
% kernelinfo.eta = @(r) r.^3;
% kernelinfo.zeta = @(r) 3*r;  kernelinfo.order = 2;
kernelinfo.eta = @(r) -r;
kernelinfo.zeta = @(r) -1./r; kernelinfo.order = 1;

ptclouds = {'raptor_head','filigree','pump','dragon','buddha','armadillo'};
cores = [2];

for j = 1:length(cores)
    core = cores(j);
    if core ~= 1
%         pool = parpool(core);
    end
    for k = 1:length(ptclouds)
        ptcloud = ptclouds{k};
        switch ptcloud
            case 'trefoil'
                load('../ptclouds/trefoil_N23064.mat');
                % load('../ptclouds/trefoil_N6144.mat');
                load('../ptclouds/trefoil_patches_M864.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'buddha'
        %         load('../ptclouds/happy_buddha.mat');
                load('../ptclouds/happy_buddha_highres.mat');
        %         load('../ptclouds/happy_buddha_highres_patches.mat');
                load('../ptclouds/happy_buddha_highres_patches_M42861.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'homer'
                load('../ptclouds/homer.mat');
                load('../ptclouds/homer_patches.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'raptor_head'
                load('../ptclouds/raptor_head.mat');
        %         load('../ptclouds/raptor_head_patches.mat');
                load('../ptclouds/raptor_head_patches_M10337.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'raptor'
                load('../ptclouds/raptor.mat');
                load('../ptclouds/raptor_patches.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'armadillo'
                load('../ptclouds/armadillo.mat');
                load('../ptclouds/armadillo_patches.mat');
                Q = eye(3);
            case 'children'
                load('../ptclouds/dancing_children.mat');
                load('../ptclouds/dancing_children_patches.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'pump'
                load('../ptclouds/pump_carter.mat');
                load('../ptclouds/pump_carter_patches.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'filigree'
                load('../ptclouds/filigree.mat');
        %         load('../ptclouds/filigree_patches.mat');
                load('../ptclouds/filigree_patches_M35130.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'gargoyle'
                load('../ptclouds/gargoyle.mat');
                load('../ptclouds/gargoyle_patches.mat');
                Q = [0 0 1;-1 0 0;0 -1 0];
            case 'bunny'
                load('../ptclouds/bunny_large.mat')
                load('../ptclouds/bunny_large_patches.mat')
                Q = [0 0 1;1 0 0;0 1 0];
            case 'dragon'
                load('../ptclouds/stanford_dragon_fullres_face_nrmls.mat'); 
                [xc,ia,ic] = unique(xc,'rows'); nrml = nrml(ia,:);
        %         load('../ptclouds/stanford_dragon_fullres_patches_N14400.mat');
                load('../ptclouds/stanford_dragon_patches_M32527.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'tori'
                load('../ptclouds/interlocked_tori_big.mat');
                load('../ptclouds/interlocked_tori_big_patches.mat');
                Q = eye(3);
            case 'hand'
                load('../ptclouds/laurent_hand.mat');
                load('../ptclouds/laurent_hand_patches.mat');
                Q = [0 0 1;1 0 0;0 1 0];
            case 'cantius_tooth'
                load('../ptclouds/cantius_tooth.mat');
                rng(6232020)
                [~,y] = kmeans(xc,200,'Replicates',5);
                Q = eye(3);
            case 'mammoth_tooth'
                load('../ptclouds/mammoth_tooth.mat');
                rng(6232020)
                [~,y] = kmeans(xc,200,'Replicates',5);
                Q = eye(3);
            case 'frog'
                load('../ptclouds/frog.mat');
                rng(6232020)
                [~,y] = kmeans(xc,100,'Replicates',5);
                Q = eye(3);
            otherwise
                load('../ptclouds/bunny_large.mat')
                load('../ptclouds/bunny_large_patches.mat')
                Q = [0 0 1;1 0 0;0 1 0];
        end
        x = xc;

        % Add noise to the normals
        % id = 1:N;
        % rng(1378271);
        % nrml(id,:) = nrml(id,:) + 0.3*randn(size(nrml(id,:)));

        % Make the normals unit vectors
        nrml = nrml./(sqrt(sum(nrml.^2,2)));

        % Rotate or shift the points
        x = x*Q';
        y = y*Q';
        nrml = nrml*Q';

        fittime = zeros(1,1);
        evaltime = fittime;

        for i = 1:1
            tic
            puminfo = cfpufit(x,nrml,y,kernelinfo,reginfo);
            fittime(i) = toc;

            tic
            [potential,X,Y,Z] = cfpuval(puminfo,mm);
            evaltime(i) = toc;            
        end
        fprintf('%s, cores = %d, fit time = %1.2f eval time = %1.2f\n\n',ptcloud,core,min(fittime),min(evaltime))
        size(X)
    end
    if core ~= 1
        delete(pool);
    end    
end