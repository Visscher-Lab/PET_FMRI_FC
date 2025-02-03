clear all
close all
clc
dataLocation = 'C:\Users\David\Documents\Tau_FC_Data\Gordon_344';
cd(dataLocation);

load('connMatrix_Gordon.mat');
load('suv_parcel_Gordon.mat');
load("Gordonp.mat");
Gordonp = sortrows(Gordonp,'Community','ascend');
p = Gordonp;


suv_gordon = squeeze(suv_pet_parcel_weightAdj_matrix);

network_names = unique(string(Gordonp.Community));
network_rValues = zeros(1,length(network_names));
for network_name = 1:length(network_names)

    network = network_names(network_name);

    withinVals = ones(size(connMatrix,1), size(connMatrix,1))* NaN; % make a matrix of NaNs
    withinVals(p.Community == network, p.Community == network) = 1; % make it NaNs where it isn't what we want
    withinValsZeroed = withinVals;
    withinValsZeroed(isnan(withinValsZeroed)) = 0;
    
    betweenVals = ones(size(connMatrix,1), size(connMatrix,1))* NaN; % make a matrix of NaNs
    betweenVals(p.Community == network, p.Community ~= network) = 1; % make it NaNs where it isn't what we want
    betweenValsZeroed = betweenVals;
    betweenValsZeroed(isnan(betweenValsZeroed)) = 0;
    
    sumValsZeroed = withinValsZeroed + betweenValsZeroed;
    
    figure;
    subplot(1,3,1)
    imagesc(withinValsZeroed)
    set(gca,'FontSize',24)
    title(['Within Comparisons, ' network ' Network'])
    
    subplot(1,3,2)
    imagesc(betweenValsZeroed)
    set(gca,'FontSize',24)
    title(['Between Comparisons, ' network ' Network'])
    
    subplot(1,3,3)
    imagesc(sumValsZeroed)
    set(gca,'FontSize',24)
    title(['All Comparisons, ' network ' Network'])
    
    subplot(1,3,3)
    imagesc(sumValsZeroed)
    set(gca,'FontSize',24)
    title(['All Comparisons, ' network ' Network'])

    % remove the diagonals
    for i = 1: size(withinVals,1)
        for j = 1: size(withinVals,2)
            if i==j
                withinVals(i,j) = NaN;
                betweenVals(i,j) = NaN;
            end
        end
    end


    for s = 1: size(connMatrix,3) % for each subject

        within(s) = mean(mean(withinVals.* connMatrix(:,:,s), "omitnan"), "omitnan");
        within_currentSubject = withinVals.* connMatrix(:,:,s);
        between(s) = mean(mean(betweenVals.* connMatrix(:,:,s), "omitnan"), "omitnan");
        between_currentSubject = betweenVals.* connMatrix(:,:,s);
        seg(s) = (within(s)-between(s))/(within(s)+between(s));
         % Sum the arrays while preserving NaNs
    summed_array = NaN(size(withinVals));
        for i = 1:333
            for j = 1:333
                if ~isnan(withinVals(i,j)) && ~isnan(betweenVals(i,j))
                    summed_array(i,j) = within_currentSubject(i,j) + between_currentSubject(i,j);
                elseif ~isnan(within_currentSubject(i,j))
                    summed_array(i,j) = within_currentSubject(i,j);
                elseif ~isnan(between_currentSubject(i,j))
                    summed_array(i,j) = between_currentSubject(i,j);
                end
            end
        end
    end
    
       
    

    network_indices = find(p.Community==network);
    
    figure;
    hold on;
    subplot(1,4,1);
    imagesc(connMatrix(:,:,s));
    caxis([0 0.5])
    
    subplot(1,4,2);
    imagesc(within_currentSubject);
    caxis([0 0.5])
    subplot(1,4,3);
    imagesc(between_currentSubject)
    caxis([0 0.5])
    
    
    subplot(1,4,4);
    imagesc(summed_array);
    caxis([0 0.5])
    

    
    colorbar(gca)
    parcel_sizes = Gordonp.SurfaceAreamm2(network_indices);
    
    suv_network = suv_gordon(network_indices,:);



    tau_burden = zeros(1,size(suv_network,2));
    for subject = 1:size(suv_network,2)
        tau_burden(subject) = sum(suv_network(:,subject) .* parcel_sizes)/sum(parcel_sizes);
    end

    figure;
    subplot(1,2,1)
    hold on
    swarmchart(ones(length(seg),1),seg,96,'Filled')
    title(['Network Segregation Values, ' network ' Network'])
    set(gca,'FontSize',24)
    ylabel('Network Segregation Value')

    subplot(1,2,2)
    hold on;
    swarmchart(ones(1,length(tau_burden)),tau_burden,96,[1 0.5 0],'filled');
    ylabel('Tau SUV')
    ylim([0 3])
    title(['Tau SUV, ' network ' Network'])
    set(gca,'FontSize',24)

    figure; 
    hold on; 
    scatter(seg, tau_burden,96,[1 0.5 0.5],'filled');
    set(gca,'FontSize',24)
    ylabel('Tau SUV')
    xlabel('Network Segregation Value')
    title(['Tau Burden vs Network Segregation Values, ' network ' Network'])

    ylim([0 3])
    
    rValue = corrcoef(seg,tau_burden);
    rValue = rValue(1,2);
    network_rValues(network_name) = rValue;
    
    save(['seg_values_' + network + '.mat'],'seg')
    save(['tau_values_' + network + '.mat'],'tau_burden')
    
end


figure;
imagesc(suv_gordon')
caxis([0 2])
ylabel('Subject')
xlabel('Region Index')
colorbar
set(gca,'FontSize',24)