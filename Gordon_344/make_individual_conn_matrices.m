clear all
close all

clc


dataLocation = 'C:\Users\David\Documents\Tau_FC_Data\Gordon_344';
individual_analysis_dir = 'C:\Users\David\Documents\Tau_FC_Data\Gordon_344\Individual_Subject_Analyses';
cd(dataLocation);

load('connMatrix_Gordon.mat');
load('suv_parcel_Gordon.mat');
load('Gordonp.mat');
Gordonp = sortrows(Gordonp,'Community','ascend');
p = Gordonp;


suv_gordon = squeeze(suv_pet_parcel_weightAdj_matrix);

network_names = unique(string(Gordonp.Community));
network_rValues = zeros(1,length(network_names));

cd(individual_analysis_dir)

for subject = 1:size(suv_gordon,2);
    if subject >= 1 && subject <=9
        subjectName = ['00' num2str(subject)];
    end
    if subject >= 10 && subject <=99
        subjectName = ['0' num2str(subject)];
    end
    if subject >= 100
        subjectName = num2str(subject);
    end
    if ~isdir(subjectName)
        mkdir(subjectName);
    end
    cd(subjectName);
    
    
end

figure;