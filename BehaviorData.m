
% load the behavior data

clear; clc;
load('/Volumes/Elham-Unifr/LongtermProject/ARCProject/Cognitive_performance_controls.mat');

SubIDs= [1;3;4;8;9;12;13;16;17;19;21;22;23;27;29;32;34;37;38;44;47;48;52;61;65;66;70;72;78;]; %these are the subject with REC and HSMT and Plast: 49 and 63 -> good subjects without plast
%%
[~,ind]=intersect(Subjects,SubIDs);

HO_SRate = HO_Pos_Res(ind)./(HO_Pos_Res(ind)+HO_Neg_Res(ind));
HO_trialn = (HO_Pos_Res(ind)+HO_Neg_Res(ind));
disp(['HO success rate = ' num2str(mean(HO_SRate)) ' +- ' num2str(std(HO_SRate))])
disp(['HO trial number = ' num2str(mean(HO_trialn)) ' +- ' num2str(std(HO_trialn))])
disp('TIME:')
disp(['HO time total= ' num2str(mean(HO_time_total)) ' +- ' num2str(std(HO_time_total))])
disp(['HO time second= ' num2str(mean(HO_time_second)) ' +- ' num2str(std(HO_time_second))])

