function [RESdata] = RawEEGtoRES(ProjectPath,SubjectInfo, EpLen, MovWin ,HF,LF)
% This function read the data of this subject in
% all conditions  and store it as ARC.RES (ARC toolbox, RES format for storing Resting EEG Structure)

% SubName: is the path to the stored ANT project in .cnet format
% written by: Elham Barzegaran, 10.2018
    
    
%%  Read the files
    SubName = fullfile(ProjectPath,[SubjectInfo.SubID '_Record' num2str(SubjectInfo.Longitude+1)]);
    files = subfiles(SubName);
    Conditions = {'rec-reofbp.cnt','hnfbp.cnt','hofbp.cnt','plastfbp.cnt','hn1fbp.cnt','hn2fbp.cnt','sc1fbp.cnt','sc2fbp.cnt','sc3fbp.cnt'};
    CondNames = { 'REC', 'REO', 'Plast', 'REC-mu', 'REC-hob', 'hob1', 'hob2', 'REC-hn', 'hn1' , 'hn2'};
    RData = cell(1,10);
    FData = cell(1,10);
    CData = cell(1,10);
    for Cond = 1:9
        display(['Reading the file:' Conditions{Cond}]);
        % Read the file and prepare events
        FileInd = find(~cellfun(@isempty,strfind(files,Conditions{Cond}))); 
        if ~isempty(FileInd)
            FileName = files{FileInd}; 
            [~,FileName,~] = fileparts(FileName); % 
            [DATA,EV,SEG] = prepare_events(fullfile(SubName,FileName));
            % average refrencing
            DATA.label{32}= 'CPz'; % I replace EOC by CPz;
            DATA.data=DATA.data(1:64,:); 
            NE = 64;
            H = eye(NE)- ((ones(NE,1)*ones(NE,1)')/(ones(NE,1)'*ones(NE,1)));
            DATA.data = H*DATA.data; % common average reference            
            % preprocess events and Calculating frequency components
            EpEv =events_preprocessing(EV,DATA, EpLen, MovWin,1);
            if (Cond==1 || Cond==4 )
                [FftData,CohData,cohf]= fftprocessing(DATA,EpEv,EpLen,0,HF,LF,1);
            else
                [FftData,CohData]= fftprocessing(DATA,EpEv,EpLen,0,HF,LF,1);
            end
            
            EpEv2 =events_preprocessing(EV,DATA, EpLen, EpLen,1);% For epoching the data, the overlap is zero
            ReData = epoching(DATA,EpEv2,EpLen);
            
            dataind = find(~cellfun(@isempty,FftData)); 
            for i= 1:numel(dataind)
                FData{dataind(i)}= cat(3,FData{dataind(i)},FftData{dataind(i)});
                RData{dataind(i)}= cat(3,RData{dataind(i)},ReData{dataind(i)});
                CData{dataind(i)}= cat(3,CData{dataind(i)},CohData{dataind(i)});
            end
        end
    end
RESdata = ARC.RES(SubjectInfo,FData,CondNames,EpLen,MovWin);
end

%% EVENT marking functions
function [DATA,EV,SEG] = prepare_events(filename)
% This function reads the data and events and prepare the events and
% segment data for furthur analysis
% INPUT: file name of data and xlsx including the events, they should have
        %the same name

% OUTPUT: DATA: EEG data
          %EV: the list of events and where they are marked in EEG time
          %seris
          % SEG: begining and the ending sample of each segment in EEG


% Elham Barzegaran, 25.11.2015
% e.barzegaran@gmail.com

% reads triggers and data from .cnt and .trg file
DATA= read_eep_cnt([filename '.cnt'],1,100);
DATA= read_eep_cnt([filename '.cnt'],1,DATA.nsample);
%C:\Users\ebarzega\Documents\My Works\long_term\Data\example\barzegaran_elham_2015-04-14_10-34-11-rec-reofbp.cnt

% reads events from .xlsx file 
% The xlsx file is made manually, since there is no function for reasing .evt file in matlab
% This file consists of following information
% the segment information with labels of 
    % Sn1 for beggining time,Sn2 for end time where n is number of segments
% the events with labels and codes as following, it is written in the first column
    % REC: 1, REO:2, Plast:3, REC-mu:4
    % REC-hob: 5, hob1: 6, hob2:7
    % REC-hn: 8, hn1: 9 , hn2:10
    % art: 11
% The second column is the time if begining of events in the recording local time hh:mm:ss.000
% The third column is the duration of the events in seconds
%
events = {'REC','REO','Plast','REC-mu','REC-hob','hob1','hob2','REC-hn','hn1','hn2','art'};
segments = {'S11','S12';'S21','S22';'S31','S32';'S41','S42'};

[~,txt,raw]=xlsread([filename '.xlsx']);
EV.time = zeros(1,3); % events time
SEG.time = zeros(1,3,2); % segments times (1-begin, 2-end)
evn = 1; % number of events
for ev = 1:size(txt,1)
    ts = regexp(txt{ev,2},':','Split'); % extract time string
    % prepare events
    if find(strcmp(events,txt(ev,1)))
        EV.time(evn,:)=str2double(ts); % event time
        EV.dur(evn,1)=cell2mat(raw(ev,3)); % events duration
        EV.label(evn,1)=txt(ev,1); % events label
        EV.code(evn,1)= find(strcmp(events,EV.label(evn,1))); % events code
        evn = evn+1;
    else if ~isempty(txt{ev})
    % prepare segments
        lt= txt(ev,1);lt= lt{1};% temporary label
        if strcmp(lt(1),'S')
            lttp= lt(end);lttp = str2num(lttp); % begin or end time of segment
            ltn = lt(2:end-1);ltn = str2num(ltn); % number of segment
            SEG.time(ltn,:,lttp) = str2double(ts);
        end
        end
    end  
end

clear ts ev txt raw lt ltn lttp evn segments events;

% convert the local time to time latency
% control for timing: 00 and 24 difference
[EV, SEG] = Preparetiming(EV,SEG);

% sort segments and calculation of segment duration 
for sg = 1: size(SEG.time,1)
    % I should add a part for sorting segments and making sure that they do not overlap
    SEG.dur(sg)= time2dur(squeeze(SEG.time(sg,:,1)),squeeze(SEG.time(sg,:,2)));
end
% check each event belongs to which segments
EV.segment = EventToSeg(EV,SEG);
% determines the time latency of events in the integrated file (all segments in one)
EV.timel = EventTimetoTimel(EV,SEG);
clear sg 
end
function [EV, SEG] = Preparetiming(EV,SEG)
% control for timing: 00 and 24 difference
TEV = EV.time;
TS = SEG.time;

%%
TSD = diff(TS,[],3);
TSDS = sign(TSD(:,1));
if sum(TSDS<0)>0
    Sind = find(TSDS<0);
    H1ind = TS (Sind,1,1);
    H2ind = TS (Sind,1,2);
    TS (Sind,1,2) = TS (Sind,1,2)+24;
    TS (Sind+1:end,1,:) = TS (Sind+1:end,1,:)+24;
    
    TEV(TEV(:,1)>=H2ind & TEV(:,1)<(H1ind-12),1) = TEV(TEV(:,1)>=H2ind & TEV(:,1)<(H1ind-12),1)+24;
end
EV.time = TEV;
SEG.time = TS;
end
function dur = time2dur(t1,t2)
% calculates the time difference of two time vectors t1 and t2 and returns
% the time difference in seconds

% Elham Barzegaran, 25.11.2015
% e.barzegaran@gmail.com

T = t2-t1;
dur = T(1)*3600+T(2)*60+T(3);
end
function EVS = EventToSeg(EV,SEG)
% This functions determines the events defined in structure EV belong to
% which segment define in structre SEG
% Input:
    %EV, events, EV.time: time of each events = ne*3 where ne is number of events and 3 is hh,mm,ss
    %SEG,Segments, SEG.time: time of each segments = ns*3 where ns is number of segments and 3 is hh,mm,ss
    
% OUtPUT:
    % EVS: ne*1 vector determines each event belong to each segment (gets values between 1 to ns)
    
    
% Elham Barzegaran, 25.11.2015
% e.barzegaran@gmail.com
    
Stime = SEG.time;
coef = [3600 60 1];
for ev = 1:size(EV.code,1)
    Etime = EV.time(ev,:);
    cmp = sum((repmat(Etime,[size(Stime,1) 1 2])-Stime).* repmat(coef,[size(Stime,1) 1 2]),2);
    EVS(ev,1)= find(squeeze(sum(sign(cmp),3))==0);
end
end
function timel= EventTimetoTimel(EV,SEG)

% Elham Barzegaran, 25.11.2015
% e.barzegaran@gmail.com


Stime = SEG.time;

for ev = 1:size(EV.code,1)
    Etime = EV.time(ev,:);
    timel(ev,1) = time2dur(squeeze(Stime(EV.segment(ev),:,1)),Etime)+sum(SEG.dur(1:EV.segment(ev)-1));
end

end
%% Preprocessing
function EpEv =events_preprocessing (EV,DATA, EpLen, MovWin,art)
% the events with labels and codes as following, it is written in the first column
    % REC: 1, REO:2, Plast:3, REC-mu:4
    % REC-hob: 5, hob1: 6, hob2:7
    % REC-hn: 8, hn1: 9 , hn2:10
    % art: 11
    
    
%% determine the label of time points
TLab = zeros(DATA.npnt,1);
if art ==1, Art=11; else Art=10;end
for ev= 1:Art
    ind = find(EV.code==ev);
    for i =1:size(ind)
        TLab(round(EV.timel(ind(i))*500:(EV.timel(ind(i))+EV.dur(ind(i)))*500))=ev; % time labeling of events
    end
end
clear ev i ind;
%% determine the epochs of analysis
EpEv= cell(10,1);
TLab(TLab==11)=0;
tlab = diff(TLab);tlab= [0; tlab];
% EpLen = 5000; % Epoch length: this can be changed ******** input parameters
% MovWin = 1000; % Moving Window: this can be changed ******** input parameters

for ev =1:10 
    ind1 = find(tlab==ev); % begining of events
    ind2 = find(tlab==-ev);% end of events
    if numel(ind1)==numel(ind2)
        epoch =[];
        for i = 1:numel(ind1)
            epoch = [epoch ind1(i):MovWin:ind2(i)-EpLen];
        end
        EpEv{ev}=epoch;
    else
        display(['Check out event #' num2str(ev)]);
    end
end
end
function RData = epoching(DATA,EpEv,EpLen)
% epochs data

for ev = 1:10
    EEG=[];
    EEGs=[];
    if ~isempty(EpEv{ev})
        EEG = zeros(size(DATA.data,1),EpLen,length(EpEv{ev}));
        Epoch = EpEv{ev};
        for ep=1:length(Epoch)
            EEG = DATA.data(:,Epoch(ep):Epoch(ep)+EpLen-1);
            EEGs(:,:,ep)=(EEG);
        end
    end
    RData{ev,1}=EEGs;

end
end
function [FData,CohData,cohf]= fftprocessing(DATA,EpEv,EpLen,Cohcom,hf,lf,param)
% calculate fft, if SurfL is 1, then surface laplacian is applied
% calculates coherence if Cohcom==1

if nargin<8
    param = 1;
end
load Sym_Elec;
Win= repmat(hann(EpLen)',[size(DATA.data,1) 1]);
f =500/2*linspace(0,1,EpLen/2+1);
HF=find(f<=hf,1,'last');
LF=find(f>=lf,1,'first');
for ev = 1:10

    FEEG=[];cm=[];
    if ~isempty(EpEv{ev})
        FEEG = zeros(size(DATA.data,1),length(LF:HF),length(EpEv{ev}));
        EEG=[];
        Epoch = EpEv{ev};
        for ep=1:length(Epoch)
            EEG = DATA.data(:,Epoch(ep):Epoch(ep)+EpLen-1);
            if param ==1 % Normalized power spepctrum
                feeg=fft((EEG.*Win)')'/EpLen;%(sum(Win(1,:)));
                PSD = 2*feeg.*conj(feeg)/(mean(Win(1,:).^2)*(500/EpLen));
                FEEG(:,:,ep) = PSD(:,LF:HF);
            elseif param ==2 % FFT 
                feeg = fft((EEG.*Win)')'/(sum(Win(1,:)));
                FEEG(:,:,ep) = abs(feeg(:,LF:HF));
            end
            
            if Cohcom==1
                for i = 1:size(ELsym,1)
                    [C F]=mscohere(EEG(ELsym(i,1),:),EEG(ELsym(i,2),:),[],[],2500,500);
                    cm(i,:,ep)= C((F>=0 & F<30));
                end
            end
        end
    end
    FData{ev,1}=FEEG;
    CohData{ev,1}= cm;
end

if Cohcom==1
    cohf= F(F>=0 & F<30);
else
    cohf=[];
end

end