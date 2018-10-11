classdef RES
    % Resting-state EEG Spectrum
    properties
        SubjectInfo
        FFTData
        EpochLen
        Overlap
    end
    
    methods
        function obj = RES(SubjectInfo, FFTData,CondNames, EpochLen, Overlap)
            if exist('SubjectInfo','var')
                obj.SubjectInfo = SubjectInfo;
            else
                obj.SubjectInfo = [];
            end
            
            if exist('FFTData','var')
                obj.FFTData = cellfun(@(x,y) ARC.FData(x,y),FFTData,CondNames); % ARRAY of FData blocks
            else
                obj.FFTData = [];
            end
            if exist('EpochLen','var')
                obj.EpochLen = EpochLen;
            else
                obj.EpochLen = [];
            end
            if exist('Overlap','var')
                obj.Overlap = Overlap;
            else
                obj.Overlap = [];
            end
            
        end
        
        function saveRES(obj,Path)
            try
                SaveName = ['FFTData_' obj.SubjectInfo.SubID '_RecordNum' num2str(obj.SubjectInfo.Longitude+1)];
            catch
                error('SubjectInfo property should have SubID and Longitude fields');
            end
            try
                if ~exist('Path','var')
                    Path = [];
                end
                save(fullfile(Path,SaveName),'obj');
            catch
            end
        end
        
        function obj = loadRES(obj,Path, subinfo)
            try
                SaveName = ['FFTData_' subinfo.SubID '_RecordNum' num2str(subinfo.Longitude+1)];
            catch
                error('SubjectInfo should have SubID and Longitude fields');
            end
                if ~exist('Path','var')
                    Path = [];
                end
                if exist(fullfile(Path,SaveName),'file')
                    load(fullfile(Path,SaveName));
                else
                    obj=0;
                end

        end
    end
end

