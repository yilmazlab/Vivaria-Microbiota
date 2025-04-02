%% Function: Set environment and get external data for AGORA2 annotaions

function [exitflag] = fSetEnvironment(options)

fn       = fieldnames(options);
exitflag = true;

% check existence of directories / make directories
fprintf('[%s] Checking paths.\n',datestr(now));
fndir = fn(~cellfun(@isempty,regexp(fn,'dir')));

for z = 1:length(fndir)
    dirName = options.(fndir{z});
    if ischar(dirName)
        fl = exist(fullfile('.',dirName),'dir');
        if ~fl
            fprintf('\t\tDirectory "%s" does not exist - creating.\n',dirName);
            [fl,msg] = mkdir(pwd,dirName);
            if ~fl
                error(msg);
            end
        end
    end
end

% check existence of data files and fetch if not available
fprintf('[%s] Checking data files.\n',datestr(now));
fndata = fn(~cellfun(@isempty,regexp(fn,'fn')));

for z = 1:length(fndata)
    fileName  = options.(fndata{z});
    
    if ischar(fileName) && ~contains(fndata{z},'Save')
        status = 'ok';
        fl = exist(fullfile(options.dirData, fileName), 'file');
                
        if ~fl
            status = 'missing';            
            sw  = strrep(fndata{z},'fn','url');
            idx = strcmp(fn,sw);
            
            if sum(idx) == 1 % url defined
                url     = options.(fn{idx});
                fprintf('\tFetching data from %s\n',url);

                if ~contains(url,'.zip') && ~contains(url,'.tar.gz')
                    outfname = websave(fullfile(options.dirData,fileName), url);
                elseif contains(url,'.zip')
                    outfname = unzip(url,options.dirData);
                elseif contains(url,'.tar.gz')
                    tmpfname = fullfile(options.dirData,'tmp.tar.gz');
                    websave(tmpfname,url);
                    outfname = untar(tmpfname,options.dirData);
                    delete(tmpfname);
                end
                    
                if ~isempty(outfname) % files were fetched
                    status = 'fetched';
                    idx = strcmp(fn, strrep(fndata{z},'fn','proc'));
                    if sum(idx) == 1 % file processing required
                        exitflag = fprocess(options, fileName);
                    end
                    
                else
                    exitflag = false;
                end
            else
                status = 'url undefined';
                exitflag = false;
            end
        end
        fprintf('\t\t[%s]\t%s \n',status,fileName);
        
    end
end

return

%% Auxiliary functions

function [exitflag] = fprocess(options, fileName)

exitflag = false;

if contains(fileName, 'metanetx') % remove commentary header
    f0 = fullfile(options.dirData, 'tmp.txt');
    f1 = fullfile(options.dirData, fileName);
    movefile(f1, f0);
    fid0 = fopen(f0, 'rt');
    fid1 = fopen(f1, 'wt');
    lastfl = true;
    lastline = '';
    
    while true
        currline = fgetl(fid0);
        if ischar(currline)
            fl = startsWith(currline,'#');
            
            if ~fl && lastfl % column identifiers
                lastline = strrep(lastline, '#', '');
                fl = false;
            end
            
            if ~fl
                fprintf(fid1, '%s\n', lastline);
            end
            lastline = currline;
            lastfl   = fl;
        else
            fprintf(fid1, '%s\n', lastline);
            break
        end
    end
    fclose(fid0);    
    fclose(fid1);
    delete(f0);
    exitflag = true;
end

return
