%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%    IMPORT OF .PHI FILE IN MATLAB WORKSPACE                              %
%                                                                         %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%INSERT THE REQUIRED DIRECTORIES:
% DIRECTORY OF THE INPUT FILE .phi IN THE SAME FOLDER IT WILL BE SAVED THE output.m file:
directoryLoad='Z:\MSc Mathematical Engineering\Third Semester\Fluid Labs\CFD\CFD2\main simulation\LowRe.phi';
directorySave='Z:\MSc Mathematical Engineering\Third Semester\Fluid Labs\CFD\CFD2\main simulation\LowRe.mat';

% INSERT YES IF THE CASE IS BODY FITTED (IS THE GRID IS CARTESIAN LEFT NO)
BFC='no'  
% INSERT THE DIRECTORY OF NON CARTESIAN FILES (THIS IS NOT MANDATORY)
directory_bfc='F:\DIDATTICA\2018-19\FL\ReferenceCases\RC1_1\xyz'; 

% DEFINE WHICH VARIABLES YOU WANT TO PROCESS
variables={'P1' 'V1' 'W1' 'DWDY' 'KE' 'EP' 'YPLS' 'STRS' 'ENUT'}; %to let process all the variables write 'all' (pay attention to time required)


% VARIABLES SAVED IN WORKSPACE
% NX, NY, NZ                number of cells in the three directions
% X_E, Y_N, Z_H             coordinates of cell's faces
% X_C, Y_C, Z_C             coordinates of cell's centers
% xn3D, yn3D, zn3D          position of the cell's nodes in case of BFC simulation  (NX+1xNY+1xNZ+1)
% xc3D, yc3D, zc3D          position of the cell's nodes in case of BFC simulation  (NXxNYxNZ)
% P1, U1, V1, W1 ...        matrix (NXxNYxNZ) of the variables
% U1_C, V1_C, W1_C          matrix (NXxNYxNZ) of velocities interpolated in cell's centers
% XX_C, YY_C, XX_E, YY_N    meshgrid 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %Procedure
%Lettura del file .phi e importazione delle componenti fondamentali.
fnameopen=[directoryLoad];
fnameopen_bfc=[directory_bfc];

newData1 = importdata(fnameopen, ' ', 2);
Parte_1=newData1.data(1:2,:);           %numero di celle
clear newData1;

formatSpec = '%s%[^\n\r]';
fileID = fopen(fnameopen,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '' , 'ReturnOnError', false);
fclose(fileID);
Parte_01=dataArray{1,1}(6:13);
clearvars  formatSpec fileID dataArray
                                                            tic
kk=1;
for i=1:8
    if i==8
        k=2;
    else k=1;
    end
    for j=1:19-k
        Parte_02(1:4)=(Parte_01{i,1}(2+4*(j-1):4*(j-1)+5));           %Nomi Variabili
        Parte_0(kk)={Parte_02};
        kk=kk+1;
    end
end
%Parte_0(strcmp(Parte_0,'nul '))=[];
clear newData1 i j k kk;
                                                             time(1)=toc;

N_col=6;
NX=Parte_1(1,1);
NY=Parte_1(1,2);
NZ=Parte_1(1,3);
txtrh_parte_2=13;
txtrh_parte_2_bfc=1;
txtrh=15+ceil(NX/N_col)+ceil(NY/N_col)+2*ceil(NZ/N_col);



                                                            tic
if BFC=='no'
    formatSpec = '%13s%13s%13s%13s%13s%s%[^\n\r]';
    fileID = fopen(fnameopen,'r');

    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,txtrh_parte_2, 'ReturnOnError', false);
    fclose(fileID);

    rowstart=ceil(NX/6)+ceil(NY/6)+2*ceil(NZ/6)+2+1;  
    for col=[1:N_col]
        rawData = dataArray{col};
        for row=rowstart:length(rawData)-2
        raw3(row-(rowstart-1),col) =str2double(rawData{row});
        end
        for row=1:length(rawData)-2
        raw2(row,col) =str2double(rawData{row});
        end
    end
    
    Parte_2 = (raw2);  %variabili salvate in ordine di variabile per piano x-y
    Parte_3 = (raw3);  %variabili salvate in ordine di variabile per piano x-y
    %%% Clear temporary variables
    clearvars phi startRow formatSpec fileID dataArray ans raw2 raw3 col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

  else
    formatSpec = '%13s%13s%13s%13s%s%[^\n\r]';
    fileID = fopen(fnameopen_bfc,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,txtrh_parte_2_bfc, 'ReturnOnError', false);
    fclose(fileID);
    %%% Convert the contents of columns containing numeric strings to numbers.
    % Replace non-numeric strings with NaN.
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));
    for col=[1,2,3,4,5]
        % Converts strings in the input cell array to numbers. Replaced non-numeric
        % strings with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1);
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData{row}, regexstr, 'names');
                numbers = result.numbers;
                
                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if any(numbers==',');
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(numbers, thousandsRegExp, 'once'));
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric strings to numbers.
                if ~invalidThousandsSeparator;
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch me
            end
        end
    end
    %%% Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells
    %%% Create output variable
    phi = cell2mat(raw);
    Parte_2=phi; %variabili salvate in ordine di variabile per piano x-y
    %%% Clear temporary variables
    clearvars phi startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
    
    
    
end

                                                        time(2)=toc;

%%
                                                          
                                                          tic
formatSpec = '%s%[^\n\r]';
fileID = fopen(fnameopen,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '' , 'ReturnOnError', false);
fclose(fileID);
Parte_00=dataArray{1,1}(txtrh-1:txtrh);
clearvars  formatSpec fileID dataArray

NVariabili = (length(find(Parte_00{1}=='T')))+(length(find(Parte_00{2}=='T')))   %numero variabili salvate nel phi
clear newData1;
Var=[Parte_0(find(Parte_00{1}=='T')-1) Parte_0(find(Parte_00{2}=='T')-1+72)];
M=NVariabili;

operazione=('acquisizione')

                                                            time(4)=toc;
                                                            
                                                            tic
if BFC=='no'
    %Acquisizione dati della geometria e costruzione dei tre vettori
    %X_E,Y_N,Z_H (NB questa parte funziona sicuramente solo SE NX,NY>6).
    
    Nrighe_x=ceil(NX/6);
    Nrighe_y=ceil(NY/6);
    Nrighe_z=ceil(NZ/6);
    
    if Nrighe_x==1
        X_E=Parte_2(1,1:NX);
    else
        for irighe_x=1:Nrighe_x-1
            X_E((1+6*(irighe_x-1)):(6*irighe_x))=Parte_2(irighe_x,1:6);
        end
        X_E((1+6*(Nrighe_x-1)):NX)=Parte_2(Nrighe_x,1:(NX-6*(Nrighe_x-1)));
    end
    
    if Nrighe_y==1
        Y_N=Parte_2(Nrighe_x+1,1:NY);
    else
        for irighe_y=1:Nrighe_y-1
            Y_N((1+6*(irighe_y-1)):(6*irighe_y))=Parte_2(Nrighe_x+irighe_y,1:6);
        end
        Y_N((1+6*(Nrighe_y-1)):NY)=Parte_2(Nrighe_x+Nrighe_y,1:(NY-6*(Nrighe_y-1)));
    end
    
    if Nrighe_z==1
        Z_H=Parte_2(Nrighe_x+Nrighe_y+1,1:NZ);
    else
        for irighe_z=1:Nrighe_z-1
            Z_H((1+6*(irighe_z-1)):(6*irighe_z))=Parte_2(Nrighe_x+Nrighe_y+irighe_z,1:6);
        end
        Z_H((1+6*(Nrighe_z-1)):NZ)=Parte_2(Nrighe_x+Nrighe_y+Nrighe_z,1:(NZ-6*(Nrighe_z-1)));
    end
    
    operazione=('scrittura coordinate facce celle')
                                                                time(5)=toc;
    %Costruzione dei tre vettori X_C, Y_C, Z_C, contenenti il centro delle
    %celle in ciascuna delle 3 direzioni.
                                                                 tic   
    X_C(1)=X_E(1)/2;
    for i=2:NX
        X_C(i)=(X_E(i)+X_E(i-1))/2;
    end
    
    Y_C(1)=Y_N(1)/2;
    for i=2:NY
        Y_C(i)=(Y_N(i)+Y_N(i-1))/2;
    end
    
    Z_C(1)=Z_H(1)/2;
    for i=2:NZ
        Z_C(i)=(Z_H(i)+Z_H(i-1))/2;
    end
    
else
    
    NXC=NX+1;
    NYC=NY+1;
    NZC=NZ+1;
    
    Nrighe=ceil((NXC*NYC)/5);
    Nrighe_x=Nrighe;
    Nrighe_y=Nrighe;
    Nrighe_z=Nrighe;
    
    Parte_2c=Parte_2';
    Parte_2c=Parte_2c(:)';
    
    for i=1:NZC
        xn(i,:)=Parte_2c((1+(Nrighe*5)*3*(i-1)):NXC*NYC+(Nrighe*5)*3*(i-1));
        yn(i,:)=Parte_2c(Nrighe*5+1+(Nrighe*5)*3*(i-1):(Nrighe*5)+NXC*NYC+(Nrighe*5)*3*(i-1));
        zn(i,:)=Parte_2c(Nrighe*5*2+1+(Nrighe*5)*3*(i-1):(Nrighe*5)*2+NXC*NYC+(Nrighe*5)*3*(i-1));
    end
    for i=1:NZC
        for j=1:NXC
            xn3D(j,:,i)=xn(i,1+NYC*(j-1):NYC*j);
            yn3D(j,:,i)=yn(i,1+NYC*(j-1):NYC*j);
            zn3D(j,:,i)=zn(i,1+NYC*(j-1):NYC*j);
        end
    end
    operazione=('scrittura coordinate facce celle')
                                                                 time(6)= toc; 
                                                                tic    
    for k=1:NZC-1
        for i=1:NXC-1
            for j=1:NYC-1
                xc3D(i,j,k)=(xn3D(i,j,k)+xn3D(i+1,j,k))/2;
                yc3D(i,j,k)=(yn3D(i,j,k)+yn3D(i,j+1,k))/2;
                zc3D(i,j,k)=(zn3D(i,j,k)+zn3D(i,j,k+1))/2;
            end
        end
    end
    operazione=('scrittura matrice dei centri delle celle')
end

                                                                time(7)=toc;
                                                                
                                                                tic

                                                                
Nrighe_var=ceil(NX*NY/6);

for i=1:length(Var)
     Var{i}(isspace(Var{i})) = [];
end

if (isequal(variables{1},'all')==0)
    for k=1:length(variables)
        for i=1:length(Var)
            if isequal(Var{i},variables{k})==1
               index(k)=i;
            end
        end
    end
else
    index=[1:length(Var)];
end




for i=index %1:length(Var)
for IZ=1:NZ    
    %Qui sta costruendo un vettore provvisorio con il contenuto della
    %pressione alla SLAB IZ    
    if Nrighe_var==1
        eval(sprintf('a%s = Parte_3((1+M*(IZ-1)+(i-1)*Nrighe_var),1:NX*NY);', Var{i}));
        
    else
        for irighe_var=1:Nrighe_var-1
            eval(sprintf('a%s((1+6*(irighe_var-1)):(6*irighe_var))=Parte_3(irighe_var+Nrighe_var*M*(IZ-1)+(i-1)*Nrighe_var,1:6);', Var{i}));
        end
        eval(sprintf('a%s((1+6*(Nrighe_var-1)):(NX*NY))=Parte_3(Nrighe_var+Nrighe_var*M*(IZ-1)+(i-1)*Nrighe_var,1:((NX*NY)-6*(Nrighe_var-1)));', Var{i}));
    end
    
    %Qui compila la matrice di pressione a partire dal vettore provvisorio.
    
    for IX=1:NX
        for IY=1:NY
            eval(sprintf('%s(IX,IY,IZ)=a%s((IX-1)*NY+IY);', Var{i}, Var{i}));
        end
    end
end
end
operazione=('scrittura matrice delle variabili')

                                                                time(8)=toc;
                                                                tic
for i=1:length(Var)
    eval(sprintf('clear a%s;', Var{i}));
end
clear directory_bfc fnameopen_bfc irighe_y NVariabili txtrh_parte_2 txtrh_parte_2_bfc Parte_2cNrighe xn yn zn i j k U_v V_v W_v p_v Parte_1 Parte_2 Parte_3 Parte_0 Parte_00 Parte_01 Parte_02 txtrh operazione Nrighe_var Nrighe_x Nrighe_y Nrighe_z N_col i irighe_var irighe_x irighe_z M N_col fnameopen

                                                                time(9)=toc;


clear ii time index rowstart 

save(directorySave)
