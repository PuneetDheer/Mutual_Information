%CODED BY: PUNEET DHEER another method for calculation
%UPDATED ON: 11-MAY-2017
%
% INPUT:
% X1= data channel 1
% Y1= data channel 2
% Ws = Window size (in sample point)
% Lp = shift the window by sample point
% ch1id,ch2id= channel id's just to track the process 
%
% OUTPUT:
% Norm_MIxy,Norm_MIyx = Normalized[0 1] mutual information both are same value(property of mutual information Norm_MIxy==Norm_MIyx)
% MIxy,MIyx = raw mutual information both are same value(property of mutual information Norm_MIxy==Norm_MIyx)
% H_e_v = Highest eigen value for common synchronization
% S_esti_MI = S-estimator based synchronization 
%%
function [Norm_MIxy,Norm_MIyx,MIxy,MIyx]= Mutual_Information(X1,Y1,Ws,Lp,ch1id,ch2id)

    % Ws=1048; %WINDOW SIZE
    % Lp=1048; %LEAVE POINT
    Lw=1;
    Z=Ws;

    %tic
    X1=X1;
    Y1=Y1;
    windows=ceil((length(X1)-Ws+1)/Lp);
    ord=3;
    t=1;
    k=2;%for joint entropy normalization
    for i=1:windows
        fprintf('T_Windows: %d_%d of Channel Pairs %d_%d \n',i,windows, ch1id,ch2id);
        X=X1(Lw:Z);
        Y=Y1(Lw:Z);
        Lw=Lw+Lp;
        Z=Z+Lp;

        OPxI = permut(X,ord,t); %feature in signal x
        unixI = unique(OPxI);
        outxI = [unixI;histc(OPxI(:)',unixI)];
        Probx = sum(outxI(2,:));  
        Probx = outxI(2,:)/sum(Probx); 
        ENx = -sum(Probx .* log(Probx)); %MARGINAL ENTROPY OF X
        %ENx=ENx/log(factorial(ord));
         
        OPyI = permut(Y,ord,t);%feature in signal y
        uniyI = unique(OPyI);
        outyI = [uniyI;histc(OPyI(:)',uniyI)];
        Proby = sum(outyI(2,:));  
        Proby = outyI(2,:)/sum(Proby); 
        ENy = -sum(Proby .* log(Proby));  %MARGINAL ENTROPY OF Y
        %ENy=ENy/log(factorial(ord));

        Hxy=[OPxI;OPyI];
        Hyx=[OPyI;OPxI];


        [unique_rows_xy,~,ind]=unique(Hxy','rows');
        Hxycounts = histc(ind,unique(ind));
        % counts = histc(ind,1:max(ind));
        Probxy=sum(Hxycounts);  
        Probxy=Hxycounts/sum(Probxy); 
        ENxy=-sum(Probxy .* log(Probxy)); %JOINT ENTROPY OF XY
        %ENxy=ENxy/(log(factorial(ord))*k);
        

        [unique_rows_yx,~,ind]=unique(Hyx','rows');
        Hyxcounts = histc(ind,unique(ind));
        Probyx=sum(Hyxcounts);  
        Probyx=Hyxcounts/sum(Probyx); 
        ENyx=-sum(Probyx .* log(Probyx)); %JOINT ENTROPY OF YX
        %ENyx=ENyx/(log(factorial(ord))*k);
        

        MIxy(i)=ENx+ENy-ENxy;%RAW
        MIyx(i)=ENx+ENy-ENyx;%RAW

        Norm_MIxy(i)=MIxy(i)/min(ENx,ENy); %[0 1]
        Norm_MIyx(i)=MIyx(i)/min(ENx,ENy); %[0 1]


    %toc
    end
