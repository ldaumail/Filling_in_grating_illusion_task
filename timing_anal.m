data = load('/Users/tonglab/Desktop/Loic/data/M017/230208/phantom_v1/sM017_phantom_sn1_rn1_date20230208T113329.mat'); 

data.ex.fixFlipT(1:168,:) = data.ex.stimFlipT;
for i =1:size(data.ex.flipTime,2)
    %diff1 = data.ex.flipTime;
    %diff(:,i) = data.ex.allFlipTimes'-data.ex.flipTime(:,i);
    diff2(:,i) = diff(data.ex.fixFlipT(:,i));
    
end

for i = 2:size(data.ex.flipTime,2)
    diff3(i) = data.ex.fixFlipT(1,i)- data.ex.fixFlipT(end,i-1);
end

data.ex.runTime

%% Task the waits for backtick at every block

data = load('/Users/tonglab/Desktop/Loic/data/Dave/Dave/session1/phantom_v1/sDave_phantom_sn1_rn1_date20230301T141627.mat');
%data.ex.fixFlipT(1:168,:) = data.ex.stimFlipT;

for i =1:size(data.ex.flipTime,2)
    %diff1 = data.ex.flipTime;
    diff1(i) = data.ex.stimFlipT(end,i)-data.ex.stimFlipT(1,i);
   
    diff2(i) = data.ex.fixFlipT(end,i)-data.ex.fixFlipT(169,i);
    
end

tr1FlipTs = diff(data.ex.fixFlipT(:,1));
tr3FlipTs = diff(data.ex.fixFlipT(:,3));


for i =1:size(ex.flipTime,2)
    %diff1 = data.ex.flipTime;
    diff1(i) = ex.stimFlipT(end,i)-ex.stimFlipT(1,i);
   
    diff2(i) = ex.fixFlipT(end,i)-ex.fixFlipT(169,i);
    
end
mean((ex.blockOnsetTime(2:end)-ex.backtickT(1:end-1)')*1000)
