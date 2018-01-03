fils=dir('*.mat');
t=0;
for i=1:length(fils)
    aa=load(fils(i).name);
    t=t+aa.time;
end
t/length(fils)