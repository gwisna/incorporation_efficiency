clear all;
clc;
%opening file
locs1 = h5read('Panel2A_Bima_R1 long normal anneal 1nM_R1 5nM_5_TIRF_50ms_15kframes_70mW_locs_render_picked_avg3_picked.hdf5', '/locs');

%input
n_picks_per_origami=12; %number of picks per origami

n_origami=(max(locs1.group)-min(locs1.group))/n_picks_per_origami; %number of origami being picked
frame=locs1.frame;
x=locs1.x;
y=locs1.y;
photons=locs1.photons;
sx=locs1.sx;
sy=locs1.sy;
bg=locs1.bg;
lpx=locs1.lpx;
lpy=locs1.lpy;
ellipticity=locs1.ellipticity;
net_gradient=locs1.net_gradient;
subgroup=locs1.group-11;

%assigning origami group to subgroup of localizations
for i=1:1:length(subgroup)
    for j=1:1:n_origami
        if mod(subgroup(i),n_origami)==j|mod(subgroup(i),n_origami)==0
            origami_group(i,1)=j;
        end
    end
end

%Counting the number of localization in each subgroup of localizations
for i=1:1:n_origami
    for j=1:1:n_picks_per_origami
        numorigami_numpick(i,j)=0;
        for k=1:1:length(subgroup)
            if subgroup(k)==i+(j-1)*n_origami
                numorigami_numpick(i,j)=numorigami_numpick(i,j)+1;
            end
        end
    end
end

%1D form of numorigami_numpick
localization=numorigami_numpick(:);
localization_sorted=sort(localization);

%determining nbins
ndata=double(n_origami*n_picks_per_origami);
std_dev=std(localization_sorted);
xbins=(1/sqrt(ndata))*std_dev*3.49;
nbins=round((max(localization_sorted)-min(localization_sorted))/xbins);
histo_plot=histogram(localization_sorted, nbins);
edges = histo_plot.BinEdges;
counts = histo_plot.BinCounts;
values = histo_plot.Values;

%determining threshold
pd=fitdist(localization_sorted,'normal');
y_pdf = pdf(pd,localization_sorted);
factor=max(counts)/max(y_pdf);
y=y_pdf*factor;
plot((localization_sorted),y)
threshold_y=max(y)/2;
index_y = find(y >= threshold_y, 1, 'first');
threshold_x=localization_sorted(index_y,1);

%determining detected docking per origami
for i=1:1:n_origami
    for j=1:1:n_picks_per_origami
         if numorigami_numpick(i,j)>=threshold_x
             numorigami_numpick_detected(i,j)=1;
         else
              numorigami_numpick_detected(i,j)=0;
         end
         detected_dockings(i,1)=sum(numorigami_numpick_detected(i,:));
    end
end

%placing into 12 binnings 
for i=1:1:n_picks_per_origami
    detecteddockings_counts(i,1)=i;
    detecteddockings_counts(i,2)=0;
    for j=1:1:n_origami
             if detected_dockings(j,1)==i
             detecteddockings_counts(i,2)=detecteddockings_counts(i,2)+1;
             end
    end
end

histogram(localization_sorted, nbins)
xlabel('Localizations') 
ylabel('Counts')
hold on
plot((localization_sorted),y)

bar(detecteddockings_counts(:,1),detecteddockings_counts(:,2))
xlabel('Detected Dockings') 
ylabel('Counts') 
average_incorporation=mean(detected_dockings)
std_dev_incorporation=std(detected_dockings)
incorporation_eff=average_incorporation/n_picks_per_origami*100
std_dev_incorporation_efficiency=std(detected_dockings)/n_picks_per_origami*100

                
                           

