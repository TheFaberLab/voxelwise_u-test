% extracts the BOLD response to a stimulation in Block design
% determines activated voxels
% tests whether the signal of an active voxel is positive
% averages the signal from the activated voxels with a positive signal
% change according to the stimulation paradigm
% for TR=1s 

function voxelwise_u_test(nr_slices, nr_scans_in_loop, nr_mrscans, nr_matrix, t_before, t_stim, t_shift, alpha, images, folder)


%nr_slices         - number of slices of the MR scans
%
%nr_scans_in_loop - number of scans in one iteration of the stimulation cycle; program works only, if this value is bigger than 5
%
%nr_mrscans       - number of MRI scans in one measurement
%
%nr_matrix        - matrix size (e.g 80x80 voxels)
%
%t_before         - time before the stimulation (in seconds)
%
%t_stim           - duration of the stimulation (in seconds)
%
%t_shift          - time shift of the paradigm to take into account the delay of the BOLD response.
%
%alpha            - From this p value, voxels are considered significant.
%                   This value is divided by the number of voxels in the selected ROI. (Bonferroni correction)
%
%images           - rs or srs data (resulting from DICOM import of SPM)
%                   imported like
%                        for k=1:nr_mrscans
%                            nr_akt_scan=1+(k-1)*nr_slices;
%                            images(k,:,:,:)=flipdim(analyze75read([path_mrscans '\' filename num2str(nr_akt_scan, '%05g') '-' num2str(nr_akt_scan, '%06g')]),1);
%                        end
%folder             - folder where the workspace will be stored

%__________________________________________________________________________
% Henriette Lambers, 2019

nr_loop=ceil(nr_mrscans/nr_scans_in_loop);
t_after=nr_scans_in_loop-t_before-t_stim;


% slicewise ROI
ROI_mask_ttest = zeros(nr_slices, nr_matrix, nr_matrix);
nr_pixel=0;
sum_tsnr_voxel = 0;
bigger=0;
smaller=0;

for slice=1:nr_slices
    
    P=squeeze(double(images(:,:,:,slice)));
    colormap(gray)
    imagesc(squeeze(P(1,:,:)));
    title(['draw a ROI in slice ' num2str(slice)], 'FontSize', 16);
    [ROI, x, y] = roipoly();
    ROI_mask(slice,:,:) = double(ROI);   
    hold on;
    plot(x, y, 'r.-', 'MarkerSize', 15);
    
    

    %%%%determine the number of pixels in the ROI & timecourse of the total ROI%%%%%%%%%
    pixelanzahl_all(slice) = sum(sum(ROI_mask(slice,:,:)));
    nr_pixel = nr_pixel + pixelanzahl_all(slice);
    
    if pixelanzahl_all(slice)>0
        for t=1:size(P,1)
        timecourse_all_sl(t,slice)=sum(sum(squeeze(P(t,:,:)).*squeeze(ROI_mask(slice,:,:))));
        end
        
        clearvars t
    end
end
for slice=1:nr_slices
    P=squeeze(double(images(:,:,:,slice)));
    for index_x=1:nr_matrix
        for index_y=1:nr_matrix
            if ROI_mask(slice,index_x,index_y)>0
                %%%%%%%%%%%%time course of the observed pixel %%%%%%%%%
                for t=1:size(P,1)
                    timecourse_pixel(t)=P(t,index_x,index_y); %time course of a pixel,within the marked ROI
                end
                mean_timec_pixel_value(index_x,index_y,slice)=mean(timecourse_pixel);
                clearvars t
                
                %%%%%%%%%%%%  signal while stimulation %%%%%%%%%
                clearvars timecourse_on
                for n=1:nr_loop
                    timecourse_on(:,n)=timecourse_pixel(t_shift+t_before+(n-1)*nr_scans_in_loop+1:t_shift+t_before+(n-1)*nr_scans_in_loop+t_stim);
                end 
                    timecourse_on=timecourse_on(:);
                 clearvars n
                %%%%%%%%%%%% signal while rest %%%%%%%%%%%%%%%%
                    timecourse_out_1=timecourse_pixel(1:t_shift+t_before)';
                    
                    clearvars timecourse_out_2
                for n=1:nr_loop-1
                    timecourse_out_2(:,n)=timecourse_pixel(t_shift+t_before+t_stim+(n-1)*nr_scans_in_loop+1:t_shift+t_before+t_stim+(n-1)*nr_scans_in_loop+t_after+t_before);
                end
                    timecourse_out_2=timecourse_out_2(:);
                    timecourse_out_3=timecourse_pixel((nr_loop-1)*nr_scans_in_loop+t_shift+t_before+t_stim+1:end)';
                    timecourse_out=[timecourse_out_1;timecourse_out_2;timecourse_out_3];
                    clearvars n
                
                
               [p_utest, h_utest] = ranksum(timecourse_on,timecourse_out,'Alpha',alpha/nr_pixel);
                p_utest_value(index_x,index_y,slice)=p_utest;
                
                
                ROI_mask_ttest(slice,index_x,index_y)=ROI_mask(slice,index_x,index_y);
                if h_utest==0
                   ROI_mask_ttest(slice,index_x,index_y)=0;
                else
                   mean_on=mean(timecourse_on);
                   mean_out=mean(timecourse_out);
                    if mean_on > mean_out
                        bigger=bigger+1;
                    else
                        smaller=smaller+1;
                        ROI_mask_ttest(slice,index_x,index_y)=0;
                    end
                end
            end   
        end
    end
    pixelanzahl_ttest(slice) = sum(sum(ROI_mask_ttest(slice,:,:)));
    if pixelanzahl_ttest(slice)>0
        for t=1:size(P,1)
        timecourse_ttest_sl(t,slice)=sum(sum(squeeze(P(t,:,:)).*squeeze(ROI_mask_ttest(slice,:,:))));
        end
        clearvars t
    end
    
end

%%%%%%%%%%%%%% determine the averaged time course for the active voxels
timecourse_ttest=sum(timecourse_ttest_sl,2);
timecourse_ttest_norm=timecourse_ttest/mean(timecourse_ttest)-1;


%%%%%%%%%%%%%% average the time course according to the stimulation
length_lastline=length(timecourse_ttest((nr_loop-1)*nr_scans_in_loop+1:end));
timecourse_ttest_av=zeros(nr_scans_in_loop,1);

if length_lastline==nr_scans_in_loop
    for k=1:nr_loop
    timecourse_ttest_av=timecourse_ttest_av+timecourse_ttest_norm((k-1)*nr_scans_in_loop+1:k*nr_scans_in_loop);
    end
    timecourse_ttest_av=timecourse_ttest_av/nr_loop;
else
    for k=1:nr_loop-1
    timecourse_ttest_av=timecourse_ttest_av+timecourse_ttest_norm((k-1)*nr_scans_in_loop+1:k*nr_scans_in_loop);
    end
    timecourse_ttest_av(1:length_lastline)=timecourse_ttest_av(1:length_lastline)+timecourse_ttest_norm(k*nr_scans_in_loop+1:end);
    timecourse_ttest_av(1:length_lastline)=timecourse_ttest_av(1:length_lastline)/nr_loop;
    timecourse_ttest_av(1+length_lastline:end)=timecourse_ttest_av(1+length_lastline:end)/(nr_loop-1);
end


%images of the averaged time courses
figure
plot(100*timecourse_ttest_av,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
XTicks = 0:5:nr_scans_in_loop;
set(gca,'box','off')
set(gca, 'TickDir', 'out')
set(gca, 'XTickMode', 'manual', 'XTick', XTicks, 'xlim', [0,nr_scans_in_loop],'FontSize',15);
set(gca, 'YTickMode', 'auto','FontSize',15);
xlabel('time / s','FontSize',15);
ylabel('BOLD-signal / %','FontSize',15)

path=strcat(folder, '\workspace_utest');
save(path)
