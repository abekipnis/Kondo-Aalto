function y=QPI_movie(xdata, dIdV_Matrix, folder, baseName, blx, bly)
    %function that based on given inputs gives you 2 movies, where times
    %axis is bias. The first one is dI/dV map changing with bias and the
    %second one is fft of the dI/dV map changing with bias. It also saves
    %the movies into folder of the .3ds file.
    
    %INPUTS:
    %xdata - bias vector
    %dI_sum - average dI/dV spectrum
    %dIdV_Matrix - matrix of dI/dV values, where the first two dimensions 
    %are spatial, defined by blx and bly, and the third one is bias axis
    %folder - folder of the original .3ds file
    %baseName - name of the original .3ds file without directory
    %blx - spatial coordinates in x direction
    %bly - spatial coordinates in y direction


    s=size(dIdV_Matrix);
    v = VideoWriter([folder, baseName(1:end-4), '_dIdV_movie.avi']);        %starting the video
    v.FrameRate=10;
    open(v)
    figure(11)                                                              %opening given figure
    set(gcf, 'Position',  [200, 200, 1000, 400])                            %setting size of the figure
    
    sorted=sort(reshape(dIdV_Matrix,[s(1)*s(2) s(3)]),1);
    upper=sorted(round(s(1)*s(2)*0.98),:);
    lower=sorted(round(s(1)*s(2)*0.02),:);
    
    dI_sum(:)=mean(mean(dIdV_Matrix,2));
        
    for i=1:s(3)
        subplot(1,2,1);                                                     %creating subfigure in the given figure
        a(1:s(1),1:s(2))=dIdV_Matrix(:,:,i);
        surf(blx, bly, a', 'Edgecolor', 'none');                            %plotting dI/dV map in the subfigure
        axis equal

        caxis([lower(i) upper(i)])
        xlim([min(blx) max(blx)])
        ylim([min(bly) max(bly)])
        view([0 0 1])
        title('dI/dV at given V')

        subplot(1,2,2);                                                     %creating other subfigure in the given figure
        ar=area(xdata,[lower(:) upper(:)], 'Edgecolor', 'none'); ar(1).FaceColor = 'none'; ar(2).FaceColor = [1 0 0]; alpha(.2)
        hold on; grid on
        plot(xdata, dI_sum, 'k', 'Linewidth', 1)                                            %plotting average dI/dV spectrum in the subfigure
        
        plot([xdata(i) xdata(i)], ylim, 'r')
        hold off
        xlim([min(xdata) max(xdata)]);
        title('average spectrum')

        c= getframe(figure(11));
        writeVideo(v,c);                                                    %writing frame into the movie
    end
    close(v)

    v = VideoWriter([folder, baseName(1:end-4), '_QPI_movie.avi']);         %starting the video
    v.FrameRate=10;
    open(v)
    figure(12)                                                              %opening given figure
    set(gcf, 'Position',  [200, 200, 1000, 400])                            %setting size of the figure
    blx_fft=blx/2/pi-median(blx/2/pi);                                      %new x coordinates for the fft
    bly_fft=bly/2/pi-median(bly/2/pi);                                      %new y coordinates for the fft
    for i=1:s(3)
        subplot(1,2,1);                                                     %creating subfigure in the given figure
        a(1:s(1),1:s(2))=dIdV_Matrix(:,:,i);
        aa=fftshift(abs(fft2(a)));
        surf(blx_fft, bly_fft, aa', 'Edgecolor', 'none');                   %plotting fft(dI/dV map) in the subfigure
        axis equal
        
        vec=sort(reshape(aa,1,s(1)*s(2)));        
        caxis([vec(round(s(1)*s(2)*0.92)) vec(round(s(1)*s(2)*0.99))])
        xlim([min(blx_fft) max(blx_fft)])
        ylim([min(bly_fft) max(bly_fft)])
        view([0 0 1])
        colormap(fliplr(gray')')
        title('FFT of dI/dV map at given V')

        subplot(1,2,2);                                                     %creating other subfigure in the given figure
        ar=area(xdata,[lower(:) upper(:)], 'Edgecolor', 'none'); ar(1).FaceColor = 'none'; ar(2).FaceColor = [1 0 0]; alpha(.2)
        hold on; grid on
        plot(xdata, dI_sum, 'k')                                            %plotting average dI/dV spectrum in the subfigure
        
        plot([xdata(i) xdata(i)], ylim, 'r')
        hold off
        xlim([min(xdata) max(xdata)]);
        title('average spectrum')
        
        c= getframe(figure(12));
        writeVideo(v,c);                                                    %writing frame into the movie
    end
    close(v)
end