function y=QPI_movie(xdata, dI_sum, dIdV_Matrix, folder, baseName, blx, bly)

    s=size(dIdV_Matrix);
    v = VideoWriter([folder, 'Analysis\' ,baseName(1:end-4), '_dIdV_movie.avi']);
    v.FrameRate=10;
    open(v)
    figure(11)
    set(gcf, 'Position',  [100, 100, 1000, 400])
    for i=1:s(3)
        subplot(1,2,1);
        a(1:s(1),1:s(2))=dIdV_Matrix(:,:,i);
        surf(blx, bly, a, 'Edgecolor', 'none');
        vec=sort(reshape(a,1,s(1)*s(2)));
        
        caxis([vec(round(s(1)*s(2)*0.05)) vec(round(s(1)*s(2)*0.95))])
        xlim([blx(1) blx(end)])
        ylim([bly(1) bly(end)])
        view([0 0 1])
        title('dI/dV at given V')

        subplot(1,2,2);
        plot(xdata, dI_sum, 'k')
        hold on
        ax=gca;
        plot([xdata(i) xdata(i)], [ax.YLim(1) ax.YLim(2)], 'r')
        hold off
        xlim([min(xdata) max(xdata)]);
%         ylim([0 2.5])
        title('average spectra')

        c= getframe(figure(11));
        writeVideo(v,c);
    end
    close(v)

    v = VideoWriter([folder, 'Analysis\', baseName(1:end-4), '_QPI_movie.avi']);
    v.FrameRate=10;
    open(v)
    figure(12)
    set(gcf, 'Position',  [100, 100, 1000, 400])
    blx_fft=blx/2/pi-median(blx/2/pi);
    bly_fft=bly/2/pi-median(bly/2/pi);
    for i=1:s(3)
        subplot(1,2,1);
        a(1:s(1),1:s(2))=dIdV_Matrix(:,:,i);
        aa=fftshift(abs(fft2(a)));
        surf(blx_fft, bly_fft, aa, 'Edgecolor', 'none');
        vec=sort(reshape(aa,1,s(1)*s(2)));
        
        caxis([vec(round(s(1)*s(2)*0.75)) vec(round(s(1)*s(2)*0.98))])
        xlim([blx_fft(1) blx_fft(end)])
        ylim([bly_fft(1) bly_fft(end)])
        view([0 0 1])
        title('fft of dI/dV')
        subplot(1,2,2);
        plot(xdata, dI_sum, 'k')
        hold on
        ax=gca;
        plot([xdata(i) xdata(i)], [ax.YLim(1) ax.YLim(2)], 'r')
        hold off
        xlim([min(xdata) max(xdata)]);
%         ylim([0 2.5])
        title('average spectra')
        
        c= getframe(figure(12));
        writeVideo(v,c);
    end
    close(v)
end