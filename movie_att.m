%Written by Rachith Aiyappa in the summer of 2018 at CRCA, Toulouse, France
%Email : rachithaiyappa96@gmail.com for any queries.

tic
%myVideo = VideoWriter('H:\rachith_jitesh_spatial\blah');
myVideo = VideoWriter('blah');
myVideo.FrameRate = 30;
open(myVideo);
for i = 2:1:size(samp_time,2)
    i
%for i = 2:1:size(time_temp,2)
    a = 1;
    hold off
    
%     while a <= N
%             %plot([x(a,i),x(a,i-1)],[y(a,i),y(a,i-1)],'k-','markersize',4)
%             plot(200*cos(phi(a,i)),200*sin(phi(a,i)))
%             axis([0 L 0 L])
%             %axis([-50 10 -50 80])
%             hold on
%             %plot(x(a,i-1),y(a,i-1),'.','markersize',10);
%             plot(x_eq(a,i-1),y_eq(a,i-1),'.','markersize',10);
%             a = a + 1;
%     end
    while a<=N
        plot(cos(phi(a,i)),sin(phi(a,i)))
        axis([0 L 0 L])
        hold on 
        quiver(x_eq(a,i-1),y_eq(a,i-1),cos(phi_eq(a,i-1)),sin(phi_eq(a,i-1)),1.5,'MaxHeadSize',1.5);
        a = a + 1;
    end
    M(i) = getframe(gcf);
    writeVideo(myVideo, M(i));
end
close(myVideo);
toc