%figure(poner numero de figura)
fig=figure(3);
set(gca,'visible','off')
filename = 'minaPC_wo_ovlpd.gif'; % Specify the output file name
%Colocar la vista en perspectiva 
%Tambien hacer click derecho y en rotate options seleccionar "fized aspect 
%ratio axes"
delay=0;
%rect=[488,342,560,420];
for i = 1:180
    camorbit(-2,0,'data',[0 0 1])
    drawnow
    frame = getframe(fig);
    temp= frame2im(frame);
    [A,map] = rgb2ind(temp,256);
    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay);
    end
end