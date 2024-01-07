function [index,xy] = fxShowFEM(Element, Node, Sigma, N)
patch('Faces',Element,'Vertices' ,Node,'FaceVertexCData' ,Sigma,'FaceColor' ,'flat' ,'EdgeColor' ,'None' );
axis off image;
% text(Node(3256,1),Node(3256,2),'1');
% text(Node(875,1),Node(875,2),'5');
% colormap(pink);
% hold on;
for i = 1:length(Element)
    xy(i,:) = mean(Node(Element(i,:),:));
%     text(xy(i,1),xy(i,2),num2str(i));
end


if nargin > 3
    [gx,gy] = ginput(N);
    for i = 1:N
        temp = [gx(i),gy(i)];
        temp2 = xy-repmat(temp,length(Element),1);
        temp3 = temp2(:,1).^2+temp2(:,2).^2
        [~,index(i)] = min(temp3);
    end
end
        