function H = snapshot(filename, column1, column2);
A = textread(filename);
A = A(4:end,:);
X = A(:,1);
Y = A(:,2);
V = A(:,3);
W = A(:,4);
Z = A(:,column1) - A(:,column2);
numpts = size(Z);
ptsPerSide = sqrt(numpts(1));
x1 = linspace(min(X),max(X),ptsPerSide);
x2 = linspace(min(Y),max(Y),ptsPerSide);

hold off;
%if numpts <= 2400
%  plot3(X,Y,Z,'.'), hold; 
%end
axis tight; 
fx = reshape(Z,ptsPerSide,ptsPerSide);
H = surf(x1,x2,fx,'FaceColor','interp');
view(157,25)
%title(filename,'FontWeight','bold','FontSize',18);
[C,I] = max(Z);
X(I);
Y(I);
Z(I);
