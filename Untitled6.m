figure
hold
mymatrix = C-A;
grade=[0.07, 0.1, 0.2, 0.3, 0.7];
 colo=['b','g','y','m','r'];
 m=max(max(max(mymatrix)));
 lev=[m*1/20, m*1/8, m*1/4, m/2, m*9/10];
 for i=1:5
     [F,V]=isosurface(mymatrix,lev(i));
     patch('Faces',F,'Vertices',V,'FaceColor',colo(i),'FaceAlpha',grade(i),'EdgeAlpha',grade(i));
 end
 view(3)
 
%  grade=[0.07, 0.1, 0.2, 0.3, 0.7];
%  colo=['b','g','y','m','r'];
%  lev=[0, 1.1, 2.5, 0.9, 1];
%  for i=1:3
%      [F,V]=isosurface(phantom,lev(i));
%      patch('Faces',F,'Vertices',V,'FaceColor',colo(i),'FaceAlpha',grade(i),'EdgeAlpha',grade(i));
%  end
%  view(3)