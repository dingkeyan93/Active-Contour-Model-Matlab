function band = isband( phi, front, width )%判断是否为窄带点
% CALCBAND Determine which points are within narrow band
%    CALCBAND( phi, front, width ) Based on the indices of the
%    front points, determine which pixels of phi are within a
%    narrow band of width 'width' of the front points. Return a
%    boolean matrix the same size as phi.
%返回一个和phi大小相同的二值矩阵，值为表示是窄带内的点
%front为所有非零点的索引值


% grab size of phi
[ m, n ] = size( phi );

% precompute width squared to save computations later
widthsq = width^2;

% retrieve indices of and total number of front points
n_front = size( front, 1 );%当第二个参数为1时返回行数，即第一维的大小

% create an boolean matrix whose value at each pixel is 0 or 1
% depending on whether that pixel is a band point or not
band = zeros( m, n );

% for each pixel in phi
for ii = 1 : m;
  for jj = 1 : n;

    % check if it is within 'width' of a front pixel
    closest_dist = inf;%inf表示正无穷大
    %判断点(ii,jj)是否在窄带内
    for k = 1 : n_front;

      dist = sum( ( front( k, : ) - [ ii, jj ] ).^2 );%front(k,:)返回(k,:)点的索引值
      if( dist < widthsq )
	band( ii, jj ) = 1;
	break;
      end;

    end;

  end;
end;
