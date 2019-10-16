function subMatrix = get_free_patches_4x(matrix,indices)

% get square subset of a matrix which is composed of 4 identical square matrices,
% given a list of indices to apply along each axis of each subset

% Written to do the matrix-subsetting part of the stress calculation for
% computing creep rates on stress-free patches in 3D

% Eric Lindsey, May 2018

n=length(indices);
part1 = matrix(1:n,1:n);
part2 = matrix(1:n,n+1:2*n);
part3 = matrix(n+1:2*n,1:n);
part4 = matrix(n+1:2*n,n+1:2*n);

sub1=part1(indices,indices);
sub2=part2(indices,indices);
sub3=part3(indices,indices);
sub4=part4(indices,indices);

subMatrix = [sub1, sub2; sub3, sub4];

end
