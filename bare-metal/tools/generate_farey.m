function farey=generate_farey(n)
%% generate Farey sequence up to denominator n
a=[0;1;1;n];
farey=NaN(2,n^2);
farey(:,1)=a(1:2);
ind=1;
while a(3)<=n
    k=floor((n+a(2))/a(4));
    a=[a(3);a(4);k*a(3)-a(1);k*a(4)-a(2)];
    ind=ind+1;
    farey(:,ind)=a(1:2);
end
farey=farey(:,1:ind);
end
%% From Wikipedia
% def farey( n, asc=True ):
%     """Python function to print the nth Farey sequence, either ascending or descending."""
%     if asc: 
%         a, b, c, d = 0, 1,  1  , n     # (*)
%     else:
%         a, b, c, d = 1, 1, n-1 , n     # (*)
%     print "%d/%d" % (a,b)
%     while (asc and c <= n) or (not asc and a > 0):
%         k = int((n + b)/d)
%         a, b, c, d = c, d, k*c - a, k*d - b
%         print "%d/%d" % (a,b)
