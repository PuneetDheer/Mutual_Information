function OPi = permut(S,ord,t)
    
    ly = length(S);
    permlist = perms(1:3);
    c(1:length(permlist))=0;
    OPi=zeros(1,ly-t*(ord-1));

     for j=1:ly-t*(ord-1)
         [a,iv]=sort(S(j:t:j+t*(ord-1)));
         for jj=1:length(permlist)
             if (abs(permlist(jj,:)-iv))==0
                 OPi(j)=jj;
             end
         end
     end 