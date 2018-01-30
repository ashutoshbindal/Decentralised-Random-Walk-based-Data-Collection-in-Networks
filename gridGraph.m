N=10;
N2=N^2;
A=zeros(N2,N2);
%{
for i=1:N2
    for j=i:N2
        if (i~=j && abs(i-j)<=N) && ((mod(abs(i-j),N)==1 && floor((i-1)/N)==floor((j-1)/N)) || mod(abs(i-j),N)==0)
            A(i,j)=1;
            A(j,i)=1;
        end
    end
end
A(1,N)=1;
A(N,1)=1;
A(1,(N-1)*N+1)=1;
A((N-1)*N+1,1)=1;

A((N-1)*N+1,N2)=1;
A(N2,(N-1)*N+1)=1;
A(N,N2)=1;
A(N2,N)=1;
%}

ll=0;
rl=0;
ul=0;
dl=N^2;
for rows=1:N
    for cols=1:N
      no=(rows-1)*N+cols;
      %get neighbors
      ln=no-1;
      rn=no+1;
      un=no-N;
      dn=no+N;
      %limits
      ll=(rows-1)*N;
      rl=ll+N;
      ll=ll+1;
      %check limits
      if (ln<ll)
          ln=rl;
      end
      if(rn>rl)
          rn=ll;
      end
      if(un<=ul)
          un=no+(N-1)*N;
      end
      if(dn>dl)
          dn=no-(N-1)*N;
      end
      A(no,[ln rn un dn])=1;
    end
end


fileID = fopen('deg5_1000.txt','w');
fprintf(fileID,'%d\n',A);
fclose(fileID);