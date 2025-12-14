function P=mmult(varargin)
P=varargin{1};
for i=2:2:length(varargin)
    interface=varargin{i};
    A=varargin{i+1};
    Psz=mshape(P,'size');
    Asz=mshape(A,'size');
    if iscell(interface)
        Aother=setdiff(1:length(Asz),interface{2});
        A=mshape(A,'p',[interface{2},Aother]);
        Asz=mshape(A,'size');
        if length(interface)>2
            nvecdim=interface{3};
        else
            nvecdim=0;
        end
        interface=interface{1};
    end
    nint=length(interface);
    [np,na]=deal(length(Psz)-nvecdim,length(Asz)-nvecdim);
    [pbefore,pafter,pvec]=deal((1:interface(1)-1),interface(end)+1:np,np+(1:nvecdim));
    Pparts=cellfun(@(s)prod(Psz(s)),{pbefore,pafter,interface,np+(1:nvecdim)});
    Aparts=cellfun(@(s)prod(Asz(s)),{1:nint,nint+1:na,na+(1:nvecdim)});
    P3=mshape(P,...
        'p',[pbefore,pafter,interface,pvec],...
        'r',{Pparts(1)*Pparts(2),Pparts(3),Pparts(4)});
    A3=mshape(A,'r',num2cell(Aparts));
    if nvecdim>0
        P3=mshape(P3,'blkdiag');
        A3=mshape(A3,'blkdiag');
    end
    P3=mshape(P3,'mat');
    A3=mshape(A3,'mat');
    posA=np-nint+(1:na-nint);
    posPafter=length(pbefore)+(1:length(pafter));
    posvec=np+na-2*nint+(1:nvecdim);
    P=mshape(P3*A3,'diagblk',[Pparts(1)*Pparts(2),Aparts(2),Pparts(4)],...
        'r',[num2cell(Psz([pbefore,pafter])),...
        num2cell(Asz(nint+1:na)),num2cell(Psz(pvec))],...
        'p',[pbefore,posA,posPafter,posvec],'struct');
end
end
