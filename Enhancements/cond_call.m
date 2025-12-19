function y=cond_call(f,args,nf,name)
st=dbstack(0);
fnames={st.name};
if any(strcmp(name,fnames))
    y=f(args{:});
else 
    y=NaN(nf);
end
end
