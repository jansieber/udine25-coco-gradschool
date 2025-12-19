function calledfrom=check_call(name)
st=dbstack(0);
fnames={st.name};
calledfrom=any(strcmp(name,fnames));
end
