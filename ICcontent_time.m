ICcont1=[];
ICcont2=[];
for l=1:size(ICArrt,2)
ICtypes = TotAgLeft(ICArrt{2,l});
ICcont1 = [ICcont1 ; ICtypes{1}];
ICcont2 = [ICcont2 ; ICtypes{2}];
end