
figure;

AVG_SURF = SurfStatReadSurf({'AVGmid_LEFT.obj','AVGmid_RIGHT.obj'});
load AAL_atlas_81924.txt;
load cth.txt;

Value=zeros(size(AAL_atlas_81924));

for i=1:81924, 
    for j=1:79, 
        if AAL_atlas_81924(i)==cth(j), 
            Value(i)=cth(j,2); 
        end; 
    end; 
end;

SurfStatView(Value, AVG_SURF, 'Cortical Surface Area');

SurfStatColormap('hot');

SurfStatColLim([0.51 1]);
