function add_path
pathPrincipal=pwd;
path(pathPrincipal,path);

pathHelp=strcat(pathPrincipal,'\Help');
path(pathHelp,path);

subfolders={'\HTM','\Images'};
for i=1:size(subfolders,2)
       subPath=strcat(pathHelp,subfolders{i});
       path(subPath,path);
end