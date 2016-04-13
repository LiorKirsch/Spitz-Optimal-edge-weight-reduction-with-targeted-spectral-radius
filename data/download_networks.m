konect_meta = textscan(fopen('konect_meta.csv'),'%q %q %q %q %f %f %f %f %f %f','delimiter',',','Headerlines',1);
konect_url = konect_meta{1};
konect_code = konect_meta{2};
konect_name = konect_meta{3};
konect_type = konect_meta{5};
konect_weight = konect_meta{6};
konect_vertices = konect_meta{8};
konect_edges = konect_meta{9};

konect_types_dict = {'undirected','directed','bipart'};
konect_type = konect_types_dict(konect_type);

if ~exist('net_mat','dir')
    mkdir('net_mat');
end

networks_with_error = false(length(konect_url),1);
for i =length(konect_url):-1:1
   network_name = konect_name{i};
 try
       if  strcmp(konect_url{i} ,'')
           fprintf('skipping %s\n', network_name);
       else
           filename = strsplit(konect_url{i},'/');
           filename = filename{end};

           just_filename = strsplit(filename,'.');
           just_filename = just_filename{1};

           matfile = sprintf('net_mat/%s.mat',just_filename);
           if exist(matfile,'file')
               fprintf('skipping, matfile exists %s\n', matfile);
           else
               if exist(just_filename,'dir')
                   fprintf('skipping directory exists %s\n', just_filename);
               else
                   fprintf('Downloading\t %s\n', network_name);
                   system(sprintf('wget %s', konect_url{i}));

                   fprintf('Extracting\t %s\n', network_name);
                   system(sprintf('tar xjf %s',filename) );
                   system(sprintf('rm %s',filename) );
               end

               net_file = dir(sprintf('%s/out.*',just_filename));
               fprintf('Parsing\t  %s\n', network_name);
               net_file = fopen(fullfile(just_filename,net_file.name));
               net_data = textscan(net_file,'%d %d %f %*[^\n]','delimiter', ' ','CommentStyle', '%');
               net_rowind  = net_data{1};
               net_colind  = net_data{2};
               net_weight  = net_data{3};
               fclose(net_file);
	       clear('net_data');

               assert(length(net_rowind) == konect_edges(i),'the parsed #edges is not the same as in the metadata');

               if isnan(net_weight(1))
                   A = sparse(double(net_rowind), double(net_colind), ones(size(net_colind)), double(konect_vertices(i)), double(konect_vertices(i)));
               else
                   A = sparse(double(net_rowind), double(net_colind), net_weight, double(konect_vertices(i)), double(konect_vertices(i)));
               end

               if strcmp(konect_type{i},'undirected')
                    A = (A +A')/2;
                    A = A - diag(diag(A));
               end
               fprintf('Saving\t %s\n', network_name);
               save(matfile, 'A', 'network_name', '-v7.3');

               system(sprintf('rm -r %s/',just_filename) );
           end
       end
  catch
      networks_with_error(i) = true;
      fprintf('something wrong with %s\n', network_name)
  end
%    pause;
end
