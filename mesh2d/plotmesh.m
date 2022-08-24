function []=plotmesh(the_coor,connect,rgb,bck)

        
        patch('faces',connect(:,:),'vertices',the_coor, ...
                            'facecolor',bck, ...
                        'edgecolor',rgb) ;
                        hold on; 
                        
                        axis image off;
                        
end
                        