function plotcell(cell_array)
    figure; 
    hold on;
    len=length(cell_array);
    cmap=hsv(len);
    for i =1:len
        plot(cell_array{i},'Color', cmap(i,:));
    end
end