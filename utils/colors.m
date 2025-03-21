blue = '#648FFF';
pink = '#DC267F';
orange = '#FE6100';
yellow = '#FFB000';
purple = '#785EF0';
black = '#000000';

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end