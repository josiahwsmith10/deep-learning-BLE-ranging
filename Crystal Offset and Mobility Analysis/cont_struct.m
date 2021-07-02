% This function will concatenate structure fields

function struct_out = cont_struct(struct_in)

for i = 1:size(struct_in,2)
    
    % Make sure to do initialization
    if i == 1
        struct_vrs = fieldnames(struct_in); % Fieldnames
        struct_out = struct_in(1); % First set of variable values
    else
        for k = 1:length(struct_vrs)
            struct_out.(struct_vrs{k}) = [struct_out.(struct_vrs{k});struct_in(i).(struct_vrs{k})];
        end
    end

end