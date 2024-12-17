


function positions = find_positions(str1, str2)
    % 确保两个字符串长度相等
    if length(str1) ~= length(str2)
        error('两个字符串长度不一致');
    end
    
    % 初始化不一致位置的列表
    positions = [];
    
    % 比较字符串中的每个字符
    for i = 1:length(str1)
        if str1(i) ~= str2(i)
            positions = [positions, i];
        end
    end
end

