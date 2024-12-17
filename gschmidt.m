function Q = gschmidt(A)
    [m, n] = size(A);
    Q = zeros(m, n);  % 存储正交化后的向量
    
    for i = 1:n
        Q(:,i) = A(:,i);  % 将原始向量集合的第i个向量赋值给Q
        
        for j = 1:i-1
            Q(:,i) = Q(:,i) - (Q(:,j)' * A(:,i)) / (Q(:,j)' * Q(:,j)) * Q(:,j);  % 正交化计算
        end
        
        % 单位化
        Q(:,i) = Q(:,i) / norm(Q(:,i));
    end
end
