classdef jacobi
    methods (Static)
        function isdom = isDD(A) 
            isdom = true;
            for r = 1:size(A,1)
                rowdom = 2 * abs(A(r,r)) >= sum(abs(A(r,:)));
                isdom = isdom && rowdom;
            end
            if isdom == 0
                disp (['Matrix A is not diagonally-dominant']);
            elseif isdom == 1
                    disp (['Matrix A is diagonally-dominant']); 
            end
        end

        function x = simultaneous(A, b)
            % Conducts Jacobi iterations to get solution to Ax = b
            arguments
                A (:,:)
                b = zeros([size(A,1) 1])
            end

            [nrow, ncol] = size(A);
            if ncol ~= nrow
                error("Try inputing a square matrix.", 2, "Does not support pivoting non-square matrices.")
            end

            A_aug = [A b];
            % A_aug = [A([3 1 4 2], :) b([3 1 4 2],:)];

            x = zeros([size(A,1) 1]);
            i = 0;
            while norm(A_aug(:,1:end-1)*x - A_aug(:,end)) > 1e-5
                i = i + 1;
                x = x + (A_aug(:,end) - A_aug(:,1:end-1)*x)./diag(A_aug(:,1:end-1))

                if norm(x) > 1e5 % if solution diverges, permute the rows randomly. 
                    disp (['Failed. Trying another row permutation.']);
                    x = zeros([size(A,1) 1]);
                    RIG = randperm(size(A_aug, 1))
                    A_aug = A_aug(RIG, :);
                end
                plot(i, norm(A_aug(:,end) - A_aug(:,1:end-1)*x), "-");
                hold on
            end
        end

        function x = seidel(A, b)
            % Conducts Seidel iterations to get solution to Ax = b
            [nrow, ncol] = size(A);
            if ncol ~= nrow
                error("Try inputing a square matrix.")
            end

            x = zeros([nrow 1]);
            
            A_aug = [A b];

            while norm(A_aug(:,1:end-1)*x - A_aug(:,end)) > 1e-5
                for i = 1:nrow
                    x(i) = x(i) + (A_aug(i,end) - A_aug(i, 1:end-1)*x)/A_aug(i,i) % updates the values of x row-wise

                    if norm(x) > 1e2 % if solution diverges, permute the rows randomly. 
                        disp (['Failed. Trying another row permutation.']);
                        x = zeros([nrow 1]);
                        RIG = randperm(size(A_aug, 1))
                        A_aug = A_aug(RIG, :);
                        eig(A_aug(:, 1:end-1))
                    end
                end
            end
        end

        function x = SOR(A, b, w)
            % Conducts Jacobi iterations to get solution to Ax = b
            arguments
                A
                b
                w = 1.2
            end

            [nrow, ncol] = size(A);
            if ncol ~= nrow
                error("Try inputing a square matrix.")
            end

            x = zeros([nrow 1]);
            
            A_aug = [A b];

            while norm(A_aug(:,1:end-1)*x - A_aug(:,end)) > 1e-5
                for i = 1:nrow
                    x(i) = x(i) + w*(A_aug(i,end) - A_aug(i, 1:end-1)*x)/A_aug(i,i)

                    if norm(x) > 1e2 % if solution diverges, permute the rows randomly. 
                        disp (['Failed. Trying another row permutation.']);
                        x = zeros([nrow 1]);
                        RIG = randperm(size(A_aug, 1))
                        A_aug = A_aug(RIG, :);
                        eig(A_aug(:, 1:end-1))
                    end
                end
            end
        end
    end
end