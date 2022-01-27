classdef gaussian
    methods (Static)
        function A_scaled = scaler(A)
            % NORMALIZES each column using highest value from each row.
            A_scaled = zeros(size(A));
            [nrow, ncol] = size(A);
            for i = 1:nrow
                [~, idx] = max(abs(A(i,:)));
                divisor = A(i,idx);
                A_scaled(i,:) = A(i,:)/divisor;
            end
        end

        function [A_pivoted, b_pivoted] = pivoter(A, b)
            % PIVOTS (row alone) the matrix without scaling

            arguments
                A (:,:)
                b = zeros([size(A,1) 1])
            end

            A_sc = gaussian.scaler(A);
            % two approaches possible to make major diagonal dominant: sort 1st col throughout, then sort 2nd col from just after last 1 in col 1, and so on.
            % the other: find the 1 in each row and use that to index the rows correctly.
            [nrow, ncol] = size(A);
            if ncol ~= nrow
                error("Try inputing a square matrix.", 2, "Does not support pivoting non-square matrices.")
            end

            max_idxs = zeros([nrow 1]);
            for i = 1:nrow
                [~, idx] = max(abs(A_sc(i,:)));
                max_idxs(i) = idx;
            end
            [~, final_idxs] = sort(max_idxs);
            A_pivoted = A(final_idxs, :);
            b_pivoted = b(final_idxs, :);
        end

        function x = jordan_eliminate(A, b)
            % ELIMINATES post pivoting
            [nrow, ncol] = size(A);
            [A_piv, b_piv] = gaussian.pivoter(A, b);

            A_aug = [A_piv b_piv]
            for i = 1:ncol
                for j = 1:nrow
                    if j==i
                        continue;
                    end
                    A_aug(j,:) = A_aug(j,:) - A_aug(j,i)*A_aug(i,:)/A_aug(i,i);
                end
            end

            x = round(A_aug(:, end)./diag(A_aug(:, 1:end-1)), 5);
        end

        function x = gauss_eliminate(A, b)
            % ELIMINATES post pivoting
            [nrow, ncol] = size(A);
            [A_piv, b_piv] = gaussian.pivoter(A, b);

            A_aug = [A_piv b_piv];
            for i = 1:ncol
                for j = i+1:nrow
                    A_aug(j,:) = A_aug(j,:) - A_aug(j,i)*A_aug(i,:)/A_aug(i,i);
                end
            end
            A_aug
            x = zeros(nrow, 1);

            for i=nrow:-1:1
                x(i) = (A_aug(i,end) - A_aug(i,1:end-1)*x)/A(i,i);
            end
        end
    end
end