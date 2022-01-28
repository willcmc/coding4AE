classdef gaussian
    methods (Static)
        function [A_scaled, b_scaled] = scaler(A, b)
            % NORMALIZES each column using highest value from each row.

            arguments
                A (:,:)
                b = zeros([size(A,1) 1])
            end

            A_scaled = zeros(size(A));
            [nrow, ncol] = size(A);

            % finding largest in row and normalizing row by it
            for i = 1:nrow
                [~, idx] = max(abs(A(i,:)));
                divisor = A(i,idx);
                A_scaled(i,:) = A(i,:)/divisor;
                b(i) = b(i)/divisor;
            end

            b_scaled = b;
        end

        function [A_pivoted, b_pivoted] = pivoter(A, b)
            % PIVOTS (row alone) the matrix without scaling, based on scaled values

            arguments
                A (:,:)
                b = zeros([size(A,1) 1])
            end

            [A_sc, b_sc] = gaussian.scaler(A,b);

            [nrow, ncol] = size(A);
            if ncol ~= nrow
                error("Try inputing a square matrix.", 2, "Does not support pivoting non-square matrices.")
            end

            A_here = [A b];

            max_idxs = zeros([ncol 1]);
            for i = 1:ncol
                [~, idx] = max(abs(A_sc(:, i))); % finding index of max abs of each column

                if abs(A_sc(i, i)) < abs(A_sc(idx, i)) % swapping rows if more dominant diagonal possibility found
                    temp = A_sc(i,:);
                    A_sc(i,:) = A_sc(idx,:);
                    A_sc(idx,:) = temp;

                    temp = A_here(i,:);
                    A_here(i,:) = A_here(idx,:);
                    A_here(idx,:) = temp;
                end
            end

            A_pivoted = A_here(:,1:end-1);
            b_pivoted = A_here(:,end);
        end

        function x = jordan_eliminate(A, b)
            % ELIMINATES post pivoting
            [nrow, ncol] = size(A);
            [A_piv, b_piv] = gaussian.pivoter(A, b);

            A_aug = [A_piv b_piv];
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

        function A_inv = jordan_inverse(A)
            % INVERTS post pivoting
            [nrow, ncol] = size(A);
            A_inv = eye(nrow);

            % gives us pivoted matrix
            A_piv = [gaussian.pivoter(A) A_inv];
            
            for i = 1:ncol
                for j = 1:nrow
                    if j==i
                        continue;
                    end
                    % eliminating every off-diagonal term and performing same operations on RHS (I)
                    A_piv(j,:) = A_piv(j,:) - A_piv(j,i)*A_piv(i,:)/A_piv(i,i);
                end
            end
            
            A_piv = A_piv./diag(A_piv(:,1:ncol));
            A_inv = A_piv(:,ncol+1:end);            


            % Restoring column order as original matrix
            rightidxs = zeros([nrow,1]);
            shouldbeI = A*A_inv;
            for i=1:ncol
                [~, idx] = max(shouldbeI(i,:));
                rightidxs(i) = idx
            end

            A_inv = A_inv(:,rightidxs)

        end

        function x = gauss_eliminate(A, b)
            % ELIMINATES post pivoting
            [nrow, ncol] = size(A);
            [A_piv, b_piv] = gaussian.pivoter(A, b);

            A_aug = [A_piv b_piv];
            for i = 1:ncol
                for j = i+1:nrow
                    % Forming upper triangular by eliminating all elements under pivots
                    A_aug(j,:) = A_aug(j,:) - A_aug(j,i)*A_aug(i,:)/A_aug(i,i);
                end
            end
            % A_aug
            x = zeros(nrow, 1);

            for i=nrow:-1:1
                x(i) = (A_aug(i,end) - A_aug(i,1:end-1)*x)/A_aug(i,i);
            end
        end
    end
end