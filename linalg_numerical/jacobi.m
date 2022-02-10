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

        function [A_dom,b_dom] = diagDom(A_aug)       
            while(1) 
                % till diagonally dominant matrix found
                if jacobi.isDD(A_aug(:,1:end-1))
                    disp (['Matrix A is diagonally-dominant']);
                    break;
                else
                    A_aug = A_aug(randperm(size(A_aug, 1)), :);
                end
            end
            A_dom = A_aug(:,1:end-1);
            b_dom = A_aug(:,end);
        end

        function x = simultaneous(A, b)
            % Conducts Jacobi iterations to get solution to Ax = b
            arguments
                A (:,:)
                b = zeros([size(A,1) 1])
            end

            A_aug = [A b];
            % A_aug = [A([3 1 4 2], :) b([3 1 4 2],:)];

            x = zeros([size(A,1) 1]);
            i = 0;
            while norm(A_aug(:,1:end-1)*x - A_aug(:,end)) > 1e-5
                i = i + 1;
                x = x + (A_aug(:,end) - A_aug(:,1:end-1)*x)./diag(A_aug(:,1:end-1))

                if norm(x) > 1e5
                    disp (['Failed. Trying another row permutation.']);
                    x = zeros([size(A,1) 1]);
                    RIG = randperm(size(A_aug, 1))
                    % RIG = [3 1 4 2];
                    A_aug = A_aug(RIG, :);
                end
                plot(i, norm(A_aug(:,end) - A_aug(:,1:end-1)*x), "-");
                hold on

                A_aug
            end
        end
    end
end