      function [ind] = index( n, j)
            % index(nx,j) returns the ind, such that s(ind) corresponds
            % to the jth stage vector
           ind = 1 + (j-1) * n : j * n; 
      end