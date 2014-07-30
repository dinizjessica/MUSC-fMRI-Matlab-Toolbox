function BC = betweenness_wei_helper(conn_matrix)

L = weight_conversion(conn_matrix, 'lengths');
BC = betweenness_wei(L);
            
end