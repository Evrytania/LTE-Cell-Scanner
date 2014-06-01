function llr = pb2llr(pbit1, pbit0, adj)
llr = log((pbit1+adj)./(pbit0+adj));
