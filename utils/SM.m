function res = SM(x1, x2)
	id1 = (x1 ~= 0);
	id2 = (x2 ~= 0);
	% res = sum(id1 & id2) / sum(id1 | id2);
    res = 100*sum(~xor(x1 ~= 0, x2 ~= 0))/numel(x1);
% 	res = sum(id1 & id2);% / numel(x1);
end