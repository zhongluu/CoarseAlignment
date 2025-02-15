function queue = mineQueue(obj, num)
	queue.head = 1;
	queue.tail = 1;
	queue.len = num + 1;
	[r,c] = size(obj);
	queue.objSizeR = r;
	queue.objSizeC = c;
	if r == 1
		queue.data = repmat(obj, queue.len, 1);
	elseif c == 1
		queue.data = repmat(obj, 1, queue.len);
	elseif r~=1&&c~=1
		queue.data = repmat(obj,1,1,queue.len);
	else
		error("can not support this object");
	end
end