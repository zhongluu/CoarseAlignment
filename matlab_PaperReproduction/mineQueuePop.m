function [obj,queue,isOk] = mineQueuePop(queue)
	isEmpty = mineQueueIsEmpty(queue);
	if isEmpty == true
		isOk = false;
		obj = [];
	end
	if queue.objSizeR == 1
		obj = queue.data(queue.tail, :);
	elseif queue.objSizeC == 1
		obj = queue.data(:, queue.tail);
	elseif queue.objSizeR ~= 1 && queue.objSizeC ~= 1
		obj = queue.data(:,:, queue.tail);
	else
		isOk = false;
	end
	queue.head = mod(queue.head + 1, queue.len);
    if queue.head == 0
        queue.head = queue.len;
    end
	isOk = true;
end