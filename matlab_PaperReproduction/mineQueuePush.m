function [queue,isOk] = mineQueuePush(queue, obj)
	isFull = mineQueueIsFull(queue);
	if isFull == true
		isOk = false;
		return
	end
	if queue.objSizeR == 1
		queue.data(queue.tail, :) = obj;
	elseif queue.objSizeC == 1
		queue.data(:, queue.tail) = obj;
	elseif queue.objSizeR ~= 1 && queue.objSizeC ~= 1
		queue.data(:,:, queue.tail) = obj;
	else
		error("object incorrect");
	end
	queue.tail = mod(queue.tail + 1, queue.len);
    if queue.tail == 0
        queue.tail = queue.len;
    end
	isOk = true;
end