function isFull = mineQueueIsFull(queue)
	if mod(queue.tail + 1, queue.len) == queue.head || queue.tail + 1 == queue.head
		isFull = true;
	else
		isFull = false;
	end
end