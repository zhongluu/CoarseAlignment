function isEmpty = mineQueueIsEmpty(Queue)
	if Queue.tail == Queue.head
		isEmpty = true;
	else
		isEmpty = false;
	end
end