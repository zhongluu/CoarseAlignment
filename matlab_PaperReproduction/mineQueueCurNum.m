function curNum = mineQueueCurNum(queue)
	curNum = mod(queue.tail - queue.head + queue.len, queue.len);
end