def mapper(char):
	char = char.lower()
	if char in "еёъыародяьбю":
		return 1
	if char in "фв":
		return 2
	return 0
	
with open("text.txt") as f:
	data = "".join([line for line in f])
	
print(sum(map(mapper, data)))
