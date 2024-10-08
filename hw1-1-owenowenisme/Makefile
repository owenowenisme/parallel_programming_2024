CC := mpicc
CFLAGS := -Wall
TARGET := hw1_1
v := 1

all: $(TARGET)

$(TARGET):f74116194_hw1_1.c
	$(CC) -o $@ $^

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f output/*.out
	rm -f $(TARGET)