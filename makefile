CC = g++
CFLAGS = -O3 -Wall -Winline -ansi -pedantic $(INCLUDES)
LIBS = 
INCLUDES =
TARGET = multigrid

SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

.PHONY : clean depend

clean: 
	@/bin/rm -f $(OBJS) 
	@/bin/rm -f $(TARGET)

depend: 
	@makedepend -- $(CFLAGS) -- $(SRC)