CC = g++
CFLAGS = -O3 -Wall -Wshadow -Winline -ansi -pedantic -g $(INCLUDES)
LIBS = 
INCLUDES =
TARGET = multigrid

SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

.PHONY : clean depend

clean: 
	@/bin/rm -f $(OBJS) 
	@/bin/rm -f $(TARGET)

depend: 
	@makedepend -- $(CFLAGS) -- $(SRC)