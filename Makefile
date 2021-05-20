# makefile
FILES = test.c mkha.c poly.c prf.c tag.c
CFLAGS = -I /usr/local/include/flint -lgmp -lflint -L. -lpbc 
CC = gcc

test: 
	$(CC) $(FILES) $(CFLAGS) -o test
.PHONY : clean
 clean:
	rm -f test