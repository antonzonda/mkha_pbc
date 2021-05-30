# makefile
FILES = test.c mkha.c poly.c prf.c tag.c
CFLAGS = -Wpedantic -Wall -Werror -Wextra -I /usr/local/include/flint -lgmp -lflint -L. -lpbc -lrelic -lrelic_s -g
CC = gcc

test: 
	$(CC) $(FILES) $(CFLAGS) -o test
.PHONY : clean
 clean:
	rm -f test