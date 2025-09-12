CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LDFLAGS = -lm

symnmf: symnmf.o symnmf.h
	$(CC) -o symnmf symnmf.o $(CFLAGS) $(LDFLAGS)

main.o: symnmf.c
	$(CC) -c symnmf.c $(CFLAGS)
