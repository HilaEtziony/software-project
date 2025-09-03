CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors

symnmf: symnmf.o symnmf.h
    $(CC) -o symnmf symnmf.o $(CFLAGS)

main.o: symnmf.c
    $(CC) -c symnmf.c $(CFLAGS)
