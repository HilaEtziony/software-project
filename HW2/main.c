#include <stdio.h>

#include <stdio.h>

void recurse(int n) {
    static int count = 0;
    count++;
    printf("n=%d, count=%d\n", n, count);
    if (n > 0)
        recurse(n - 1);
}

int main() {
    recurse(3);
    return 0;
}
