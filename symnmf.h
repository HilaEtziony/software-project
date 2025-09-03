struct cord {
    double value;
    struct cord *next;
};

struct vector {
    struct vector *next;
    struct cord *cords;
};

void symnmf_algo(int k, struct cord** H, struct vector* W);

void free_vectors(struct vector *head);

void free_cords(struct cord *head);
