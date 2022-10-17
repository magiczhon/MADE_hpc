typedef enum { false, true } bool;

void seedThreads();
bool circle(float x, float y, float r);
float parallel_compute_pi(int *dots_in_circle, const int all_dots);
float compute_pi(int *dots_in_circle, const int all_dots);
