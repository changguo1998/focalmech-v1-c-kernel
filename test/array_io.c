#include <stdio.h>
#include <stdlib.h>

int main() {
    float* arr = NULL;
    arr = (float*)malloc(5 * sizeof(float));
    FILE *fp = fopen("arrayio.bin", "r");
    fread(arr, sizeof(float), 5, fp);
    fclose(fp);
    for (int i = 0; i < 5; i++)
        printf("arr[%d] = %f\n", i, arr[i]);
    free(arr);
    return 0;
}
