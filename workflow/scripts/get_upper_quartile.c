
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define DYN_ARRAY_INIT 50
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/*
#upper quartile computed on n. tpx for each position of the dbd
    # Only way is keeping count for each position of dbd.
    #upper quartile stability computed on all stability values of all tpx of the dbd
    # No need to track every position of the dbd, stability values can be kept in a single array,
    # but, they must be repeated for the number of positions, or the information must be tracked in some way
    #Es: (Value, Count). Sort by value, iterate once to get Sum(Count),iterate second time to get quartile :)
*/

typedef struct {
    int b; int e; float stability;
} TPX;
typedef struct {
    int b; int e;
} DBD;
typedef struct {
    float stability;
    int count;
} Stability;
typedef struct {
    int b; int e; int len;
    int *pos; 
    Stability *stability;
    int stab_array_len; int stab_array_max;
} DBD_PROCESSED;

int compare_ints(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
}
int compare_stability(const void *x, const void *y){
    Stability *a = (Stability*)x;
    Stability *b = (Stability*)y;
    return (a->stability - b->stability > 0) ? -1 : ((a->stability - b->stability < 0) ? 1 : 0);
}

DBD binary_search(int b, int e, DBD_PROCESSED *dbd_processed, int num_dbd){
    int first=-1; int last=-1;
    int L = 0; int H = num_dbd;
    while(L <= H){
        int M = (L+H)/2;
        if ((dbd_processed[M].b<=b && dbd_processed[M].e > b) || (dbd_processed[M].b<e && dbd_processed[M].e >= b)){
            first = M; last = M+1;
            while(first>=1 && ((dbd_processed[first-1].b<=b && dbd_processed[first-1].e > b) || (dbd_processed[first-1].b<e && dbd_processed[first-1].e >= b))){
                first--;
            }
            while(last<num_dbd && ((dbd_processed[last].b<=b && dbd_processed[last].e > b) || (dbd_processed[last].b<e && dbd_processed[last].e >= b))){
                last++;
            }
            break;
        } else if (dbd_processed[M].b >= e){
            H = M-1;
        } else {
            L = M+1;
        }
    }
    DBD to_return;
    to_return.b = first; to_return.e = last;
    return to_return;
}

void get_upper_quartile(TPX* tpx, int num_tpx, DBD* dbd, int num_dbd){
    //Initialize dbd_processed
    DBD_PROCESSED *dbd_processed = (DBD_PROCESSED*)malloc(num_dbd*sizeof(DBD_PROCESSED));
    for (int i=0; i<num_dbd; i++){
        dbd_processed[i].b = dbd[i].b;
        dbd_processed[i].e = dbd[i].e;
        dbd_processed[i].len = dbd[i].e-dbd[i].b;
        dbd_processed[i].pos = (int*)malloc((dbd[i].e-dbd[i].b)*sizeof(int));
        dbd_processed[i].stability = (Stability*)malloc((DYN_ARRAY_INIT)*sizeof(Stability));
        memset(dbd_processed[i].pos, 0, (dbd[i].e-dbd[i].b)*sizeof(int));
        dbd_processed[i].stab_array_len = 0;
        dbd_processed[i].stab_array_max = DYN_ARRAY_INIT;
    }
    //Iterate over tpx
    for (int i=0; i<num_tpx; i++){
        //Get range of arrays in the vector of DBDS
        DBD range = binary_search(tpx[i].b, tpx[i].e, dbd_processed, num_dbd);
        int b = range.b; int e = range.e;
        while(b<e){
            //DBD d = dbd_processed[b];
            int b_dbd = MAX(tpx[i].b - dbd_processed[b].b,0);
            int e_dbd = MIN(MAX(tpx[i].e -  dbd_processed[b].b,0),  dbd_processed[b].e- dbd_processed[b].b);
            if (dbd_processed[b].stab_array_len==dbd_processed[b].stab_array_max){
                dbd_processed[b].stab_array_max*=2;
                dbd_processed[b].stability = (Stability*)realloc(dbd_processed[b].stability, dbd_processed[b].stab_array_max*sizeof(Stability));
            }
            dbd_processed[b].stability[dbd_processed[b].stab_array_len].stability = tpx[i].stability;
            dbd_processed[b].stability[dbd_processed[b].stab_array_len].count = e_dbd-b_dbd;
            dbd_processed[b].stab_array_len++;
            while(b_dbd<e_dbd){
                dbd_processed[b].pos[b_dbd]++;
                b_dbd++;
            }
            b++;
        }
    }
    //Iterate over dbd_processed and compute quartiles
    for (int i=0; i<num_dbd; i++){
        int quartile=0; float stability_quartile=0;
        //Get quartile for pos
        if (dbd_processed[i].len>0){
            qsort(dbd_processed[i].pos, dbd_processed[i].len, sizeof(int), compare_ints);
            quartile = dbd_processed[i].pos[dbd_processed[i].len/4*3];
        }
        if ( dbd_processed[i].stab_array_len>0){
            //Get quartile for stability - more complicated, must sort both stability and stability_count
            qsort(dbd_processed[i].stability, dbd_processed[i].stab_array_len, sizeof(Stability), compare_stability);
            int count = 0;
            for (int j=0; j<dbd_processed[i].stab_array_len; j++){
                count+=dbd_processed[i].stability[j].count;
            }
            int threshold = count/4;
            count = 0;
            for (int j=0; j<dbd_processed[i].stab_array_len; j++){
                count += dbd_processed[i].stability[j].count;
                if (count>=threshold){
                    stability_quartile = dbd_processed[i].stability[j].stability;
                    break;
                }
            }
        }
        printf("%d;%f\n", quartile, stability_quartile);
    }
}

DBD * parse_dbd(char *dbd_file, int *num_dbd){
    FILE *file = fopen(dbd_file, "r");
    if (!file) {
        exit(1);
    }
    int size = 0; int i = 0;
    int col1, col2;
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        size++;
    }
    if (size == 0){
        return NULL;
    }
    DBD *dbd = (DBD*)malloc(size*sizeof(DBD));
    rewind(file);
    while (fgets(line, sizeof(line), file)) {
        // Parse the two integers separated by \t
        if (sscanf(line, "%*s\t%d\t%d", &col1, &col2) != 2) {
            fprintf(stderr, "Error parsing line: %s", line);
            continue;
        }
        dbd[i].b = col1; dbd[i].e = col2;
        i++;
    }
    fclose(file);
    *num_dbd = size;
    return dbd;
}
TPX * parse_tpx(char *tpx_file, int *num_tpx){
    FILE *file = fopen(tpx_file, "r");
    if (!file) {
        exit(1);
    }
    int size = 0; int i = 0;
    int col1, col2; float col3;
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        size++;
    }
    if (size == 0){
        return NULL;
    }
    TPX *tpx = (TPX*)malloc(size*sizeof(TPX));
    rewind(file);
    while (fgets(line, sizeof(line), file)) {
        // Parse the two integers separated by \t
        if (sscanf(line, "%d\t%d\t%f", &col1, &col2, &col3) != 3) {
            fprintf(stderr, "Error parsing line: %s", line);
            continue;
        }
        tpx[i].b = col1; tpx[i].e = col2; tpx[i].stability = col3;
        i++;
    }
    fclose(file);
    *num_tpx = size;
    return tpx;
}

int main(int argc, char *argv[]){
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <dbd_file> <tpx_file>\n", argv[0]);
        return EXIT_FAILURE;
    }
    char *dbd_file = argv[1];
    char *tpx_file = argv[2];
    int num_dbd, num_tpx;
    DBD *dbd = parse_dbd(dbd_file, &num_dbd);
    TPX *tpx = parse_tpx(tpx_file, &num_tpx);
    if (dbd == NULL || tpx == NULL){
        return 0;
    }
    get_upper_quartile(tpx, num_tpx, dbd, num_dbd);
    return 0;
}