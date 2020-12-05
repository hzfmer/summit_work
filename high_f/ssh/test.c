#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>


int main(int argc, char *argv []) {
    int nprocs, rank;

    float std = 0.1;
    int nx = 19440, ny = 14904, nz = 1250;
    int px = 1, py = 1;
    int nz_ssh = nz, nz_tap = nz * 3 / 4;
    int nvar = 3;
    float vpmax = 8000., vsmax=5000.;
    char fh_ssh[1024], fh_hom[1024], fh_het[1024];


    int c = 0;
    static struct option long_options[] = {
        {"std",     required_argument,   0,  'a'},
        {"nx",      required_argument,   0,  'x'},
        {"ny",      required_argument,   0,  'y'},
        {"nz",      required_argument,   0,  'z'},
        {"nz_ssh",  required_argument,   0,  'b'},
        {"nz_tap",  optional_argument,   0,  'c'},
        {"vpmax",   required_argument,   0,  'p'},
        {"vsmax",   required_argument,   0,  's'},
        {"fh_ssh",  required_argument,   0,  'd'},
        {"fh_hom",  required_argument,   0,  'e'},
        {"fh_het",  required_argument,   0,  'f'},
        {"px",      required_argument,   0,  'g'},
        {"py",      required_argument,   0,  'i'},
        {0, 0, 0, 0}
    };
    while (1) {
        c = getopt_long_only(argc, argv, "a:x:y:z:b:c::p:s:d:e:f:g:i:",
                long_options, NULL);
        if (c == -1) break;
        switch (c) {
            case 'a':
                std = atof(optarg); break;
            case 'x':
                nx = atoi(optarg); break;
            case 'y':
                ny = atoi(optarg); break;
            case 'z':
                nz = atoi(optarg); break;
            case 'b':
                printf("nz_ssh: optind=%d, optarg=%s, argv[optind]=%s\n",
                            optind, optarg, argv[optind]);
                nz_ssh = atoi(optarg); break;
            case 'c':
                printf("nz_tap: optind=%d, optarg=%s, argv[optind]=%s\n", 
                            optind, optarg, argv[optind]);
                if (optarg == NULL && argv[optind] != NULL
                        && argv[optind][0] != '-') {
                    nz_tap = atoi(argv[optind]);
                    optind++;
                } else if (optarg != NULL)  {
                    nz_tap = atoi(optarg);
                } else {
                    nz_tap = nz_ssh + 1;
                }
                //printf("fh_hom=%s\n",fh_hom);
                break;
            case 'p':
                vpmax = atof(optarg); break;
            case 's':
                vsmax = atof(optarg); break;
            case 'd':
                strcpy(fh_ssh, optarg); break;
            case 'e':
                strcpy(fh_hom, optarg); break;
            case 'f':
                strcpy(fh_het, optarg); break;
            case 'g':
                px = atoi(optarg); break;
            case 'i':
                py = atoi(optarg); break;
            default:
                printf("Usage %s\nOptions:\n--std\t--nx\t--ny\t--nz\t--nz_ssh\t--nz_tap\t--vpmax\t--vsmax\n", argv[0]);
                printf("nz_tap controls how deep (# of layers) to start tapering, default is no tapering.\n");
                abort();
        }
    }

    if (optind < argc) {
        printf("Non-option ARGV-elements: ");
        while (optind < argc) 
          printf("%s\t", argv[optind++]);
        printf("\n");
    }

    
    int size[3] = {nx, ny, nz};
    int part[3] = {px, py};
    int err = 0;
    char mpiErrStr[1024];
    int mpiErrStrLen;
        
    fprintf(stdout, "nx=%d, ny=%d, nz=%d, std=%f, nz_ssh=%d, nz_tap=%d, vpmax=%f, vsmax=%f\n", nx, ny, nz, std, nz_ssh, nz_tap, vpmax, vsmax);
    fprintf(stdout, "fh_ssh=%s\n, fh_hom=%s\n, fh_het=%s\n", fh_ssh, fh_hom, fh_het);
    fflush(stdout);
    return err;
}

