#ifndef CONFIG_H
#define CONFIG_H

#define BLAS_UNDERSCORE
#define LAPACK_UNDERSCORE
#define LARPACK(routine) LARPACK_ ## routine

#define LARPACK_CROSSOVER 8

#endif /* CONFIG_H */