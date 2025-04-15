#include <gperftools/profiler.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <thread>
#include <vector>
#include <cstring>

#include <omp.h>

#ifndef PROFILE
#define PROFILE 0
#endif

float tdiff(struct timeval *start, struct timeval *end) {
    return (end->tv_sec - start->tv_sec) +
           1e-6 * (end->tv_usec - start->tv_usec);
}

struct alignas(64) Planet {
    double x;
    double y;
    double vx;
    double vy;
};

struct PlanetProp {
    double mass;
};

using Double2d = std::vector<std::vector<double>>;

unsigned long long seed = 100;

unsigned long long randomU64() {
    seed ^= (seed << 21);
    seed ^= (seed >> 35);
    seed ^= (seed << 4);
    return seed;
}

double randomDouble() {
    unsigned long long next = randomU64();
    next >>= (64 - 26);
    unsigned long long next2 = randomU64();
    next2 >>= (64 - 26);
    return ((next << 27) + next2) / (double)(1LL << 53);
}

int nplanets;
int timesteps;
constexpr double dt = 0.001;

int main(int argc, const char **argv) {
    if (argc < 2) {
        printf("Usage: %s <nplanets> <timesteps>\n", argv[0]);
        return 1;
    }
    nplanets = atoi(argv[1]);
    timesteps = atoi(argv[2]);

    std::vector<Planet> planets(nplanets);
    std::vector<PlanetProp> props(nplanets);

    for (int i = 0; i < nplanets; i++) {
        props[i].mass = randomDouble() * 10 + 0.2;
        planets[i].x = (randomDouble() - 0.5) * 100 * pow(1 + nplanets, 0.4);
        planets[i].y = (randomDouble() - 0.5) * 100 * pow(1 + nplanets, 0.4);
        planets[i].vx = randomDouble() * 5 - 2.5;
        planets[i].vy = randomDouble() * 5 - 2.5;
    }

    Double2d masses6(nplanets, std::vector<double>(nplanets));
    for (int i = 0; i < nplanets; i++) {
        for (int j = 0; j < nplanets; j++) {
            double mass2 = props[i].mass * props[j].mass;
            masses6[i][j] = mass2 * mass2 * mass2;
        }
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);
#if PROFILE
    ProfilerStart("my_profile.prof");
#endif

    int nthreads = 1;
    if (nplanets >= 500) {
        // Parallel algorithm inspired by
        // https://www.cs.usask.ca/~spiteri/CMPT851/notes/nBody.pdf
        int nthreads = std::thread::hardware_concurrency();
        double *vx = (double *)malloc(nthreads * nplanets * sizeof(double));
        double *vy = (double *)malloc(nthreads * nplanets * sizeof(double));
        for (int i = 0; i < timesteps; i++) {
            #pragma omp parallel
            {
                int offset = omp_get_thread_num() * nplanets;
                std::memset(vx + offset, 0, nplanets * sizeof(double));
                std::memset(vy + offset, 0, nplanets * sizeof(double));
            }

            #pragma omp parallel
            {
                int base = omp_get_thread_num() * nplanets;

                #pragma omp for schedule(dynamic)
                for (int i = 0; i < nplanets; i++) {
                    for (int j = i + 1; j < nplanets; j++) {
                        double dx = planets[j].x - planets[i].x;
                        double dy = planets[j].y - planets[i].y;
                        double distSqr = dx * dx + dy * dy + 0.0001;
                        double distSqrt = sqrt(distSqr);
                        double invDist3 =
                            dt * masses6[i][j] / (distSqr * distSqrt);
                        double xInvDist3 = dx * invDist3;
                        double yInvDist3 = dy * invDist3;
                        vx[base + i] += xInvDist3;
                        vy[base + i] += yInvDist3;
                        vx[base + j] -= xInvDist3;
                        vy[base + j] -= yInvDist3;
                    }
                }
            }

            #pragma omp parallel for schedule(static)
            for (int i = 0; i < nplanets; ++i) {
                double vxCum = 0.0;
                double vyCum = 0.0;
                for (int base = 0; base < nthreads * nplanets;
                     base += nplanets) {
                    vxCum += vx[base + i];
                    vyCum += vy[base + i];
                }
                planets[i].vx += vxCum;
                planets[i].vy += vyCum;
                planets[i].x += dt * planets[i].vx;
                planets[i].y += dt * planets[i].vy;
            }
        }
    } else {
        std::vector<Planet> nextPlanets(nplanets);

        for (int i = 0; i < timesteps; i++) {
            for (int i = 0; i < nplanets; i++) {
                nextPlanets[i].x = planets[i].x;
                nextPlanets[i].y = planets[i].y;
                nextPlanets[i].vx = planets[i].vx;
                nextPlanets[i].vy = planets[i].vy;
            }

            for (int i = 0; i < nplanets; i++) {
                for (int j = i + 1; j < nplanets; j++) {
                    double dx = planets[j].x - planets[i].x;
                    double dy = planets[j].y - planets[i].y;
                    double distSqr = dx * dx + dy * dy + 0.0001;
                    double distSqrt = sqrt(distSqr);
                    double invDist3 = dt * masses6[i][j] / (distSqr * distSqrt);
                    double xInvDist3 = dx * invDist3;
                    double yInvDist3 = dy * invDist3;
                    nextPlanets[i].vx += xInvDist3;
                    nextPlanets[i].vy += yInvDist3;
                    nextPlanets[j].vx -= xInvDist3;
                    nextPlanets[j].vy -= yInvDist3;
                }
                nextPlanets[i].x += dt * nextPlanets[i].vx;
                nextPlanets[i].y += dt * nextPlanets[i].vy;
            }

            nextPlanets.swap(planets);
        }
    }
#if PROFILE
    ProfilerStop();
#endif
    gettimeofday(&end, NULL);

    printf("Total time to run simulation %0.6f seconds with %d threads\n",
           tdiff(&start, &end), nthreads);
    for (int i = 0; i < nplanets; ++i) {
        printf("%f,%f\n", planets[i].x, planets[i].y);
    }
    return 0;
}
