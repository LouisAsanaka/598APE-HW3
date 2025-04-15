#include <gperftools/profiler.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <vector>

#ifndef PROFILE
#define PROFILE 0
#endif

float tdiff(struct timeval *start, struct timeval *end) {
    return (end->tv_sec - start->tv_sec) +
           1e-6 * (end->tv_usec - start->tv_usec);
}

struct Planet {
    double x;
    double y;
};

struct PlanetProp {
    double mass;
    double vx;
    double vy;
};

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

void inline __attribute__((always_inline))
next(std::vector<Planet> &planets, std::vector<Planet> &nextPlanets,
     std::vector<PlanetProp> &props,
     const std::vector<std::vector<double>> &masses6) {
    for (int i = 0; i < nplanets; i++) {
        nextPlanets[i].x = planets[i].x;
        nextPlanets[i].y = planets[i].y;
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
            props[i].vx += xInvDist3;
            props[i].vy += yInvDist3;
            props[j].vx -= xInvDist3;
            props[j].vy -= yInvDist3;
        }
        nextPlanets[i].x += dt * props[i].vx;
        nextPlanets[i].y += dt * props[i].vy;
    }
}

int main(int argc, const char **argv) {
    if (argc < 2) {
        printf("Usage: %s <nplanets> <timesteps>\n", argv[0]);
        return 1;
    }
    nplanets = atoi(argv[1]);
    timesteps = atoi(argv[2]);

    std::vector<Planet> planets(nplanets);
    std::vector<Planet> nextPlanets(nplanets);
    std::vector<PlanetProp> props(nplanets);
    std::vector<std::vector<double>> masses6(nplanets,
                                             std::vector<double>(nplanets));

    for (int i = 0; i < nplanets; i++) {
        props[i].mass = randomDouble() * 10 + 0.2;
        planets[i].x = (randomDouble() - 0.5) * 100 * pow(1 + nplanets, 0.4);
        planets[i].y = (randomDouble() - 0.5) * 100 * pow(1 + nplanets, 0.4);
        props[i].vx = randomDouble() * 5 - 2.5;
        props[i].vy = randomDouble() * 5 - 2.5;
    }

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
    for (int i = 0; i < timesteps; i++) {
        next(planets, nextPlanets, props, masses6);
        planets.swap(nextPlanets);
    }
#if PROFILE
    ProfilerStop();
#endif
    gettimeofday(&end, NULL);

    printf("Total time to run simulation %0.6f seconds\n", tdiff(&start, &end));
    for (int i = 0; i < nplanets; ++i) {
        printf("%f,%f\n", planets[i].x, planets[i].y);
    }
    return 0;
}
