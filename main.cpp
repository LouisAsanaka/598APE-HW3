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
    double mass;
    double x;
    double y;
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
next(std::vector<Planet> &planets, std::vector<Planet> &nextPlanets) {
    for (int i = 0; i < nplanets; i++) {
        nextPlanets[i].vx = planets[i].vx;
        nextPlanets[i].vy = planets[i].vy;
        nextPlanets[i].mass = planets[i].mass;
        nextPlanets[i].x = planets[i].x;
        nextPlanets[i].y = planets[i].y;
    }

    for (int i = 0; i < nplanets; i++) {
        for (int j = 0; j < nplanets; j++) {
            double dx = planets[j].x - planets[i].x;
            double dy = planets[j].y - planets[i].y;
            double distSqr = dx * dx + dy * dy + 0.0001;
            double invDist = planets[i].mass * planets[j].mass / sqrt(distSqr);
            double invDist3 = invDist * invDist * invDist;
            nextPlanets[i].vx += dt * dx * invDist3;
            nextPlanets[i].vy += dt * dy * invDist3;
        }
        nextPlanets[i].x += dt * nextPlanets[i].vx;
        nextPlanets[i].y += dt * nextPlanets[i].vy;
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

    for (int i = 0; i < nplanets; i++) {
        planets[i].mass = randomDouble() * 10 + 0.2;
        planets[i].x = (randomDouble() - 0.5) * 100 * pow(1 + nplanets, 0.4);
        planets[i].y = (randomDouble() - 0.5) * 100 * pow(1 + nplanets, 0.4);
        planets[i].vx = randomDouble() * 5 - 2.5;
        planets[i].vy = randomDouble() * 5 - 2.5;
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);
#if PROFILE
    ProfilerStart("my_profile.prof");
#endif
    for (int i = 0; i < timesteps; i++) {
        next(planets, nextPlanets);
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
