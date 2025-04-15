FUNC := g++
FLAGS := -O3 -ffast-math -march=native -lm -g -Werror -lprofiler -fopenmp

# Enable profiling via a macro
ifndef PROFILE
	PROFILE := 0
endif
FLAGS += -DPROFILE=$(PROFILE)

all:
	$(FUNC) ./main.cpp -o ./main.exe $(FLAGS)

clean:
	rm -f ./*.exe