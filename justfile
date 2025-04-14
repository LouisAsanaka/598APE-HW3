alias b := build
alias c := clean

build profile='0':
    make -j PROFILE={{profile}}

clean:
    make clean

view-profile:
	~/go/bin/pprof -http "0.0.0.0:8080" ./main.exe ./my_profile.prof

test bodies timesteps: build
    ./test.sh {{bodies}} {{timesteps}}
