CPPFLAGS=
LDFLAGS=-lstdc++ -lm
ROOTFLAGS=$(shell root-config --cflags)
ROOTLDFLAGS=$(shell root-config --glibs)

spots2Dij: spots2Dij.cpp
	g++ -O3 -o $@ $^ ${LDFLAGS}

readDij: readDij.cpp
	g++ -O3 ${CPPFLAGS} -o $@ $^ ${LDFLAGS}

phsp2root: phsp2root.cpp
	g++ -O3 ${CPPFLAGS} ${ROOTFLAGS} -o $@ $^ ${LDFLAGS} ${ROOTLDFLAGS}

combineTopasPhsp: combineTopasPhsp.cpp
	g++ ${CPPFLAGS} -o $@ $^ ${LDFLAGS}

comparePhsp: comparePhsp.cpp
	clang++ -O3 ${CPPFLAGS} ${ROOTFLAGS} -o $@ $^ ${LDFLAGS} ${ROOTLDFLAGS}

gpmc2xio: gpmc2xio.cpp
	g++ -O3 ${CPPFLAGS} -o $@ $^ ${LDFLAGS}

gpmc2xioLET: gpmc2xioLET.cpp
	g++ -O3 ${CPPFLAGS} -o $@ $^ ${LDFLAGS}

gpmc2xioLETd: gpmc2xioLETd.cpp
	g++ ${CPPFLAGS} -o $@ $^ ${LDFLAGS}

sumCubes: sumCubes.cpp
	g++ ${CPPFLAGS} -o $@ $^ ${LDFLAGS}

applyThreshold: applyThreshold.cpp
	g++ -O3 -o $@ $^ ${LDFLAGS}

let2letd: let2letd.cpp
	g++ -O3 ${CPPFLAGS} -o $@ $^ ${LDFLAGS}

# clean:
# 	rm -f phsp2root
# 	# rm -f geometrygPMC2XiO
# 	# rm -f pdg2fakeIAEA
# 	# rm -f normSumLET

.PHONY: clean
